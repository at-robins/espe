use std::{
    ffi::OsString,
    hash::{Hash, Hasher},
    io::{BufWriter, Write},
    path::Path,
    process::{Child, Command, Output, Stdio},
};

use actix_web::web;
use chrono::NaiveDateTime;
use diesel::{ExpressionMethods, QueryDsl, RunQueryDsl};
use twox_hash::XxHash64;

use crate::{
    application::{
        config::Configuration,
        database::DatabaseManager,
        error::{SeqError, SeqErrorType},
    },
    model::{
        db::experiment_execution::{ExecutionStatus, ExperimentExecution},
        exchange::experiment_pipeline::{
            ExperimentPipelineBlueprint, ExperimentPipelineStepBlueprint,
        },
        internal::pipeline_blueprint::PipelineStepBlueprint,
    },
};

use super::pipeline_service::LoadedPipelines;

/// Builds the specifc pipeline step container at the specified context.
///
/// # Parameters
///
/// * `step` - the [`PipelineStepBlueprint`] to build the container for
/// * `context` - the context directory contianing the pipeline
pub fn build_pipeline_step<P: AsRef<Path>, T: AsRef<str>>(
    step: &PipelineStepBlueprint,
    pipeline_id: T,
    context: P,
) -> Result<Child, SeqError> {
    let mut pipeline_step_path = context.as_ref().to_path_buf();
    pipeline_step_path.push("container");
    pipeline_step_path.push(step.container());
    let build_arg: OsString = "build".into();
    let name_spec: OsString = "-t".into();
    let name_arg: OsString = format_container_name(pipeline_id, step.id()).into();

    let child = Command::new("docker")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .args([
            build_arg.as_os_str(),
            name_spec.as_os_str(),
            name_arg.as_os_str(),
            pipeline_step_path.as_os_str(),
        ])
        .spawn()?;
    Ok(child)
}

/// Runs the specifc pipeline step replacing all its previous output.
///
/// # Parameters
///
/// * `pipeline_id` - the ID of the containing [`PipelineBlueprint`]
/// * `step` - the [`PipelineStepBlueprint`] to run
/// * `variables` - the variables that were specified for step execution
/// * `experiment_id` - the ID of the experiment
/// * `app_cofig` - the app [`Configuration`]
pub fn run_pipeline_step<T: AsRef<str>>(
    pipeline_id: T,
    step: &ExperimentPipelineStepBlueprint,
    experiment_id: i32,
    app_config: web::Data<Configuration>,
) -> Result<Child, SeqError> {
    let experiment_id = experiment_id.to_string();
    // Create the output directory.
    let output_path = app_config.experiment_step_path(&experiment_id, step.id());
    // Clear the output folder if the step has been run before.
    if output_path.exists() {
        std::fs::remove_dir_all(&output_path)?;
    }
    // Then create the output directory.
    std::fs::create_dir_all(&output_path)?;
    // Set basic arguments.
    let mut arguments: Vec<OsString> = vec![
        "run".into(),
        "--name".into(),
        format_container_name(&pipeline_id, step.id()).into(),
        "--rm".into(),
    ];
    arguments.extend(pipeline_step_mount(output_path, "/output", false));
    // Set initial sample input mount.
    arguments.extend(pipeline_step_mount(
        app_config.experiment_input_path(&experiment_id),
        "/input/samples",
        true,
    ));
    // Set input mounts / dependencies.
    for dependency_id in step.dependencies() {
        arguments.extend(pipeline_step_mount(
            app_config.experiment_step_path(&experiment_id, dependency_id),
            format!("/input/steps/{}", dependency_id),
            true,
        ));
    }
    // Set global mounts.
    step.variables()
        .iter()
        .filter(|var_instance| var_instance.is_global_data_reference())
        // Filter out variables without values.
        .filter_map(|var_instance| var_instance.value().as_ref().map(|value| (var_instance.id(), value)))
        .for_each(|(global_var_id, global_var_value)| {
            arguments.extend(pipeline_step_mount(
                app_config.global_data_path(global_var_value),
                format!("/input/globals/{}", global_var_id),
                true,
            ));
        });

    // Set other variables.
    step.variables()
        .iter()
        .filter(|var_instance| !var_instance.is_global_data_reference())
        // Filter out variables without values.
        .filter_map(|var_instance| var_instance.value().as_ref().map(|value| (var_instance.id(), value)))
        .for_each(|(other_var_id, other_var_value)| {
            arguments.push("--env".into());
            arguments.push(format!("{}={}", other_var_id, other_var_value).into());
        });

    // Set container to run.
    arguments.push(format_container_name(pipeline_id, step.id()).into());

    let output = Command::new("docker")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .args(arguments)
        .spawn()?;
    Ok(output)
}

/// Creates an [`OsString`] for a container mount.
///
/// # Parameters
///
/// * `source` - the path to the local directory
/// * `target` - the path to the bound container directory
fn pipeline_step_mount<P: AsRef<Path>, Q: AsRef<Path>>(
    source: P,
    target: Q,
    readonly: bool,
) -> Vec<OsString> {
    let mut mount_args: Vec<OsString> = vec!["--mount".into()];
    let mut mount_arg: OsString = "type=bind,source=".into();
    mount_arg.push(source.as_ref());
    mount_arg.push(OsString::from(",target="));
    mount_arg.push(target.as_ref());
    if readonly {
        mount_arg.push(OsString::from(",readonly"));
    }
    mount_args.push(mount_arg);
    mount_args
}

/// Converts the specified name to a valid container name.
///
/// # Parameters
///
/// * `name` - the name to convert
fn format_container_name<T: AsRef<str>, R: AsRef<str>>(
    pipeline_id: T,
    pipeline_step_id: R,
) -> String {
    let name = format!("{}{}", pipeline_id.as_ref(), pipeline_step_id.as_ref());
    let mut hasher = XxHash64::with_seed(154);
    name.hash(&mut hasher);
    hasher.finish().to_string()
}

pub struct ContainerHandler {
    config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    loaded_pipelines: web::Data<LoadedPipelines>,
    executed_step: Option<ExperimentExecution>,
    build_process: Option<Child>,
    run_process: Option<Child>,
}

impl ContainerHandler {
    /// Creates a new container handler with the specified [`Configuration`].
    ///
    /// # Parameters
    ///
    /// * `config` - the application's configuration
    /// * `loaded_pipelines` - the pipelines loaded by the application
    pub fn new(
        config: web::Data<Configuration>,
        database_manager: web::Data<DatabaseManager>,
        loaded_pipelines: web::Data<LoadedPipelines>,
    ) -> Self {
        Self {
            config,
            database_manager,
            loaded_pipelines,
            executed_step: None,
            build_process: None,
            run_process: None,
        }
    }

    /// Returns ```true``` if the handler is currently executing a pipeline step.
    pub fn is_running(&self) -> bool {
        self.executed_step.is_some()
    }

    /// Kills all processes attached to this container handler.
    fn kill(&mut self) -> Result<(), SeqError> {
        if self.is_running() {
            Self::kill_process(&self.executed_step, &mut self.build_process)?;
            Self::kill_process(&self.executed_step, &mut self.run_process)?;
        }
        self.reset();
        Ok(())
    }

    /// Kills the specified process.
    ///
    /// # Parameters
    ///
    /// * `step` - the corresponding execution step
    /// * `process` - the process to kill
    fn kill_process(
        step: &Option<ExperimentExecution>,
        process: &mut Option<Child>,
    ) -> Result<(), SeqError> {
        if let Some(ref mut current_process) = process {
            if let Err(error) = current_process.kill() {
                return Err(SeqError::new(
                    "Unkillable process",
                    SeqErrorType::InternalServerError,
                    format!(
                        "Cannot kill {:?} with process ID {}: {}",
                        step,
                        current_process.id(),
                        error
                    ),
                    "Process cannot be terminated.",
                ));
            }
        }
        Ok(())
    }

    /// Resets the internal state of the handler.
    fn reset(&mut self) {
        self.build_process = None;
        self.run_process = None;
        self.executed_step = None;
    }

    /// Starts the specified execution.
    ///
    /// # Parameters
    ///
    /// * `step` - the [`ExperimentExecution`] step to start
    pub fn start(&mut self, step: ExperimentExecution) -> Result<(), SeqError> {
        if self.is_running() {
            return Err(SeqError::new(
                "Process already running",
                SeqErrorType::InternalServerError,
                format!("Cannot start {:?} while is still running {:?}", step, self.executed_step),
                "Process is already running.",
            ));
        }
        self.reset();
        let mut connection = self.database_manager.database_connection()?;
        connection.immediate_transaction(|connection| {
            let running_status: String = ExecutionStatus::Running.into();
            let end_time: Option<NaiveDateTime> = None;
            diesel::update(crate::schema::experiment_execution::table.find(step.id))
                .set((
                    crate::schema::experiment_execution::execution_status.eq(running_status),
                    crate::schema::experiment_execution::start_time
                        .eq(Some(chrono::Utc::now().naive_local())),
                    crate::schema::experiment_execution::end_time.eq(end_time),
                ))
                .execute(connection)
        })?;
        log::info!("Starting {:?}", &step);
        self.executed_step = Some(step);
        Ok(())
    }

    /// Fails the current pipeline step and all other pipeline steps belonging to the experiment.
    /// Returns an error if no step is currently running.
    pub fn fail(&mut self) -> Result<(), SeqError> {
        if !self.is_running() {
            return Err(SeqError::new(
                "No process running",
                SeqErrorType::InternalServerError,
                "No process is running and can thus not be failed.",
                "No process is currently running.",
            ));
        }
        if let Some(step) = &self.executed_step {
            let mut connection = self.database_manager.database_connection()?;
            let experiment_id = step.experiment_id;
            self.kill()?;
            connection.immediate_transaction(|connection| {
                ExperimentExecution::update_scheduled_status_by_experiment(
                    experiment_id,
                    ExecutionStatus::Failed,
                    connection,
                )
            })?;
        } else {
            return Err(SeqError::new(
                "No process running",
                SeqErrorType::InternalServerError,
                "No process is specified and can thus not be failed.",
                "No process is currently running.",
            ));
        }
        Ok(())
    }

    /// Aborts the the pipeline step belonging to the specified experiment
    /// if currently executed.
    /// 
    /// # Parameters
    /// 
    /// * `experiment_id` - the ID of the experiment to abort 
    pub fn abort(&mut self, experiment_id: i32) -> Result<(), SeqError> {
        if let Some(step) = &self.executed_step {
            if experiment_id == step.experiment_id {
                self.kill()?;
            }
        }
        Ok(())
    }

    /// Returns the [`ProcessStatus`] of the specified process.
    ///
    /// # Parameters
    ///
    /// * `process` - the process to get the status for
    fn get_process_status(process: &mut Option<Child>) -> Result<ProcessStatus, SeqError> {
        if let Some(process) = process {
            if process.try_wait()?.is_some() {
                Ok(ProcessStatus::Finished)
            } else {
                Ok(ProcessStatus::Running)
            }
        } else {
            Ok(ProcessStatus::NotStarted)
        }
    }

    /// Returns the current [`ExperimentExecution`] or an error if none is set.
    fn get_executed_step(&self) -> Result<&ExperimentExecution, SeqError> {
        (&self.executed_step).as_ref().ok_or_else(|| {
            SeqError::new(
                "Invalid execution state",
                SeqErrorType::InternalServerError,
                "Expected an execution step to be set.",
                "Execution state is invalid.",
            )
        })
    }

    /// Logs the output and returns an error if the exit status was unsuccessful.
    ///
    /// # Parameters
    ///
    /// * `output` - the output to parse and log
    /// * `build` - ```true``` if the build output is parsed, ```false``` if the run output is parsed
    fn parse_output(&self, output: Output, build: bool) -> Result<(), SeqError> {
        if let Some(step) = &self.executed_step {
            let mut log_path = self
                .config
                .experiment_logs_path(step.experiment_id.to_string());
            let process_type = if build { "build" } else { "run" };
            std::fs::create_dir_all(&log_path)?;
            log_path.push(format!(
                "{}_{}.log",
                format_container_name(&step.pipeline_id, &step.pipeline_step_id),
                &process_type
            ));
            let log_file = std::fs::OpenOptions::new()
                .create(true)
                .write(true)
                .append(false)
                .truncate(true)
                .open(log_path)?;
            let mut buffered_writer = BufWriter::new(log_file);
            buffered_writer.write_all("[[ STDOUT ]]\n".as_bytes())?;
            buffered_writer.write_all(&output.stdout)?;
            buffered_writer.write_all("\n\n[[ STDERR ]]\n".as_bytes())?;
            buffered_writer.write_all(&output.stderr)?;
            buffered_writer.write_all("\n\n[[ EXIT STATUS ]]\n".as_bytes())?;
            buffered_writer.write_all(output.status.to_string().as_bytes())?;
            if output.status.success() {
                Ok(())
            } else {
                Err(SeqError::new(
                    "Non-zero exit status",
                    SeqErrorType::InternalServerError,
                    format!(
                        "The {} execution step {:?} failed with exit status {:?}.",
                        process_type, &self.executed_step, output.status
                    ),
                    "The container exit status was not zero.",
                ))
            }
        } else {
            Err(SeqError::new(
                "Invalid execution state",
                SeqErrorType::InternalServerError,
                "Cannot parse an output if no execution step is set.",
                "Execution state is invalid.",
            ))
        }
    }

    /// Updates the handler and returns ```true``` if ready for processing another step.
    pub fn update(&mut self) -> Result<bool, SeqError> {
        match self.update_inner() {
            Ok(status) => Ok(status),
            Err(err) => {
                self.fail()?;
                Err(err)
            },
        }
    }

    /// Updates the handler and returns ```true``` if ready for processing another step.
    fn update_inner(&mut self) -> Result<bool, SeqError> {
        // Ready to handle new containers if none is actively running.
        if !self.is_running() {
            return Ok(true);
        }
        let step_db_id: i32 = self.get_executed_step()?.id;
        // Check for run processes.
        match Self::get_process_status(&mut self.run_process)? {
            ProcessStatus::Finished => {
                if let Some(run) = self.run_process.take() {
                    // Handle output.
                    self.parse_output(run.wait_with_output()?, false)?;
                    // Sets the status to finished.
                    let mut connection = self.database_manager.database_connection()?;
                    connection.immediate_transaction(|connection| {
                        let finished_status: String = ExecutionStatus::Finished.into();
                        diesel::update(crate::schema::experiment_execution::table.find(step_db_id))
                            .set((
                                crate::schema::experiment_execution::execution_status
                                    .eq(finished_status),
                                crate::schema::experiment_execution::end_time
                                    .eq(Some(chrono::Utc::now().naive_local())),
                            ))
                            .execute(connection)
                    })?;
                    log::info!("Finished {:?}", &self.executed_step);
                    // Reset the internal state if finished successfully.
                    self.reset();
                    Ok(true)
                } else {
                    Err(SeqError::new(
                        "Invalid container state",
                        SeqErrorType::InternalServerError,
                        format!("The run process status of {:?} indicated a finsihed process, but none was found.", &self.executed_step),
                        "The execution state is invalid.",
                    ))
                }
            },
            // Checks for build processes if the run process has not yet been started.
            ProcessStatus::NotStarted => match Self::get_process_status(&mut self.build_process)? {
                ProcessStatus::Finished => {
                    if let Some(build) = self.build_process.take() {
                        // Handle output.
                        self.parse_output(build.wait_with_output()?, true)?;
                        // Start the subsequent run process.
                        self.start_run_process()?;
                        Ok(false)
                    } else {
                        Err(SeqError::new(
                            "Invalid container state",
                            SeqErrorType::InternalServerError,
                            format!("The build process status of {:?} indicated a finsihed process, but none was found.", &self.executed_step),
                            "The execution state is invalid.",
                        ))
                    }
                },
                // Starts the build process if not yet done.
                ProcessStatus::NotStarted => {
                    self.start_build_process()?;
                    Ok(false)
                },
                ProcessStatus::Running => Ok(false),
            },
            ProcessStatus::Running => Ok(false),
        }
    }

    /// Starts the build process.
    fn start_build_process(&mut self) -> Result<(), SeqError> {
        let step = self.get_executed_step()?;
        if let Some(pipeline) = self.loaded_pipelines.get(&step.pipeline_id) {
            if let Some(step_blueprint) = pipeline
                .pipeline()
                .steps()
                .iter()
                .find(|s| s.id() == &step.pipeline_step_id)
            {
                log::info!("Building {:?}", &step);
                self.build_process = Some(build_pipeline_step(
                    step_blueprint,
                    &step.pipeline_id,
                    pipeline.context(),
                )?);
                Ok(())
            } else {
                Err(SeqError::new(
                    "Invalid container state",
                    SeqErrorType::InternalServerError,
                    format!("The exectued pipeline step {:?} is not loaded.", &self.executed_step),
                    "The pipeline step is not defined.",
                ))
            }
        } else {
            Err(SeqError::new(
                "Invalid container state",
                SeqErrorType::InternalServerError,
                format!(
                    "The pipeline including exectued pipeline step {:?} is not loaded.",
                    &self.executed_step
                ),
                "The pipeline is not defined.",
            ))
        }
    }

    /// Starts the run process.
    fn start_run_process(&mut self) -> Result<(), SeqError> {
        let step = self.get_executed_step()?;
        if let Some(pipeline) = self.loaded_pipelines.get(&step.pipeline_id) {
            let mut connection = self.database_manager.database_connection()?;
            let values = crate::model::db::pipeline_step_variable::PipelineStepVariable::get_values_by_experiment_and_pipeline(step.experiment_id, pipeline.pipeline().id(), &mut connection)?;
            let pipeline = ExperimentPipelineBlueprint::from_internal(pipeline.pipeline(), values);
            if let Some(step_blueprint) = pipeline
                .steps()
                .iter()
                .find(|s| s.id() == &step.pipeline_step_id)
            {
                log::info!("Running {:?}", &step);
                self.run_process = Some(run_pipeline_step(
                    &step.pipeline_id,
                    step_blueprint,
                    step.experiment_id,
                    web::Data::clone(&self.config),
                )?);
                Ok(())
            } else {
                Err(SeqError::new(
                    "Invalid container state",
                    SeqErrorType::InternalServerError,
                    format!("The exectued pipeline step {:?} is not loaded.", &self.executed_step),
                    "The pipeline step is not defined.",
                ))
            }
        } else {
            Err(SeqError::new(
                "Invalid container state",
                SeqErrorType::InternalServerError,
                format!(
                    "The pipeline including exectued pipeline step {:?} is not loaded.",
                    &self.executed_step
                ),
                "The pipeline is not defined.",
            ))
        }
    }
}

/// The status of a container process.
enum ProcessStatus {
    Finished,
    NotStarted,
    Running,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pipeline_step_mount() {
        let input_arg = pipeline_step_mount(
            "/app/context/experiments/experiment_id/steps/step_id",
            "/input/steps/step_id",
            true,
        );
        let correct_input_arg: Vec<OsString> =
            vec!["--mount".into(), "type=bind,source=/app/context/experiments/experiment_id/steps/step_id,target=/input/steps/step_id,readonly"
                .into()];
        assert_eq!(input_arg, correct_input_arg);
        let output_arg = pipeline_step_mount(
            "/app/context/experiments/experiment_id/steps/step_id",
            "/output",
            false,
        );
        let correct_output_arg: Vec<OsString> = vec![
            "--mount".into(),
            "type=bind,source=/app/context/experiments/experiment_id/steps/step_id,target=/output"
                .into(),
        ];
        assert_eq!(output_arg, correct_output_arg);
    }
}
