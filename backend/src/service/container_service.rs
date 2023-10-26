use std::{
    collections::HashMap,
    ffi::OsString,
    io::{BufWriter, Write},
    path::Path,
    process::{Child, Command, Output},
};

use actix_web::web;
use chrono::NaiveDateTime;
use diesel::{ExpressionMethods, QueryDsl, RunQueryDsl};

use crate::{
    application::{
        config::{Configuration, LogOutputType, LogProcessType},
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

/// The container environment variable specifying all mounts.
const CONTAINER_ENV_MOUNT: &str = "MOUNT_PATHS";

/// Builds the specifc pipeline step container at the specified context.
///
/// # Parameters
///
/// * `step` - the [`PipelineStepBlueprint`] to build the container for
/// * `pipeline_id` - the ID of the containing [`PipelineBlueprint`]
/// * `context` - the context directory contianing the pipeline
/// * `experiment_id` - the ID of the experiment
/// * `app_cofig` - the app [`Configuration`]
pub fn build_pipeline_step<P: AsRef<Path>, T: AsRef<str>>(
    step: &PipelineStepBlueprint,
    pipeline_id: T,
    context: P,
    experiment_id: i32,
    app_config: web::Data<Configuration>,
) -> Result<Child, SeqError> {
    let mut pipeline_step_path = context.as_ref().to_path_buf();
    pipeline_step_path.push("container");
    pipeline_step_path.push(step.container());
    let build_arg: OsString = "build".into();
    let name_spec: OsString = "-t".into();
    let name_arg: OsString = format_container_name(&pipeline_id, step.id()).into();
    let progress_arg: OsString = "--progress=plain".into();

    // Create log directory.
    let logs_path = app_config.experiment_logs_path(experiment_id.to_string());
    std::fs::create_dir_all(&logs_path)?;
    // Open stdout log file.
    let log_path_stdout = app_config.experiment_log_path(
        experiment_id.to_string(),
        &pipeline_id,
        step.id(),
        LogProcessType::Build,
        LogOutputType::StdOut,
    );
    let log_file_stdout = std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .append(false)
        .truncate(true)
        .open(log_path_stdout)?;
    // Open stderr log file.
    let log_path_stderr = app_config.experiment_log_path(
        experiment_id.to_string(),
        &pipeline_id,
        step.id(),
        LogProcessType::Build,
        LogOutputType::StdErr,
    );
    let log_file_stderr = std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .append(false)
        .truncate(true)
        .open(log_path_stderr)?;

    let child = Command::new("docker")
        .stdout(log_file_stdout)
        .stderr(log_file_stderr)
        .args([
            build_arg.as_os_str(),
            name_spec.as_os_str(),
            name_arg.as_os_str(),
            progress_arg.as_os_str(),
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
        "/input/base",
        true,
    ));
    // Set input mounts / dependencies.
    let mut mount_map_dependencies = serde_json::Map::new();
    for dependency_id in step.dependencies() {
        let target = format!("/input/steps/{}", Configuration::hash_string(dependency_id));
        mount_map_dependencies
            .insert(dependency_id.to_string(), serde_json::Value::String(target.clone()));
        arguments.extend(pipeline_step_mount(
            app_config.experiment_step_path(&experiment_id, dependency_id),
            target,
            true,
        ));
    }
    // Set global mounts.
    let mut mount_map_globals = serde_json::Map::new();
    step.variables()
        .iter()
        .filter(|var_instance| var_instance.is_global_data_reference())
        // Filter out variables without values.
        .filter_map(|var_instance| var_instance.value().as_ref().map(|value| (var_instance.id(), value)))
        .for_each(|(global_var_id, global_var_value)| {
            let target = format!("/input/globals/{}", Configuration::hash_string(global_var_id));
            mount_map_globals.insert(
                global_var_id.to_string(),
                serde_json::Value::String(target.clone()),
            );
            arguments.extend(pipeline_step_mount(
                app_config.global_data_path(global_var_value),
                target,
                true,
            ));
        });

    // Set mount envrionment variable.
    let mut mount_map_top = serde_json::Map::new();
    mount_map_top.insert("input".to_string(), serde_json::Value::String("/input/base".to_string()));
    mount_map_top.insert("output".to_string(), serde_json::Value::String("/output".to_string()));
    mount_map_top
        .insert("dependencies".to_string(), serde_json::Value::Object(mount_map_dependencies));
    mount_map_top.insert("globals".to_string(), serde_json::Value::Object(mount_map_globals));
    let mount_paths = serde_json::Value::Object(mount_map_top);
    arguments.push("--env".into());
    arguments.push(format!("{}={}", CONTAINER_ENV_MOUNT, mount_paths.to_string()).into());

    // Set other variables.
    step.variables()
        .iter()
        .filter(|var_instance| !var_instance.is_global_data_reference())
        // Filter out variables without values.
        .filter_map(|var_instance| var_instance.value().as_ref().map(|value| (var_instance.id(), value)))
        .for_each(|(other_var_id, other_var_value)| {
            if other_var_id != CONTAINER_ENV_MOUNT {
                arguments.push("--env".into());
                arguments.push(format!("{}={}", other_var_id, other_var_value).into());
            } else {
                log::warn!("Pipeline {} step {} tried to overwrite the reserved environment variable {} with value {}.", pipeline_id.as_ref(), step.id(), other_var_id, other_var_value);
            }
        });

    // Set container to run.
    arguments.push(format_container_name(&pipeline_id, step.id()).into());

    // Create log directory.
    let logs_path = app_config.experiment_logs_path(experiment_id.to_string());
    std::fs::create_dir_all(&logs_path)?;
    // Open stdout log file.
    let log_path_stdout = app_config.experiment_log_path(
        experiment_id.to_string(),
        &pipeline_id,
        step.id(),
        LogProcessType::Run,
        LogOutputType::StdOut,
    );
    let log_file_stdout = std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .append(false)
        .truncate(true)
        .open(log_path_stdout)?;
    // Open stderr log file.
    let log_path_stderr = app_config.experiment_log_path(
        experiment_id.to_string(),
        &pipeline_id,
        step.id(),
        LogProcessType::Run,
        LogOutputType::StdErr,
    );
    let log_file_stderr = std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .append(false)
        .truncate(true)
        .open(log_path_stderr)?;

    let output = Command::new("docker")
        .stdout(log_file_stdout)
        .stderr(log_file_stderr)
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
    Configuration::hash_string(name)
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
    fn parse_output(&self, output: Output, process_type: LogProcessType) -> Result<(), SeqError> {
        if let Some(step) = &self.executed_step {
            let logs_path = self
                .config
                .experiment_logs_path(step.experiment_id.to_string());
            std::fs::create_dir_all(&logs_path)?;
            let log_path = self.config.experiment_log_path(
                step.experiment_id.to_string(),
                &step.pipeline_id,
                &step.pipeline_step_id,
                process_type,
                LogOutputType::ExitCode,
            );
            let log_file = std::fs::OpenOptions::new()
                .create(true)
                .write(true)
                .append(false)
                .truncate(true)
                .open(log_path)?;
            let mut buffered_writer = BufWriter::new(log_file);
            let exit_code = output
                .status
                .code()
                .map(|code| code.to_string())
                .unwrap_or("Terminated by signal".to_string());
            buffered_writer.write_all(exit_code.as_bytes())?;
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
                    self.parse_output(run.wait_with_output()?, LogProcessType::Run)?;
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
                        self.parse_output(build.wait_with_output()?, LogProcessType::Build)?;
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

    /// Starts the run process.
    fn start_run_process(&mut self) -> Result<(), SeqError> {
        let step = self.get_executed_step()?;
        if let Some(pipeline) = self.loaded_pipelines.get(&step.pipeline_id) {
            let mut connection = self.database_manager.database_connection()?;
            let values = crate::model::db::pipeline_step_variable::PipelineStepVariable::get_values_by_experiment_and_pipeline(step.experiment_id, pipeline.pipeline().id(), &mut connection)?;
            // The stati of the pipeline steps should be None at this point so an empty map is supplied instead of loading them from the database.
            let pipeline = ExperimentPipelineBlueprint::from_internal(
                pipeline.pipeline(),
                values,
                HashMap::new(),
            );
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
