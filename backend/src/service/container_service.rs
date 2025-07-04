use std::{
    collections::HashMap,
    ffi::OsString,
    io::{BufWriter, Write},
    path::Path,
    process::{Child, Command, Output, Stdio},
    sync::Arc,
    thread,
};

use actix_web::web;
use chrono::NaiveDateTime;
use diesel::{ExpressionMethods, QueryDsl, RunQueryDsl, SqliteConnection};

use crate::{
    application::{
        config::{ApplicationMode, Configuration, LogOutputType, LogProcessType},
        database::DatabaseManager,
        error::{SeqError, SeqErrorType},
    },
    model::{
        db::{
            experiment_execution::{ExecutionStatus, ExperimentExecution},
            pipeline_build_register::PipelineBuildRegister,
        },
        exchange::experiment_pipeline::{
            ExperimentPipelineBlueprint, ExperimentPipelineStepBlueprint,
        },
        internal::pipeline_blueprint::{ContextualisedPipelineBlueprint, PipelineStepBlueprint},
    },
};

use super::{
    experiment_service::{open_step_log, prepare_context_for_build, prepare_context_for_run},
    pipeline_service::LoadedPipelines,
};

/// The container environment variable specifying all mounts.
const CONTAINER_ENV_MOUNT: &str = "MOUNT_PATHS";
/// The timeout for checking if a container exists.
const CONTAINER_QUERY_TIMEOUT: std::time::Duration = std::time::Duration::from_secs(2);

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
    prepare_context_for_build(
        &pipeline_id,
        step.id(),
        experiment_id,
        web::Data::clone(&app_config),
    )?;
    let mut pipeline_step_path = context.as_ref().to_path_buf();
    pipeline_step_path.push("container");
    pipeline_step_path.push(step.container());
    let build_arg: OsString = "build".into();
    let name_spec: OsString = "-t".into();
    let name_arg: OsString = format_container_name(&pipeline_id, step.id()).into();
    let progress_arg: OsString = "--progress=plain".into();

    let child = Command::new("docker")
        .stdout(open_step_log(
            &pipeline_id,
            step.id(),
            experiment_id,
            LogProcessType::Build,
            LogOutputType::StdOut,
            web::Data::clone(&app_config),
        )?)
        .stderr(open_step_log(
            &pipeline_id,
            step.id(),
            experiment_id,
            LogProcessType::Build,
            LogOutputType::StdErr,
            app_config,
        )?)
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
pub fn run_pipeline_step(
    pipeline: &ExperimentPipelineBlueprint,
    step: &ExperimentPipelineStepBlueprint,
    experiment_id: i32,
    app_config: web::Data<Configuration>,
) -> Result<Child, SeqError> {
    prepare_context_for_run(
        pipeline.id(),
        step.id(),
        experiment_id,
        web::Data::clone(&app_config),
    )?;
    let experiment_id_string = experiment_id.to_string();
    let output_path = app_config.experiment_step_path(&experiment_id_string, step.id());
    let mut arguments: Vec<OsString> = vec![
        "run".into(),
        "--name".into(),
        format_container_name(pipeline.id(), step.id()).into(),
        "--rm".into(),
    ];

    arguments.extend(pipeline_step_mount(output_path, "/output", false));
    // Set initial sample input mount.
    arguments.extend(pipeline_step_mount(
        app_config.experiment_input_path(&experiment_id_string),
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
            app_config.experiment_step_path(&experiment_id_string, dependency_id),
            target,
            true,
        ));
    }
    // Set global mounts.
    let mut mount_map_globals = serde_json::Map::new();
    for (global_var_id, global_var_value) in pipeline.global_variables()
        .iter()
        .chain(step.variables().iter())
        .filter(|var_instance| var_instance.is_global_data_reference())
        // Filter out variables without values.
        .filter_map(|var_instance| var_instance.value().as_ref().map(|value| (var_instance.id(), value)))
    {
        let target = format!("/input/globals/{}", Configuration::hash_string(global_var_id));
        mount_map_globals
            .insert(global_var_id.to_string(), serde_json::Value::String(target.clone()));
        let global_data_path = app_config.global_data_path(global_var_value);
        arguments.extend(pipeline_step_mount(&global_data_path, target, true));
        // Create global data directory in case the reposiory is empty.
        std::fs::create_dir_all(&global_data_path)?;
    }

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
    pipeline.global_variables()
        .iter()
        .chain(step.variables().iter())
        .filter(|var_instance| !var_instance.is_global_data_reference())
        // Filter out variables without values.
        .filter_map(|var_instance| var_instance.value().as_ref().map(|value| (var_instance.id(), value)))
        .for_each(|(other_var_id, other_var_value)| {
            if other_var_id != CONTAINER_ENV_MOUNT {
                arguments.push("--env".into());
                arguments.push(format!("{}={}", other_var_id, other_var_value).into());
            } else {
                log::warn!("Pipeline {} step {} tried to overwrite the reserved environment variable {} with value {}.", pipeline.id(), step.id(), other_var_id, other_var_value);
            }
        });

    // Set container to run.
    arguments.push(format_container_name(pipeline.id(), step.id()).into());

    let output = Command::new("docker")
        .stdout(open_step_log(
            pipeline.id(),
            step.id(),
            experiment_id,
            LogProcessType::Run,
            LogOutputType::StdOut,
            web::Data::clone(&app_config),
        )?)
        .stderr(open_step_log(
            pipeline.id(),
            step.id(),
            experiment_id,
            LogProcessType::Run,
            LogOutputType::StdErr,
            app_config,
        )?)
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

/// Returns `true` if the container should be built / re-built.
///
/// # Parameters
///
/// * `pipeline_id` - the ID of the pipeline
/// * `pipeline_step_id` - the ID of the pipeline step
/// * `pipeline_version` - the version of the pipeline
/// * `connection` - the database connection
fn should_build<
    PipelineIdIdType: Into<String>,
    PipelineStepIdIdType: Into<String>,
    PipelineVersionType: Into<String>,
>(
    pipeline_id: PipelineIdIdType,
    pipeline_step_id: PipelineStepIdIdType,
    pipeline_version: PipelineVersionType,
    connection: &mut SqliteConnection,
    app_config: web::Data<Configuration>,
) -> Result<bool, SeqError> {
    // Always re-builds the container if the application runs
    // in development mode.
    if app_config.mode() == ApplicationMode::Development {
        log::trace!("Running in development mode and skipping the build cache check.");
        return Ok(true);
    }

    let pipeline_id = pipeline_id.into();
    let pipeline_step_id = pipeline_step_id.into();
    let pipeline_version = pipeline_version.into();
    // If the container is registered in the database,
    // check if it exists.
    if PipelineBuildRegister::is_built(
        &pipeline_id,
        &pipeline_step_id,
        &pipeline_version,
        connection,
    )? {
        log::info!("Container image {} for pipeline {} step {} version {} is present in the database. Checking physical image.",
                    format_container_name(&pipeline_id, &pipeline_step_id),
                    &pipeline_id,
                    &pipeline_step_id,
                    &pipeline_version
                );
        let image_arg: OsString = "image".into();
        let inspect_arg: OsString = "inspect".into();
        let name_arg: OsString = format_container_name(&pipeline_id, &pipeline_step_id).into();
        let mut child = Command::new("docker")
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .args([
                image_arg.as_os_str(),
                inspect_arg.as_os_str(),
                name_arg.as_os_str(),
            ])
            .spawn()?;
        let start = std::time::Instant::now();
        // Blocking here and waiting for the command to finish is fine, as this should
        // typically be very fast and the subsequent execution time is way longer than
        // the maximum timeout.
        while std::time::Instant::now() - start <= CONTAINER_QUERY_TIMEOUT
            && child.try_wait()?.is_none()
        {
            thread::sleep(std::time::Duration::from_millis(100))
        }
        match child.try_wait()? {
            Some(status) => {
                if status.success() {
                    // Re-building the container can only be skipped if
                    // its entered into the database and if the
                    // image physically exists.
                    return Ok(false);
                } else {
                    log::warn!("Database inconsistency. Container image {} for pipeline {} step {} version {} is present in the database, but not as physical image.",
                        format_container_name(&pipeline_id, &pipeline_step_id),
                        &pipeline_id,
                        &pipeline_step_id,
                        &pipeline_version
                    );
                }
            },
            // Kills the process if still running.
            None => {
                log::warn!("Checking existance of container image {} for pipeline {} step {} version {} timed out. The process is killed.",
                    format_container_name(&pipeline_id, &pipeline_step_id),
                    &pipeline_id,
                    &pipeline_step_id,
                    &pipeline_version
                );
                child.kill()?;
            },
        }
    }
    Ok(true)
}

pub struct ContainerHandler {
    config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    loaded_pipelines: web::Data<LoadedPipelines>,
    executed_step: Option<ExperimentExecution>,
    build_process: Option<Child>,
    build_version: Option<String>,
    run_process: Option<Child>,
    stop_processes: Vec<Child>,
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
            build_version: None,
            run_process: None,
            stop_processes: Vec::new(),
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
            self.stop_container()?;
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

    /// Tries to stop the container associated with the current execution step.
    fn stop_container(&mut self) -> Result<(), SeqError> {
        log::info!("Trying to stop container {:?}.", &self.executed_step);
        if self.run_process.is_some() {
            if let Some(step) = &self.executed_step {
                let stop_arg: OsString = "stop".into();
                let name_arg: OsString =
                    format_container_name(&(step.pipeline_id), &(step.pipeline_step_id)).into();
                self.stop_processes.push(
                    Command::new("docker")
                        .stdout(Stdio::piped())
                        .stderr(Stdio::piped())
                        .args([stop_arg.as_os_str(), name_arg.as_os_str()])
                        .spawn()?,
                );
            }
        } else {
            log::warn!(
                "Cannot stop container {:?} as there is no according run process.",
                &self.executed_step
            );
        }
        Ok(())
    }

    /// Resets the internal state of the handler.
    fn reset(&mut self) {
        self.build_process = None;
        self.build_version = None;
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

    /// Returns the current [`ContextualisedPipelineBlueprint`] or an error
    /// if no [`ExperimentExecution`] is set or the pipeline cannot be found.
    fn get_pipeline_blueprint(&self) -> Result<Arc<ContextualisedPipelineBlueprint>, SeqError> {
        let step = self.get_executed_step()?;
        self.loaded_pipelines.get(&step.pipeline_id).ok_or_else(|| {
            SeqError::new(
                "Invalid container state",
                SeqErrorType::InternalServerError,
                format!(
                    "The pipeline including exectued pipeline step {:?} is not loaded.",
                    &self.executed_step
                ),
                "The pipeline is not defined.",
            )
        })
    }

    /// Loads the pipeline version for the step to be built / executed
    /// and fixes it.
    fn load_pipeline_version(&mut self) -> Result<&str, SeqError> {
        // Loads the pipeline version and fixes it for the current step exection / build.
        if self.build_version.is_none() {
            let pipeline_version = self
                .get_pipeline_blueprint()
                .map_err(|err| err.chain("The pipeline for finding out the pipeline step version could not be loaded."))?
                .pipeline()
                .version()
                .to_string();
            self.build_version = Some(pipeline_version);
        }
        Ok(self
            .build_version
            .as_ref()
            .expect("The pipeline step build version must be set at this point.")
            .as_str())
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

    /// Updates and logs process that should stop containers.
    fn update_stop_processes(&mut self) {
        let mut i = 0;
        while i < self.stop_processes.len() {
            if self.stop_processes[i].try_wait().is_ok() {
                match self.stop_processes.remove(i).wait_with_output() {
                    Ok(output) => {
                        if output.status.success() {
                            log::info!("Container successfully stopped.");
                        } else {
                            log::error!(
                                "Stopping a container failed with exit code {}",
                                output.status
                            );
                        }
                        log::info!(
                            "Stopped container StdOut: {}",
                            String::from_utf8_lossy(&output.stdout)
                        );
                        log::info!(
                            "Stopped container StdErr: {}",
                            String::from_utf8_lossy(&output.stderr)
                        );
                    },
                    Err(err) => {
                        // Just log the error as we cannot
                        // do anything about it.
                        log::error!(
                            "Getting the output of a container stop process failed with error: {}",
                            err
                        );
                    },
                }
            } else {
                i += 1;
            }
            if !self.stop_processes.is_empty() {
                log::info!("{} containers are currently being stopped.", self.stop_processes.len());
            }
        }
    }

    /// Updates the handler and returns ```true``` if ready for processing another step.
    fn update_inner(&mut self) -> Result<bool, SeqError> {
        // First handle active stop processes.
        self.update_stop_processes();
        // Ready to handle new containers if none is actively running.
        if !self.is_running() {
            return Ok(true);
        }
        let step_db_id: i32 = self.get_executed_step()?.id;
        let pipeline_id: String = self.get_executed_step()?.pipeline_id.to_string();
        let pipeline_step_id: String = self.get_executed_step()?.pipeline_step_id.to_string();
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
            ProcessStatus::NotStarted => {
                let mut connection = self.database_manager.database_connection()?;
                let pipeline_version = self
                    .load_pipeline_version()
                    .map(|version| version.to_string())
                    .map_err(|err| err.chain(""))?;
                match Self::get_process_status(&mut self.build_process)? {
                    ProcessStatus::Finished => {
                        if let Some(build) = self.build_process.take() {
                            // Handle output.
                            self.parse_output(build.wait_with_output()?, LogProcessType::Build)?;
                            // Start the subsequent run process.
                            self.start_run_process()?;
                            // Saves the successful build to the database.
                            PipelineBuildRegister::set_built(
                                pipeline_id,
                                pipeline_step_id,
                                pipeline_version,
                                &mut connection,
                            )?;
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
                        // Starts the build process if a (re-)build should be carried.
                        // Otherwise directly starts the run process.
                        if should_build(
                            &pipeline_id,
                            &pipeline_step_id,
                            &pipeline_version,
                            &mut connection,
                            web::Data::clone(&self.config),
                        )? {
                            self.start_build_process()?;
                        } else {
                            log::info!("The container for pipeline {} step {} version {} has already been built. Skipping build step.",&pipeline_id, &pipeline_step_id, pipeline_version);
                            self.start_run_process()?;
                        }
                        Ok(false)
                    },
                    ProcessStatus::Running => Ok(false),
                }
            },
            ProcessStatus::Running => Ok(false),
        }
    }

    /// Starts the build process.
    fn start_build_process(&mut self) -> Result<(), SeqError> {
        let step = self.get_executed_step()?;
        let pipeline = self.get_pipeline_blueprint()?;
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
    }

    /// Starts the run process.
    fn start_run_process(&mut self) -> Result<(), SeqError> {
        let step = self.get_executed_step()?;
        let pipeline = self.get_pipeline_blueprint()?;
        let mut connection = self.database_manager.database_connection()?;
        let values_global = crate::model::db::pipeline_global_variable::PipelineGlobalVariable::get_values_by_experiment_and_pipeline(step.experiment_id, pipeline.pipeline().id(), &mut connection)?;
        let values_step = crate::model::db::pipeline_step_variable::PipelineStepVariable::get_values_by_experiment_and_pipeline(step.experiment_id, pipeline.pipeline().id(), &mut connection)?;
        // The stati of the pipeline steps should be None at this point so an empty map is supplied instead of loading them from the database.
        let pipeline = ExperimentPipelineBlueprint::from_internal(
            pipeline.pipeline(),
            values_global,
            values_step,
            HashMap::new(),
        );
        if let Some(step_blueprint) = pipeline
            .steps()
            .iter()
            .find(|s| s.id() == &step.pipeline_step_id)
        {
            log::info!("Running {:?}", &step);
            self.run_process = Some(run_pipeline_step(
                &pipeline,
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
