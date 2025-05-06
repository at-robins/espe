//! The `config` module provides configuration setting for application parameters.

/// The context for UUID generation.
const UUID_CONTEXT: Context = Context::new(0);
/// The node ID for UUID generation.
const UUID_NODE_ID: &[u8; 6] = &[12, 221, 33, 14, 35, 16];
/// The context path where temporary data are stored.
const PATH_TEMPORARY: &str = "tmp";
/// The context path where temporary upload files are stored.
const PATH_TEMPORARY_UPLOAD: &str = "upload";
/// The context path where temporary download files are stored.
const PATH_TEMPORARY_DOWNLOAD: &str = "download";
/// The context path where data related to specific experiments or samples is stored.
const PATH_FILES_EXPERIMENTS: &str = "experiments";
/// The file inside each pipeline folder defining the pipeline.
pub const PIPELINE_DEFINITION_FILE: &str = "pipeline.json";
/// The sub-folder where pipeline step output of an experiment is stored.
pub const PATH_FILES_EXPERIMENTS_STEPS: &str = "steps";
/// The sub-folder where experiment input is stored.
pub const PATH_FILES_EXPERIMENTS_INPUT: &str = "input";
/// The sub-folder where all logs are stored.
pub const PATH_FILES_EXPERIMENTS_LOGS: &str = "logs";
/// The folder where global data is stored.
pub const PATH_FILES_GLOBAL_DATA: &str = "globals";
/// The folder where pipeline attachments are stored.
pub const PATH_PIPELINE_ATTACHMENTS: &str = "attachments";

use std::{
    collections::HashSet,
    fmt::Display,
    hash::{Hash, Hasher},
    path::PathBuf,
    time::SystemTime,
};

use getset::{CopyGetters, Getters};
use serde::{Deserialize, Serialize};
use twox_hash::XxHash64;
use uuid::{
    v1::{Context, Timestamp},
    Uuid,
};

use super::{
    environment::{
        CONTEXT_FOLDER, DATABASE_URL, LOG_LEVEL, MODE, PIPELINE_FOLDER, SERVER_ADDRESS, SERVER_PORT,
    },
    error::{SeqError, SeqErrorType},
};

#[derive(Debug, Getters, CopyGetters, PartialEq, Serialize, Deserialize)]
/// A configuration that defines basic parameters of the application.
pub struct Configuration {
    /// The path to the application database file.
    #[getset(get = "pub")]
    database_url: String,
    /// The logging level.
    #[getset(get = "pub")]
    log_level: String,
    /// The address of the server.
    #[getset(get = "pub")]
    server_address: String,
    /// The port of the server.
    #[getset(get = "pub")]
    server_port: String,
    /// The folder where all context relevant data is stored.
    #[getset(get = "pub")]
    context_folder: PathBuf,
    /// The folder where all pipeline definitions are stored.
    #[getset(get = "pub")]
    pipeline_folder: PathBuf,
    /// The application run mode.
    #[getset(get_copy = "pub")]
    mode: ApplicationMode,
}

impl Configuration {
    /// Creates a new configuration with the specified parameters.
    ///
    /// # Parameters
    ///
    /// * `database_url` - the URL  / URI of the database  
    /// * `log_level` - the logging level
    /// * `server_address` - the address of the server
    /// * `server_port` - the port of the server
    /// * `context_folder` - the folder, in which all context related resources are stored
    /// * `pipeline_folder` - the folder, in which all pipeline definitions are stored
    /// * `mode` - the [`ApplicationMode`] the application should be running in
    pub fn new<
        DatabaseUrlType: Into<String>,
        LogLevelType: Into<String>,
        ServerAddressType: Into<String>,
        ServerPortType: Into<String>,
        ContextFolderType: Into<PathBuf>,
        PipelineFolderType: Into<PathBuf>,
    >(
        database_url: DatabaseUrlType,
        log_level: LogLevelType,
        server_address: ServerAddressType,
        server_port: ServerPortType,
        context_folder: ContextFolderType,
        pipeline_folder: PipelineFolderType,
        mode: ApplicationMode,
    ) -> Self {
        Self {
            database_url: database_url.into(),
            log_level: log_level.into(),
            server_address: server_address.into(),
            server_port: server_port.into(),
            context_folder: context_folder.into(),
            pipeline_folder: pipeline_folder.into(),
            mode,
        }
    }

    /// Creates a new configuration if all enviroment variables are setup correctly.
    pub fn create_from_environment() -> Result<Self, SeqError> {
        Ok(Self::new(
            Self::get_environment_variable(DATABASE_URL)?,
            Self::get_environment_variable(LOG_LEVEL)?,
            Self::get_environment_variable(SERVER_ADDRESS)?,
            Self::get_environment_variable(SERVER_PORT)?,
            Self::get_environment_variable(CONTEXT_FOLDER)?,
            Self::get_environment_variable(PIPELINE_FOLDER)?,
            Self::get_environment_variable(MODE)
                .and_then(|env_value| ApplicationMode::try_from(env_value.as_str()))?,
        ))
    }

    /// Retrieves an environment variable by name and returns an error in case of it not being set or being invalid.
    ///
    /// # Parameters
    ///
    /// * `environment_variable` - the environment variable to retrieve
    fn get_environment_variable(environment_variable: &str) -> Result<String, SeqError> {
        std::env::var(environment_variable).map_err(|error| {
            SeqError::from_var_error(error, environment_variable).chain(format!(
                "Error while loading environment variable \"{}\".",
                environment_variable
            ))
        })
    }

    /// The context path where temporary data are stored.
    pub fn temporary_path(&self) -> PathBuf {
        self.context_folder().join(PATH_TEMPORARY)
    }

    /// The context path where temporary upload files are stored.
    pub fn temporary_upload_path(&self) -> PathBuf {
        let mut path: PathBuf = self.temporary_path();
        path.push(PATH_TEMPORARY_UPLOAD);
        path
    }

    /// The context path where temporary download files are stored.
    pub fn temporary_download_path(&self) -> PathBuf {
        let mut path: PathBuf = self.temporary_path();
        path.push(PATH_TEMPORARY_DOWNLOAD);
        path
    }

    /// The context path where a specific temporary download file is stored.
    ///
    /// # Parameters
    ///
    /// * `file_id` - the ID of the temporary file
    pub fn temporary_download_file_path<T: Into<String>>(&self, file_id: T) -> PathBuf {
        let mut path: PathBuf = self.temporary_download_path();
        path.push(file_id.into());
        path
    }

    /// The context path where all global data is stored.
    pub fn globals_path(&self) -> PathBuf {
        let mut path: PathBuf = self.context_folder().clone();
        path.push(PATH_FILES_GLOBAL_DATA);
        path
    }

    /// The context path where a specific pipeline is located.
    ///
    /// # Parameters
    ///
    /// * `pipeline_directory` - the directory the pipeline is located at
    pub fn pipeline_path<T: Into<String>>(&self, pipeline_directory: T) -> PathBuf {
        self.pipeline_folder.join(pipeline_directory.into())
    }

    /// The context path where a specific attachment of the specified pipeline is stored.
    ///
    /// # Parameters
    ///
    /// * `pipeline_directory` - the directory the pipeline is located at
    /// * `attachment_name` - the file name of the attachment
    pub fn pipeline_attachment_path<S: Into<String>, T: Into<String>>(
        &self,
        pipeline_directory: T,
        attachment_name: S,
    ) -> PathBuf {
        let mut path = self.pipeline_path(pipeline_directory);
        path.push(PATH_PIPELINE_ATTACHMENTS);
        path.push(attachment_name.into());
        path
    }

    /// The context path where the specified global data is stored.
    ///
    /// # Parameters
    ///
    /// * `global_id` - the ID of the global data
    pub fn global_data_path<P: AsRef<str>>(&self, global_id: P) -> PathBuf {
        let mut path: PathBuf = self.globals_path();
        path.push(global_id.as_ref());
        path
    }

    /// The context path where data related to specific experiments or samples is stored.
    pub fn experiments_path(&self) -> PathBuf {
        let mut path: PathBuf = self.context_folder().clone();
        path.push(PATH_FILES_EXPERIMENTS);
        path
    }

    /// The context path where data related to the specified experiment is stored.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment
    pub fn experiment_path<P: AsRef<str>>(&self, experiment_id: P) -> PathBuf {
        let mut path: PathBuf = self.experiments_path();
        path.push(experiment_id.as_ref());
        path
    }

    /// The context path where data related to the pipeline input
    /// of the specified experiment is stored.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment
    pub fn experiment_input_path<P: AsRef<str>>(&self, experiment_id: P) -> PathBuf {
        let mut path: PathBuf = self.experiment_path(experiment_id);
        path.push(PATH_FILES_EXPERIMENTS_INPUT);
        path
    }

    /// The context path where data related to the pipeline execution steps
    /// of the specified experiment is stored.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment
    pub fn experiment_steps_path<P: AsRef<str>>(&self, experiment_id: P) -> PathBuf {
        let mut path: PathBuf = self.experiment_path(experiment_id);
        path.push(PATH_FILES_EXPERIMENTS_STEPS);
        path
    }

    /// The context path where data related to the pipeline execution step
    /// of the specified experiment is stored.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment
    /// * `step_id` - the ID of the step
    pub fn experiment_step_path<P: AsRef<str>, Q: AsRef<str>>(
        &self,
        experiment_id: P,
        step_id: Q,
    ) -> PathBuf {
        let mut path: PathBuf = self.experiment_steps_path(experiment_id);
        // Hashing the ID prevents invalid characters in file paths.
        path.push(Self::hash_string(step_id));
        path
    }

    /// The context path where data related to the pipeline logs
    /// of a specified experiment is stored.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment
    pub fn experiment_logs_path<P: AsRef<str>>(&self, experiment_id: P) -> PathBuf {
        let mut path: PathBuf = self.experiment_path(experiment_id);
        path.push(PATH_FILES_EXPERIMENTS_LOGS);
        path
    }

    /// The context path where a specific pipeline log file is stored.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment
    /// * `pipeline_id` - the ID of the pipeline
    /// * `step_id` - the ID of the pipeline step
    /// * `process_type` - the type of process to log
    /// * `output_type` - the output type of the process to log
    pub fn experiment_log_path<P: AsRef<str>, Q: AsRef<str>, R: AsRef<str>>(
        &self,
        experiment_id: P,
        pipeline_id: Q,
        step_id: R,
        process_type: LogProcessType,
        output_type: LogOutputType,
    ) -> PathBuf {
        let mut path: PathBuf = self.experiment_logs_path(experiment_id);
        path.push(format!(
            "{}_{}_{}.log",
            Self::hash_string(format!("{}{}", pipeline_id.as_ref(), step_id.as_ref())),
            process_type,
            output_type,
        ));
        path
    }

    /// The context paths of pipeline log files of the specified types.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment
    /// * `pipeline_id` - the ID of the pipeline
    /// * `step_id` - the ID of the pipeline step
    /// * `log_types` - the [`LogProcessType`]s to get the log files for
    pub fn experiment_log_paths<P: AsRef<str>, Q: AsRef<str>, R: AsRef<str>>(
        &self,
        experiment_id: P,
        pipeline_id: Q,
        step_id: R,
        log_types: &[LogProcessType],
    ) -> HashSet<PathBuf> {
        let log_sub_types = [
            LogOutputType::StdOut,
            LogOutputType::StdErr,
            LogOutputType::ExitCode,
        ];
        let mut paths = HashSet::with_capacity(log_types.len() * log_sub_types.len());
        for process_type in log_types {
            for output_type in log_sub_types {
                paths.insert(self.experiment_log_path(
                    experiment_id.as_ref(),
                    pipeline_id.as_ref(),
                    step_id.as_ref(),
                    *process_type,
                    output_type,
                ));
            }
        }
        paths
    }

    /// Generates a V1 UUID.
    pub fn generate_uuid() -> Uuid {
        let now = SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .expect("The UNIX epoch must be the earliest possible time point.");
        let timestamp = Timestamp::from_unix(UUID_CONTEXT, now.as_secs(), now.subsec_nanos());
        Uuid::new_v1(timestamp, UUID_NODE_ID)
    }

    /// Returns the full server address including port information.
    pub fn server_address_and_port(&self) -> String {
        format!("{}:{}", self.server_address(), self.server_port())
    }

    /// Hashes the specified string and returns the resulting number as a string.
    /// The hash is constant for the same string on repeated uses.
    ///
    /// # Parameters
    ///
    /// * `value` - the string to hash
    pub fn hash_string<T: AsRef<str>>(value: T) -> String {
        let mut hasher = XxHash64::with_seed(154);
        value.as_ref().hash(&mut hasher);
        hasher.finish().to_string()
    }
}

#[derive(Debug, Clone, Copy)]
/// The process types of log files.
pub enum LogProcessType {
    /// The build process.
    Build,
    // The run process.
    Run,
}

impl Display for LogProcessType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                LogProcessType::Build => "build",
                LogProcessType::Run => "run",
            }
        )
    }
}

#[derive(Debug, Clone, Copy)]
/// The output types of log files.
pub enum LogOutputType {
    /// Standard output stream.
    StdOut,
    /// Standard error stream.
    StdErr,
    /// The exit code.
    ExitCode,
}

impl Display for LogOutputType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                LogOutputType::StdOut => "stdout",
                LogOutputType::StdErr => "stderr",
                LogOutputType::ExitCode => "exitcode",
            }
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
/// The application run mode.
pub enum ApplicationMode {
    /// Release mode optimised for performance.
    Release,
    /// Development mode.
    Development,
}

impl Display for ApplicationMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                ApplicationMode::Release => "release",
                ApplicationMode::Development => "development",
            }
        )
    }
}

impl TryFrom<&str> for ApplicationMode {
    type Error = SeqError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        if value == Self::Release.to_string() {
            Ok(Self::Release)
        } else if value == Self::Development.to_string() {
            Ok(Self::Development)
        } else {
            Err(SeqError::new(
                "Application mode error",
                SeqErrorType::InternalServerError,
                format!("\"{}\" is not a valid application run mode.", value),
                "Invalid application run mode.",
            ))
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_server_address_and_port() {
        let config =
            Configuration::new("", "", "127.0.0.1", "8080", "", "", ApplicationMode::Release);
        assert_eq!(&config.server_address_and_port(), "127.0.0.1:8080");
    }

    #[test]
    fn test_temporary_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        let path: PathBuf = "./application/context/tmp".into();
        assert_eq!(config.temporary_path(), path);
    }

    #[test]
    fn test_temporary_upload_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        let path: PathBuf = "./application/context/tmp/upload".into();
        assert_eq!(config.temporary_upload_path(), path);
    }

    #[test]
    fn test_temporary_download_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        let path: PathBuf = "./application/context/tmp/download".into();
        assert_eq!(config.temporary_download_path(), path);
    }

    #[test]
    fn test_temporary_download_file_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        let id = "01234567-89ab-cdef-0123-456789abcdef";
        let path: PathBuf = format!("./application/context/tmp/download/{}", id).into();
        assert_eq!(config.temporary_download_file_path(id), path);
    }

    #[test]
    fn test_globals_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        let path: PathBuf = "./application/context/globals".into();
        assert_eq!(config.globals_path(), path);
    }

    #[test]
    fn test_pipeline_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "./application/pipelines",
            ApplicationMode::Release,
        );
        let path: PathBuf = "./application/pipelines/test".into();
        assert_eq!(config.pipeline_path("test"), path);
    }

    #[test]
    fn test_pipeline_attachment_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "./application/pipelines",
            ApplicationMode::Release,
        );
        let path: PathBuf = "./application/pipelines/test/attachments/test.txt".into();
        assert_eq!(config.pipeline_attachment_path("test", "test.txt"), path);
    }

    #[test]
    fn test_global_data_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        let path: PathBuf = "./application/context/globals/global_id".into();
        assert_eq!(config.global_data_path("global_id"), path);
    }

    #[test]
    fn test_experiments_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        let path: PathBuf = "./application/context/experiments".into();
        assert_eq!(config.experiments_path(), path);
    }

    #[test]
    fn test_experiment_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        let path: PathBuf = "./application/context/experiments/test_id".into();
        assert_eq!(config.experiment_path("test_id"), path);
    }

    #[test]
    fn test_experiment_steps_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        let path: PathBuf = "./application/context/experiments/test_id/steps".into();
        assert_eq!(config.experiment_steps_path("test_id"), path);
    }

    #[test]
    fn test_experiment_step_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        // Hash of step_id.
        let path: PathBuf =
            "./application/context/experiments/experiment_id/steps/4363919453614495606".into();
        assert_eq!(config.experiment_step_path("experiment_id", "step_id"), path);
    }

    #[test]
    fn test_experiment_logs_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        // Hash of step_id.
        let path: PathBuf = "./application/context/experiments/experiment_id/logs".into();
        assert_eq!(config.experiment_logs_path("experiment_id"), path);
    }

    #[test]
    fn test_experiment_log_path() {
        let config = Configuration::new(
            "",
            "",
            "",
            "",
            "./application/context",
            "",
            ApplicationMode::Release,
        );
        // Hash of step_id.
        let path: PathBuf =
            "./application/context/experiments/experiment_id/logs/13269802908832430007_build_stderr.log"
                .into();
        assert_eq!(
            config.experiment_log_path(
                "experiment_id",
                "pipeline_id",
                "step_id",
                LogProcessType::Build,
                LogOutputType::StdErr
            ),
            path
        );
    }

    #[test]
    fn test_hash_string() {
        let random_string = "39012rtuj132-0t1jp41-9/n\n\t@#$%^&*()|}{\"?>¡ªº£€˚„";
        let hash = Configuration::hash_string(random_string);
        let allowed_characters = vec!['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'];
        assert!(hash.len() > 0);
        // u64 max
        assert!(hash.len() <= 20);
        for character in hash.chars() {
            assert!(allowed_characters.contains(&character));
        }
    }
}
