//! The `config` module provides configuration setting for application parameters.

/// The context for UUID generation.
const UUID_CONTEXT: Context = Context::new(0);
/// The node ID for UUID generation.
const UUID_NODE_ID: &[u8; 6] = &[12, 221, 33, 14, 35, 16];
/// The context path where temporary files are stored.
const PATH_FILES_TEMPORARY: &str = "tmp/files";
/// The context path where data related to specific experiments or samples is stored.
const PATH_FILES_EXPERIMENTS: &str = "experiments";
/// The file name of the initially submitted sample before processing.
pub const PATH_FILES_EXPERIMENT_INITIAL_FASTQ: &str = "00_initial.fastq.gz";
/// The file inside each pipeline folder defining the pipeline.
pub const PIPELINE_DEFINITION_FILE: &str = "pipeline.json";
/// The sub-folder where pipeline step output is stored.
pub const PATH_FILES_EXPERIMENTS_STEPS: &str = "steps";
/// The sub-folder where initial pipeline input samples are stored.
pub const PATH_FILES_EXPERIMENTS_SAMPLES: &str = "samples";
/// The folder where global data is stored.
pub const PATH_FILES_GLOBAL_DATA: &str = "globals";

use std::{path::PathBuf, time::SystemTime};

use diesel::{Connection, SqliteConnection};
use getset::Getters;
use serde::{Deserialize, Serialize};
use uuid::{
    v1::{Context, Timestamp},
    Uuid,
};

use super::{
    environment::{
        CONTEXT_FOLDER, DATABASE_URL, LOG_LEVEL, PIPELINE_FOLDER, SERVER_ADDRESS, SERVER_PORT,
    },
    error::SeqError,
};

#[derive(Debug, Getters, PartialEq, Serialize, Deserialize)]
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
    ) -> Self {
        Self {
            database_url: database_url.into(),
            log_level: log_level.into(),
            server_address: server_address.into(),
            server_port: server_port.into(),
            context_folder: context_folder.into(),
            pipeline_folder: pipeline_folder.into(),
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
        ))
    }

    /// Retrieves an environment variable by name and returns an error in case of it not being set or being invalid.
    ///
    /// # Parameters
    ///
    /// * `environment_variable` - the environment variable to retrieve
    fn get_environment_variable(environment_variable: &str) -> Result<String, SeqError> {
        std::env::var(environment_variable)
            .map_err(|error| SeqError::from_var_error(error, environment_variable))
    }

    /// Returns a connection to the database if possible.
    pub fn database_connection(&self) -> Result<SqliteConnection, SeqError> {
        let connection = SqliteConnection::establish(self.database_url())?;
        connection.execute("PRAGMA foreign_keys = ON;")?;
        Ok(connection)
    }

    /// The context path where temporary files are stored.
    pub fn temporary_file_path(&self) -> PathBuf {
        let mut path: PathBuf = self.context_folder().clone();
        path.push(PATH_FILES_TEMPORARY);
        path
    }

    /// The context path where all global data is stored.
    pub fn globals_path(&self) -> PathBuf {
        let mut path: PathBuf = self.context_folder().clone();
        path.push(PATH_FILES_GLOBAL_DATA);
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
        path.push(step_id.as_ref());
        path
    }

    /// The context path where data related to the initial pipeline input samples
    /// of a specified experiment is stored.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment
    pub fn experiment_samples_path<P: AsRef<str>>(&self, experiment_id: P) -> PathBuf {
        let mut path: PathBuf = self.experiment_path(experiment_id);
        path.push(PATH_FILES_EXPERIMENTS_SAMPLES);
        path
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
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_server_address_and_port() {
        let config = Configuration::new("", "", "127.0.0.1", "8080", "", "");
        assert_eq!(&config.server_address_and_port(), "127.0.0.1:8080");
    }

    #[test]
    fn test_temporary_file_path() {
        let config = Configuration::new("", "", "", "", "./application/context", "");
        let path: PathBuf = "./application/context/tmp/files".into();
        assert_eq!(config.temporary_file_path(), path);
    }

    #[test]
    fn test_globals_path() {
        let config = Configuration::new("", "", "", "", "./application/context", "");
        let path: PathBuf = "./application/context/globals".into();
        assert_eq!(config.globals_path(), path);
    }

    #[test]
    fn test_global_data_path() {
        let config = Configuration::new("", "", "", "", "./application/context", "");
        let path: PathBuf = "./application/context/globals/global_id".into();
        assert_eq!(config.global_data_path("global_id"), path);
    }

    #[test]
    fn test_experiments_path() {
        let config = Configuration::new("", "", "", "", "./application/context", "");
        let path: PathBuf = "./application/context/experiments".into();
        assert_eq!(config.experiments_path(), path);
    }

    #[test]
    fn test_experiment_path() {
        let config = Configuration::new("", "", "", "", "./application/context", "");
        let path: PathBuf = "./application/context/experiments/test_id".into();
        assert_eq!(config.experiment_path("test_id"), path);
    }

    #[test]
    fn test_experiment_steps_path() {
        let config = Configuration::new("", "", "", "", "./application/context", "");
        let path: PathBuf = "./application/context/experiments/test_id/steps".into();
        assert_eq!(config.experiment_steps_path("test_id"), path);
    }

    #[test]
    fn test_experiment_step_path() {
        let config = Configuration::new("", "", "", "", "./application/context", "");
        let path: PathBuf = "./application/context/experiments/experiment_id/steps/step_id".into();
        assert_eq!(config.experiment_step_path("experiment_id", "step_id"), path);
    }

    #[test]
    fn test_experiment_samples_path() {
        let config = Configuration::new("", "", "", "", "./application/context", "");
        let path: PathBuf = "./application/context/experiments/test_id/samples".into();
        assert_eq!(config.experiment_samples_path("test_id"), path);
    }
}
