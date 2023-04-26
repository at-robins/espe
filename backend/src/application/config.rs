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

use std::time::SystemTime;

use diesel::{Connection, SqliteConnection};
use getset::Getters;
use serde::{Deserialize, Serialize};
use uuid::{
    v1::{Context, Timestamp},
    Uuid,
};

use super::{
    environment::{CONTEXT_FOLDER, DATABASE_URL, LOG_LEVEL, SERVER_ADDRESS, SERVER_PORT, PIPELINE_FOLDER},
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
    context_folder: String,
    /// The folder where all pipeline definitions are stored.
    #[getset(get = "pub")]
    pipeline_folder: String,
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
        ContextFolderType: Into<String>,
        PipelineFolderType: Into<String>,
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
        Ok(Self {
            database_url: Self::get_environment_variable(DATABASE_URL)?,
            log_level: Self::get_environment_variable(LOG_LEVEL)?,
            server_address: Self::get_environment_variable(SERVER_ADDRESS)?,
            server_port: Self::get_environment_variable(SERVER_PORT)?,
            context_folder: Self::get_environment_variable(CONTEXT_FOLDER)?,
            pipeline_folder: Self::get_environment_variable(PIPELINE_FOLDER)?,
        })
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
    pub fn temporary_file_path(&self) -> String {
        format!("{}/{}", self.context_folder(), PATH_FILES_TEMPORARY)
    }

    /// The context path where data related to specific experiments or samples is stored.
    pub fn experiment_path(&self) -> String {
        format!("{}/{}", self.context_folder(), PATH_FILES_EXPERIMENTS)
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
