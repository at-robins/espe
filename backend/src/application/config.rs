//! The `config` module provides configuration setting for application parameters.

/// The folder in which all information relevant for application configuration is stored.
const DEFAULT_FOLDER_APPLICATION: &str = "application";
/// The name of the default configuration file.
const DEFAULT_FILE_APPLICATION_CONFIGURATION: &str = "configuration";
/// The file extension of a JSON.
const EXTENSION_JSON: &str = "json";
/// The default server address.
const DEFAULT_CONFIG_SERVER_ADDRESS: &str = "127.0.0.1";
/// The default server port.
const DEFAULT_CONFIG_SERVER_PORT: &str = "8080";
/// The context for UUID generation.
const UUID_CONTEXT: Context = Context::new(0);
/// The node ID for UUID generation.
const UUID_NODE_ID: &[u8; 6] = &[12, 21, 33, 4, 35, 116];

use std::{
    fs::{File, OpenOptions},
    path::PathBuf,
    time::SystemTime,
};

use serde::{Deserialize, Serialize};
use uuid::{
    v1::{Context, Timestamp},
    Uuid,
};

use super::error::SeqError;

#[derive(Debug, PartialEq, Serialize, Deserialize)]
/// A configuration that defines basic parameters of the application.
pub struct Configuration {
    server_address: Option<String>,
    server_port: Option<String>,
}

impl Configuration {
    /// Tries to load a previously stored configuration from the default location.
    pub fn load_from_file() -> Result<Self, SeqError> {
        let file = File::open(Configuration::application_configuration_file_path())?;
        Ok(serde_json::from_reader(file)?)
    }

    /// Tries to save the current configuration to the default file location.
    pub fn save_to_file(&self) -> Result<(), SeqError> {
        let path = Configuration::application_configuration_file_path();
        if let Some(parent_path) = path.parent() {
            std::fs::create_dir_all(parent_path)?;
        }
        let file = OpenOptions::new()
            .create(true)
            .append(false)
            .write(true)
            .open(path)?;
        Ok(serde_json::to_writer_pretty(file, self)?)
    }

    /// Checks if the configuration exists as a physical file.
    pub fn exists() -> bool {
        Configuration::application_configuration_file_path().exists()
    }

    /// Generates a V1 UUID.
    pub fn generate_uuid() -> Uuid {
        let now = SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .expect("The UNIX epoch must be the earliest possible time point.");
        let timestamp = Timestamp::from_unix(UUID_CONTEXT, now.as_secs(), now.subsec_nanos());
        Uuid::new_v1(timestamp, UUID_NODE_ID).expect("The node ID must be of length 6.")
    }

    /// The path to the application configuration folder.
    pub fn application_configuration_folder_path() -> PathBuf {
        PathBuf::from(DEFAULT_FOLDER_APPLICATION)
    }

    /// The path to the application configuration file.
    pub fn application_configuration_file_path() -> PathBuf {
        let mut path = Configuration::application_configuration_folder_path();
        path.push(DEFAULT_FILE_APPLICATION_CONFIGURATION);
        path.set_extension(EXTENSION_JSON);
        path
    }

    /// Returns the address of the server.
    pub fn server_address(&self) -> String {
        self.server_address
            .clone()
            .unwrap_or(DEFAULT_CONFIG_SERVER_ADDRESS.to_string())
    }

    /// Returns the port of the server.
    pub fn server_port(&self) -> String {
        self.server_port
            .clone()
            .unwrap_or(DEFAULT_CONFIG_SERVER_PORT.to_string())
    }

    /// Returns the full server address including port information.
    pub fn server_address_and_port(&self) -> String {
        format!("{}:{}", self.server_address(), self.server_port())
    }
}

impl Default for Configuration {
    fn default() -> Self {
        Self {
            server_address: None,
            server_port: None,
        }
    }
}

#[cfg(test)]
mod test;
