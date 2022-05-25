//! The `config` module provides configuration setting for application parameters.

/// The context for UUID generation.
const UUID_CONTEXT: Context = Context::new(0);
/// The node ID for UUID generation.
const UUID_NODE_ID: &[u8; 6] = &[12, 221, 33, 14, 35, 16];

use std::{
    time::SystemTime,
};

use getset::Getters;
use serde::{Deserialize, Serialize};
use uuid::{
    v1::{Context, Timestamp},
    Uuid,
};

use super::{environment::{DATABASE_URL, SERVER_PORT, SERVER_ADDRESS}, error::SeqError};

#[derive(Debug, Getters, PartialEq, Serialize, Deserialize)]
/// A configuration that defines basic parameters of the application.
pub struct Configuration {
    /// The path to the application database file.
    #[getset(get = "pub")]
    database_url: String,
    /// The address of the server.
    #[getset(get = "pub")]
    server_address: String,
    /// The port of the server.
    #[getset(get = "pub")]
    server_port: String,
}

impl Configuration {
    /// Creates a new configuration if all enviroment variables are setup correctly. 
    pub fn new() -> Result<Self, SeqError> {
        Ok(Self {
            database_url: std::env::var(DATABASE_URL)?,
            server_address: std::env::var(SERVER_ADDRESS)?,
            server_port: std::env::var(SERVER_PORT)?,
        })
    }

    /// Generates a V1 UUID.
    pub fn generate_uuid() -> Uuid {
        let now = SystemTime::now()
            .duration_since(SystemTime::UNIX_EPOCH)
            .expect("The UNIX epoch must be the earliest possible time point.");
        let timestamp = Timestamp::from_unix(UUID_CONTEXT, now.as_secs(), now.subsec_nanos());
        Uuid::new_v1(timestamp, UUID_NODE_ID).expect("The node ID must be of length 6.")
    }

    /// Returns the full server address including port information.
    pub fn server_address_and_port(&self) -> String {
        format!("{}:{}", self.server_address(), self.server_port())
    }
}
