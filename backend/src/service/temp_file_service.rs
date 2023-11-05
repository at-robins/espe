use std::{fs::DirEntry, time::SystemTime};

use actix_web::web;

use crate::application::{
    config::Configuration,
    error::{SeqError, SeqErrorType},
};

/// The maximum age in seconds a temporary download file is kept before being deleted.
const MAX_AGE_TEMPORARY_DOWNLOAD: u64 = 60 * 60 * 24;

/// A manager for temporary files, e.g. automatic deletion.
pub struct TemporaryFileManager {
    config: web::Data<Configuration>,
}

impl TemporaryFileManager {
    /// Creates a new `TemporaryFileManager`.
    ///
    /// # Parameters
    ///
    /// * `config` - the application's [`Configuration`]
    pub fn new(config: web::Data<Configuration>) -> Self {
        Self { config }
    }

    /// Updates the file manager.
    pub fn update(&self) -> Result<(), SeqError> {
        log::info!("Managing temporary files...");
        let temp_download_path = self.config.temporary_download_path();

        // Checks all temporary download files / folders and collects errors.
        let errors: Vec<std::io::Error> = std::fs::read_dir(temp_download_path)?
            // Log errors and return processable files / folders.
            // Continues despite errors to at least clean up all the data that can be cleaned up. 
            .filter_map(Self::filter_log)
            // Filters old temporary data.
            .filter_map(|entry| {
                Self::filter_log(Self::dir_entry_created(&entry))
                    .and_then(|created| Self::filter_log(created.elapsed()))
                    .filter(|lifetime| lifetime.as_secs() > MAX_AGE_TEMPORARY_DOWNLOAD)
                    .map(|_| entry.path())
            })
            // Tries to delete old data and collects errors.
            .filter_map(|entry| {
                match if entry.is_dir() {
                    log::info!("Deleting temporary download directory {}.", entry.display());
                    std::fs::remove_dir_all(entry)
                } else {
                    log::info!("Deleting temporary download file {}.", entry.display());
                    std::fs::remove_file(entry)
                } {
                    Ok(_) => None,
                    Err(err) => Some(err),
                }
            })
            .collect();

        if !errors.is_empty() {
            // Returns errors if present.
            let combined_error = errors.into_iter().fold(String::new(), |mut acc, error| {
                acc.push_str(&error.to_string());
                acc.push('\n');
                acc
            });
            Err(SeqError::new(
                    "std::io::Error",
                    SeqErrorType::InternalServerError,
                    combined_error,
                    "An unforseen error occured during temporary file management. Please consult the logs.",
                ))
        } else {
            // Returns if successful.
            log::info!("Done managing temporary files.");
            Ok(())
        }
    }

    /// Gets the metadata of the directory entry.
    ///
    /// # Parameters
    ///
    /// * `entry` - the diretory entry
    fn dir_entry_created(entry: &DirEntry) -> Result<SystemTime, std::io::Error> {
        entry.metadata()?.created()
    }

    /// Filters and logs errors.
    ///
    /// # Parameters
    ///
    /// * `value` - the value to filter and log
    fn filter_log<T, E: Into<SeqError>>(value: Result<T, E>) -> Option<T> {
        match value {
            Ok(inner) => Some(inner),
            Err(err) => {
                // Automatically logs the error.
                let _err: SeqError = err.into();
                None
            },
        }
    }
}
