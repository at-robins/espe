use std::{
    io::Write,
    path::{Path, PathBuf},
    sync::Arc,
};

use actix_multipart::Multipart;
use actix_web::web;
use futures_util::TryStreamExt;
use serde::de::DeserializeOwned;
use uuid::Uuid;

use crate::application::{
    config::Configuration,
    error::{SeqError, SeqErrorType},
};

/// A form that can be uploaded to the backend.
pub trait UploadForm: std::fmt::Debug + DeserializeOwned {
    /// Validates the form.
    fn validate(&self) -> Result<(), String>;
}

const MAX_MULTIPART_FORM_SIZE: usize = 524_288;

pub async fn parse_multipart_file<T: UploadForm, P: AsRef<Path>>(
    mut payload: Multipart,
    temporary_file_path: P,
) -> Result<(T, Arc<Path>), SeqError> {
    // Validate that all information is provided.
    let mut is_file_provided = false;
    let mut upload_info: Option<T> = None;

    let temp_file_path: Arc<Path> = temporary_file_path.as_ref().into();

    // Iterate over the multipart stream and save the file.
    while let Some(mut field) = payload.try_next().await? {
        match field.name() {
            Some("form") => {
                let mut body = web::BytesMut::new();
                while let Some(chunk) = field.try_next().await? {
                    let current_length = body.len() + chunk.len();
                    if (current_length) > MAX_MULTIPART_FORM_SIZE {
                        return Err(SeqError::new(
                        "Multipart overflow", 
                        SeqErrorType::BadRequestError,
                        format!("The maximum length for multipart form data is {} bytes, but the current chunk adds up to {} bytes.", MAX_MULTIPART_FORM_SIZE, current_length), 
                        "The multipart body is too large."));
                    } else {
                        body.extend_from_slice(&chunk);
                    }
                }
                upload_info = Some(serde_json::from_slice::<T>(&body)?);
            },
            Some("file") => {
                // Write the file to disk.
                let file_path_ref = Arc::clone(&temp_file_path);
                let mut file = web::block(|| std::fs::File::create(file_path_ref)).await??;
                while let Some(chunk) = field.try_next().await? {
                    file = web::block(move || file.write_all(&chunk).map(|_| file)).await??;
                }
                is_file_provided = true;
            },
            name => log::warn!("Unknown content name option: {:?}", name),
        }
    }

    let upload_info = upload_info.ok_or_else(|| {
        SeqError::new(
            "Missing multipart data",
            SeqErrorType::BadRequestError,
            "The multipart upload did not provide form data.",
            "The upload requires both a file and the according form data.",
        )
    })?;

    if !is_file_provided {
        return Err(SeqError::new(
            "Missing multipart data",
            SeqErrorType::BadRequestError,
            format!("The multipart upload {:?} did not provide a file.", upload_info),
            "The upload requires both a file and the according form data.",
        ));
    }

    // Validate the uploaded form data.
    upload_info.validate().map_err(|e| {
        SeqError::new(
            "Form invalid",
            SeqErrorType::BadRequestError,
            format!("Validation of form data {:?} failed with error: {}", upload_info, e),
            "The provided form data is invalid.",
        )
    })?;

    Ok((upload_info, temp_file_path))
}

/// Deletes the temporary file associated with the specified UUID.
///
/// # Parameters
///
/// * `uuid` - the UUID of the temporary file
/// * `app_config` - the app [`Configuration`]
pub fn delete_temporary_file(uuid: Uuid, app_config: Arc<Configuration>) -> Result<(), SeqError> {
    let mut file_path: PathBuf = app_config.temporary_upload_path();
    file_path.push(uuid.to_string());
    if file_path.exists() {
        std::fs::remove_file(&file_path)?;
        log::info!("Deleted temporary file {}.", file_path.to_string_lossy());
    } else {
        log::warn!("Tried to delete non existing temporary file {:?}.", file_path)
    }
    Ok(())
}

/// Creates a temporary file path and the associated folder structure.
/// Returns the path and [`Uuid`] of the file.
///
/// # Parameters
///
/// * `app_config` - the app [`Configuration`]
pub fn create_temporary_file(app_config: Arc<Configuration>) -> Result<(PathBuf, Uuid), SeqError> {
    // Generate a new UUID for the temporary file.
    let uuid = Configuration::generate_uuid();

    // Create the temporary folder.
    let mut temp_file_path: PathBuf = app_config.temporary_upload_path();
    std::fs::create_dir_all(&temp_file_path)?;
    temp_file_path.push(uuid.to_string());
    Ok((temp_file_path, uuid))
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use serial_test::serial;

    use crate::{application::config::Configuration, test_utility::TestContext};

    use super::*;

    // The tests need to be serial to prevent the same context UUID to be issued
    // to different tests at the same time.

    #[test]
    #[serial]
    fn test_create_temporary_file() {
        let context = TestContext::new();
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let app_config: Arc<Configuration> = Arc::new((&context).into());
        let (path, id) = create_temporary_file(Arc::clone(&app_config)).unwrap();
        assert!(path.ends_with(id.to_string()));
        assert!(path.starts_with(app_config.temporary_upload_path()));
        assert!(!path.exists());
    }

    #[test]
    #[serial]
    fn test_delete_temporary_file() {
        let context = TestContext::new();
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let app_config: Arc<Configuration> = Arc::new((&context).into());
        let (path, id) = create_temporary_file(Arc::clone(&app_config)).unwrap();
        std::fs::File::create(&path).unwrap();
        assert!(path.exists());
        assert!(path.is_file());
        delete_temporary_file(id, Arc::clone(&app_config)).unwrap();
        assert!(!path.exists());
        // Deleting a non existant file should not throw an error.
        delete_temporary_file(id, app_config).unwrap();
        assert!(!path.exists());
    }
}
