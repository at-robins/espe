use std::{
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::{
    application::{
        config::Configuration,
        database::DatabaseManager,
        error::{SeqError, SeqErrorType},
    },
    model::{
        db::{
            experiment::Experiment, experiment_execution::ExperimentExecution,
            global_data::GlobalData,
        },
        exchange::file_path::{FileDetails, FilePath},
    },
    service::{
        download_service::{ArchiveStream, DownloadTrackerManager},
        experiment_service::is_experiment_locked_err,
        global_data_service::is_global_data_locked_err,
        multipart_service::{
            create_temporary_file, delete_temporary_file, parse_multipart_file, UploadForm,
        },
        pipeline_service::LoadedPipelines,
    },
};
use actix_files::NamedFile;
use actix_multipart::Multipart;
use actix_web::{http::header, web, HttpResponse};
use diesel::SqliteConnection;
use futures::StreamExt;
use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
/// The category or type of a file system related request.
pub enum FileRequestCategory {
    /// A request related to global data repositories.
    Globals,
    /// A request related to experiments.
    Experiments,
}

impl FileRequestCategory {
    /// Returns the base file path depending on the request category.
    ///
    /// # Parameters
    /// * `app_config` - the application's [`Configuration`]
    /// * `id` - the ID corresponding to the request category entity
    pub fn base_path(&self, app_config: web::Data<Configuration>, id: i32) -> PathBuf {
        match self {
            FileRequestCategory::Globals => app_config.global_data_path(id.to_string()),
            FileRequestCategory::Experiments => app_config.experiment_input_path(id.to_string()),
        }
    }

    /// Returns if an entity of the request category exists or an error otherwise.
    ///
    /// # Parameters
    /// * `id` - the ID corresponding to the request category entity
    /// * `connection` - a database connection
    pub fn entity_exists(
        &self,
        id: i32,
        connection: &mut SqliteConnection,
    ) -> Result<(), SeqError> {
        match self {
            FileRequestCategory::Globals => GlobalData::exists_err(id, connection),
            FileRequestCategory::Experiments => Experiment::exists_err(id, connection),
        }
    }

    /// Returns an error if the corresponding entity is locked.
    ///
    /// # Parameters
    /// * `id` - the ID corresponding to the request category entity
    /// * `pipelines` - all currently loaded pipelines
    /// * `connection` - a database connection
    pub fn is_locked_err(
        &self,
        id: i32,
        pipelines: web::Data<LoadedPipelines>,
        download_tracker: web::Data<DownloadTrackerManager>,
        connection: &mut SqliteConnection,
    ) -> Result<(), SeqError> {
        match self {
            FileRequestCategory::Globals => {
                is_global_data_locked_err(id, pipelines, download_tracker, connection).map_err(
                    |err| {
                        err.chain(format!(
                            "The file request could not be processed as \
                            global data repository {} is locked.",
                            id
                        ))
                    },
                )
            },
            FileRequestCategory::Experiments => {
                is_experiment_locked_err(id, download_tracker, connection).map_err(|err| {
                    err.chain(format!(
                        "The file request could not be \
                        processed as experiment {} is locked.",
                        id
                    ))
                })
            },
        }
    }
}

pub async fn get_files(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    params: web::Path<(FileRequestCategory, i32)>,
) -> Result<HttpResponse, SeqError> {
    let (category, id) = params.into_inner();
    let mut connection = database_manager.database_connection().map_err(|err| {
        err.chain(format!(
            "Getting files for category {:?} ID {} failed \
            as no database connection could be established.",
            category, id
        ))
    })?;

    category.entity_exists(id, &mut connection).map_err(|err| {
        err.chain(format!(
            "No files for category {:?} ID {} can be retrieved as it does not exist.",
            category, id
        ))
    })?;

    let mut all_files: Vec<FileDetails> = Vec::new();
    let data_path = category.base_path(app_config, id);
    if data_path.exists() {
        for entry in walkdir::WalkDir::new(&data_path) {
            let entry = entry.map_err(|err| {
                SeqError::from(err).chain(format!(
                    "Directory entry in {} could not be retrieved \
                    while getting files for category {:?} ID {}.",
                    data_path.display(),
                    category,
                    id
                ))
            })?;
            let relative_path = entry
                .path()
                .strip_prefix(&data_path)
                .expect("The base directory must be a parent of the contained components.");
            let (components, invalid_path) =
                relative_path
                    .components()
                    .fold((Vec::new(), false), |mut acc, comp| {
                        if let Some(valid_component) = comp.as_os_str().to_str() {
                            acc.0.push(valid_component.to_string());
                            acc
                        } else {
                            (acc.0, true)
                        }
                    });
            if invalid_path {
                return Err(SeqError::new(
                    "Invalid path",
                    SeqErrorType::InternalServerError,
                    format!("The path {} contains invalid characters.", entry.path().display()),
                    "A path is invalid.",
                ));
            }
            // Excludes the root folder, which has an empty component vector.
            if !components.is_empty() {
                all_files.push(FileDetails {
                    path_components: components,
                    is_file: entry.path().is_file(),
                });
            }
        }
        all_files.sort();
    }
    Ok(HttpResponse::Ok().json(all_files))
}

pub async fn delete_files_by_path(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    pipelines: web::Data<LoadedPipelines>,
    download_tracker: web::Data<DownloadTrackerManager>,
    params: web::Path<(FileRequestCategory, i32)>,
    path: web::Json<FilePath>,
) -> Result<HttpResponse, SeqError> {
    let (category, id) = params.into_inner();
    let delete_info = path.into_inner();
    let delete_path = delete_info.file_path();
    let mut connection = database_manager.database_connection().map_err(|err| {
        err.chain(format!(
            "Deleting files for category {:?} ID {} failed \
            as no database connection could be established.",
            category, id
        ))
    })?;

    category.entity_exists(id, &mut connection).map_err(|err| {
        err.chain(format!(
            "No files for category {:?} ID {} can be deleted as it does not exist.",
            category, id
        ))
    })?;
    category
        .is_locked_err(id, pipelines, download_tracker, &mut connection)
        .map_err(|err| {
            err.chain(format!(
                "The file deletion for path {} request failed as {:?}/{} is locked.",
                delete_path.display(),
                category,
                id
            ))
        })?;
    if let Err(validation_error) = delete_info.validate() {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            format!(
                "The file path {:?} is invalid and deletion faild with error: {}.",
                delete_info, validation_error
            ),
            "The requested path is invalid.",
        ));
    }

    let data_path = category.base_path(app_config, id);
    let full_path: PathBuf = data_path.join(&delete_path);
    if !full_path.exists() {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::NotFoundError,
            format!("The path {} does not exist.", full_path.display()),
            "The resource does not exist.",
        ));
    } else if full_path.is_dir() {
        std::fs::remove_dir_all(&full_path).map_err(|err| {
            SeqError::from(err).chain(format!(
                "Deleting directory {} for category {:?} ID {} failed.",
                full_path.display(),
                category,
                id
            ))
        })?;
    } else {
        std::fs::remove_file(&full_path).map_err(|err| {
            SeqError::from(err).chain(format!(
                "Deleting file {} for category {:?} ID {} failed.",
                full_path.display(),
                category,
                id
            ))
        })?;
    }
    Ok(HttpResponse::Ok().finish())
}

/// This function is needed as an unhandled payload, which might be caused
/// by a preceeding error, causes the connection / request to hang.
/// Dropping the payload strangely causes the next request to hang,
/// so the payload is handled in case of an error as long as this bug exists.
async fn mock_handle_payload_on_error(payload: &mut Multipart) {
    while let Some(_) = payload.next().await {}
}

/// Necessary to mitigate the bug described in [`mock_handle_payload_on_error`].
async fn setup_post_add_file(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    pipelines: web::Data<LoadedPipelines>,
    download_tracker: web::Data<DownloadTrackerManager>,
    params: web::Path<(FileRequestCategory, i32)>,
) -> Result<
    (
        FileRequestCategory,
        i32,
        PathBuf,
        uuid::Uuid,
        diesel::r2d2::PooledConnection<diesel::r2d2::ConnectionManager<SqliteConnection>>,
    ),
    SeqError,
> {
    let (category, id) = params.into_inner();

    let (temporary_file_path, temporary_file_id) = create_temporary_file(Arc::clone(&app_config))
        .map_err(|err| {
        err.chain(format!(
            "Temporary file creation for file upload failed for {:?}/{}.",
            category, id
        ))
    })?;

    let mut connection = database_manager.database_connection().map_err(|err| {
        SeqError::from(err).chain(format!(
            "No database connection could be obtained when trying to add a file to {:?}/{}.",
            category, id
        ))
    })?;
    category
        .is_locked_err(id, pipelines, download_tracker, &mut connection)
        .map_err(|err| {
            err.chain(format!("The file upload request failed as {:?}/{} is locked.", category, id))
        })?;

    Ok((category, id, temporary_file_path, temporary_file_id, connection))
}

pub async fn post_add_file(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    pipelines: web::Data<LoadedPipelines>,
    download_tracker: web::Data<DownloadTrackerManager>,
    params: web::Path<(FileRequestCategory, i32)>,
    mut payload: Multipart,
) -> Result<HttpResponse, SeqError> {
    match setup_post_add_file(
        web::Data::clone(&app_config),
        database_manager,
        pipelines,
        download_tracker,
        params,
    )
    .await
    {
        Ok((category, id, temporary_file_path, temporary_file_id, mut connection)) => {
            persist_multipart(
                payload,
                id,
                category,
                temporary_file_path.as_path(),
                web::Data::clone(&app_config),
                &mut connection,
            )
            .await
            .map_err(|error| {
                // Delete temporary file on error.
                if delete_temporary_file(temporary_file_id, Arc::clone(&app_config)).is_err() {
                    log::error!(
                        "Failed to delete temporary file {} upon error {}.",
                        temporary_file_path.display(),
                        error
                    );
                }
                error
            })?;

            Ok(HttpResponse::Created().finish())
        },
        Err(err) => {
            mock_handle_payload_on_error(&mut payload).await;
            Err(err.chain("File upload validation failed before the payload was handled."))
        },
    }
}

pub async fn post_add_folder(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    pipelines: web::Data<LoadedPipelines>,
    download_tracker: web::Data<DownloadTrackerManager>,
    params: web::Path<(FileRequestCategory, i32)>,
    upload_info: web::Json<FilePath>,
) -> Result<HttpResponse, SeqError> {
    let (category, id) = params.into_inner();
    let upload_info = upload_info.into_inner();
    // Error if no path was specified.
    if upload_info.path_components.is_empty() {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            "The specified path is empty.",
            "The specified path is invalid.",
        ));
    }
    // Error if path is invalid.
    if let Err(validation_error) = upload_info.validate() {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            format!(
                "The file path {:?} is invalid and upload faild with error: {}.",
                upload_info, validation_error
            ),
            "The requested path is invalid.",
        ));
    }

    let mut connection = database_manager.database_connection().map_err(|err| {
        err.chain(format!(
            "Adding folder for category {:?} ID {} failed \
            as no database connection could be established.",
            category, id
        ))
    })?;

    category.entity_exists(id, &mut connection).map_err(|err| {
        err.chain(format!(
            "No folder for category {:?} ID {} can be added as it does not exist.",
            category, id
        ))
    })?;

    // Validate the according entity is not locked.
    category
        .is_locked_err(id, pipelines, download_tracker, &mut connection)
        .map_err(|err| {
            err.chain(format!(
                "The folder creation for {:?} request failed as {:?}/{} is locked.",
                upload_info, category, id
            ))
        })?;

    // Validate that the file path is not already existant.
    let full_path = category
        .base_path(app_config, id)
        .join(upload_info.file_path());
    if full_path.exists() {
        return Err(SeqError::new(
            "Conflicting request",
            SeqErrorType::Conflict,
            format!(
                "File at path {} in category {:?} does already exist.",
                upload_info.file_path().display(),
                category
            ),
            "The resource does already exist.",
        ));
    } else {
        std::fs::create_dir_all(&full_path).map_err(|err| {
            SeqError::from(err).chain(format!(
                "Creation of directory {} for category {:?} ID {} failed.",
                full_path.display(),
                category,
                id
            ))
        })?;
    }

    Ok(HttpResponse::Created().finish())
}

async fn persist_multipart<P: AsRef<Path>>(
    payload: Multipart,
    id: i32,
    category: FileRequestCategory,
    temporary_file_path: P,
    app_config: web::Data<Configuration>,
    connection: &mut SqliteConnection,
) -> Result<(), SeqError> {
    let (upload_info, temp_file_path) =
        parse_multipart_file::<FilePath, P>(payload, temporary_file_path).await?;

    // Error if no path was specified.
    if upload_info.path_components.is_empty() {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            "The specified path is empty.",
            "The specified path is invalid.",
        ));
    }

    // Validate the existance of the entity.
    category.entity_exists(id, connection).map_err(|err| {
        err.chain(format!(
            "Category {:?} ID {} does not exist. Persisting multipart data failed.",
            category, id
        ))
    })?;

    // Validate that the file path is not already existent.
    let full_path = category
        .base_path(web::Data::clone(&app_config), id)
        .join(upload_info.file_path());
    if full_path.exists() {
        return Err(SeqError::new(
            "Conflicting request",
            SeqErrorType::Conflict,
            format!(
                "File at path {} in category {:?} does already exist.",
                upload_info.file_path().display(),
                category
            ),
            "The resource does already exist.",
        ));
    }

    // Create folder and copy file to destination.
    temp_file_to_data_file(&full_path, &temp_file_path).map_err(|err| {
        err.chain(format!(
            "Moving temporary multipart file {} to \n
            final location {} failed for entitiy {:?}/{}.",
            temp_file_path.display(),
            full_path.display(),
            category,
            id
        ))
    })?;

    Ok(())
}

pub async fn get_experiment_archive_step_results(
    database_manager: web::Data<DatabaseManager>,
    app_config: web::Data<Configuration>,
    download_tracker_manager: web::Data<DownloadTrackerManager>,
    path_variables: web::Path<(i32, String)>,
) -> Result<HttpResponse, SeqError> {
    let (experiment_id, step_hash) = path_variables.into_inner();
    let mut connection = database_manager.database_connection().map_err(|err| {
        err.chain(format!(
            "Archiving step results for experiment {} step {} failed \
            as no database connection could be established.",
            experiment_id, step_hash
        ))
    })?;
    Experiment::exists_err(experiment_id, &mut connection).map_err(|err| {
        err.chain(format!(
            "Archiving step results for experiment {} step {} failed \
            as the experiment does not exist.",
            experiment_id, step_hash
        ))
    })?;

    // Map hash back to pipeline and step from the executions stored in the database.
    let (pipeline_id, step_id) =
        ExperimentExecution::get_by_experiment(experiment_id, &mut connection)
            .map_err(|err| {
                SeqError::from(err).chain(format!(
                    "Archiving step results for experiment {} step {} failed \
                    as the experiment execution could not be retrieved.",
                    experiment_id, step_hash
                ))
            })?
            .into_iter()
            .find(|execution| {
                step_hash
                    == Configuration::hash_pipeline_step_id(
                        &execution.pipeline_id,
                        &execution.pipeline_step_id,
                    )
            })
            .map(|execution| (execution.pipeline_id, execution.pipeline_step_id))
            .ok_or(SeqError::new(
                "Archive directory not found",
                SeqErrorType::NotFoundError,
                format!(
                    "Hash {} could not be mapped back into \
                    a pipeline step ID in experiment {}, \
                    indicating a missing output directory.",
                    step_hash, experiment_id
                ),
                "The directory to archive and download is not present.",
            ))?;
    let _download_tracker = download_tracker_manager.track_experiment_output_download_step(
        experiment_id,
        &pipeline_id,
        &step_id,
    );
    // // Defines the source path.
    let source = app_config.experiment_step_path(experiment_id.to_string(), &pipeline_id, &step_id);

    let file_name = format!("{}.zip", sanitize_filename::sanitize(&step_id));
    let archive_stream = ArchiveStream::new(source).map_err(|err| {
        err.chain(format!(
            "Archive stream generation failed for experiment {} step {} - {}.",
            experiment_id, pipeline_id, step_id
        ))
    })?;

    //Return the archive as stream.
    let content_disposition = header::ContentDisposition {
        disposition: header::DispositionType::Attachment,
        parameters: vec![header::DispositionParam::FilenameExt(
            header::ExtendedValue {
                charset: header::Charset::Ext(String::from("UTF-8")),
                language_tag: None,
                value: file_name.into_bytes(),
            },
        )],
    };

    Ok(HttpResponse::Ok()
        .content_type("application/zip")
        .insert_header((header::CONTENT_DISPOSITION, content_disposition.to_string()))
        .streaming(archive_stream))
}

pub async fn get_pipeline_attachment(
    app_config: web::Data<Configuration>,
    info: web::Path<(String, String)>,
) -> Result<NamedFile, SeqError> {
    let (pipeline_directory, attachment_name) = info.into_inner();

    // If the attachment name or pipeline folder contain a path seperator
    // return an error as this allows attacks by using relative components.
    if PathBuf::from(&pipeline_directory).iter().count() != 1 {
        return Err(SeqError::new(
            "Bad request",
            SeqErrorType::BadRequestError,
            format!(
                "The pipeline directory {} is invalid and is probalbly supposed to compromise the system.",
                pipeline_directory
            ),
            "Invalid pipeline ID.",
        ));
    }
    if PathBuf::from(&attachment_name).iter().count() != 1 {
        return Err(SeqError::new(
            "Bad request",
            SeqErrorType::BadRequestError,
            format!(
                "The attachment name {} is invalid and is probalbly supposed to compromise the system.",
                attachment_name
            ),
            "Invalid attachment name.",
        ));
    }

    let path = app_config.pipeline_attachment_path(&pipeline_directory, &attachment_name);

    if !path.exists() {
        return Err(SeqError::new(
            "Not found",
            SeqErrorType::NotFoundError,
            format!("Attachment at path {} does not exist.", path.display()),
            "File not found.",
        ));
    }

    Ok(NamedFile::open(&path).map_err(|err| {
        SeqError::from(err).chain(format!(
            "Getting pipeline attachement {} from directory {}\
            failed as it could not be opened.",
            attachment_name, pipeline_directory
        ))
    })?)
}

fn temp_file_to_data_file<P: AsRef<Path>, Q: AsRef<Path>>(
    file_path: P,
    temp_file_path: Q,
) -> Result<(), SeqError> {
    let file_path = file_path.as_ref();
    log::info!(
        "Saving temporary file {} to {}.",
        temp_file_path.as_ref().display(),
        file_path.display()
    );
    std::fs::create_dir_all(&file_path.parent().ok_or_else(|| {
        SeqError::new(
            "Invalid file path",
            SeqErrorType::BadRequestError,
            format!("The file path {} is not a valid path.", file_path.display()),
            "The file path is invalid.",
        )
    })?)
    .map_err(|err| {
        SeqError::from(err).chain(format!(
            "Moving temporary file {} to {} failed \
            as the directory structure could not be created.",
            temp_file_path.as_ref().display(),
            file_path.display()
        ))
    })?;
    std::fs::rename(&temp_file_path, &file_path).map_err(|err| {
        SeqError::from(err).chain(format!(
            "Moving temporary file {} to {} failed \
            as the file could not be moved.",
            temp_file_path.as_ref().display(),
            file_path.display()
        ))
    })?;
    Ok(())
}

#[cfg(test)]
mod global_data_tests;

#[cfg(test)]
mod experiment_tests;
