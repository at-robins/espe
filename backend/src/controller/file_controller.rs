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
        db::{experiment::Experiment, global_data::GlobalData},
        exchange::file_path::{FileDetails, FilePath},
        internal::archive::ArchiveMetadata,
    },
    service::multipart_service::{
        create_temporary_file, delete_temporary_file, parse_multipart_file, UploadForm,
    },
};
use actix_files::NamedFile;
use actix_multipart::Multipart;
use actix_web::{web, HttpResponse};
use diesel::SqliteConnection;
use serde::{Deserialize, Serialize};
use zip_extensions::zip_create_from_directory;

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
    /// * `connection` - a databse connection
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
}

pub async fn get_files(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    params: web::Path<(FileRequestCategory, i32)>,
) -> Result<HttpResponse, SeqError> {
    let (category, id) = params.into_inner();
    let mut connection = database_manager.database_connection()?;
    category.entity_exists(id, &mut connection)?;

    let mut all_files: Vec<FileDetails> = Vec::new();
    let data_path = category.base_path(app_config, id);
    if data_path.exists() {
        for entry in walkdir::WalkDir::new(&data_path) {
            let entry = entry?;
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
    params: web::Path<(FileRequestCategory, i32)>,
    path: web::Json<FilePath>,
) -> Result<HttpResponse, SeqError> {
    let (category, id) = params.into_inner();
    let delete_info = path.into_inner();
    let delete_path = delete_info.file_path();
    let mut connection = database_manager.database_connection()?;
    category.entity_exists(id, &mut connection)?;
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
        std::fs::remove_dir_all(full_path)?;
    } else {
        std::fs::remove_file(full_path)?;
    }
    Ok(HttpResponse::Ok().finish())
}

pub async fn post_add_file(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
    params: web::Path<(FileRequestCategory, i32)>,
    payload: Multipart,
) -> Result<HttpResponse, SeqError> {
    let (category, id) = params.into_inner();

    let (temporary_file_path, temporary_file_id) = create_temporary_file(Arc::clone(&app_config))?;

    persist_multipart(
        payload,
        id,
        category,
        temporary_file_path.as_path(),
        web::Data::clone(&app_config),
        web::Data::clone(&database_manager),
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
}

pub async fn post_add_folder(
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
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

    let mut connection = database_manager.database_connection()?;

    // Validate the existance of the entity.
    category.entity_exists(id, &mut connection)?;

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
        std::fs::create_dir_all(full_path)?;
    }

    Ok(HttpResponse::Created().finish())
}

async fn persist_multipart<P: AsRef<Path>>(
    payload: Multipart,
    id: i32,
    category: FileRequestCategory,
    temporary_file_path: P,
    app_config: web::Data<Configuration>,
    database_manager: web::Data<DatabaseManager>,
) -> Result<(), SeqError> {
    let mut connection = database_manager.database_connection()?;

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
    category.entity_exists(id, &mut connection)?;

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
    temp_file_to_data_file(full_path, temp_file_path)?;

    Ok(())
}

pub async fn post_experiment_archive_step_results(
    database_manager: web::Data<DatabaseManager>,
    app_config: web::Data<Configuration>,
    experiment_id: web::Path<i32>,
    step_id: web::Json<String>,
) -> Result<String, SeqError> {
    let experiment_id: i32 = experiment_id.into_inner();
    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(experiment_id, &mut connection)?;

    // Sets up the required information.
    let step_id: String = step_id.into_inner();
    let archive_id = Configuration::generate_uuid();
    let archive_meta = ArchiveMetadata::new(format!("{}.zip", &step_id));

    // Defines source and target paths.
    let source = app_config.experiment_step_path(experiment_id.to_string(), &step_id);
    let target = app_config.temporary_download_file_path(archive_id);
    let target_meta = ArchiveMetadata::metadata_path(&target);

    // Creates the parent directory if necessary.
    if let Some(target_parent) = target.parent() {
        std::fs::create_dir_all(target_parent)?;
    }
    // Creates the archive.
    zip_create_from_directory(&target, &source).map_err(|err| {
        SeqError::new(
            "Archiving error",
            SeqErrorType::InternalServerError,
            format!("Creation of a downloadable archive for experiment {} ({}) from {} to {} failed with error: {}", experiment_id, step_id, source.display(), target.display(), err),
            "Downloadable archive could not be created.",
        )
    })?;
    // Creates the archive metadata.
    serde_json::to_writer(std::fs::File::create(target_meta)?, &archive_meta)?;

    //Return the archive ID.
    Ok(archive_id.to_string())
}

pub async fn get_experiment_download_step_results(
    database_manager: web::Data<DatabaseManager>,
    app_config: web::Data<Configuration>,
    info: web::Path<(i32, String)>,
) -> Result<NamedFile, SeqError> {
    let (experiment_id, archive_id) = info.into_inner();

    // If the archive ID contains a path seperator return an error
    // as it is not a valid ID and allows attacks by using relative
    // components.
    if archive_id.contains(std::path::MAIN_SEPARATOR) {
        return Err(SeqError::new(
            "Bad request",
            SeqErrorType::BadRequestError,
            format!(
                "The archive id {} is invalid and is probalbly supposed to compromise the system.",
                archive_id
            ),
            "Invalid archive ID.",
        ));
    }

    let mut connection = database_manager.database_connection()?;
    Experiment::exists_err(experiment_id, &mut connection)?;

    let archive_path = app_config.temporary_download_file_path(archive_id);
    if !archive_path.exists() {
        return Err(SeqError::new(
            "Not found",
            SeqErrorType::NotFoundError,
            format!("Archive file at path {} does not exist.", archive_path.display()),
            "File not found.",
        ));
    }

    let archive_meta_path = ArchiveMetadata::metadata_path(&archive_path);
    if !archive_meta_path.exists() {
        return Err(SeqError::new(
            "Not found",
            SeqErrorType::NotFoundError,
            format!(
                "Archive metadata file at path {} does not exist.",
                archive_meta_path.display()
            ),
            "File not found.",
        ));
    }

    let archive_meta: ArchiveMetadata =
        serde_json::from_reader(std::fs::File::open(&archive_meta_path)?)?;

    Ok(NamedFile::from_file(std::fs::File::open(&archive_path)?, archive_meta.file_name())?)
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
    })?)?;
    std::fs::rename(temp_file_path, file_path)?;
    Ok(())
}

#[cfg(test)]
mod global_data_tests;

#[cfg(test)]
mod experiment_tests;
