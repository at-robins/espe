use std::{
    path::{Path, PathBuf},
    sync::Arc,
};

use crate::{
    application::{
        config::Configuration,
        error::{SeqError, SeqErrorType},
    },
    diesel::{ExpressionMethods, QueryDsl, RunQueryDsl},
    model::{
        db::global_data::{GlobalData, NewGlobalData},
        exchange::global_data_details::{
            GlobalDataDetails, GlobalDataFileDetails, GlobalDataFilePath,
        },
    },
    service::multipart_service::{
        create_temporary_file, delete_temporary_file, parse_multipart_file,
    },
};
use actix_multipart::Multipart;
use actix_web::{web, HttpRequest, HttpResponse};

/// The maximum length a global data title is allowed to have.
const MAXIMUM_TITLE_LENGTH: usize = 512;

pub async fn create_global_data(
    request: HttpRequest,
    name: actix_web::web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let name: String = name.into_inner();
    validate_global_data_name(&name)?;
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    log::info!("Creating global data repository with name {}.", &name);
    let new_record = NewGlobalData::new(name, None);
    let inserted_id: i32 = diesel::insert_into(crate::schema::global_data::table)
        .values(&new_record)
        .returning(crate::schema::global_data::id)
        .get_result(&mut connection)?;
    // Return the ID of the created global data.
    Ok(HttpResponse::Created().json(inserted_id))
}

pub async fn delete_global_data(
    request: HttpRequest,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    GlobalData::exists_err(id, &mut connection)?;
    log::info!("Deleting global data repository with ID {}.", id);
    // Remove all files belonging to the global data repository.
    let global_path = app_config.global_data_path(id.to_string());
    if global_path.exists() {
        std::fs::remove_dir_all(global_path)?;
    }
    // Delete the repository from the database.
    connection.immediate_transaction(|connection| {
        diesel::delete(crate::schema::global_data::table)
            .filter(crate::schema::global_data::id.eq(id))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn get_global_data(
    request: HttpRequest,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    GlobalData::exists_err(id, &mut connection)?;
    let global_repo_details: GlobalDataDetails = crate::schema::global_data::table
        .find(id)
        .first::<GlobalData>(&mut connection)?
        .into();
    Ok(HttpResponse::Ok().json(global_repo_details))
}

pub async fn get_global_data_files(
    request: HttpRequest,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    GlobalData::exists_err(id, &mut connection)?;

    let mut global_repo_files: Vec<GlobalDataFileDetails> = Vec::new();
    let global_data_path = app_config.global_data_path(id.to_string());
    if global_data_path.exists() {
        for entry in walkdir::WalkDir::new(&global_data_path) {
            let entry = entry?;
            let relative_path = entry
                .path()
                .strip_prefix(&global_data_path)
                .expect("The global data directory must be a parent of the contained components.");
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
                global_repo_files.push(GlobalDataFileDetails {
                    path_components: components,
                    is_file: entry.path().is_file(),
                });
            }
        }
        global_repo_files.sort();
    }
    Ok(HttpResponse::Ok().json(global_repo_files))
}

pub async fn patch_global_data_name(
    request: HttpRequest,
    id: web::Path<i32>,
    new_name: web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_name = new_name.into_inner();
    validate_global_data_name(&new_name)?;
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    GlobalData::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::global_data::table.find(id))
            .set(crate::schema::global_data::global_data_name.eq(new_name))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_global_data_comment(
    request: HttpRequest,
    id: web::Path<i32>,
    new_comment: web::Json<Option<String>>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Sanitise the HTML.
    let new_comment = new_comment.into_inner().map(|inner| ammonia::clean(&inner));
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    GlobalData::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::global_data::table.find(id))
            .set(crate::schema::global_data::comment.eq(new_comment))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn list_global_data(
    request: HttpRequest,
) -> Result<web::Json<Vec<GlobalDataDetails>>, SeqError> {
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    let global_repos: Vec<GlobalDataDetails> = GlobalData::get_all(&mut connection)?
        .into_iter()
        .map(|val| val.into())
        .collect();
    Ok(web::Json(global_repos))
}

pub async fn delete_global_data_files_by_path(
    request: HttpRequest,
    id: web::Path<i32>,
    path: web::Json<GlobalDataFilePath>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let delete_info = path.into_inner();
    let delete_path = delete_info.file_path();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    GlobalData::exists_err(id, &mut connection)?;

    let global_data_path: PathBuf = app_config.global_data_path(id.to_string());
    let full_path: PathBuf = global_data_path.join(&delete_path);
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

pub async fn post_global_data_add_file(
    request: HttpRequest,
    id: web::Path<i32>,
    payload: Multipart,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");

    let (temporary_file_path, temporary_file_id) = create_temporary_file(Arc::clone(&app_config))?;

    persist_multipart(payload, id, temporary_file_path.as_path(), Arc::clone(app_config))
        .await
        .map_err(|error| {
            // Delete temporary file on error.
            if delete_temporary_file(temporary_file_id, Arc::clone(app_config)).is_err() {
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

pub async fn post_global_data_add_folder(
    request: HttpRequest,
    id: web::Path<i32>,
    upload_info: web::Json<GlobalDataFilePath>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let upload_info = upload_info.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;

    // Validate the existance of the global data repository.
    GlobalData::exists_err(id, &mut connection)?;

    // Validate that the file path is not already existant.
    let full_path = app_config
        .global_data_path(id.to_string())
        .join(upload_info.file_path());
    if full_path.exists() {
        return Err(SeqError::new(
            "Conflicting request",
            SeqErrorType::Conflict,
            format!(
                "Global data file at path {} does already exist.",
                upload_info.file_path().display()
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
    global_data_id: i32,
    temporary_file_path: P,
    app_config: Arc<Configuration>,
) -> Result<(), SeqError> {
    let mut connection = app_config.database_connection()?;

    let (upload_info, temp_file_path) =
        parse_multipart_file::<GlobalDataFilePath, P>(payload, temporary_file_path).await?;

    // Validate the existance of the global data repository.
    GlobalData::exists_err(global_data_id, &mut connection)?;

    // Validate that the file path is not already existant.
    let full_path = app_config
        .global_data_path(global_data_id.to_string())
        .join(upload_info.file_path());
    if full_path.exists() {
        return Err(SeqError::new(
            "Conflicting request",
            SeqErrorType::Conflict,
            format!(
                "Global data file at path {} does already exist.",
                upload_info.file_path().display()
            ),
            "The resource does already exist.",
        ));
    }

    // Create folder and copy file to destination.
    temp_file_to_global_data(global_data_id, upload_info.file_path(), temp_file_path, app_config)?;

    Ok(())
}

fn temp_file_to_global_data<P: AsRef<Path>, Q: AsRef<Path>>(
    global_data_id: i32,
    file_path: P,
    temp_file_path: Q,
    app_config: Arc<Configuration>,
) -> Result<(), SeqError> {
    let mut final_file_path: PathBuf = app_config.global_data_path(global_data_id.to_string());
    final_file_path.push(file_path);
    log::info!(
        "Saving temporary file {} to {}.",
        temp_file_path.as_ref().display(),
        final_file_path.display()
    );
    std::fs::create_dir_all(&final_file_path.parent().ok_or_else(|| {
        SeqError::new(
            "Invalid file path",
            SeqErrorType::BadRequestError,
            format!("The file path {} is not a valid path.", final_file_path.display()),
            "The file path is invalid.",
        )
    })?)?;
    std::fs::rename(temp_file_path, final_file_path)?;
    Ok(())
}

fn validate_global_data_name<T: AsRef<str>>(name: T) -> Result<(), SeqError> {
    let name: &str = name.as_ref();
    if name.is_empty() {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            "The title may not be empty.",
            "The title is invalid.",
        ));
    }
    if name.len() > MAXIMUM_TITLE_LENGTH {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            format!("The title {} exceeds the limit of {} characters.", name, MAXIMUM_TITLE_LENGTH),
            "The title is invalid.",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
