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
        db::global_data::{GlobalData, NewGlobalData, NewGlobalDataFile},
        exchange::{
            global_data_details::{GlobalDataDetails, GlobalDataFileDetails},
            global_data_file_upload::GlobalDataFileUpload,
        },
    },
    service::multipart_service::{
        create_temporary_file, delete_temporary_file, parse_multipart_file,
    },
};
use actix_multipart::Multipart;
use actix_web::{web, HttpRequest, HttpResponse};
use diesel::dsl::exists;

pub async fn create_global_data(
    request: HttpRequest,
    name: actix_web::web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let name: String = name.into_inner();
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
    let global_repo_files: Vec<GlobalDataFileDetails> = GlobalData::files(id, &mut connection)?
        .into_iter()
        .map(|value| value.into())
        .collect();
    Ok(HttpResponse::Ok().json(global_repo_files))
}

pub async fn patch_global_data_name(
    request: HttpRequest,
    id: web::Path<i32>,
    new_name: web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_name = new_name.into_inner();
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
    let new_comment = new_comment.into_inner();
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

    let inserted_id =
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

    // Return the ID of the created attachment.
    Ok(HttpResponse::Created().json(inserted_id))
}

async fn persist_multipart<P: AsRef<Path>>(
    payload: Multipart,
    global_data_id: i32,
    temporary_file_path: P,
    app_config: Arc<Configuration>,
) -> Result<i32, SeqError> {
    let mut connection = app_config.database_connection()?;

    let (upload_info, temp_file_path) =
        parse_multipart_file::<GlobalDataFileUpload, P>(payload, temporary_file_path).await?;

    // Validate the existance of the global data repository.
    GlobalData::exists_err(global_data_id, &mut connection)?;

    // Write to database.
    let new_global_data_file: NewGlobalDataFile =
        NewGlobalDataFile::new(global_data_id, upload_info.file_path_as_string());
    let inserted_id = diesel::insert_into(crate::schema::global_data_file::table)
        .values(&new_global_data_file)
        .returning(crate::schema::global_data_file::id)
        .get_result(&mut connection)?;

    // Create folder and copy file to destination.
    temp_file_to_global_data(global_data_id, upload_info.file_path(), temp_file_path, app_config).map_err(|error| {
        // Roll back database if there is an error while moving the temporary file.
        if let Err(e) = diesel::delete(
            crate::schema::global_data_file::dsl::global_data_file
                .filter(crate::schema::global_data_file::id.eq(inserted_id)),
        )
        .execute(&mut connection)
        {
            log::error!(
                "Roll back of database after insertion of global data file {} failed with error: {}.",
                inserted_id,
                e
            );
        }
        error
    })?;

    // Return the ID of the created attachment.
    Ok(inserted_id)
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
