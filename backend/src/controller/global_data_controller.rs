use std::sync::Arc;

use crate::{
    application::{
        config::Configuration,
        error::{SeqError, SeqErrorType},
    },
    diesel::{ExpressionMethods, QueryDsl, RunQueryDsl},
    model::{
        db::global_data::{GlobalData, NewGlobalData},
        exchange::global_data_details::GlobalDataDetails,
    },
};
use actix_web::{web, HttpRequest, HttpResponse};

/// The maximum length a global data name is allowed to have.
const MAXIMUM_NAME_LENGTH: usize = 512;

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

fn validate_global_data_name<T: AsRef<str>>(name: T) -> Result<(), SeqError> {
    let name: &str = name.as_ref();
    if name.is_empty() {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            "The name may not be empty.",
            "The name is invalid.",
        ));
    }
    if name.len() > MAXIMUM_NAME_LENGTH {
        return Err(SeqError::new(
            "Invalid request",
            SeqErrorType::BadRequestError,
            format!("The name {} exceeds the limit of {} characters.", name, MAXIMUM_NAME_LENGTH),
            "The name is invalid.",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
