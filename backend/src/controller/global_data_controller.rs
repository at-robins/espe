use std::sync::Arc;

use crate::{
    application::{config::Configuration, error::SeqError},
    diesel::{ExpressionMethods, QueryDsl, RunQueryDsl},
    model::{
        db::global_data::{GlobalData, NewGlobalData},
        exchange::global_data_details::GlobalDataDetails,
    },
    service::validation_service::{validate_comment, validate_entity_name},
};
use actix_web::{web, HttpRequest, HttpResponse};

pub async fn create_global_data(
    request: HttpRequest,
    name: actix_web::web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let name: String = name.into_inner();
    validate_entity_name(&name)?;
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    log::info!("Creating global data repository with name {}.", &name);
    let new_record = NewGlobalData::new(name);
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
    validate_entity_name(&new_name)?;
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
    // Sanitise the HTML and validate.
    let new_comment = new_comment.into_inner().map(|inner| ammonia::clean(&inner));
    if let Some(inner) = &new_comment {
        validate_comment(inner)?;
    }
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

#[cfg(test)]
mod tests;
