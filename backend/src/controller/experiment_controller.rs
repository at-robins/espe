use std::sync::Arc;

use crate::{
    application::{config::Configuration, error::SeqError, error::SeqErrorType},
    diesel::RunQueryDsl,
    model::{db::experiment::{NewExperiment, Experiment}, exchange::experiment_details::ExperimentDetails},
};
use actix_web::{HttpRequest, HttpResponse, web};
use diesel::{ExpressionMethods, QueryDsl};

/// The maximum length an experiment name is allowed to have.
const MAXIMUM_NAME_LENGTH: usize = 512;

pub async fn create_experiment(
    request: HttpRequest,
    name: actix_web::web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let name: String = name.into_inner();
    validate_experiment_name(&name)?;
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    log::info!("Creating experiment with name {}.", &name);
    let new_record = NewExperiment::new(name);
    let inserted_id: i32 = diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .returning(crate::schema::experiment::id)
        .get_result(&mut connection)?;
    // Return the ID of the created experiment.
    Ok(HttpResponse::Created().json(inserted_id))
}

pub async fn delete_experiment(
    request: HttpRequest,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    log::info!("Deleting experiment with ID {}.", id);
    // Remove all files belonging to the experiment.
    let experiment_path = app_config.experiment_path(id.to_string());
    if experiment_path.exists() {
        std::fs::remove_dir_all(experiment_path)?;
    }
    // Delete the experiment from the database.
    connection.immediate_transaction(|connection| {
        diesel::delete(crate::schema::experiment::table)
            .filter(crate::schema::experiment::id.eq(id))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn get_experiment(
    request: HttpRequest,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    let experiment_details: ExperimentDetails = crate::schema::experiment::table
        .find(id)
        .first::<Experiment>(&mut connection)?
        .into();
    Ok(HttpResponse::Ok().json(experiment_details))
}

pub async fn patch_experiment_name(
    request: HttpRequest,
    id: web::Path<i32>,
    new_name: web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_name = new_name.into_inner();
    validate_experiment_name(&new_name)?;
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::experiment_name.eq(new_name))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_experiment_comment(
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
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::comment.eq(new_comment))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn list_experiment(
    request: HttpRequest,
) -> Result<web::Json<Vec<ExperimentDetails>>, SeqError> {
    // Retrieve the app config.
    let app_config = request
        .app_data::<Arc<Configuration>>()
        .expect("The configuration must be accessible.");
    let mut connection = app_config.database_connection()?;
    let experiments: Vec<ExperimentDetails> = Experiment::get_all(&mut connection)?
        .into_iter()
        .map(|val| val.into())
        .collect();
    Ok(web::Json(experiments))
}

fn validate_experiment_name<T: AsRef<str>>(name: T) -> Result<(), SeqError> {
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
