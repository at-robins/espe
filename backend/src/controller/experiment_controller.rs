use std::sync::Arc;

use crate::{
    application::{config::Configuration, error::SeqError, error::SeqErrorType},
    diesel::RunQueryDsl,
    model::db::experiment::NewExperiment,
};
use actix_web::{HttpRequest, HttpResponse};

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
    log::info!("Creating experiment repository with name {}.", &name);
    let new_record = NewExperiment::new(name);
    let inserted_id: i32 = diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .returning(crate::schema::experiment::id)
        .get_result(&mut connection)?;
    // Return the ID of the created experiment.
    Ok(HttpResponse::Created().json(inserted_id))
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
