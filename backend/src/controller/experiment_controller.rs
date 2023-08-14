use crate::{
    application::{
        config::Configuration,
        error::{SeqError, SeqErrorType},
    },
    diesel::RunQueryDsl,
    model::{
        db::experiment::{Experiment, NewExperiment},
        exchange::{
            experiment_details::ExperimentDetails, experiment_pipeline::ExperimentPipelineBlueprint,
        },
    },
    service::{
        pipeline_service::LoadedPipelines,
        validation_service::{validate_comment, validate_entity_name, validate_mail},
    },
};
use actix_web::{web, HttpResponse};
use diesel::{ExpressionMethods, QueryDsl};

pub async fn create_experiment(
    app_config: web::Data<Configuration>,
    name: actix_web::web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let name: String = name.into_inner();
    validate_entity_name(&name)?;
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
    app_config: web::Data<Configuration>,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Retrieve the app config.
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
    app_config: web::Data<Configuration>,
    id: web::Path<i32>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    let experiment_details: ExperimentDetails = crate::schema::experiment::table
        .find(id)
        .first::<Experiment>(&mut connection)?
        .into();
    Ok(HttpResponse::Ok().json(experiment_details))
}

pub async fn patch_experiment_name(
    app_config: web::Data<Configuration>,
    id: web::Path<i32>,
    new_name: web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_name = new_name.into_inner();
    validate_entity_name(&new_name)?;
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::experiment_name.eq(new_name))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_experiment_mail(
    app_config: web::Data<Configuration>,
    id: web::Path<i32>,
    new_name: web::Json<String>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_mail = new_name.into_inner();
    validate_mail(&new_mail)?;
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::mail.eq(new_mail))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_experiment_comment(
    app_config: web::Data<Configuration>,
    id: web::Path<i32>,
    new_comment: web::Json<Option<String>>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    // Sanitise the HTML and validate.
    let new_comment = new_comment.into_inner().map(|inner| ammonia::clean(&inner));
    if let Some(inner) = &new_comment {
        validate_comment(inner)?;
    }
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::comment.eq(new_comment))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn patch_experiment_pipeline(
    app_config: web::Data<Configuration>,
    pipelines: web::Data<LoadedPipelines>,
    id: web::Path<i32>,
    new_pipeline: web::Json<Option<String>>,
) -> Result<HttpResponse, SeqError> {
    let id: i32 = id.into_inner();
    let new_pipeline = new_pipeline.into_inner();
    if let Some(new_pipeline_id) = &new_pipeline {
        if !pipelines.is_loaded(new_pipeline_id) {
            return Err(SeqError::new(
                "Not Found",
                SeqErrorType::NotFoundError,
                format!("No pipeline with ID {} is currently loaded.", new_pipeline_id),
                "The pipeline ID is invalid.",
            ));
        }
    }
    let mut connection = app_config.database_connection()?;
    Experiment::exists_err(id, &mut connection)?;
    connection.immediate_transaction(|connection| {
        diesel::update(crate::schema::experiment::table.find(id))
            .set(crate::schema::experiment::pipeline_id.eq(new_pipeline))
            .execute(connection)
    })?;
    Ok(HttpResponse::Ok().finish())
}

pub async fn list_experiment(
    app_config: web::Data<Configuration>,
) -> Result<web::Json<Vec<ExperimentDetails>>, SeqError> {
    let mut connection = app_config.database_connection()?;
    let experiments: Vec<ExperimentDetails> = Experiment::get_all(&mut connection)?
        .into_iter()
        .map(|val| val.into())
        .collect();
    Ok(web::Json(experiments))
}

pub async fn get_experiment_pipelines(
    app_config: web::Data<Configuration>,
    pipelines: web::Data<LoadedPipelines>,
    id: web::Path<i32>,
) -> Result<web::Json<Vec<ExperimentPipelineBlueprint>>, SeqError> {
    let experiment_id: i32 = id.into_inner();
    let mut connection = app_config.database_connection()?;
    let mut experiment_pipelines = Vec::new();
    for pipeline in pipelines.pipelines() {
        let values = crate::model::db::pipeline_step_variable::PipelineStepVariable::get_values_by_experiment_and_pipeline(experiment_id, pipeline.pipeline().id(), &mut connection)?;
        experiment_pipelines
            .push(ExperimentPipelineBlueprint::from_internal(pipeline.pipeline(), values));
    }
    Ok(web::Json(experiment_pipelines))
}

#[cfg(test)]
mod tests;
