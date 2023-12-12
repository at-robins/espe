use actix_web::{web, HttpResponse};

use crate::{
    application::{
        config::Configuration,
        error::{SeqError, SeqErrorType},
    },
    model::internal::pipeline_blueprint::PipelineBlueprint,
    service::pipeline_service::LoadedPipelines,
};

/// Return all pipeline blueprints that are currently loaded.
pub async fn get_pipeline_blueprints(
    pipelines: web::Data<LoadedPipelines>,
) -> Result<web::Json<Vec<PipelineBlueprint>>, SeqError> {
    Ok(web::Json(
        pipelines
            .pipelines()
            .iter()
            .map(|pipeline| pipeline.pipeline().clone())
            .collect::<Vec<PipelineBlueprint>>(),
    ))
}

/// Return all pipeline blueprints that are currently loaded.
pub async fn get_pipeline_blueprint(
    pipelines: web::Data<LoadedPipelines>,
    id: web::Json<String>,
) -> Result<web::Json<PipelineBlueprint>, SeqError> {
    let id = id.into_inner();
    if let Some(pipeline) = pipelines.get(&id) {
        Ok(web::Json(pipeline.pipeline().clone()))
    } else {
        Err(SeqError::new(
            "Not Fount",
            SeqErrorType::NotFoundError,
            format!("No pipeline with ID {} is loaded.", id),
            "Invalid ID.",
        ))
    }
}

/// Update the currently loaded pipeline blueprints.
pub async fn patch_pipeline_blueprints(
    pipelines: web::Data<LoadedPipelines>,
    app_config: web::Data<Configuration>,
) -> Result<HttpResponse, SeqError> {
    pipelines.update_loaded_pipelines(app_config)?;
    Ok(HttpResponse::Ok().finish())
}
