use actix_web::{web, HttpResponse, Responder};
use chrono::Utc;
use rand::{thread_rng, Rng};

use crate::{
    application::{
        config::Configuration,
        error::{SeqError, SeqErrorType},
    },
    model::{
        exchange::{
            pipeline_blueprint_details::PipelineBlueprintDetails,
            pipeline_step_details::PipelineStepDetails,
        },
        internal::step::PipelineStepStatus,
    },
    service::pipeline_service::LoadedPipelines,
};

pub async fn get_pipeline_instance(wrapped_id: web::Path<u64>) -> Result<impl Responder, SeqError> {
    let _id = wrapped_id.into_inner();
    let dummy_response: Vec<PipelineStepDetails> = vec![
        PipelineStepDetails {
            id: 0,
            name: "FASTQ-QC".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 1,
            name: "Trimming".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 2,
            name: "FASTQ-QC".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 3,
            name: "Alignment".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 4,
            name: "SAM-BAM-conversion".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 5,
            name: "BAM indexing".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 6,
            name: "BAM sorting".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 7,
            name: "Remove duplicates".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 8,
            name: "Remove mitochondrial reads".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 9,
            name: "Remove blacklisted regions".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 10,
            name: "Shift reads".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 11,
            name: "Post-alignment QC".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
        PipelineStepDetails {
            id: 12,
            name: "Peak calling".to_string(),
            status: random_status().to_string(),
            creation_time: Utc::now(),
        },
    ];
    Ok(web::Json(dummy_response))
}

/// Return all pipeline blueprints that are currently loaded.
pub async fn get_pipeline_blueprints(
    pipelines: web::Data<LoadedPipelines>,
) -> Result<web::Json<Vec<PipelineBlueprintDetails>>, SeqError> {
    Ok(web::Json(
        pipelines
            .pipelines()
            .iter()
            .map(|pipeline| PipelineBlueprintDetails::from(pipeline))
            .collect::<Vec<PipelineBlueprintDetails>>(),
    ))
}

/// Return all pipeline blueprints that are currently loaded.
pub async fn get_pipeline_blueprint(
    pipelines: web::Data<LoadedPipelines>,
    id: web::Json<String>,
) -> Result<web::Json<PipelineBlueprintDetails>, SeqError> {
    let id = id.into_inner();
    if let Some(pipeline) = pipelines.get(&id) {
        Ok(web::Json(PipelineBlueprintDetails::from(&pipeline)))
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

fn random_status() -> PipelineStepStatus {
    match thread_rng().gen_range(0..4) {
        0 => PipelineStepStatus::Failed,
        1 => PipelineStepStatus::Pending,
        2 => PipelineStepStatus::Running,
        3 => PipelineStepStatus::Success,
        _ => PipelineStepStatus::Failed,
    }
}
