use std::sync::Arc;

use actix_web::{web, HttpRequest, Responder};
use chrono::Utc;
use diesel::{RunQueryDsl, OptionalExtension};
use rand::{thread_rng, Rng};

use crate::{
    application::{
        config::Configuration,
        error::{InternalError, SeqError},
    },
    model::{
        exchange::pipeline_step_details::PipelineStepDetails, internal::step::PipelineStepStatus, db::pipeline::Pipeline,
    },
};

pub async fn get_pipeline_instance(wrapped_id: web::Path<u64>) -> Result<impl Responder, SeqError> {
    let id = wrapped_id.into_inner();
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

pub async fn get_pipeline_blueprints(request: HttpRequest) -> Result<impl Responder, SeqError> {
    // Retrieve the app config.
    let app_config = SeqError::log_error(request.app_data::<Arc<Configuration>>().ok_or_else(|| {
        SeqError::InternalServerError(InternalError::new(
            "Configuration",
            "The server configuration could not be accessed.",
            "Missing configuration.",
        ))
    }))?;
    let connection = SeqError::log_error(app_config.database_connection())?;
    let pipelines = SeqError::log_error(crate::schema::pipeline::table.load::<Pipeline>(&connection))?;
    Ok(web::Json(pipelines.iter().map(|pipeline| pipeline.id).collect::<Vec<i32>>()))
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
