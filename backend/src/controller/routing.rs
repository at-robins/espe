use actix_files::{Files, NamedFile};
use actix_web::web::{self, ServiceConfig};
use log::error;

use crate::application::error::SeqError;

use super::{
    experiment_controller::{
        create_experiment, delete_experiment, get_experiment, get_experiment_locked,
        get_experiment_execution_status, get_experiment_pipeline_run, get_experiment_pipelines,
        list_experiment, patch_experiment_comment, patch_experiment_mail, patch_experiment_name,
        patch_experiment_pipeline, post_execute_experiment, post_execute_experiment_step,
        post_experiment_execution_abort, post_experiment_execution_reset,
        post_experiment_pipeline_global_variable, post_experiment_pipeline_step_variable,
    },
    file_controller::{
        delete_files_by_path, get_experiment_download_step_results, get_files,
        get_pipeline_attachment, post_add_file, post_add_folder,
        post_experiment_archive_step_results,
    },
    global_data_controller::{
        create_global_data, delete_global_data, get_global_data, list_global_data,
        patch_global_data_comment, patch_global_data_name,
    },
    log_controller::get_experiment_step_logs,
    pipeline_controller::{
        get_pipeline_blueprint, get_pipeline_blueprints, patch_pipeline_blueprints,
    },
};

/// Serve the entry point to the single page web app.
async fn index() -> Result<NamedFile, SeqError> {
    NamedFile::open("./static_dist/index.html").map_err(|error| {
        let internal_error: SeqError = error.into();
        error!("{}", internal_error);
        internal_error
    })
}

/// Configures the routing.
pub fn routing_config(cfg: &mut ServiceConfig) {
    cfg

    // UI redirect
    .route("/", web::get().to(index))
    .route("/ui", web::get().to(index))
    .route("/ui/{rest:.*}", web::get().to(index))
    // Pipelines
    .route("/api/pipelines/attachments/{pipeline}/{attachment}", web::get().to(get_pipeline_attachment))
    .route("/api/pipelines/blueprint", web::get().to(get_pipeline_blueprint))
    .route("/api/pipelines/blueprints", web::get().to(get_pipeline_blueprints))
    .route("/api/pipelines/blueprints", web::patch().to(patch_pipeline_blueprints))
    // Experiments
    .route("/api/experiments", web::get().to(list_experiment))
    .route("/api/experiments", web::post().to(create_experiment))
    .route("/api/experiments/{id}", web::delete().to(delete_experiment))
    .route("/api/experiments/{id}", web::get().to(get_experiment))
    .route("/api/experiments/{id}", web::post().to(post_execute_experiment))
    .route("/api/experiments/{id}/abort", web::post().to(post_experiment_execution_abort))
    .route("/api/experiments/{id}/archive", web::post().to(post_experiment_archive_step_results))
    .route("/api/experiments/{id}/comment", web::patch().to(patch_experiment_comment))
    .route("/api/experiments/{id}/download/{archive}", web::get().to(get_experiment_download_step_results))
        // This method is only POST to support the JSON message body.
    .route("/api/experiments/{id}/logs", web::post().to(get_experiment_step_logs))
    .route("/api/experiments/{id}/locked", web::get().to(get_experiment_locked))
    .route("/api/experiments/{id}/mail", web::patch().to(patch_experiment_mail))
    .route("/api/experiments/{id}/name", web::patch().to(patch_experiment_name))
    .route("/api/experiments/{id}/pipeline", web::patch().to(patch_experiment_pipeline))
    .route("/api/experiments/{id}/pipelines", web::get().to(get_experiment_pipelines))
    .route("/api/experiments/{id}/rerun", web::post().to(post_execute_experiment_step))
    .route("/api/experiments/{id}/reset", web::post().to(post_experiment_execution_reset))
    .route("/api/experiments/{id}/run", web::get().to(get_experiment_pipeline_run))
    .route("/api/experiments/{id}/status", web::get().to(get_experiment_execution_status))
    .route("/api/experiments/{id}/variable/global", web::post().to(post_experiment_pipeline_global_variable))
    .route("/api/experiments/{id}/variable/step", web::post().to(post_experiment_pipeline_step_variable))
    // Global data repositories
    .route("/api/globals", web::get().to(list_global_data))
    .route("/api/globals", web::post().to(create_global_data))
    .route("/api/globals/{id}", web::get().to(get_global_data))
    .route("/api/globals/{id}", web::delete().to(delete_global_data))
    .route("/api/globals/{id}/comment", web::patch().to(patch_global_data_comment))
    .route("/api/globals/{id}/name", web::patch().to(patch_global_data_name))
    // Files
    .route("/api/files/{category}/{id}", web::get().to(get_files))
    .route("/api/files/{category}/{id}", web::post().to(post_add_file))
    .route("/api/files/{category}/{id}", web::delete().to(delete_files_by_path))
    .route("/api/folders/{category}/{id}", web::post().to(post_add_folder))
    // Registers static frontend resources. Needs to be last to not overwrite other routes.
    .service(Files::new("/", "./static_dist").show_files_listing());
}
