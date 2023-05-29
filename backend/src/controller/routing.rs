use actix_files::{Files, NamedFile};
use actix_web::web::{self, ServiceConfig};
use log::error;

use crate::application::error::SeqError;

use super::{
    global_data_controller::{
        create_global_data, delete_global_data, get_global_data, get_global_data_files,
        list_global_data, patch_global_data_comment, patch_global_data_name,
        post_global_data_add_file, delete_global_data_files_by_path, post_global_data_add_folder,
    },
    pipeline_controller::{get_pipeline_blueprints, get_pipeline_instance},
    sample_controller::upload_sample,
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

    .route("/api/pipeline/instance/{id}", web::get().to(get_pipeline_instance))
    .route("/api/pipeline/blueprint", web::get().to(get_pipeline_blueprints))
    .route("/api/experiment", web::post().to(upload_sample))
    .route("/api/globals", web::get().to(list_global_data))
    .route("/api/globals", web::post().to(create_global_data))
    .route("/api/globals/{id}", web::get().to(get_global_data))
    .route("/api/globals/{id}", web::delete().to(delete_global_data))
    .route("/api/globals/{id}/comment", web::patch().to(patch_global_data_comment))
    .route("/api/globals/{id}/name", web::patch().to(patch_global_data_name))
    .route("/api/globals/{id}/files", web::get().to(get_global_data_files))
    .route("/api/globals/{id}/files", web::post().to(post_global_data_add_file))
    .route("/api/globals/{id}/files", web::delete().to(delete_global_data_files_by_path))
    .route("/api/globals/{id}/folders", web::post().to(post_global_data_add_folder))

    // Registers static frontend resources. Needs to be last to not overwrite other routes.
    .service(Files::new("/", "./static_dist").show_files_listing());
}
