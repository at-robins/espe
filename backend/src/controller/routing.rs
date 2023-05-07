use actix_files::{Files, NamedFile};
use actix_web::web::{self, ServiceConfig};
use log::error;

use crate::application::error::SeqError;

use super::{pipeline_controller::{get_pipeline_instance, get_pipeline_blueprints}, sample_controller::upload_sample, global_data_controller::{create_global_data, list_global_data, delete_global_data}};

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
    .route("/api/globals/{id}", web::delete().to(delete_global_data))

    // Registers static frontend resources. Needs to be last to not overwrite other routes.
    .service(Files::new("/", "./static_dist").show_files_listing());
}
