use actix_files::{Files, NamedFile};
use actix_web::web::{self, ServiceConfig};
use log::error;

use crate::application::error::SeqError;

use super::pipeline_controller::get_pipeline;

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

    .route("/api/pipeline/{id}", web::get().to(get_pipeline))

    // Registers static frontend resources. Needs to be last to not overwrite other routes.
    .service(Files::new("/", "./static_dist").show_files_listing());
}
