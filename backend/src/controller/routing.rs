use actix_files::{Files, NamedFile};
use actix_web::web::{self, ServiceConfig};

async fn index() -> actix_web::Result<NamedFile> {
    Ok(NamedFile::open("./static_dist/index.html")?)
}

pub fn routing_config(cfg: &mut ServiceConfig) {
    cfg

    // UI redirect
    .route("/", web::get().to(index))
    .route("/ui", web::get().to(index))
    .route("/ui/{rest:.*}", web::get().to(index))

    // Registers static frontend resources. Needs to be last to not overwrite other routes.
    .service(Files::new("/", "./static_dist").show_files_listing());
}
