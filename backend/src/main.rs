use std::sync::Arc;

use actix_web::{middleware, App, HttpServer};
use application::{config::Configuration, error::SeqError};
use controller::routing::routing_config;

#[actix_web::main]
async fn main() -> Result<(), SeqError> {
    let app_config = Arc::new(init_config());
    let app_config_internal = Arc::clone(&app_config);
    Ok(HttpServer::new(move || {
        App::new()
            .wrap(middleware::Logger::default())
            .app_data(Arc::clone(&app_config_internal))
            .configure(routing_config)
    })
    .bind(app_config.server_address_and_port())?
    .run()
    .await?)
}

fn init_config() -> Configuration {
    let config = Configuration::load_from_file()
        .map_err(|config_loading_error| {
            println!("Configuration could not be loaded using default: {:?}", config_loading_error);
            config_loading_error
        })
        .unwrap_or_default();
    if !Configuration::exists() {
        if let Err(config_save_error) = config.save_to_file() {
            println!("Configuration could not be saved: {:?}", config_save_error);
        }
    }
    config
}

mod application;
mod controller;
