use crate::{application::config::Configuration, controller::routing::routing_config};
use actix_web::{
    body::MessageBody,
    dev::{ServiceFactory, ServiceRequest, ServiceResponse},
    middleware, App, Error,
};
use diesel::{Connection, SqliteConnection};
use dotenv::dotenv;
use log::{error, warn};
use std::{path::PathBuf, sync::Arc};
use uuid::Uuid;

/**
 * Creates a fully configured app for testing purposes.
 */
pub fn create_test_app(
    database_url: String,
) -> App<
    impl ServiceFactory<
        ServiceRequest,
        Response = ServiceResponse<impl MessageBody>,
        Config = (),
        InitError = (),
        Error = Error,
    >,
> {
    dotenv().unwrap();
    env_logger::try_init_from_env(env_logger::Env::new().filter("debug")).ok();
    App::new()
        .wrap(middleware::Logger::default())
        .app_data(Arc::clone(&Arc::new(
            Configuration::new(
                database_url,
                "info".to_string(),
                "127.0.0.1".to_string(),
                "8080".to_string(),
            ),
        )))
        .configure(routing_config)
}

pub fn create_test_database() -> TestDatabaseContext {
    TestDatabaseContext::new()
}

pub struct TestDatabaseContext {
    id: Uuid,
}

impl TestDatabaseContext {
    pub fn new() -> TestDatabaseContext {
        let id = Configuration::generate_uuid();
        let context = TestDatabaseContext { id };
        diesel_migrations::run_pending_migrations(&context.get_connection()).unwrap();
        context
    }

    pub fn path_to_test_database(&self) -> PathBuf {
        let mut file_path: PathBuf = "../testing_resources/databases/".into();
        file_path.push(self.id.to_string());
        file_path
    }

    pub fn get_connection(&self) -> SqliteConnection {
        let connection = SqliteConnection::establish(
            &self
                .path_to_test_database()
                .into_os_string()
                .into_string()
                .unwrap(),
        )
        .unwrap();
        connection.execute("PRAGMA foreign_keys = ON;").unwrap();
        connection
    }
}

impl Drop for TestDatabaseContext {
    fn drop(&mut self) {
        let file_path = self.path_to_test_database();
        if file_path.exists() {
            if let Err(e) = std::fs::remove_file(file_path) {
                error!("Dropping test database context {} failed with error: {}", self.id, e);
            }
        } else {
            warn!("Tried to delete non existing temporary file {:?}.", file_path)
        }
    }
}
