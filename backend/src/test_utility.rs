use crate::{application::config::Configuration, controller::routing::routing_config};
use actix_web::{
    body::MessageBody,
    dev::{ServiceFactory, ServiceRequest, ServiceResponse},
    middleware, App, Error,
};
use diesel::{Connection, SqliteConnection};
use dotenv::dotenv;
use log::{error, warn};
use std::{path::{PathBuf}, sync::Arc};
use uuid::Uuid;

/**
 * Creates a fully configured app for testing purposes
 * for testing against the specified test database.
 *
 * # Parameters
 *
 * * `database_url` - the URL / URI of the test database
 */
pub fn create_test_app<T: Into<String>>(
    database_url: T,
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
        .app_data(Arc::clone(&Arc::new(Configuration::new(
            database_url,
            "info",
            "127.0.0.1",
            "8080",
        ))))
        .configure(routing_config)
}

/**
 * A test context that provides clean test resources (e.g. test databases) and according
 * initialisation and cleanup on a per test basis.
 */
pub struct TestContext {
    id: Uuid,
}

impl TestContext {
    /**
     * Creates a new `TestContext`.
     */
    pub fn new() -> TestContext {
        let id = Configuration::generate_uuid();
        let context = TestContext { id };
        diesel_migrations::run_pending_migrations(&context.get_connection()).unwrap();
        context
    }

    /**
     * Returns the URL / URI of the test database.
     */
    pub fn database_url(&self) -> String {
        format!("../testing_resources/databases/{}", self.id)
    }

    /**
     * Opens a connection to the test database.
     */
    pub fn get_connection(&self) -> SqliteConnection {
        let connection = SqliteConnection::establish(&self.database_url()).unwrap();
        connection.execute("PRAGMA foreign_keys = ON;").unwrap();
        connection
    }
}

impl Drop for TestContext {
    fn drop(&mut self) {
        let file_path: PathBuf = self.database_url().into();
        if file_path.exists() {
            if let Err(e) = std::fs::remove_file(file_path) {
                error!("Dropping test database context {} failed with error: {}", self.id, e);
            }
        } else {
            warn!("Tried to delete non existing temporary file {:?}.", file_path)
        }
    }
}
