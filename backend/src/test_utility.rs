use crate::{
    application::{
        config::{ApplicationMode, Configuration},
        database::DatabaseManager,
        environment::LOG_LEVEL,
    },
    controller::routing::routing_config,
    model::{
        db::{
            experiment::Experiment,
            experiment_execution::{ExecutionStatus, NewExperimentExecution},
        },
        internal::archive::ArchiveMetadata,
    },
    service::{execution_service::ExecutionScheduler, pipeline_service::LoadedPipelines},
};
use actix_web::{
    body::MessageBody,
    dev::{ServiceFactory, ServiceRequest, ServiceResponse},
    middleware, web, App, Error,
};
use diesel::{connection::SimpleConnection, Connection, RunQueryDsl, SqliteConnection};
use diesel_migrations::MigrationHarness;
use dotenv::dotenv;
use parking_lot::Mutex;
use std::path::PathBuf;
use uuid::Uuid;

/// The test resource path.
pub const TEST_RESOURCES_PATH: &str = "../testing_resources";

/// Creates a fully configured app for testing purposes
/// for testing against the specified test database.
///
/// # Parameters
///
/// * `database_url` - the URL / URI of the test database
pub fn create_test_app(
    context: &TestContext,
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
    env_logger::try_init_from_env(env_logger::Env::new().filter(LOG_LEVEL)).ok();
    let app_config = &web::Data::<Configuration>::new(context.into());
    let database_manager =
        web::Data::new(DatabaseManager::new(web::Data::clone(&app_config)).unwrap());
    let loaded_pipelines =
        web::Data::new(LoadedPipelines::new(web::Data::clone(&app_config)).unwrap());
    App::new()
        .wrap(middleware::Logger::default())
        .app_data(web::Data::clone(&app_config))
        .app_data(web::Data::clone(&loaded_pipelines))
        .app_data(web::Data::clone(&database_manager))
        .app_data(web::Data::new(Mutex::new(ExecutionScheduler::new(
            web::Data::clone(&app_config),
            web::Data::clone(&database_manager),
            web::Data::clone(&loaded_pipelines),
        ))))
        .configure(routing_config)
}

/// A test context that provides clean test resources (e.g. test databases) and according
/// initialisation and cleanup on a per test basis.
pub struct TestContext {
    id: Uuid,
    pipeline_folder_override: Option<String>,
}

impl TestContext {
    /// Creates a new `TestContext`.
    pub fn new() -> TestContext {
        let id = Uuid::new_v4();
        let context = TestContext {
            id,
            pipeline_folder_override: None,
        };
        std::fs::create_dir_all(context.context_folder()).unwrap();
        std::fs::create_dir_all(context.pipeline_folder()).unwrap();
        let con = &mut context.get_connection();
        con.run_pending_migrations(
            diesel_migrations::FileBasedMigrations::find_migrations_directory().unwrap(),
        )
        .unwrap();
        context
    }

    /// Returns the URL / URI of the test database.
    pub fn database_url(&self) -> String {
        format!("{}/{}.db", self.context_folder(), self.id)
    }

    /// Returns the context folder that stores all context information.
    pub fn context_folder(&self) -> String {
        format!("{}/{}", TEST_RESOURCES_PATH, self.id)
    }

    /// Returns the pipeline definition folder.
    pub fn pipeline_folder(&self) -> String {
        if let Some(custom_folder) = &self.pipeline_folder_override {
            custom_folder.clone()
        } else {
            format!("{}/{}_pipelines", self.context_folder(), self.id)
        }
    }

    /// Overrides the default testing pipeline folder with a custom one.
    ///
    /// # Parameters
    ///
    /// * `pipeline_folder` - the custom testing pipeline folder
    pub fn set_pipeline_folder<T: Into<String>>(&mut self, pipeline_folder: T) {
        self.pipeline_folder_override = Some(pipeline_folder.into());
    }

    /// Opens a connection to the test database.
    pub fn get_connection(&self) -> SqliteConnection {
        let mut connection = SqliteConnection::establish(&self.database_url()).unwrap();
        connection
            .batch_execute(
                "PRAGMA foreign_keys = ON;
            PRAGMA journal_mode = WAL;
            PRAGMA synchronous = NORMAL;
            PRAGMA busy_timeout = 10000;",
            )
            .unwrap();
        connection
    }
}

impl Drop for TestContext {
    fn drop(&mut self) {
        let context_path: PathBuf = self.context_folder().into();
        if context_path.exists() {
            if let Err(e) = std::fs::remove_dir_all(&context_path) {
                log::error!("Dropping test context {} failed with error: {}", self.id, e);
            } else {
                log::info!("Removed test context {}.", context_path.display());
            }
        } else {
            log::warn!(
                "Tried to delete non existing temporary context folder {}.",
                context_path.display()
            );
        }
    }
}

impl From<&TestContext> for Configuration {
    fn from(context: &TestContext) -> Self {
        Configuration::new(
            context.database_url(),
            "info",
            "127.0.0.1",
            "8080",
            context.context_folder(),
            context.pipeline_folder(),
            ApplicationMode::Release,
        )
    }
}

pub const DEFAULT_EXPERIMENT_ID: i32 = 42;
pub const DEFAULT_PIPELINE_ID: &str = "testing_pipeline";
pub const DEFAULT_PIPELINE_STEP_ID: &str = "fastqc";
pub const DEFAULT_ARCHIVE_ID: &str = "42";
// let global_variable_id = "global_number";
// let step_variable_id = "number";
// let new_variable_value = "42";

/// Creates a default dummy experiment for testing.
///
/// # Parameters
///
/// * `connection` - a connection to the test database
pub fn create_default_experiment(connection: &mut SqliteConnection) {
    let new_record = Experiment {
        id: DEFAULT_EXPERIMENT_ID,
        experiment_name: "Dummy record".to_string(),
        comment: Some("A comment".to_string()),
        mail: Some("a.b@c.de".to_string()),
        pipeline_id: Some(DEFAULT_PIPELINE_ID.to_string()),
        creation_time: chrono::Utc::now().naive_local(),
    };
    diesel::insert_into(crate::schema::experiment::table)
        .values(&new_record)
        .execute(connection)
        .unwrap();
}

/// Creates a default dummy archive for testing.
///
/// # Parameters
///
/// * `connection` - a connection to the test database
pub fn create_default_temporary_download_archive(context: &TestContext) {
    let test_config = Configuration::from(context);
    std::fs::create_dir_all(test_config.temporary_download_path()).unwrap();
    let archive_path = test_config.temporary_download_file_path(DEFAULT_ARCHIVE_ID);
    std::fs::File::create_new(&archive_path).unwrap();
    let archive_metadata = ArchiveMetadata::new(format!("{}.zip", DEFAULT_ARCHIVE_ID));
    let archive_metadata_path = ArchiveMetadata::metadata_path(&archive_path);
    serde_json::to_writer(
        std::fs::File::create_new(archive_metadata_path).unwrap(),
        &archive_metadata,
    )
    .unwrap();
}

/// Creates an execution step with the specified status for the default dummy experiment for testing.
///
/// # Parameters
///
/// * `connection` - a connection to the test database
/// * `status` - the [`ExecutionStatus`] of the step execution
pub fn create_default_experiment_execution(
    connection: &mut SqliteConnection,
    status: ExecutionStatus,
) {
    let new_execution = [NewExperimentExecution::new_with_status(
        DEFAULT_EXPERIMENT_ID,
        DEFAULT_PIPELINE_ID,
        DEFAULT_PIPELINE_STEP_ID,
        status,
    )];
    diesel::insert_into(crate::schema::experiment_execution::table)
        .values(&new_execution)
        .execute(connection)
        .unwrap();
}
