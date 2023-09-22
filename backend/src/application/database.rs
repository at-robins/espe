use actix_web::web;
use diesel::{
    connection::SimpleConnection,
    r2d2::{ConnectionManager, Pool, PooledConnection},
    SqliteConnection,
};

use super::{
    config::Configuration,
    error::{SeqError, SeqErrorType, DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE},
};

#[derive(Debug)]
pub struct CustomConnection {}

impl diesel::r2d2::CustomizeConnection<SqliteConnection, diesel::r2d2::Error> for CustomConnection {
    fn on_acquire(&self, conn: &mut SqliteConnection) -> Result<(), diesel::r2d2::Error> {
        (|| {
            conn.batch_execute("PRAGMA journal_mode = WAL; PRAGMA synchronous = NORMAL; PRAGMA foreign_keys = ON; PRAGMA busy_timeout = 10000;")?;
            Ok(())
        })()
        .map_err(diesel::r2d2::Error::QueryError)
    }
}

pub struct DatabaseManager {
    config: web::Data<Configuration>,
    /// The database connection pool.
    db_connection_pool: Pool<ConnectionManager<SqliteConnection>>,
}

impl DatabaseManager {
    pub fn new(config: web::Data<Configuration>) -> Result<Self, SeqError> {
        let pool = Pool::builder()
            .max_size(16)
            .connection_customizer(Box::new(CustomConnection {}))
            .build(ConnectionManager::<SqliteConnection>::new(config.database_url()))
            .map_err(|err| {
                SeqError::new(
                    "r2d2::Error",
                    SeqErrorType::InternalServerError,
                    err,
                    DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
                )
            })?;
        Ok(Self {
            config,
            db_connection_pool: pool,
        })
    }

    /// Returns a connection to the database if possible.
    pub fn database_connection(
        &self,
    ) -> Result<PooledConnection<ConnectionManager<SqliteConnection>>, SeqError> {
        Ok(self.db_connection_pool.get().map_err(|err| {
            SeqError::new(
                "r2d2::Error",
                SeqErrorType::InternalServerError,
                err,
                DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
            )
        })?)
    }
}
