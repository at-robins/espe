use crate::{
    application::error::{SeqError, SeqErrorType},
    schema::experiment::{self},
};
use chrono::{NaiveDateTime, Utc};
use diesel::{
    ExpressionMethods, Identifiable, Insertable, QueryDsl, Queryable, RunQueryDsl, SqliteConnection,
};
use getset::Getters;

#[derive(Identifiable, Queryable, Insertable, PartialEq, Debug)]
#[diesel(table_name = experiment)]
/// A queryable experiment database entry.
pub struct Experiment {
    pub id: i32,
    pub experiment_name: String,
    pub mail: Option<String>,
    pub pipeline_id: Option<String>,
    pub comment: Option<String>,
    pub creation_time: NaiveDateTime,
}

impl Experiment {
    /// Returns `true` if the entity with the specified ID exists and `false` otherwise.
    ///
    /// # Parameters
    ///
    /// * `id` - the entity ID
    /// * `connection` - the database connection
    pub fn exists(
        id: i32,
        connection: &mut SqliteConnection,
    ) -> Result<bool, diesel::result::Error> {
        diesel::select(diesel::dsl::exists(
            crate::schema::experiment::table.filter(crate::schema::experiment::id.eq(id)),
        ))
        .get_result(connection)
    }

    /// Returns [`Ok`] if the entity with the specified ID exists
    /// and a `NotFound` [`Err`] if not present.
    ///
    /// # Parameters
    ///
    /// * `id` - the entity ID
    /// * `connection` - the database connection
    pub fn exists_err(id: i32, connection: &mut SqliteConnection) -> Result<(), SeqError> {
        if Self::exists(id, connection)? {
            Ok(())
        } else {
            Err(SeqError::new(
                "Invalid request",
                SeqErrorType::NotFoundError,
                format!("Experiment with ID {} does not exist.", id),
                "The entity does not exist.",
            ))
        }
    }

    /// Returns the entity with the specified ID.
    ///
    /// # Parameters
    ///
    /// * `id` - the entity ID
    /// * `connection` - the database connection
    pub fn get(
        id: i32,
        connection: &mut SqliteConnection,
    ) -> Result<Experiment, diesel::result::Error> {
        crate::schema::experiment::table
            .find(id)
            .first::<Experiment>(connection)
    }

    /// Returns all entities.
    ///
    /// # Parameters
    ///
    /// * `connection` - the database connection
    pub fn get_all(
        connection: &mut SqliteConnection,
    ) -> Result<Vec<Experiment>, diesel::result::Error> {
        crate::schema::experiment::table.load::<Experiment>(connection)
    }
}

#[derive(Insertable, PartialEq, Debug, Getters)]
#[diesel(table_name = experiment)]
/// A new experiment database record.
pub struct NewExperiment {
    #[getset(get = "pub")]
    experiment_name: String,
    #[getset(get = "pub")]
    mail: Option<String>,
    #[getset(get = "pub")]
    pipeline_id: Option<String>,
    #[getset(get = "pub")]
    comment: Option<String>,
    #[getset(get = "pub")]
    creation_time: NaiveDateTime,
}

impl NewExperiment {
    /// Creates a new experiment record for insertion into the database.
    ///
    /// # Parameters
    ///
    /// * `name` - the experiment's name
    pub fn new(name: String) -> Self {
        Self {
            experiment_name: name,
            mail: None,
            pipeline_id: None,
            comment: None,
            creation_time: Utc::now().naive_utc(),
        }
    }
}
