use crate::application::error::{SeqError, SeqErrorType};
use crate::schema::global_data::{self};
use chrono::{NaiveDateTime, Utc};
use diesel::{
    ExpressionMethods, Identifiable, Insertable, QueryDsl, Queryable,
    RunQueryDsl, SqliteConnection,
};
use getset::{CopyGetters, Getters};

#[derive(Identifiable, Queryable, Insertable, PartialEq, Debug)]
#[diesel(table_name = global_data)]
/// A queryable global data database entry.
pub struct GlobalData {
    pub id: i32,
    pub global_data_name: String,
    pub comment: Option<String>,
    pub creation_time: NaiveDateTime,
}

impl GlobalData {
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
            crate::schema::global_data::table.filter(crate::schema::global_data::id.eq(id)),
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
                format!("Global data with ID {} does not exist.", id),
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
    ) -> Result<GlobalData, diesel::result::Error> {
        crate::schema::global_data::table
            .find(id)
            .first::<GlobalData>(connection)
    }

    /// Returns all entities.
    ///
    /// # Parameters
    ///
    /// * `connection` - the database connection
    pub fn get_all(
        connection: &mut SqliteConnection,
    ) -> Result<Vec<GlobalData>, diesel::result::Error> {
        crate::schema::global_data::table.load::<GlobalData>(connection)
    }
}

#[derive(Insertable, PartialEq, Debug, Getters, CopyGetters)]
#[diesel(table_name = global_data)]
/// A new global data database record.
pub struct NewGlobalData {
    #[getset(get = "pub")]
    global_data_name: String,
    #[getset(get = "pub")]
    comment: Option<String>,
    #[getset(get = "pub")]
    creation_time: NaiveDateTime,
}

impl NewGlobalData {
    /// Creates a new global data record for insertion into the database.
    ///
    /// # Parameters
    ///
    /// * `name` - the global data's name
    /// * `comment` - the optional comment describing the global data repository
    pub fn new<S: Into<String>, T: Into<Option<String>>>(name: S, comment: T) -> Self {
        Self {
            global_data_name: name.into(),
            comment: comment.into(),
            creation_time: Utc::now().naive_utc(),
        }
    }
}
