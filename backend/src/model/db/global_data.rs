use crate::application::error::{SeqError, SeqErrorType};
use crate::model::exchange::global_data_file_upload::DATABASE_PATH_COMPONENT_SEPARATOR;
use crate::schema::global_data::{self};
use crate::schema::global_data_file::{self};
use chrono::{NaiveDateTime, Utc};
use diesel::{
    BelongingToDsl, BoolExpressionMethods, ExpressionMethods, Identifiable, Insertable, QueryDsl,
    Queryable, RunQueryDsl, SqliteConnection, TextExpressionMethods,
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

    /// Returns `true` if the file path for the [`GlobalData`] with the specified ID exists and `false` otherwise.
    ///
    /// # Parameters
    ///
    /// * `global_data_id` - the entity ID
    /// * `path` - the file path
    /// * `connection` - the database connection
    pub fn exists_path<P: AsRef<str>>(
        global_data_id: i32,
        path: P,
        connection: &mut SqliteConnection,
    ) -> Result<bool, diesel::result::Error> {
        diesel::select(diesel::dsl::exists(
            crate::schema::global_data_file::table.filter(
                crate::schema::global_data_file::global_data_id
                    .eq(global_data_id)
                    .and(crate::schema::global_data_file::file_path.eq(path.as_ref())),
            ),
        ))
        .get_result(connection)
    }

    /// Deletes all files of the [`GlobalData`] with the specifed ID
    /// that start with the specified path prefix.
    ///
    /// # Parameters
    ///
    /// * `global_data_id` - the entity ID
    /// * `path_prefix` - the file path prefix
    /// * `connection` - the database connection
    pub fn delete_files_by_path_prefix<P: AsRef<str>>(
        global_data_id: i32,
        path_prefix: P,
        connection: &mut SqliteConnection,
    ) -> Result<usize, diesel::result::Error> {
        diesel::delete(
            crate::schema::global_data_file::table.filter(
                crate::schema::global_data_file::global_data_id
                    .eq(global_data_id)
                    .and(
                        crate::schema::global_data_file::file_path
                            .like(format!("{}{}%", path_prefix.as_ref(), DATABASE_PATH_COMPONENT_SEPARATOR))
                            .or(crate::schema::global_data_file::file_path
                                .eq(path_prefix.as_ref())),
                    ),
            ),
        )
        .execute(connection)
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

    /// Returns all [`GlobalDataFile`]s belonging to the entity with the specified ID.
    ///
    /// # Parameters
    ///
    /// * `id` - the entity ID
    /// * `connection` - the database connection
    pub fn files(
        id: i32,
        connection: &mut SqliteConnection,
    ) -> Result<Vec<GlobalDataFile>, diesel::result::Error> {
        GlobalDataFile::belonging_to(&Self::get(id, connection)?).load(connection)
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

#[derive(Identifiable, Queryable, Insertable, Associations, PartialEq, Debug)]
#[diesel(belongs_to(GlobalData, foreign_key = global_data_id))]
#[diesel(table_name = global_data_file)]
/// A queryable global data associated file database entry.
pub struct GlobalDataFile {
    pub id: i32,
    pub global_data_id: i32,
    pub file_path: String,
    pub creation_time: NaiveDateTime,
}

#[derive(Insertable, PartialEq, Debug, Getters, CopyGetters)]
#[diesel(table_name = global_data_file)]
/// A new global data associated file database record.
pub struct NewGlobalDataFile {
    #[getset(get = "pub")]
    global_data_id: i32,
    #[getset(get = "pub")]
    file_path: String,
    #[getset(get = "pub")]
    creation_time: NaiveDateTime,
}

impl NewGlobalDataFile {
    /// Creates a new global data file record for insertion into the database.
    ///
    /// # Parameters
    ///
    /// * `global_data_id` - the global data's ID
    /// * `file_path` - the relative file path of where to store the files
    pub fn new<T: Into<String>>(global_data_id: i32, file_path: T) -> Self {
        Self {
            global_data_id,
            file_path: file_path.into(),
            creation_time: Utc::now().naive_utc(),
        }
    }
}
