use crate::schema::global_data::{self};
use crate::schema::global_data_file::{self};
use chrono::{NaiveDateTime, Utc};
use diesel::{Identifiable, Insertable, Queryable};
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
