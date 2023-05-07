use crate::schema::global_data::{self};
use chrono::{NaiveDateTime, Utc};
use diesel::{Identifiable, Insertable, Queryable};
use getset::{CopyGetters, Getters};

#[derive(Identifiable, Queryable, Insertable, PartialEq, Debug)]
#[table_name = "global_data"]
/// A queryable global data database entry.
pub struct GlobalData {
    pub id: i32,
    pub global_data_name: String,
    pub comment: Option<String>,
    pub creation_time: NaiveDateTime,
}

#[derive(Insertable, PartialEq, Debug, Getters, CopyGetters)]
#[table_name = "global_data"]
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
