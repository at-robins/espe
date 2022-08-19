use crate::schema::pipeline;
use chrono::{NaiveDateTime, Utc};
use diesel::{Identifiable, Queryable};
use getset::Getters;

#[derive(Identifiable, Queryable, PartialEq, Debug)]
#[table_name = "pipeline"]
/// A pipeline database entry.
pub struct Pipeline {
    pub id: i32,
    pub pipeline_name: String,
    pub comment: String,
    pub creation_time: NaiveDateTime,
}

#[derive(Insertable, PartialEq, Debug, Getters)]
#[table_name = "pipeline"]
/// A new experiment database record.
pub struct NewPipeline {
    #[getset(get = "pub")]
    pipeline_name: String,
    #[getset(get = "pub")]
    comment: String,
    #[getset(get = "pub")]
    creation_time: NaiveDateTime,
}

impl NewPipeline {
    /// Creates a new pipeline record for insertion into the database.
    ///
    /// # Parameters
    ///
    /// * `name` - the pipeline's name
    /// * `comment` - the comment describing the pipeline
    pub fn new(name: String, comment: String) -> Self {
        Self {
            pipeline_name: name,
            comment,
            creation_time: Utc::now().naive_utc(),
        }
    }
}