use chrono::{DateTime, Utc};
use diesel::{Identifiable, Queryable};
use crate::schema::pipeline;

#[derive(Identifiable, Queryable, PartialEq, Debug)]
#[table_name = "pipeline"]
/// A pipeline database entry.
pub struct Pipeline {
    pub id: i64,
    pub pipeline_name: String,
    pub comment: String,
    pub creation_time: DateTime<Utc>,
}
