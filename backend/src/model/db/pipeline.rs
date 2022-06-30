use crate::schema::pipeline;
use chrono::NaiveDateTime;
use diesel::{Identifiable, Queryable};

#[derive(Identifiable, Queryable, PartialEq, Debug)]
#[table_name = "pipeline"]
/// A pipeline database entry.
pub struct Pipeline {
    pub id: i32,
    pub pipeline_name: String,
    pub comment: String,
    pub creation_time: NaiveDateTime,
}
