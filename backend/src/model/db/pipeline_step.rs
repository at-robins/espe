use crate::schema::pipeline_step;
use chrono::NaiveDateTime;
use diesel::{Identifiable, Queryable};

use super::pipeline::Pipeline;

#[derive(Identifiable, Queryable, Associations, PartialEq, Debug)]
#[belongs_to(Pipeline, foreign_key = "pipeline_id")]
#[table_name = "pipeline_step"]
/// A pipeline database entry.
pub struct PipelineStep {
    pub id: i32,
    pub pipeline_id: i32,
    pub execution_type: String,
    pub execution_configuration: String,
    pub ordering: i32,
    pub creation_time: NaiveDateTime,
}
