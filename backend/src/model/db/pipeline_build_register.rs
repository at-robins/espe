use crate::schema::pipeline_build_register::{self};
use chrono::{NaiveDateTime, Utc};
use diesel::{
    BoolExpressionMethods, ExpressionMethods, Identifiable, Insertable, QueryDsl, Queryable,
    RunQueryDsl, SqliteConnection,
};
use getset::Getters;

#[derive(Identifiable, Queryable, Insertable, PartialEq, Debug)]
#[diesel(table_name = pipeline_build_register)]
/// A queryable pipeline / container build register database entry.
pub struct PipelineBuildRegister {
    pub id: i32,
    pub pipeline_id: String,
    pub pipeline_step_id: String,
    pub pipeline_version: String,
    pub creation_time: NaiveDateTime,
}

impl PipelineBuildRegister {
    /// Returns `true` if the pipeline step has already been built with the specific version.
    ///
    /// # Parameters
    ///
    /// * `pipeline_id` - the ID of the pipeline
    /// * `pipeline_step_id` - the ID of the pipeline step
    /// * `pipeline_version` - the version of the pipeline
    /// * `connection` - the database connection
    pub fn is_built<
        PipelineIdIdType: Into<String>,
        PipelineStepIdIdType: Into<String>,
        PipelineVersionType: Into<String>,
    >(
        pipeline_id: PipelineIdIdType,
        pipeline_step_id: PipelineStepIdIdType,
        pipeline_version: PipelineVersionType,
        connection: &mut SqliteConnection,
    ) -> Result<bool, diesel::result::Error> {
        diesel::select(diesel::dsl::exists(
            crate::schema::pipeline_build_register::table.filter(
                crate::schema::pipeline_build_register::pipeline_id
                    .eq(pipeline_id.into())
                    .and(
                        crate::schema::pipeline_build_register::pipeline_step_id
                            .eq(pipeline_step_id.into()),
                    )
                    .and(
                        crate::schema::pipeline_build_register::pipeline_version
                            .eq(pipeline_version.into()),
                    ),
            ),
        ))
        .get_result(connection)
    }

    /// Inserts the specifed build record into the database.
    ///
    /// # Parameters
    ///
    /// * `pipeline_id` - the ID of the pipeline
    /// * `pipeline_step_id` - the ID of the pipeline step
    /// * `pipeline_version` - the version of the pipeline
    /// * `connection` - the database connection
    pub fn set_built<
        PipelineIdIdType: Into<String>,
        PipelineStepIdIdType: Into<String>,
        PipelineVersionType: Into<String>,
    >(
        pipeline_id: PipelineIdIdType,
        pipeline_step_id: PipelineStepIdIdType,
        pipeline_version: PipelineVersionType,
        connection: &mut SqliteConnection,
    ) -> Result<(), diesel::result::Error> {
        let build_record =
            NewPipelineBuildRegister::new(pipeline_id, pipeline_step_id, pipeline_version);
        diesel::insert_into(crate::schema::pipeline_build_register::table)
            .values(&build_record)
            .on_conflict((
                crate::schema::pipeline_build_register::pipeline_id,
                crate::schema::pipeline_build_register::pipeline_step_id,
            ))
            .do_update()
            .set(&build_record)
            .execute(connection)?;
        Ok(())
    }
}

#[derive(Insertable, PartialEq, Debug, Getters, AsChangeset)]
#[diesel(table_name = pipeline_build_register)]
/// A new pipeline / container build register database record.
pub struct NewPipelineBuildRegister {
    pipeline_id: String,
    pipeline_step_id: String,
    pipeline_version: String,
    creation_time: NaiveDateTime,
}

impl NewPipelineBuildRegister {
    /// Creates a new pipeline build register record for insertion into the database.
    ///
    /// # Parameters
    ///
    /// * `pipeline_id` - the ID of the pipeline
    /// * `pipeline_step_id` - the ID of the pipeline step
    /// * `pipeline_version` - the version of the pipeline
    pub fn new<
        PipelineIdIdType: Into<String>,
        PipelineStepIdIdType: Into<String>,
        PipelineVersionType: Into<String>,
    >(
        pipeline_id: PipelineIdIdType,
        pipeline_step_id: PipelineStepIdIdType,
        pipeline_version: PipelineVersionType,
    ) -> Self {
        Self {
            pipeline_id: pipeline_id.into(),
            pipeline_step_id: pipeline_step_id.into(),
            pipeline_version: pipeline_version.into(),
            creation_time: Utc::now().naive_utc(),
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::test_utility::TestContext;

    use super::*;

    #[test]
    fn test_is_built() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        let pipeline_id = "test id";
        let pipeline_step_id = "test step id";
        let pipeline_version = "1.0.0";
        let new_record =
            NewPipelineBuildRegister::new(pipeline_id, pipeline_step_id, pipeline_version);
        diesel::insert_into(crate::schema::pipeline_build_register::table)
            .values(&new_record)
            .execute(&mut connection)
            .unwrap();
        assert!(PipelineBuildRegister::is_built(
            pipeline_id,
            pipeline_step_id,
            pipeline_version,
            &mut connection
        )
        .unwrap());
        assert!(!PipelineBuildRegister::is_built(
            pipeline_id.to_owned() + "a",
            pipeline_step_id,
            pipeline_version,
            &mut connection
        )
        .unwrap());
        assert!(!PipelineBuildRegister::is_built(
            pipeline_id,
            pipeline_step_id.to_owned() + "a",
            pipeline_version,
            &mut connection
        )
        .unwrap());
        assert!(!PipelineBuildRegister::is_built(
            pipeline_id,
            pipeline_step_id,
            pipeline_version.to_owned() + "a",
            &mut connection
        )
        .unwrap());
    }

    #[test]
    fn test_set_built() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        let pipeline_id = "test id";
        let pipeline_step_id = "test step id";
        let pipeline_version = "1.0.0";
        let pipeline_version_updated = "1.0.1";
        // Insert a dummy record.
        let dummy_record = PipelineBuildRegister {
            id: 0,
            pipeline_id: "dummy id".to_string(),
            pipeline_step_id: "dummy step id".to_string(),
            pipeline_version: "0.4.2".to_string(),
            creation_time: chrono::Utc::now().naive_local(),
        };
        diesel::insert_into(crate::schema::pipeline_build_register::table)
            .values(&dummy_record)
            .execute(&mut connection)
            .unwrap();
        // Inserts build record.
        PipelineBuildRegister::set_built(
            pipeline_id,
            pipeline_step_id,
            pipeline_version,
            &mut connection,
        )
        .unwrap();
        let all_records = crate::schema::pipeline_build_register::table
            .load::<PipelineBuildRegister>(&mut connection)
            .unwrap();
        assert_eq!(all_records.len(), 2);
        assert_eq!(all_records[0], dummy_record);
        assert!(PipelineBuildRegister::is_built(
            pipeline_id,
            pipeline_step_id,
            pipeline_version,
            &mut connection
        )
        .unwrap());
        assert!(!PipelineBuildRegister::is_built(
            pipeline_id,
            pipeline_step_id,
            pipeline_version_updated,
            &mut connection
        )
        .unwrap());
        // Updates the inserted build record.
        PipelineBuildRegister::set_built(
            pipeline_id,
            pipeline_step_id,
            pipeline_version_updated,
            &mut connection,
        )
        .unwrap();
        let all_records_updated = crate::schema::pipeline_build_register::table
            .load::<PipelineBuildRegister>(&mut connection)
            .unwrap();
        assert_eq!(all_records_updated.len(), 2);
        assert_eq!(all_records_updated[0], dummy_record);
        assert!(!PipelineBuildRegister::is_built(
            pipeline_id,
            pipeline_step_id,
            pipeline_version,
            &mut connection
        )
        .unwrap());
        assert!(PipelineBuildRegister::is_built(
            pipeline_id,
            pipeline_step_id,
            pipeline_version_updated,
            &mut connection
        )
        .unwrap());
    }
}
