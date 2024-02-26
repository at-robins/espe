use std::collections::HashMap;

use crate::schema::pipeline_global_variable::{self};
use chrono::{NaiveDateTime, Utc};
use diesel::{
    BoolExpressionMethods, ExpressionMethods, Identifiable, Insertable, OptionalExtension,
    QueryDsl, Queryable, RunQueryDsl, SelectableHelper, SqliteConnection,
};
use getset::Getters;

#[derive(Identifiable, Queryable, Insertable, Selectable, PartialEq, Debug)]
#[diesel(belongs_to(Experiment, foreign_key = experiment_id))]
#[diesel(table_name = pipeline_global_variable)]
/// A queryable global pipeline variable database entry.
pub struct PipelineGlobalVariable {
    pub id: i32,
    pub experiment_id: i32,
    pub pipeline_id: String,
    pub variable_id: String,
    pub variable_value: Option<String>,
    pub creation_time: NaiveDateTime,
}

impl PipelineGlobalVariable {
    /// Returns the variable with the specified IDs if present in the database.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment the variable belongs to
    /// * `pipeline_id` - the ID of the pipeline the variable belongs to
    /// * `variable_id` - the ID of the variable
    /// * `connection` - the database connection
    pub fn get<T: Into<String>, R: Into<String>>(
        experiment_id: i32,
        pipeline_id: T,
        variable_id: R,
        connection: &mut SqliteConnection,
    ) -> Result<Option<PipelineGlobalVariable>, diesel::result::Error> {
        crate::schema::pipeline_global_variable::table
            .filter(
                crate::schema::pipeline_global_variable::experiment_id
                    .eq(experiment_id)
                    .and(
                        crate::schema::pipeline_global_variable::pipeline_id.eq(pipeline_id.into()),
                    )
                    .and(
                        crate::schema::pipeline_global_variable::variable_id.eq(variable_id.into()),
                    ),
            )
            .first(connection)
            .optional()
    }

    /// Returns all variables belonging to the specified experiment and pipeline.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment for which to load variables
    /// * `pipeline_id` - the ID of the pipeline for which to load variables
    /// * `connection` - the database connection
    pub fn get_by_experiment_and_pipeline<T: Into<String>>(
        experiment_id: i32,
        pipeline_id: T,
        connection: &mut SqliteConnection,
    ) -> Result<Vec<PipelineGlobalVariable>, diesel::result::Error> {
        crate::schema::pipeline_global_variable::table
            .filter(
                crate::schema::pipeline_global_variable::experiment_id
                    .eq(experiment_id)
                    .and(
                        crate::schema::pipeline_global_variable::pipeline_id.eq(pipeline_id.into()),
                    ),
            )
            .select(PipelineGlobalVariable::as_select())
            .load(connection)
    }

    /// Returns all variable values belonging to the specified experiment and pipeline.
    /// The keys of the returned map are the variable IDs.
    /// The values of the map are the string representations of the variable values.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment for which to load variables
    /// * `pipeline_id` - the ID of the pipeline for which to load variables
    /// * `connection` - the database connection
    pub fn get_values_by_experiment_and_pipeline<T: Into<String>>(
        experiment_id: i32,
        pipeline_id: T,
        connection: &mut SqliteConnection,
    ) -> Result<HashMap<String, String>, diesel::result::Error> {
        let ps_variables =
            Self::get_by_experiment_and_pipeline(experiment_id, pipeline_id, connection)?;
        let mut variable_map = HashMap::with_capacity(ps_variables.len());
        for ps_variable in ps_variables {
            if let Some(ps_value) = ps_variable.variable_value {
                variable_map.insert(ps_variable.variable_id, ps_value);
            }
        }
        Ok(variable_map)
    }
}

#[derive(Insertable, PartialEq, Debug, Getters)]
#[diesel(table_name = pipeline_global_variable)]
/// A new global pipeline variable database record.
pub struct NewPipelineGlobalVariable {
    #[getset(get = "pub")]
    pub experiment_id: i32,
    #[getset(get = "pub")]
    pub pipeline_id: String,
    #[getset(get = "pub")]
    pub variable_id: String,
    #[getset(get = "pub")]
    pub variable_value: Option<String>,
    #[getset(get = "pub")]
    pub creation_time: NaiveDateTime,
}

impl NewPipelineGlobalVariable {
    /// Creates a new global pipeline variable record for insertion into the database.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment the variable belongs to
    /// * `pipeline_id` - the ID of the pipeline the variable belongs to
    /// * `variable_id` - the id of the variable
    /// * `variable_value` - the value of the variable
    pub fn new<Q: Into<String>, S: Into<String>, T: Into<Option<String>>>(
        experiment_id: i32,
        pipeline_id: Q,
        variable_id: S,
        variable_value: T,
    ) -> Self {
        Self {
            experiment_id,
            pipeline_id: pipeline_id.into(),
            variable_id: variable_id.into(),
            variable_value: variable_value.into(),
            creation_time: Utc::now().naive_utc(),
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::{model::db::experiment::NewExperiment, test_utility::TestContext};

    use super::*;

    #[test]
    fn test_get_by_experiment_and_pipeline() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        // Create dummy experiments.
        let experiment_0 = NewExperiment::new("0".to_string());
        let experiment_1 = NewExperiment::new("1".to_string());
        diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_0)
            .execute(&mut connection)
            .unwrap();
        let experiment_id: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_1)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        // Dummy variable setup.
        let pipeline_id = "Dummy pipeline";
        assert!(PipelineGlobalVariable::get_by_experiment_and_pipeline(
            experiment_id,
            pipeline_id,
            &mut connection
        )
        .unwrap()
        .is_empty());
        let number_of_records = 42;
        let new_records: Vec<PipelineGlobalVariable> = (0..number_of_records)
            .map(|id| PipelineGlobalVariable {
                id,
                experiment_id,
                pipeline_id: pipeline_id.to_string(),
                variable_id: id.to_string(),
                variable_value: Some(id.to_string()),
                creation_time: chrono::Utc::now().naive_local(),
            })
            .collect();
        diesel::insert_into(crate::schema::pipeline_global_variable::table)
            .values(&new_records)
            .execute(&mut connection)
            .unwrap();
        assert_eq!(
            new_records,
            PipelineGlobalVariable::get_by_experiment_and_pipeline(
                experiment_id,
                pipeline_id,
                &mut connection
            )
            .unwrap()
        );
        // After deletion of the experiment the variables should be deleted as well.
        diesel::delete(crate::schema::experiment::table)
            .filter(crate::schema::experiment::id.eq(experiment_id))
            .execute(&mut connection)
            .unwrap();
        assert!(PipelineGlobalVariable::get_by_experiment_and_pipeline(
            experiment_id,
            pipeline_id,
            &mut connection
        )
        .unwrap()
        .is_empty());
    }

    #[test]
    fn test_get_values_by_experiment_and_pipeline() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        // Create dummy experiments.
        let experiment_0 = NewExperiment::new("0".to_string());
        let experiment_1 = NewExperiment::new("1".to_string());
        diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_0)
            .execute(&mut connection)
            .unwrap();
        let experiment_id: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_1)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        // Dummy variable setup.
        let pipeline_id = "Dummy pipeline";
        assert!(PipelineGlobalVariable::get_by_experiment_and_pipeline(
            experiment_id,
            pipeline_id,
            &mut connection
        )
        .unwrap()
        .is_empty());
        let number_of_records = 42;
        let mut expected_map = HashMap::new();
        let new_records: Vec<PipelineGlobalVariable> = (0..number_of_records)
            .map(|id| {
                expected_map.insert(id.to_string(), id.to_string());
                PipelineGlobalVariable {
                    id,
                    experiment_id,
                    pipeline_id: pipeline_id.to_string(),
                    variable_id: id.to_string(),
                    variable_value: Some(id.to_string()),
                    creation_time: chrono::Utc::now().naive_local(),
                }
            })
            .collect();
        diesel::insert_into(crate::schema::pipeline_global_variable::table)
            .values(&new_records)
            .execute(&mut connection)
            .unwrap();
        assert_eq!(
            expected_map,
            PipelineGlobalVariable::get_values_by_experiment_and_pipeline(
                experiment_id,
                pipeline_id,
                &mut connection
            )
            .unwrap()
        );
    }

    #[test]
    fn test_get() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        // Create dummy experiments.
        let experiment = NewExperiment::new("0".to_string());
        let experiment_id: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        // Dummy variable setup.
        let pipeline_variable_row_id: i32 = 42;
        let pipeline_id = "Dummy pipeline";
        let variable_id = "Dummy variable";
        let mut pipeline_variable = PipelineGlobalVariable {
            id: pipeline_variable_row_id,
            experiment_id,
            pipeline_id: pipeline_id.to_string(),
            variable_id: variable_id.to_string(),
            variable_value: Some("Dummy value".to_string()),
            creation_time: chrono::Utc::now().naive_local(),
        };
        assert!(PipelineGlobalVariable::get(
            experiment_id,
            pipeline_id,
            variable_id,
            &mut connection
        )
        .unwrap()
        .is_none());
        // Setting the variable
        diesel::insert_into(crate::schema::pipeline_global_variable::table)
            .values(&pipeline_variable)
            .execute(&mut connection)
            .unwrap();
        assert_eq!(
            &PipelineGlobalVariable::get(
                experiment_id,
                pipeline_id,
                variable_id,
                &mut connection
            )
            .unwrap()
            .unwrap(),
            &pipeline_variable
        );
        // Clearing the variable value.
        pipeline_variable.variable_value = None;
        diesel::update(
            crate::schema::pipeline_global_variable::table
                .filter(crate::schema::pipeline_global_variable::id.eq(pipeline_variable.id)),
        )
        .set(
            crate::schema::pipeline_global_variable::variable_value
                .eq(pipeline_variable.variable_value.clone()),
        )
        .execute(&mut connection)
        .unwrap();
        assert_eq!(
            &PipelineGlobalVariable::get(
                experiment_id,
                pipeline_id,
                variable_id,
                &mut connection
            )
            .unwrap()
            .unwrap(),
            &pipeline_variable
        );
    }
}
