use std::{borrow::Borrow, collections::HashMap};

use crate::{
    application::error::{SeqError, SeqErrorType},
    model::internal::pipeline_blueprint::PipelineBlueprint,
    schema::pipeline_global_variable,
};
use chrono::{NaiveDateTime, Utc};
use diesel::{
    BoolExpressionMethods, ExpressionMethods, Identifiable, Insertable, OptionalExtension,
    QueryDsl, Queryable, RunQueryDsl, SelectableHelper, SqliteConnection,
};
use getset::Getters;

#[derive(Identifiable, Queryable, Insertable, Selectable, PartialEq, Debug, Clone)]
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

    /// Returns all variable entries with the specified variable ID from the specified pipeline.
    ///
    /// # Parameters
    ///
    /// * `pipeline_id` - the ID of the pipeline the variable belongs to
    /// * `variable_id` - the ID of the variable
    /// * `connection` - the database connection
    pub fn get_by_pipeline_and_variable_id<T: Into<String>, R: Into<String>>(
        pipeline_id: T,
        variable_id: R,
        connection: &mut SqliteConnection,
    ) -> Result<Vec<PipelineGlobalVariable>, diesel::result::Error> {
        crate::schema::pipeline_global_variable::table
            .filter(
                crate::schema::pipeline_global_variable::variable_id
                    .eq(variable_id.into())
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

    /// Returns an error if required pipeline step variables have not been set.
    ///
    /// # Parameters
    ///
    /// * `step` - the step to validate variables for
    /// * `experiment_id` - the ID of the experiment the step belongs to
    /// * `connection` - the database connection
    pub fn validate_global_variables<T: Borrow<PipelineBlueprint>>(
        pipeline: T,
        experiment_id: i32,
        connection: &mut SqliteConnection,
    ) -> Result<(), SeqError> {
        let pipeline: &PipelineBlueprint = pipeline.borrow();
        let experiment_variables =
            Self::get_values_by_experiment_and_pipeline(experiment_id, pipeline.id(), connection)?;
        for variable in pipeline.global_variables() {
            if variable.required().unwrap_or(false) {
                // Error if required variables are not set.
                if !experiment_variables.contains_key(variable.id()) {
                    return Err(SeqError::new(
                                    "Invalid run",
                                    SeqErrorType::BadRequestError,
                                    format!("The experiment {} is missing the required global variable with pipeline id {} and variable id {}.", experiment_id, pipeline.id(), variable.id()),
                                    "The requested run parameters are invalid.",
                                ));
                }
            }
        }
        Ok(())
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
            &PipelineGlobalVariable::get(experiment_id, pipeline_id, variable_id, &mut connection)
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
            &PipelineGlobalVariable::get(experiment_id, pipeline_id, variable_id, &mut connection)
                .unwrap()
                .unwrap(),
            &pipeline_variable
        );
    }

    #[test]
    fn test_get_by_pipeline_and_variable_id() {
        // Use a reference to the context, so the context is not dropped early
        // and messes up test context folder deletion.
        let context = TestContext::new();
        let mut connection = context.get_connection();
        // Create dummy experiments.
        let experiment_0 = NewExperiment::new("0".to_string());
        let experiment_1 = NewExperiment::new("1".to_string());
        let experiment_id_0: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_0)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        let experiment_id_1: i32 = diesel::insert_into(crate::schema::experiment::table)
            .values(&experiment_1)
            .returning(crate::schema::experiment::id)
            .get_result(&mut connection)
            .unwrap();
        // Dummy variable setup.
        let pipeline_id = "Dummy pipeline";
        assert!(PipelineGlobalVariable::get_by_experiment_and_pipeline(
            experiment_id_1,
            pipeline_id,
            &mut connection
        )
        .unwrap()
        .is_empty());
        let number_of_records = 42;
        let new_records_exp_0: Vec<PipelineGlobalVariable> = (0..number_of_records)
            .map(|id| PipelineGlobalVariable {
                id,
                experiment_id: experiment_id_0,
                pipeline_id: pipeline_id.to_string(),
                variable_id: id.to_string(),
                variable_value: Some(id.to_string()),
                creation_time: chrono::Utc::now().naive_local(),
            })
            .collect();
        let new_records_exp_1: Vec<PipelineGlobalVariable> = new_records_exp_0
            .iter()
            .map(|glob_var_0| {
                let mut glob_var_1 = glob_var_0.clone();
                glob_var_1.experiment_id = experiment_id_1;
                glob_var_1.id = glob_var_1.id + number_of_records;
                glob_var_1
            })
            .collect();
        diesel::insert_into(crate::schema::pipeline_global_variable::table)
            .values(&new_records_exp_0)
            .execute(&mut connection)
            .unwrap();
        diesel::insert_into(crate::schema::pipeline_global_variable::table)
            .values(&new_records_exp_1)
            .execute(&mut connection)
            .unwrap();
        let query_variable_id = 3;
        assert_eq!(
            vec![
                new_records_exp_0[query_variable_id].clone(),
                new_records_exp_1[query_variable_id].clone()
            ],
            PipelineGlobalVariable::get_by_pipeline_and_variable_id(
                pipeline_id,
                query_variable_id.to_string(),
                &mut connection
            )
            .unwrap()
        );
    }
}
