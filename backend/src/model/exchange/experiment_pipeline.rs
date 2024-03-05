use std::{borrow::Borrow, collections::HashMap};

use getset::Getters;
use serde::{Deserialize, Serialize};

use crate::model::internal::pipeline_blueprint::{
    PipelineBlueprint, PipelineStepBlueprint, PipelineStepVariable, PipelineStepVariableCategory,
};

/// The definition of a pipeline.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct ExperimentPipelineBlueprint {
    /// The unique ID of the pipeline.
    #[getset(get = "pub")]
    id: String,
    /// The name for display.
    #[getset(get = "pub")]
    name: String,
    /// The version of the pipeline.
    #[getset(get = "pub")]
    version: String,
    /// A description of the pipeline.
    #[getset(get = "pub")]
    description: String,
    /// The global variables that can be specified for the pipeline.
    #[getset(get = "pub")]
    global_variables: Vec<ExperimentPipelineStepVariable>,
    /// The [`PipelineStepBlueprint`] that make up the pipeline.
    #[getset(get = "pub")]
    steps: Vec<ExperimentPipelineStepBlueprint>,
}

impl ExperimentPipelineBlueprint {
    /// Creates an experiment related pipeline with variable values from
    /// a [`PipelineBlueprint`].
    ///
    /// # Parameters
    ///
    /// * `pipeline` - the pipeline to convert
    /// * `values_global` - a map of global variable values, where the keys are the variable ID
    /// * `values_step` - a map of step variable values, , where the keys are a concatenation of the
    /// pipeline step ID and variable ID
    /// * `stati` - a map of stati, where the keys are the pipeline step ID
    pub fn from_internal<
        PipelineType: Borrow<PipelineBlueprint>,
        ValueGlobalMapType: Borrow<HashMap<String, String>>,
        ValueStepMapType: Borrow<HashMap<String, String>>,
        StatusMapType: Borrow<HashMap<String, String>>,
    >(
        pipeline: PipelineType,
        values_global: ValueGlobalMapType,
        values_step: ValueStepMapType,
        stati: StatusMapType,
    ) -> Self {
        let global_variables = pipeline
            .borrow()
            .global_variables()
            .iter()
            .map(|v| {
                ExperimentPipelineStepVariable::from_internal(
                    v,
                    values_global.borrow().get(v.id()).map(|s| s.clone()),
                )
            })
            .collect();
        let steps = pipeline
            .borrow()
            .steps()
            .iter()
            .map(|s| {
                ExperimentPipelineStepBlueprint::from_internal(
                    s,
                    values_step.borrow(),
                    stati.borrow(),
                )
            })
            .collect();
        Self {
            id: pipeline.borrow().id().clone(),
            name: pipeline.borrow().name().clone(),
            version: pipeline.borrow().version().clone(),
            description: pipeline.borrow().description().clone(),
            global_variables,
            steps,
        }
    }
}

/// The definition of a pipeline step.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct ExperimentPipelineStepBlueprint {
    /// The unique ID of the pipeline step.
    #[getset(get = "pub")]
    id: String,
    /// The name for display.
    #[getset(get = "pub")]
    name: String,
    /// A description of the pipeline step.
    #[getset(get = "pub")]
    description: String,
    /// The container used to execute the pipeline step.
    #[getset(get = "pub")]
    container: String,
    /// The IDs of pipeline steps this step depends on.
    #[getset(get = "pub")]
    dependencies: Vec<String>,
    /// The variables that can be specified for the pipeline step.
    #[getset(get = "pub")]
    variables: Vec<ExperimentPipelineStepVariable>,
    /// The execution status of the pipeline step.
    #[getset(get = "pub")]
    status: Option<String>,
}

impl ExperimentPipelineStepBlueprint {
    /// Creates an experiment related pipeline step with variable values from
    /// a [`PipelineStepBlueprint`].
    ///
    /// # Parameters
    ///
    /// * `pipeline_step` - the step to convert
    /// * `values` - a map of variable values, where the keys are a concatenation of the
    /// pipeline step ID and variable ID
    /// * `stati` - a map of stati, where the keys are the pipeline step ID
    pub fn from_internal<
        StepType: Borrow<PipelineStepBlueprint>,
        ValueMapType: Borrow<HashMap<String, String>>,
        StatusMapType: Borrow<HashMap<String, String>>,
    >(
        pipeline_step: StepType,
        values: ValueMapType,
        stati: StatusMapType,
    ) -> Self {
        let variables = pipeline_step
            .borrow()
            .variables()
            .iter()
            .map(|v| {
                ExperimentPipelineStepVariable::from_internal(
                    v,
                    values
                        .borrow()
                        .get(&format!("{}{}", pipeline_step.borrow().id(), v.id()))
                        .map(|s| s.clone()),
                )
            })
            .collect();
        let status: Option<String> = stati
            .borrow()
            .get(pipeline_step.borrow().id())
            .map(|s| s.clone());
        Self {
            id: pipeline_step.borrow().id().clone(),
            name: pipeline_step.borrow().name().clone(),
            description: pipeline_step.borrow().description().clone(),
            container: pipeline_step.borrow().container().clone(),
            dependencies: pipeline_step.borrow().dependencies().clone(),
            variables,
            status,
        }
    }
}

/// The definition of a pipeline step variable associated with
/// an experiment containing an optional value.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct ExperimentPipelineStepVariable {
    /// The unique ID of the variable.
    #[getset(get = "pub")]
    id: String,
    /// The name for display.
    #[getset(get = "pub")]
    name: String,
    /// A description of the variable.
    #[getset(get = "pub")]
    description: String,
    /// The type of variable.
    #[getset(get = "pub")]
    category: PipelineStepVariableCategory,
    /// If the variable must be set by its instance.
    #[getset(get = "pub")]
    required: Option<bool>,
    /// The value of the variable.
    #[getset(get = "pub")]
    value: Option<String>,
}

impl ExperimentPipelineStepVariable {
    /// Creates an experiment related pipeline step variable with a value from
    /// a [`PipelineStepVariable`].
    ///
    /// # Parameters
    ///
    /// * `pipeline_step_variable` - the variable to convert
    /// * `value` - the value of the variable
    pub fn from_internal<T: Borrow<PipelineStepVariable>>(
        pipeline_step_variable: T,
        value: Option<String>,
    ) -> Self {
        Self {
            id: pipeline_step_variable.borrow().id().clone(),
            name: pipeline_step_variable.borrow().name().clone(),
            description: pipeline_step_variable.borrow().description().clone(),
            category: pipeline_step_variable.borrow().category().clone(),
            required: pipeline_step_variable.borrow().required().clone(),
            value,
        }
    }

    /// Returns `true` if the variable is an instance of a reference to global data.
    pub fn is_global_data_reference(&self) -> bool {
        self.category.eq(&PipelineStepVariableCategory::Global)
    }
}

#[cfg(test)]
mod tests {

    use crate::model::internal::step::PipelineStepStatus;

    use super::*;

    #[test]
    fn test_variable_is_global_data_reference_false() {
        let pipeline_step_variable: PipelineStepVariable = serde_json::from_str(
            "
            {
                \"id\": \"option\",
                \"name\": \"Option\",
                \"description\": \"An option dropdown.\",
                \"category\": {
                    \"tag\": \"Option\",
                    \"content\": [
                        {
                            \"name\": \"Option 1\",
                            \"value\": \"option1\"
                        },
                        {
                            \"name\": \"Option 2\",
                            \"value\": \"option2\"
                        }
                    ]
                }
            }
            ",
        )
        .unwrap();
        let value = Some("dummy".to_string());
        let experiment_variable =
            ExperimentPipelineStepVariable::from_internal(&pipeline_step_variable, value.clone());
        assert!(!experiment_variable.is_global_data_reference());
    }

    #[test]
    fn test_variable_is_global_data_reference_true() {
        let pipeline_step_variable: PipelineStepVariable = serde_json::from_str(
            "
            {
                \"id\": \"global\",
                \"name\": \"Global\",
                \"description\": \"A global reference.\",
                \"category\": {
                    \"tag\": \"Global\"
                }
            }
            ",
        )
        .unwrap();
        let value = Some("dummy".to_string());
        let experiment_variable =
            ExperimentPipelineStepVariable::from_internal(&pipeline_step_variable, value.clone());
        assert!(experiment_variable.is_global_data_reference());
    }

    #[test]
    fn test_variable_from_internal() {
        let pipeline_step_variable: PipelineStepVariable = serde_json::from_str(
            "
            {
                \"id\": \"option\",
                \"name\": \"Option\",
                \"description\": \"An option dropdown.\",
                \"category\": {
                    \"tag\": \"Option\",
                    \"content\": [
                        {
                            \"name\": \"Option 1\",
                            \"value\": \"option1\"
                        },
                        {
                            \"name\": \"Option 2\",
                            \"value\": \"option2\"
                        }
                    ]
                }
            }
            ",
        )
        .unwrap();
        let value = Some("dummy".to_string());
        let experiment_variable =
            ExperimentPipelineStepVariable::from_internal(&pipeline_step_variable, value.clone());
        assert_eq!(experiment_variable.id(), pipeline_step_variable.id());
        assert_eq!(experiment_variable.name(), pipeline_step_variable.name());
        assert_eq!(experiment_variable.description(), pipeline_step_variable.description());
        assert_eq!(experiment_variable.category(), pipeline_step_variable.category());
        assert_eq!(experiment_variable.required(), pipeline_step_variable.required());
        assert_eq!(experiment_variable.value(), &value);
    }

    #[test]
    fn test_step_from_internal() {
        let pipeline_step: PipelineStepBlueprint = serde_json::from_str(
            "
            {
                \"id\": \"fastqc\",
                \"name\": \"FastQC\",
                \"description\": \"Performs a quality control.\",
                \"container\": \"fastqc\",
                \"dependencies\": [\"123\", \"456\"],
                \"variables\": [
                    {
                        \"id\": \"bool\",
                        \"name\": \"Boolean\",
                        \"description\": \"A boolean checkbox.\",
                        \"category\": {
                            \"tag\": \"Boolean\"
                        },
                        \"required\": true
                    },
                    {
                        \"id\": \"global\",
                        \"name\": \"Global\",
                        \"description\": \"A global data reference.\",
                        \"category\": {
                            \"tag\": \"Global\"
                        },
                        \"required\": false
                    }
                ]
            }
            ",
        )
        .unwrap();
        let mut values = HashMap::new();
        let value_bool = "true".to_string();
        let value_global = "global ID".to_string();

        values.insert("fastqcglobal".to_string(), value_global.clone());
        values.insert("fastqcbool".to_string(), value_bool.clone());

        let mut stati = HashMap::new();
        stati.insert("fastqc".to_string(), PipelineStepStatus::Pending.to_string());

        let experiment_step =
            ExperimentPipelineStepBlueprint::from_internal(&pipeline_step, values, stati);
        assert_eq!(experiment_step.id(), pipeline_step.id());
        assert_eq!(experiment_step.name(), pipeline_step.name());
        assert_eq!(experiment_step.description(), pipeline_step.description());
        assert_eq!(experiment_step.container(), pipeline_step.container());
        assert_eq!(experiment_step.dependencies(), pipeline_step.dependencies());
        assert_eq!(experiment_step.variables().len(), pipeline_step.variables().len());
        assert_eq!(experiment_step.status(), &Some(PipelineStepStatus::Pending.to_string()));

        let experiment_vars = experiment_step.variables();
        let pipeline_vars = pipeline_step.variables();
        for i in 0..experiment_vars.len() {
            assert_eq!(experiment_vars[i].id(), pipeline_vars[i].id());
            assert_eq!(experiment_vars[i].name(), pipeline_vars[i].name());
            assert_eq!(experiment_vars[i].description(), pipeline_vars[i].description());
            assert_eq!(experiment_vars[i].category(), pipeline_vars[i].category());
            assert_eq!(experiment_vars[i].required(), pipeline_vars[i].required());
            if experiment_vars[i].id() == "bool" {
                assert_eq!(experiment_vars[i].value(), &Some(value_bool.clone()));
            } else {
                assert_eq!(experiment_vars[i].value(), &Some(value_global.clone()));
            }
        }
    }

    #[test]
    fn test_pipeline_from_internal() {
        let pipeline: PipelineBlueprint = serde_json::from_str(
            "
            {
                \"id\": \"testing_pipeline\",
                \"name\": \"Testing pipeline\",
                \"version\": \"1.0.0\",
                \"description\": \"This pipeline is for testing purposes.\",
                \"global_variables\": [
                    {
                        \"id\": \"global_bool\",
                        \"name\": \"Global boolean\",
                        \"description\": \"A global boolean checkbox.\",
                        \"category\": {
                            \"tag\": \"Boolean\"
                        },
                        \"required\": true
                    }
                ],
                \"steps\": [
                    {
                        \"id\": \"fastqc1\",
                        \"name\": \"FastQC\",
                        \"description\": \"Performs a quality control.\",
                        \"container\": \"fastqc\",
                        \"dependencies\": [\"123\", \"456\"],
                        \"variables\": [
                            {
                                \"id\": \"bool\",
                                \"name\": \"Boolean\",
                                \"description\": \"A boolean checkbox.\",
                                \"category\": {
                                    \"tag\": \"Boolean\"
                                },
                                \"required\": true
                            },
                            {
                                \"id\": \"global\",
                                \"name\": \"Global\",
                                \"description\": \"A global data reference.\",
                                \"category\": {
                                    \"tag\": \"Global\"
                                },
                                \"required\": false
                            }
                        ]
                    },
                    {
                        \"id\": \"fastqc2\",
                        \"name\": \"FastQC\",
                        \"description\": \"Performs a quality control.\",
                        \"container\": \"fastqc\",
                        \"dependencies\": [\"123\", \"456\"],
                        \"variables\": [
                            {
                                \"id\": \"bool\",
                                \"name\": \"Boolean\",
                                \"description\": \"A boolean checkbox.\",
                                \"category\": {
                                    \"tag\": \"Boolean\"
                                },
                                \"required\": true
                            },
                            {
                                \"id\": \"global\",
                                \"name\": \"Global\",
                                \"description\": \"A global data reference.\",
                                \"category\": {
                                    \"tag\": \"Global\"
                                },
                                \"required\": false
                            }
                        ]
                    } 
                ]
            }        
            ",
        )
        .unwrap();

        let mut values_global = HashMap::new();
        values_global.insert("global_bool".to_string(), "true".to_string());

        let mut values_step = HashMap::new();
        values_step.insert("fastqc1global".to_string(), "01".to_string());
        values_step.insert("fastqc1bool".to_string(), "00".to_string());
        values_step.insert("fastqc2bool".to_string(), "10".to_string());
        values_step.insert("fastqc2global".to_string(), "11".to_string());

        let mut stati = HashMap::new();
        stati.insert("fastqc2".to_string(), PipelineStepStatus::Failed.to_string());

        let experiment_pipeline = ExperimentPipelineBlueprint::from_internal(
            &pipeline,
            values_global,
            values_step,
            stati,
        );
        assert_eq!(experiment_pipeline.id(), pipeline.id());
        assert_eq!(experiment_pipeline.name(), pipeline.name());
        assert_eq!(experiment_pipeline.version(), pipeline.version());
        assert_eq!(experiment_pipeline.description(), pipeline.description());
        assert_eq!(experiment_pipeline.global_variables().len(), pipeline.global_variables().len());
        assert_eq!(experiment_pipeline.steps().len(), pipeline.steps().len());

        let experiment_globals = experiment_pipeline.global_variables();
        let pipeline_globals = pipeline.global_variables();
        for j in 0..experiment_globals.len() {
            assert_eq!(experiment_globals[j].id(), pipeline_globals[j].id());
            assert_eq!(experiment_globals[j].name(), pipeline_globals[j].name());
            assert_eq!(experiment_globals[j].description(), pipeline_globals[j].description());
            assert_eq!(experiment_globals[j].category(), pipeline_globals[j].category());
            assert_eq!(experiment_globals[j].required(), pipeline_globals[j].required());
            assert_eq!(experiment_globals[j].value(), &Some("true".to_string()));
        }

        let experiment_steps = experiment_pipeline.steps();
        let pipeline_steps = pipeline.steps();
        for i in 0..experiment_steps.len() {
            assert_eq!(experiment_steps[i].id(), pipeline_steps[i].id());
            assert_eq!(experiment_steps[i].name(), pipeline_steps[i].name());
            assert_eq!(experiment_steps[i].description(), pipeline_steps[i].description());
            assert_eq!(experiment_steps[i].container(), pipeline_steps[i].container());
            assert_eq!(experiment_steps[i].dependencies(), pipeline_steps[i].dependencies());
            assert_eq!(experiment_steps[i].variables().len(), pipeline_steps[i].variables().len());

            if experiment_steps[i].id() == "fastqc2" {
                assert_eq!(
                    experiment_steps[i].status(),
                    &Some(PipelineStepStatus::Failed.to_string())
                );
            } else {
                assert_eq!(experiment_steps[i].status(), &None);
            }

            let experiment_vars = experiment_steps[i].variables();
            let pipeline_vars = pipeline_steps[i].variables();
            for j in 0..experiment_vars.len() {
                assert_eq!(experiment_vars[j].id(), pipeline_vars[j].id());
                assert_eq!(experiment_vars[j].name(), pipeline_vars[j].name());
                assert_eq!(experiment_vars[j].description(), pipeline_vars[j].description());
                assert_eq!(experiment_vars[j].category(), pipeline_vars[j].category());
                assert_eq!(experiment_vars[j].required(), pipeline_vars[j].required());
                let expected_value = format!("{}{}", i, j);
                assert_eq!(experiment_vars[j].value(), &Some(expected_value));
            }
        }
    }
}
