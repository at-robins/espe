use std::{borrow::Borrow, collections::HashMap};

use chrono::NaiveDateTime;
use getset::Getters;
use serde::{Deserialize, Serialize};

use crate::model::{
    db::experiment_execution::ExperimentExecution,
    internal::pipeline_blueprint::{
        PipelineBlueprint, PipelineStepBlueprint, PipelineStepVariable,
        PipelineStepVariableCategory,
    },
};

/// The definition of a pipeline.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct ExperimentPipelineBlueprint {
    /// The unique raw ID of the pipeline.
    #[getset(get = "pub")]
    id: String,
    /// The unique sanitised ID of the pipeline step which can safely be used as file name.
    #[getset(get = "pub")]
    sanitised_id: String,
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
    /// * `metadata` - a map of step metadata, where the keys are the pipeline step ID
    pub fn from_internal<
        PipelineType: Borrow<PipelineBlueprint>,
        ValueGlobalMapType: Borrow<HashMap<String, String>>,
        ValueStepMapType: Borrow<HashMap<String, String>>,
        MetadataMapType: Borrow<HashMap<String, ExperimentPipelineStepBlueprintMetadata>>,
    >(
        pipeline: PipelineType,
        values_global: ValueGlobalMapType,
        values_step: ValueStepMapType,
        metadata: MetadataMapType,
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
                    metadata.borrow().get(s.id()),
                )
            })
            .collect();
        Self {
            id: pipeline.borrow().id().clone(),
            sanitised_id: pipeline.borrow().sanitised_id().clone(),
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
    /// The unique raw ID of the pipeline step.
    #[getset(get = "pub")]
    id: String,
    /// The unique sanitised ID of the pipeline step which can safely be used as file name.
    #[getset(get = "pub")]
    sanitised_id: String,
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
    /// The timestamp for start of the execution.
    #[getset(get = "pub")]
    time_start: Option<String>,
    /// The timestamp for end of the execution.
    #[getset(get = "pub")]
    time_end: Option<String>,
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
    /// * `metadata` - optional metadata for the pipeline step
    pub fn from_internal<
        StepType: Borrow<PipelineStepBlueprint>,
        ValueMapType: Borrow<HashMap<String, String>>,
        MetadataType: Borrow<ExperimentPipelineStepBlueprintMetadata>,
    >(
        pipeline_step: StepType,
        values: ValueMapType,
        metadata: Option<MetadataType>,
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
        let (status, time_start, time_end): (Option<String>, Option<String>, Option<String>) =
            if let Some(step_metadata) = metadata {
                let step_metadata = step_metadata.borrow();
                (
                    step_metadata.status().clone(),
                    step_metadata
                        .time_start()
                        .map(|start_time| start_time.to_string()),
                    step_metadata
                        .time_end()
                        .map(|end_time| end_time.to_string()),
                )
            } else {
                (None, None, None)
            };
        Self {
            id: pipeline_step.borrow().id().clone(),
            sanitised_id: pipeline_step.borrow().sanitised_id().clone(),
            name: pipeline_step.borrow().name().clone(),
            description: pipeline_step.borrow().description().clone(),
            container: pipeline_step.borrow().container().clone(),
            dependencies: pipeline_step.borrow().dependencies().clone(),
            variables,
            status,
            time_start,
            time_end,
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

/// The optional metadata attached to a pipeline step.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct ExperimentPipelineStepBlueprintMetadata {
    /// The execution status of the pipeline step.
    #[getset(get = "pub")]
    status: Option<String>,
    /// The timestamp for start of the execution.
    #[getset(get = "pub")]
    time_start: Option<NaiveDateTime>,
    /// The timestamp for end of the execution.
    #[getset(get = "pub")]
    time_end: Option<NaiveDateTime>,
}

impl ExperimentPipelineStepBlueprintMetadata {
    /// Extracts the metadata for a specific pipeline step from a vector of tagged metadata.
    /// The elements matching the specified pipeline ID are removed from the input vector.
    /// Ordering of the elements in the source vector is not preserved.
    ///
    /// # Returns
    ///
    /// A map where the keys are the respective pipeline step ID and the value holds the
    /// corresponding metadata.
    pub fn extract_metadata_map<T: AsRef<str>>(
        pipeline_id: T,
        tagged_metadata: &mut Vec<ExperimentPipelineStepBlueprintTaggedMetadata>,
    ) -> HashMap<String, Self> {
        let mut metadata_map = HashMap::new();
        for metadata_index in (0..tagged_metadata.len()).rev() {
            if tagged_metadata[metadata_index].pipeline_id() == pipeline_id.as_ref() {
                let step_metadata = tagged_metadata.swap_remove(metadata_index);
                metadata_map.insert(step_metadata.step_id, step_metadata.metadata);
            }
        }
        metadata_map
    }
}

impl From<ExperimentPipelineStepBlueprintTaggedMetadata>
    for ExperimentPipelineStepBlueprintMetadata
{
    fn from(tagged_metadata: ExperimentPipelineStepBlueprintTaggedMetadata) -> Self {
        tagged_metadata.metadata
    }
}

/// The optional metadata attached to a pipeline step
/// tagged with the respective origin IDs.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct ExperimentPipelineStepBlueprintTaggedMetadata {
    /// The pipeline ID the step ID belongs to.
    #[getset(get = "pub")]
    pipeline_id: String,
    /// The pipeline step ID the metadata belongs to.
    #[getset(get = "pub")]
    step_id: String,
    /// The pipeline step metadata.
    #[getset(get = "pub")]
    metadata: ExperimentPipelineStepBlueprintMetadata,
}

impl From<ExperimentExecution> for ExperimentPipelineStepBlueprintTaggedMetadata {
    fn from(execution: ExperimentExecution) -> Self {
        ExperimentPipelineStepBlueprintTaggedMetadata {
            pipeline_id: execution.pipeline_id,
            step_id: execution.pipeline_step_id,
            metadata: ExperimentPipelineStepBlueprintMetadata {
                status: Some(execution.execution_status),
                time_start: execution.start_time,
                time_end: execution.end_time,
            },
        }
    }
}

#[cfg(test)]
mod tests {

    use chrono::{NaiveDate, NaiveTime};

    use crate::model::db::experiment_execution::ExecutionStatus;

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

        let time_start = NaiveDateTime::new(
            NaiveDate::from_ymd_opt(1871, 1, 1).unwrap(),
            NaiveTime::from_hms_opt(12, 12, 12).unwrap(),
        );
        let time_end = NaiveDateTime::new(
            NaiveDate::from_ymd_opt(1871, 1, 1).unwrap(),
            NaiveTime::from_hms_opt(16, 16, 16).unwrap(),
        );
        let metadata = ExperimentPipelineStepBlueprintMetadata {
            status: Some(ExecutionStatus::Waiting.to_string()),
            time_start: Some(time_start.clone()),
            time_end: Some(time_end.clone()),
        };

        let experiment_step =
            ExperimentPipelineStepBlueprint::from_internal(&pipeline_step, values, Some(metadata));
        assert_eq!(experiment_step.id(), pipeline_step.id());
        assert_eq!(experiment_step.sanitised_id(), pipeline_step.sanitised_id());
        assert_eq!(experiment_step.name(), pipeline_step.name());
        assert_eq!(experiment_step.description(), pipeline_step.description());
        assert_eq!(experiment_step.container(), pipeline_step.container());
        assert_eq!(experiment_step.dependencies(), pipeline_step.dependencies());
        assert_eq!(experiment_step.variables().len(), pipeline_step.variables().len());
        assert_eq!(experiment_step.status(), &Some(ExecutionStatus::Waiting.to_string()));
        assert_eq!(experiment_step.time_start(), &Some(time_start.to_string()));
        assert_eq!(experiment_step.time_end(), &Some(time_end.to_string()));

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

        let time_start = NaiveDateTime::new(
            NaiveDate::from_ymd_opt(1871, 1, 1).unwrap(),
            NaiveTime::from_hms_opt(12, 12, 12).unwrap(),
        );
        let time_end = NaiveDateTime::new(
            NaiveDate::from_ymd_opt(1871, 1, 1).unwrap(),
            NaiveTime::from_hms_opt(16, 16, 16).unwrap(),
        );

        let mut metadata = HashMap::new();
        metadata.insert(
            "fastqc2".to_string(),
            ExperimentPipelineStepBlueprintMetadata {
                status: Some(ExecutionStatus::Failed.to_string()),
                time_start: Some(time_start.clone()),
                time_end: Some(time_end.clone()),
            },
        );

        let experiment_pipeline = ExperimentPipelineBlueprint::from_internal(
            &pipeline,
            values_global,
            values_step,
            metadata,
        );
        assert_eq!(experiment_pipeline.id(), pipeline.id());
        assert_eq!(experiment_pipeline.sanitised_id(), pipeline.sanitised_id());
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
            assert_eq!(experiment_steps[i].sanitised_id(), pipeline_steps[i].sanitised_id());
            assert_eq!(experiment_steps[i].name(), pipeline_steps[i].name());
            assert_eq!(experiment_steps[i].description(), pipeline_steps[i].description());
            assert_eq!(experiment_steps[i].container(), pipeline_steps[i].container());
            assert_eq!(experiment_steps[i].dependencies(), pipeline_steps[i].dependencies());
            assert_eq!(experiment_steps[i].variables().len(), pipeline_steps[i].variables().len());

            if experiment_steps[i].id() == "fastqc2" {
                assert_eq!(
                    experiment_steps[i].status(),
                    &Some(ExecutionStatus::Failed.to_string())
                );
                assert_eq!(experiment_steps[i].time_start(), &Some(time_start.to_string()));
                assert_eq!(experiment_steps[i].time_end(), &Some(time_end.to_string()));
            } else {
                assert_eq!(experiment_steps[i].status(), &None);
                assert_eq!(experiment_steps[i].time_start(), &None);
                assert_eq!(experiment_steps[i].time_end(), &None);
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

    #[test]
    fn test_tagged_metadata_from_experiment_execution() {
        let time_start = Some(NaiveDateTime::new(
            NaiveDate::from_ymd_opt(1871, 1, 1).unwrap(),
            NaiveTime::from_hms_opt(12, 12, 12).unwrap(),
        ));
        let time_end = None;
        let pipeline_id = "test_pipeline_id";
        let pipeline_step_id = "test_step_id";
        let status = ExecutionStatus::Waiting.to_string();

        let expected_metadata = ExperimentPipelineStepBlueprintMetadata {
            status: Some(status.clone()),
            time_start: time_start.clone(),
            time_end: time_end.clone(),
        };

        let execution = ExperimentExecution {
            id: 42,
            experiment_id: 21,
            pipeline_id: pipeline_id.to_string(),
            pipeline_step_id: pipeline_step_id.to_string(),
            execution_status: status,
            start_time: time_start,
            end_time: time_end,
            creation_time: time_start.clone().unwrap(),
        };
        let tagged_metadata = ExperimentPipelineStepBlueprintTaggedMetadata::from(execution);
        assert_eq!(tagged_metadata.pipeline_id(), pipeline_id);
        assert_eq!(tagged_metadata.step_id(), pipeline_step_id);
        assert_eq!(tagged_metadata.metadata(), &expected_metadata);
    }

    #[test]
    fn test_metadata_from_tagged_metadata() {
        let time_start = NaiveDateTime::new(
            NaiveDate::from_ymd_opt(1871, 1, 1).unwrap(),
            NaiveTime::from_hms_opt(12, 12, 12).unwrap(),
        );
        let time_end = NaiveDateTime::new(
            NaiveDate::from_ymd_opt(1871, 1, 1).unwrap(),
            NaiveTime::from_hms_opt(16, 16, 16).unwrap(),
        );
        let metadata = ExperimentPipelineStepBlueprintMetadata {
            status: Some(ExecutionStatus::Waiting.to_string()),
            time_start: Some(time_start.clone()),
            time_end: Some(time_end.clone()),
        };

        let tagged_metadata = ExperimentPipelineStepBlueprintTaggedMetadata {
            pipeline_id: "test_pipeline_id".to_string(),
            step_id: "test_step_id".to_string(),
            metadata: metadata.clone(),
        };
        assert_eq!(tagged_metadata.metadata(), &metadata)
    }

    #[test]
    fn test_metadata_extract_metadata_map() {
        // let time_start = Some(NaiveDateTime::new(
        //     NaiveDate::from_ymd_opt(1871, 1, 1).unwrap(),
        //     NaiveTime::from_hms_opt(12, 12, 12).unwrap(),
        // ));
        // let time_end = None;
        let pipeline_id_0 = "test_pipeline_id_0";
        let pipeline_id_1 = "test_pipeline_id_1";
        let status = ExecutionStatus::Waiting.to_string();
        let number_of_items: usize = 42;

        let mut tagged_metadata_0 = Vec::new();
        let mut tagged_metadata_1 = Vec::new();

        for step_id in 0..number_of_items {
            let dummy_metadata = ExperimentPipelineStepBlueprintMetadata {
                status: Some(status.clone()),
                time_start: Some(NaiveDateTime::new(
                    NaiveDate::from_ymd_opt(1871 + (number_of_items as i32), 1, 1).unwrap(),
                    NaiveTime::from_hms_opt(12, 12, 12).unwrap(),
                )),
                time_end: None,
            };
            tagged_metadata_0.push(ExperimentPipelineStepBlueprintTaggedMetadata {
                pipeline_id: pipeline_id_0.to_string(),
                step_id: step_id.to_string(),
                metadata: dummy_metadata.clone(),
            });
            tagged_metadata_1.push(ExperimentPipelineStepBlueprintTaggedMetadata {
                pipeline_id: pipeline_id_1.to_string(),
                step_id: step_id.to_string(),
                metadata: dummy_metadata,
            });
        }

        let mut tagged_metadata = tagged_metadata_0.clone();
        tagged_metadata.extend(tagged_metadata_1.clone());
        assert_eq!(tagged_metadata.len(), number_of_items * 2);
        // Removes and returns first set of pipeline associated metadata.
        let metadata_map_0 = ExperimentPipelineStepBlueprintMetadata::extract_metadata_map(
            pipeline_id_0,
            &mut tagged_metadata,
        );
        assert_eq!(tagged_metadata.len(), number_of_items);
        assert_eq!(metadata_map_0.len(), number_of_items);
        for metadata in tagged_metadata_0 {
            assert_eq!(metadata_map_0.get(metadata.step_id()).unwrap(), metadata.metadata());
        }
        // Removes and returns second set of pipeline associated metadata.
        let metadata_map_1 = ExperimentPipelineStepBlueprintMetadata::extract_metadata_map(
            pipeline_id_1,
            &mut tagged_metadata,
        );
        assert_eq!(tagged_metadata.len(), 0);
        assert_eq!(metadata_map_1.len(), number_of_items);
        for metadata in tagged_metadata_1 {
            assert_eq!(metadata_map_1.get(metadata.step_id()).unwrap(), metadata.metadata());
        }
    }
}
