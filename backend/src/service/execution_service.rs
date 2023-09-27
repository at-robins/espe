use std::collections::{HashMap, HashSet};

use actix_web::web;

use crate::{
    application::{config::Configuration, error::SeqError, database::DatabaseManager},
    model::{
        db::experiment_execution::{ExecutionStatus, ExperimentExecution},
        internal::pipeline_blueprint::PipelineBlueprint,
    },
};

use super::{container_service::ContainerHandler, pipeline_service::LoadedPipelines};

/// A scheduler for execution of pipeline steps.
pub struct ExecutionScheduler {
    database_manager: web::Data<DatabaseManager>,
    loaded_pipelines: web::Data<LoadedPipelines>,
    handler: ContainerHandler,
}

impl ExecutionScheduler {
    pub fn new(
        config: web::Data<Configuration>,
        database_manager: web::Data<DatabaseManager>,
        loaded_pipelines: web::Data<LoadedPipelines>,
    ) -> Self {
        Self {
            database_manager: web::Data::clone(&database_manager),
            loaded_pipelines: web::Data::clone(&loaded_pipelines),
            handler: ContainerHandler::new(config, database_manager, loaded_pipelines),
        }
    }

    pub fn update_pipeline_execution(&mut self) -> Result<(), SeqError> {
        if self.handler.update()? {
            let mut connection = self.database_manager.database_connection()?;
            // First check if there are execution steps that have the running status, but are not currently executed.
            // This might happen if the application closed unexpectedly.
            let mut running =
                ExperimentExecution::get_by_status(ExecutionStatus::Running, &mut connection)?;
            if !running.is_empty() {
                log::info!("Found inactive execution entries with a running status: {:?}", running);
                self.handler.start(running.remove(0))?;
            } else {
                if let Some(next_step) = self.get_next_execution_step()? {
                    self.handler.start(next_step)?;
                }
            }
        }
        Ok(())
    }

    /// Returns the next execution step to be processed if any.
    fn get_next_execution_step(&self) -> Result<Option<ExperimentExecution>, SeqError> {
        let mut connection = self.database_manager.database_connection()?;
        let mut dependency_map_blueprint: HashMap<PipelineStepKey, Vec<PipelineStepKey>> =
            HashMap::new();
        let contextualised_pipelines = self.loaded_pipelines.pipelines();
        // The pipelines are only collected into a variable to prevent references of the dependency map being dropped.
        let _pipelines: Vec<&PipelineBlueprint> = contextualised_pipelines
            .iter()
            .map(|pipeline| pipeline.pipeline())
            .map(|pipeline| {
                // This happens inside the mapping for efficiency reasons.
                for step in pipeline.steps() {
                    dependency_map_blueprint.insert(
                        PipelineStepKey {
                            pipeline_id: pipeline.id(),
                            pipeline_step_id: step.id(),
                        },
                        step.dependencies()
                            .iter()
                            .map(|dependency| PipelineStepKey {
                                pipeline_id: pipeline.id(),
                                pipeline_step_id: dependency,
                            })
                            .collect(),
                    );
                }
                pipeline
            })
            .collect();
        let waiting =
            ExperimentExecution::get_by_status(ExecutionStatus::Waiting, &mut connection)?;
        let finished =
            ExperimentExecution::get_by_status(ExecutionStatus::Finished, &mut connection)?;
        let finished_step_ids: HashSet<ExperimentPipelineStepKey> = finished
            .iter()
            .map(|step| ExperimentPipelineStepKey {
                experiment_id: step.experiment_id,
                pipeline_id: &step.pipeline_id,
                pipeline_step_id: &step.pipeline_step_id,
            })
            .collect();
        let next: Option<ExperimentExecution> = waiting.into_iter().find(|step| {
            // Gets all dependencies of the current step.
            if let Some(dependencies) = dependency_map_blueprint.get(&PipelineStepKey {
                pipeline_id: &step.pipeline_id,
                pipeline_step_id: &step.pipeline_step_id,
            }) {
                dependencies
                    .iter()
                    // Maps the dependencies to experiment specific dependencies related to the tested execution step.
                    .map(|dependency| ExperimentPipelineStepKey {
                        experiment_id: step.experiment_id,
                        pipeline_id: dependency.pipeline_id,
                        pipeline_step_id: dependency.pipeline_step_id,
                    })
                    // Tests if all dependencies did already complete successfully.
                    .all(|dependency| finished_step_ids.contains(&dependency))
            } else {
                log::warn!(
                    "The pipeline step {:?} is not present in the pipeline blueprints.",
                    step
                );
                false
            }
        });
        Ok(next)
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
struct PipelineStepKey<'a> {
    pub pipeline_id: &'a String,
    pub pipeline_step_id: &'a String,
}

#[derive(Debug, PartialEq, Eq, Hash)]
struct ExperimentPipelineStepKey<'a> {
    pub experiment_id: i32,
    pub pipeline_id: &'a String,
    pub pipeline_step_id: &'a String,
}
