use std::{
    collections::HashMap,
    sync::{Arc, Weak},
};

use parking_lot::Mutex;

/// Manages the tracking of download events.
#[derive(Debug)]
pub struct DownloadTrackerManager {
    experiment_map: Arc<Mutex<HashMap<(i32, Option<(String, String)>), u64>>>,
}

impl DownloadTrackerManager {
    /// Creates a new manager to track download events.
    pub fn new() -> Self {
        Self {
            experiment_map: Arc::new(Mutex::new(HashMap::new())),
        }
    }

    /// Registers a single download event for a whole experiment output.
    /// Returns a [`DownloadTracker`] that when dropped unregisteres the event.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment to register
    pub fn track_experiment_output_download_experiment(
        &self,
        experiment_id: i32,
    ) -> DownloadTracker {
        let mut map_lock = self.experiment_map.lock();
        let key = (experiment_id, None);
        let new_value = map_lock.get(&key).map(|value| *value).unwrap_or(0) + 1;
        map_lock.insert(key.clone(), new_value);
        DownloadTracker {
            key_experiment: experiment_id,
            key_pipeline_step: None,
            map_reference: Arc::downgrade(&self.experiment_map),
        }
    }

    /// Registers a single download event for a specific experiment output step.
    /// Returns a [`DownloadTracker`] that when dropped unregisteres the event.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment to register
    /// * `pipeline_id` - the ID of the pipeline to register
    /// * `pipeline_step_id` - the ID of the pipeline step to register
    pub fn track_experiment_output_download_step<S: Into<String>, T: Into<String>>(
        &self,
        experiment_id: i32,
        pipeline_id: S,
        pipeline_step_id: T,
    ) -> DownloadTracker {
        let mut map_lock = self.experiment_map.lock();
        let pipeline_id = pipeline_id.into();
        let pipeline_step_id = pipeline_step_id.into();
        // Tracks the experiment.
        let key_experiment = (experiment_id, None);
        let new_value_experiment = map_lock
            .get(&key_experiment)
            .map(|value| *value)
            .unwrap_or(0)
            + 1;
        map_lock.insert(key_experiment.clone(), new_value_experiment);
        // Tracks the step.
        let key_step = (experiment_id, Some((pipeline_id.clone(), pipeline_step_id.clone())));
        let new_value_step = map_lock.get(&key_step).map(|value| *value).unwrap_or(0) + 1;
        map_lock.insert(key_step.clone(), new_value_step);
        // Returns the tracker.
        DownloadTracker {
            key_experiment: experiment_id,
            key_pipeline_step: Some((pipeline_id, pipeline_step_id)),
            map_reference: Arc::downgrade(&self.experiment_map),
        }
    }

    /// Returns `true` if one or more downloads are currently tracked for
    /// the specified experiment.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment to check
    pub fn is_experiment_output_download_experiment_tracked(&self, experiment_id: i32) -> bool {
        self.experiment_map
            .lock()
            .get(&(experiment_id, None))
            .map(|value| *value > 0)
            .unwrap_or(false)
    }

    /// Returns `true` if one or more downloads are currently tracked for
    /// the specified experiment step.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment to check
    /// * `pipeline_id` - the ID of the pipeline to check
    /// * `pipeline_step_id` - the ID of the pipeline step to check
    pub fn is_experiment_output_download_step_tracked<S: Into<String>, T: Into<String>>(
        &self,
        experiment_id: i32,
        pipeline_id: S,
        pipeline_step_id: T,
    ) -> bool {
        self.experiment_map
            .lock()
            .get(&(experiment_id, Some((pipeline_id.into(), pipeline_step_id.into()))))
            .map(|value| *value > 0)
            .unwrap_or(false)
    }
}

/// A download tracker that represents a tracked download.
/// When this element is dropped, the event will be unregistered.
#[derive(Debug)]
pub struct DownloadTracker {
    key_experiment: i32,
    key_pipeline_step: Option<(String, String)>,
    map_reference: Weak<Mutex<HashMap<(i32, Option<(String, String)>), u64>>>,
}

impl Drop for DownloadTracker {
    fn drop(&mut self) {
        if let Some(map) = self.map_reference.upgrade() {
            let mut map_lock = map.lock();
            if self.key_pipeline_step.is_some() {
                Self::end_download_tracking(
                    &mut map_lock,
                    (self.key_experiment, self.key_pipeline_step.clone()),
                );
            }
            Self::end_download_tracking(&mut map_lock, (self.key_experiment, None));
        }
    }
}

impl DownloadTracker {
    /// Ends the download tracking for the specified key.
    ///
    /// # Parameters
    ///
    /// * `map_lock` - the lock holding the map
    /// * `key` - the key of the download
    fn end_download_tracking(
        map_lock: &mut parking_lot::lock_api::MutexGuard<
            '_,
            parking_lot::RawMutex,
            HashMap<(i32, Option<(String, String)>), u64>,
        >,
        key: (i32, Option<(String, String)>),
    ) {
        let value = map_lock.get(&key).map(|value| *value).unwrap_or(0);
        if value > 1 {
            map_lock.insert(key, value - 1);
        } else {
            // Keeps the map reasonable in size even with a lot of experiments.
            map_lock.remove(&key);
            map_lock.shrink_to_fit();
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::test_utility::{
        DEFAULT_EXPERIMENT_ID, DEFAULT_PIPELINE_ID, DEFAULT_PIPELINE_STEP_ID,
    };

    use super::*;

    #[test]
    fn test_track_experiment_output_download_experiment() {
        let manager = DownloadTrackerManager::new();
        assert!(!manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        let _tracker_00 =
            manager.track_experiment_output_download_experiment(DEFAULT_EXPERIMENT_ID + 1);
        assert!(!manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        let tracker_01 = manager.track_experiment_output_download_experiment(DEFAULT_EXPERIMENT_ID);
        let tracker_02 = manager.track_experiment_output_download_experiment(DEFAULT_EXPERIMENT_ID);
        assert!(manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        std::mem::drop(tracker_01);
        assert!(manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        std::mem::drop(tracker_02);
        assert!(!manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
    }

    #[test]
    fn test_track_experiment_output_download_step() {
        // Tests default state before tracking.
        let manager = DownloadTrackerManager::new();
        assert!(!manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        assert!(!manager.is_experiment_output_download_step_tracked(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID
        ));
        // Adds some background.
        let _tracker_experiment_00 =
            manager.track_experiment_output_download_experiment(DEFAULT_EXPERIMENT_ID + 1);
        let _tracker_step_00 = manager.track_experiment_output_download_step(
            DEFAULT_EXPERIMENT_ID + 1,
            format!("{}a", DEFAULT_PIPELINE_ID),
            format!("{}a", DEFAULT_PIPELINE_STEP_ID),
        );
        assert!(!manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        assert!(!manager.is_experiment_output_download_step_tracked(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID
        ));
        // Adds a step tracker.
        let tracker_step_01 = manager.track_experiment_output_download_step(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
        );
        assert!(manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        assert!(manager.is_experiment_output_download_step_tracked(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID
        ));
        // Adds an additional step and an experiment tracker.
        let tracker_step_02 = manager.track_experiment_output_download_step(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID,
        );
        let tracker_exp_01 =
            manager.track_experiment_output_download_experiment(DEFAULT_EXPERIMENT_ID);
        assert!(manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        assert!(manager.is_experiment_output_download_step_tracked(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID
        ));
        // Drops the first tracker.
        std::mem::drop(tracker_step_01);
        assert!(manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        assert!(manager.is_experiment_output_download_step_tracked(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID
        ));
        // Drops the second tracker.
        std::mem::drop(tracker_step_02);
        assert!(manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        assert!(!manager.is_experiment_output_download_step_tracked(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID
        ));
        // Drops the experiment tracker.
        std::mem::drop(tracker_exp_01);
        assert!(!manager.is_experiment_output_download_experiment_tracked(DEFAULT_EXPERIMENT_ID));
        assert!(!manager.is_experiment_output_download_step_tracked(
            DEFAULT_EXPERIMENT_ID,
            DEFAULT_PIPELINE_ID,
            DEFAULT_PIPELINE_STEP_ID
        ));
    }
}
