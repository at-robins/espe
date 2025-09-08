use std::{
    collections::HashMap,
    sync::{Arc, Weak},
};

use parking_lot::Mutex;

/// Manages the tracking of download events.
#[derive(Debug)]
pub struct DownloadTrackerManager {
    experiment_map: Arc<Mutex<HashMap<i32, u64>>>,
}

impl DownloadTrackerManager {
    /// Creates a new manager to track download events.
    pub fn new() -> Self {
        Self {
            experiment_map: Arc::new(Mutex::new(HashMap::new())),
        }
    }

    /// Registers a single download event for experiment output.
    /// Returns a [`DownloadTracker`] that when dropped unregisteres the event.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment to register
    pub fn track_experiment_output_download(&self, experiment_id: i32) -> DownloadTracker {
        let mut map_lock = self.experiment_map.lock();
        let new_value = map_lock
            .get(&experiment_id)
            .map(|value| *value)
            .unwrap_or(0)
            + 1;
        map_lock.insert(experiment_id, new_value);
        DownloadTracker {
            key: experiment_id,
            map_reference: Arc::downgrade(&self.experiment_map),
        }
    }

    /// Returns `true` if one or more downloads are currently tracked for
    /// the specified experiment.
    ///
    /// # Parameters
    ///
    /// * `experiment_id` - the ID of the experiment to check
    pub fn is_experiment_output_download_tracked(&self, experiment_id: i32) -> bool {
        self.experiment_map
            .lock()
            .get(&experiment_id)
            .map(|value| *value > 0)
            .unwrap_or(false)
    }
}

/// A download tracker that represents a tracked download.
/// When this element is dropped, the event will be unregistered.
#[derive(Debug)]
pub struct DownloadTracker {
    key: i32,
    map_reference: Weak<Mutex<HashMap<i32, u64>>>,
}

impl Drop for DownloadTracker {
    fn drop(&mut self) {
        if let Some(map) = self.map_reference.upgrade() {
            let mut map_lock = map.lock();
            let value = map_lock.get(&self.key).map(|value| *value).unwrap_or(0);
            if value > 1 {
                map_lock.insert(self.key, value - 1);
            } else {
                // Keeps the map reasonable in size even with a lot of experiments.
                map_lock.remove(&self.key);
                map_lock.shrink_to_fit();
            }
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::test_utility::DEFAULT_EXPERIMENT_ID;

    use super::*;

    #[test]
    fn test_track_experiment_output_download() {
        let manager = DownloadTrackerManager::new();
        assert!(!manager.is_experiment_output_download_tracked(DEFAULT_EXPERIMENT_ID));
        let _tracker_00 = manager.track_experiment_output_download(DEFAULT_EXPERIMENT_ID + 1);
        assert!(!manager.is_experiment_output_download_tracked(DEFAULT_EXPERIMENT_ID));
        let tracker_01 = manager.track_experiment_output_download(DEFAULT_EXPERIMENT_ID);
        let tracker_02 = manager.track_experiment_output_download(DEFAULT_EXPERIMENT_ID);
        assert!(manager.is_experiment_output_download_tracked(DEFAULT_EXPERIMENT_ID));
        std::mem::drop(tracker_01);
        assert!(manager.is_experiment_output_download_tracked(DEFAULT_EXPERIMENT_ID));
        std::mem::drop(tracker_02);
        assert!(!manager.is_experiment_output_download_tracked(DEFAULT_EXPERIMENT_ID));
    }
}
