use std::{
    cell::RefCell,
    collections::HashMap,
    fs::File,
    io::{Cursor, Read, Write},
    path::{Path, PathBuf},
    rc::Rc,
    sync::{Arc, Weak},
    task::Poll,
};

use actix_web::web::Bytes;
use parking_lot::Mutex;
use zip::{write::StreamWriter, ZipWriter};

use crate::application::error::{
    SeqError, SeqErrorType, DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
};

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

const ARCHIVE_STREAM_READ_BUFFER_SIZE: usize = 512 * 1024;
const ARCHIVE_STREAM_WRITE_BUFFER_SIZE: usize = 4 * ARCHIVE_STREAM_READ_BUFFER_SIZE;

#[derive(Clone)]
/// A sharable cursor allowing interior mutability of the archive writer.
struct ArchiveCursor {
    inner: Rc<RefCell<Cursor<[u8; ARCHIVE_STREAM_WRITE_BUFFER_SIZE]>>>,
}

impl ArchiveCursor {
    /// Creates a new sharable [`ArchiveCursor`].
    fn new() -> Self {
        Self {
            inner: Rc::new(RefCell::new(Cursor::new([0; ARCHIVE_STREAM_WRITE_BUFFER_SIZE]))),
        }
    }

    /// Returns the content up to the current position.
    fn get_content(&self) -> Bytes {
        let pos = self.inner.borrow().position() as usize;
        if pos == 0 {
            Bytes::new()
        } else {
            Bytes::copy_from_slice(&self.inner.borrow().get_ref()[0..pos])
        }
    }

    /// Resets the current position to the start of the cursor.
    fn reset(&self) {
        self.inner.borrow_mut().set_position(0);
    }

    /// Gets the current content and then resets the cursor.
    fn take_content(&self) -> Bytes {
        let content = self.get_content();
        self.reset();
        content
    }
}

impl std::io::Write for ArchiveCursor {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.inner.borrow_mut().write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.borrow_mut().flush()
    }
}

/// An archive creation stream.
pub struct ArchiveStream {
    source: PathBuf,
    path_queue: Vec<PathBuf>,
    writer: Option<ZipWriter<StreamWriter<ArchiveCursor>>>,
    cursor: ArchiveCursor,
    options: zip::write::SimpleFileOptions,
    processing_file: Option<File>,
    buffer_read: Box<[u8; ARCHIVE_STREAM_READ_BUFFER_SIZE]>,
}

impl ArchiveStream {
    /// Creates a new [`ArchiveStream`] from the specified source path.
    ///
    /// # Parameters
    ///
    /// * `source` - the path to archive
    pub fn new<T: AsRef<Path>>(source: T) -> Result<Self, SeqError> {
        let source = std::fs::canonicalize(&source).map_err(|err| {
            SeqError::from(err)
                .chain(format!("Failed to resolve path {}.", source.as_ref().display()))
        })?;
        let cursor = ArchiveCursor::new();
        let writer = Some(ZipWriter::new_stream(cursor.clone()));
        let options = zip::write::SimpleFileOptions::default()
            .large_file(true)
            .compression_method(zip::CompressionMethod::Stored)
            .compression_level(None);
        let buffer_read = Box::new([0; ARCHIVE_STREAM_READ_BUFFER_SIZE]);
        if source.is_file() {
            Ok(Self {
                source: source
                    .parent()
                    .ok_or(SeqError::new(
                        "Archive stream generation error",
                        SeqErrorType::InternalServerError,
                        format!(
                            "The specified file \"{}\"does not have a parent.",
                            source.display()
                        ),
                        DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
                    ))?
                    .to_path_buf(),
                path_queue: vec![source],
                writer,
                cursor,
                options,
                processing_file: None,
                buffer_read,
            })
        } else {
            let path_queue = Self::directory_paths(&source).map_err(|err| {
                err.chain(format!(
                    "Failed to get directory entries for archive generation of \"{}\"",
                    source.display()
                ))
            })?;
            Ok(Self {
                source,
                path_queue,
                writer,
                cursor,
                options,
                processing_file: None,
                buffer_read,
            })
        }
    }

    /// Get the relative path of a path queue element.
    ///
    /// # Panics
    ///
    /// If the supplied path is not a child of the source path.
    fn relative_path<'a, T: AsRef<Path>>(&self, path: &'a T) -> &'a Path {
        path.as_ref().strip_prefix(&self.source).expect(&format!(
            "The specified path \"{}\" must be a child of the source path \"{}\"!",
            path.as_ref().display(),
            &self.source.display()
        ))
    }

    /// Gets all the entries inside the specified directory as paths.
    ///
    /// # Parameters
    ///
    /// * `path` - the directory to get the entry paths from
    fn directory_paths<T: AsRef<Path>>(path: T) -> Result<Vec<PathBuf>, SeqError> {
        match std::fs::read_dir(&path) {
            Ok(dir) => {
                let mut valid_paths = Vec::new();
                for entry in dir {
                    match entry {
                        Ok(valid_entry) => valid_paths.push(valid_entry.path()),
                        Err(err) => {
                            return Err(SeqError::from(err).chain(format!(
                                "Failed to retrieve directory entry in \"{}\".",
                                path.as_ref().display()
                            )))
                        },
                    }
                }
                Ok(valid_paths)
            },
            Err(err) => Err(SeqError::from(err)
                .chain(format!("Failed to read directory \"{}\".", path.as_ref().display()))),
        }
    }

    /// Reads a chunk of the currently open file into the archive cursor.
    /// If the end of the file is reached, it closes the open file.
    ///
    /// # Errors
    ///
    /// Returns errors if the file could not be read or the archive write failed.
    fn read_file_into_archive_cursor(&mut self) -> Result<(), SeqError> {
        if let Some(file) = self.processing_file.as_mut() {
            match file.read(self.buffer_read.as_mut()) {
                Ok(bytes_read) => {
                    if bytes_read == 0 {
                        // The end of the file was reached.
                        self.processing_file = None;
                    } else {
                        self.writer
                            .as_mut()
                            .ok_or(SeqError::new(
                                "Archive writer closed",
                                SeqErrorType::InternalServerError,
                                "A file read was called after the archive reader was closed.",
                                DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
                            ))?
                            .write_all(&self.buffer_read[0..bytes_read])
                            .map_err(|err| {
                                SeqError::from(err)
                                    .chain("Failed to write file read into the archive stream.")
                            })?;
                    }
                    Ok(())
                },
                Err(err) => {
                    Err(SeqError::from(err)
                        .chain("Failed to read bytes from file while archiving."))
                },
            }
        } else {
            Err(SeqError::new(
                "No file open",
                SeqErrorType::InternalServerError,
                "The archive stream file reader was called wihtout an open file.",
                DEFAULT_INTERNAL_SERVER_ERROR_EXTERNAL_MESSAGE,
            ))
        }
    }
}

impl futures::Stream for ArchiveStream {
    type Item = Result<Bytes, SeqError>;

    fn poll_next(
        self: std::pin::Pin<&mut Self>,
        _: &mut std::task::Context<'_>,
    ) -> Poll<Option<Self::Item>> {
        let archive_stream = self.get_mut();
        if archive_stream.writer.is_none() {
            // All contents have been written.
            Poll::Ready(None)
        } else if archive_stream.path_queue.is_empty() {
            // All files and directories have been processed.
            // Finialises the archive.
            // The unwrap works, as the contents of the writer have been tested above.
            let finished_writer = archive_stream.writer.take().unwrap();
            match finished_writer.finish() {
                Ok(_) => Poll::Ready(Some(Ok(archive_stream.cursor.take_content()))),
                Err(err) => {
                    Poll::Ready(Some(Err(SeqError::from(err).chain("Failed to finish archive."))))
                },
            }
        } else {
            if archive_stream
                .path_queue
                .last()
                .expect("The path stack must have been verified to not be empty here!")
                .is_dir()
            {
                let directory_path_absolute = archive_stream.path_queue.pop().unwrap();
                let directory_path_relative =
                    archive_stream.relative_path(&directory_path_absolute);
                if let Err(err) = archive_stream
                    .writer
                    .as_mut()
                    .unwrap()
                    .add_directory_from_path(directory_path_relative, archive_stream.options)
                {
                    return Poll::Ready(Some(Err(SeqError::from(err).chain(format!(
                        "Failed to add directory \"{}\" to archive.",
                        directory_path_absolute.display()
                    )))));
                }
                match ArchiveStream::directory_paths(directory_path_absolute) {
                    Ok(directories) => archive_stream.path_queue.extend(directories),
                    Err(err) => {
                        return Poll::Ready(Some(Err(err.chain(format!(
                        "Failed to load directory structure during archive generation of \"{}\".",
                        archive_stream.source.display()
                    )))))
                    },
                }
            } else {
                // Should be a file.
                // If a file is currently read contiunue processing the file.
                if archive_stream.processing_file.is_some() {
                    if let Err(err) = archive_stream.read_file_into_archive_cursor() {
                        return Poll::Ready(Some(Err(err.chain(format!(
                            "Failed to read file \"{:?}\" during creation of archive \"{}\".",
                            archive_stream.path_queue.last(),
                            archive_stream.source.display()
                        )))));
                    }
                    // If the file was finished remove it from the processing list.
                    if archive_stream.processing_file.is_none() {
                        archive_stream.path_queue.pop();
                    }
                } else {
                    // A new file will be opened and processed.
                    let current_file_path_absolute = archive_stream.path_queue.last().unwrap();
                    let current_file_path_relative =
                        archive_stream.relative_path(current_file_path_absolute);
                    if let Err(err) = archive_stream
                        .writer
                        .as_mut()
                        .unwrap()
                        .start_file_from_path(current_file_path_relative, archive_stream.options)
                    {
                        return Poll::Ready(Some(Err(SeqError::from(err).chain(format!(
                            "Failed to start file \"{:?}\" during creation of archive \"{}\".",
                            current_file_path_absolute.display(),
                            archive_stream.source.display()
                        )))));
                    }
                    match std::fs::File::open(current_file_path_absolute) {
                        Ok(file) => archive_stream.processing_file = Some(file),
                        Err(err) => {
                            return Poll::Ready(Some(Err(SeqError::from(err).chain(format!(
                                "Failed to open file \"{:?}\" during creation of archive \"{}\".",
                                current_file_path_absolute.display(),
                                archive_stream.source.display()
                            )))))
                        },
                    }
                }
            }
            Poll::Ready(Some(Ok(archive_stream.cursor.take_content())))
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
