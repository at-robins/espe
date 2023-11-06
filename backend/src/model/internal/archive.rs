use std::path::PathBuf;

use getset::Getters;
use serde::{Deserialize, Serialize};

/// The file extension of archive metadata files.
const ARCHIVE_METADATA_FILE_EXTENSION: &str = "meta";

/// The metadata of a downloadable archive.
#[derive(Debug, Clone, Getters, PartialEq, Serialize, Deserialize)]
pub struct ArchiveMetadata {
    /// The file name of the archive.
    #[getset(get = "pub")]
    file_name: String,
}

impl ArchiveMetadata {
    /// Creates a new archive metadata struct.
    ///
    /// # Parameters
    ///
    /// * `file_name` - the name of the archive file
    pub fn new<T: AsRef<str>>(file_name: T) -> Self {
        Self {
            file_name: sanitize_filename::sanitize(file_name),
        }
    }

    /// Returns the archive metadata path for the specified archive file.
    ///
    /// # Parameters
    ///
    /// * `archive_path` - the path of the archive file
    pub fn metadata_path<T: Into<PathBuf>>(archive_path: T) -> PathBuf {
        let mut archive_meta_path: PathBuf = archive_path.into();
        archive_meta_path.set_extension(ARCHIVE_METADATA_FILE_EXTENSION);
        archive_meta_path
    }
}

#[cfg(test)]
mod tests {

    use crate::application::config::Configuration;

    use super::*;

    #[test]
    fn test_archive_metadata_file_name() {
        let archive_name = "some_name";
        let metadata = ArchiveMetadata::new(archive_name);
        assert_eq!(metadata.file_name(), archive_name);
    }

    #[test]
    fn test_archive_metadata_file_name_invalid_characters() {
        let archive_name = "some/*?name";
        let metadata = ArchiveMetadata::new(archive_name);
        assert_eq!(metadata.file_name(), "somename");
    }

    #[test]
    fn test_metadata_path() {
        let config = Configuration::new("", "", "", "", "./application/context", "");
        let archive_id = "01234567-89ab-cdef-0123-456789abcdef";
        let archive_path = config.temporary_download_file_path(archive_id);
        let expected_meta_path: PathBuf =
            "./application/context/tmp/download/01234567-89ab-cdef-0123-456789abcdef.meta".into();
        assert_eq!(ArchiveMetadata::metadata_path(archive_path), expected_meta_path);
    }
}
