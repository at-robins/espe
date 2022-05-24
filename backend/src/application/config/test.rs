use super::*;

#[test]
/// Tests if the `application_configuration_folder_path` function returns the correct path.
fn test_application_configuration_folder_path() {
    let p = Configuration::application_configuration_folder_path();
    assert_eq!(p.to_str().unwrap(), "application");
    assert!(p.is_relative());
}

#[test]
/// Tests if the `application_configuration_file_path` function returns the correct path.
fn test_application_configuration_file_path() {
    let p = Configuration::application_configuration_file_path();
    assert_eq!(p.to_str().unwrap(), "application/configuration.json");
    assert!(p.is_relative());
}
