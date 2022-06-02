use std::process::Command;

fn main() {
    // Only build the frontend in release mode to prevent "cargo check" from being blocked.
    if std::env::var("PROFILE").map_or(false, |profile| profile.to_lowercase() == "release") {
        let status_install = Command::new("npm")
            .arg("install")
            .current_dir("../frontend/")
            .status()
            .expect("Failed to execute npm install.");
        if !status_install.success() {
            panic!("{:?}", status_install);
        }
        let status_build = Command::new("npm")
            .arg("run")
            .arg("build")
            .current_dir("../frontend/")
            .status()
            .expect("Failed to execute npm run build.");
        if !status_build.success() {
            panic!("{:?}", status_build);
        }
        println!("cargo:rerun-if-changed=../fronted");
    }
}
