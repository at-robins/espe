use std::process::Command;

fn main() {
    // Only build the frontend in release mode to prevent "cargo check" from being blocked.
    if std::env::var("PROFILE").map_or(false, |profile| profile.to_lowercase() == "release") {
        let status_install = npm_command()
            .arg("install")
            .current_dir("../frontend/")
            .status()
            .expect("Failed to execute npm install.");
        if !status_install.success() {
            panic!("{:?}", status_install);
        }
        let status_build = npm_command()
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

    /// Returns the command to run npm.
    /// This is necessary since the [`Command`] does not recognise ```npm``` on Windows
    /// from the ```PATH``` since its missing ```.exe``` extension.
    fn npm_command() -> Command {
        if cfg!(target_os = "windows") {
            let mut command = Command::new("cmd");
            command.args(["/C", "npm"]);
            command
        } else {
            Command::new("npm")
        }
    }
}
