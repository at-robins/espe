use std::process::Command;

fn main() {
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
