# MriSequencingPipeline

## Design philosophy

This sequencing pipeline is intended as a simple, streamlined, well documented and tested solution for common sequencing approaches.
Since transformation and analysis of sequencing data is typically a time and resource consuming task,
the main focus is not placed on sclability of a centralised processing pipeline for e.g. a core facility.
Rather a Docker image is provided to allow easy setup of the pipeline within a specific workgroup.
Each image is thereby self-contained, so all processed data and associated information (database etc.) is 
stored within the respective container. While this severely limits the scalability it allows for an easy setup, documentation 
and improves general usability.
For the same reasons execution of pipelines themselves is not paralellised. Instead each pipeline step uses paralellisation whenever possible.

## Architecture

### Frontend
Vue.js with the Quasar framework is used for easy implementation of the frontend.

### Backend
An actix-Rust-backend is used as a statically typed, fast and reliable backend. 

### Database
SQLite is used as it provides a simple, self-contained and fast database.
Diesel is used for automated database migration.

## Usage
TODO: Insert once the Docker image is created.

## Enviroment variables
Also see ```backend/.env```.

| Variable | Description |
| --- | --- |
| DATABASE_URL | the database location |
| SERVER_ADDRESS | the address of the server |
| SERVER_PORT | the port of the server |

## Build (for developers)
Prerequisites:
```bash
# Install Rust.
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# Install Node.
curl -fsSL https://deb.nodesource.com/setup_lts.x | sudo bash -
sudo apt update
sudo apt install nodejs
# Install the requirements for the openssl crate.
sudo apt install libssl-dev
```

Conventional release build:
```bash
cd backend
cargo build --release
```

Building for development:
```bash
cd backend
cargo build
cd ../frontend
npm install
npm run build
```

## Run (for developers)
Prerequisites:
```bash
# Install SQLite and the Diesel dependency.
sudo apt install sqlite3 libsqlite3-dev
# Install Diesel.
cargo install diesel_cli --no-default-features --features sqlite
# Initialise the database.
cd backend
mkdir application # without the sub-folder Diesel will fail to create the database
diesel setup
```

Running a release build:
```bash
cd backend
cargo run --release
```

Running for development:
```bash
cd backend
cargo run
cd ../frontend
npm run build:dev
```