# Encapsulated Sequencing Pipeline Engine

## Design philosophy

This encapsulated sequencing pipeline engine (ESPE) is intended as a simple, streamlined, well documented and tested solution for common sequencing data analysis setups.
Since transformation and analysis of sequencing data is typically a time and resource consuming task,
the main focus is not placed on scalability of a centralised processing pipeline for e.g. a core facility.
Rather an engine is provided to allow easy setup of the pipelines of choice within a specific workgroup.
The software itself is self-contained, so all processed data and associated information (database etc.) is stored locally with the software.
While this severely limits the scalability it allows for an easy setup, documentation and improves general usability.
For the same reasons execution of pipelines themselves is not paralellised.
Instead each pipeline step uses paralellisation whenever possible.
Additonally, each pipeline step can be configured as desired by the user and runs inside a container to prevent data and package conflicts, allow for running arbitrary pipeline steps without the need to install additional software and mitigate difficulties during setup.

## Architecture

### Frontend

Vue.js with the Quasar framework is used for easy implementation of the frontend.

### Backend

An actix-Rust-backend is used as a statically typed, fast and reliable backend.

### Database

SQLite is used as it provides a simple, self-contained and fast database.
Diesel is used for automated database migration.

### Containerisation

Docker was selected as tool for running and managing containers due to its ease of use and widespread adoption.

## Usage

The easiest way to run the application is by using Docker.

1. Install [Docker](https://docs.docker.com/engine/install/).
2. Build the container with the following command (requires root privileges): 
```bash
docker build -t "espe:dev" "https://github.com/at-robins/espe.git#master:container"
```
3. Create an ESPE run directory for the application (ensure enough storage space is available) with the following folder structure:
```bash
espe/context
espe/database
```
4. Start the container with the following command replacing the ```/path/to/espe``` 
with the path to the ESPE run directory (requires root privileges):
```bash
docker run \
   -e SERVER_ADDRESS=0.0.0.0 \
   --mount type=bind,source=/path/to/espe/context,target=/srv/espe/application/context \
   --mount type=bind,source=/path/to/espe/database,target=/srv/espe/application/database \
   -v /var/run/docker.sock:/var/run/docker.sock \
   -p 80:8080 \
   -d \
   espe:dev
```
6. Open a browser and enter [localhost](http://localhost).

## Enviroment variables

Also see `backend/.env`.

| Variable        | Description                                                                 |
| --------------- | --------------------------------------------------------------------------- |
| CONTEXT_FOLDER  | the context folder of the application where all relevant data is stored     |
| DATABASE_URL    | the database location                                                       |
| LOG_LEVEL       | the minimum log level of the application (`debug`, `info`, `warn`, `error`) |
| PIPELINE_FOLDER | the folder storing pipeline definitions                                     |
| SERVER_ADDRESS  | the address of the server                                                   |
| SERVER_PORT     | the port of the server                                                      |

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
# The version requirement for SQLite is >= 3.35.
sudo apt install sqlite3 libsqlite3-dev
# Install docker for containerisation of the pipeline.
sudo apt install docker.io
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

## Testing

Backend testing:

```bash
cd backend
cargo test
```

Frontend testing:

```bash
cd frontend
npm run test:unit
```
