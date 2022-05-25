# MriSequencingPipeline

## Design philosophy

This sequencing pipeline is intended as a simple, streamlined, well documented and tested solution for common sequencing approaches.
Since transformation and analysis of sequencing data is typically a time and resource consuming task,
the main focus is not placed on sclability of a centralised processing pipeline for e.g. a core facility.
Rather a Docker image is provided to allow easy setup of the pipeline within a specific workgroup.
Each image is thereby self-contained, so all processed data and associated information (database etc.) is 
stored within the respective container. While this severely limits the scalability it allows for an easy setup, documentation 
and improves general usability.

## Architecture

### Frontend
Vue.js with the Quasar framework is used for easy implementation of the frontend.

### Backend
An actix-Rust-backend is used as a statically typed, fast and reliable backend. 

### Database
SQLite is used as it provides a simple, self-contained and fast database.
Diesel is used for autmoated database migration.

## Enviroment variables
| Variable | Description |
| --- | --- |
| DATABASE_URL | the database location |
| SERVER_ADDRESS | the address of the server |
| SERVER_PORT | the port of the server |