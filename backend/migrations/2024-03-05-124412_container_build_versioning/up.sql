CREATE TABLE pipeline_build_register (
    id INTEGER PRIMARY KEY NOT NULL,
    pipeline_id TEXT NOT NULL,
    pipeline_step_id TEXT NOT NULL,
    pipeline_version TEXT NOT NULL,
    creation_time DATETIME NOT NULL,
    UNIQUE(pipeline_id, pipeline_step_id)
);