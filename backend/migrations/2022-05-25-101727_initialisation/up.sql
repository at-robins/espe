CREATE TABLE experiment (
    id INTEGER PRIMARY KEY NOT NULL,
    experiment_name TEXT NOT NULL,
    mail TEXT NOT NULL,
    pipeline_id INTEGER,
    comment TEXT NOT NULL,
    creation_time TEXT NOT NULL,
    FOREIGN KEY (pipeline_id) REFERENCES pipeline (id) ON UPDATE CASCADE ON DELETE
    SET NULL
);
CREATE TABLE pipeline (
    id INTEGER PRIMARY KEY NOT NULL,
    pipeline_name TEXT NOT NULL,
    comment TEXT NOT NULL,
    creation_time TEXT NOT NULL
);
CREATE TABLE pipeline_step (
    id INTEGER PRIMARY KEY NOT NULL,
    pipeline_id INTEGER,
    execution_type TEXT NOT NULL,
    execution_configuration TEXT NOT NULL,
    ordering INTEGER NOT NULL,
    creation_time TEXT NOT NULL,
        FOREIGN KEY (pipeline_id) REFERENCES pipeline (id) ON UPDATE CASCADE ON DELETE CASCADE
)