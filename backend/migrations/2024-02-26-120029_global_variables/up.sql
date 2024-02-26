CREATE TABLE pipeline_global_variable (
    id INTEGER PRIMARY KEY NOT NULL,
    experiment_id INTEGER NOT NULL,
    pipeline_id TEXT NOT NULL,
    variable_id TEXT NOT NULL,
    variable_value TEXT,
    creation_time DATETIME NOT NULL,
    FOREIGN KEY (experiment_id) REFERENCES experiment (id) ON UPDATE CASCADE ON DELETE CASCADE
);