export type PipelineStepDetail = {
  id: number;
  name: string;
  status: PipelineStepStatus;
  creationTime: string;
};

export enum PipelineStepStatus {
  Running = "Running",
  Pending = "Pending",
  Success = "Success",
  Failed = "Failed",
}

export type ExperimentUpload = {
  name: string;
  mail: string | undefined;
  comment: string | undefined;
  pipelineId: string;
};

export type ExperimentDetails = {
  id: number;
  name: string;
  pipelineId: string | null | undefined;
  mail: string | null | undefined;
  comment: string | null | undefined;
  creationTime: string;
};

export type GlobalDataDetails = {
  id: number;
  name: string;
  comment: string | null | undefined;
  creationTime: string;
};

export type FileDetails = {
  pathComponents: string[];
  isFile: boolean;
};

export type FilePath = {
  pathComponents: string[];
};

export type FileTreeNode = {
  id: string;
  label: string;
  children: FileTreeNode[];
  parents: string[];
  isFile: boolean;
  isUploaded: boolean;
  error: string | null;
};

export type ErrorResponse = {
  code: number;
  uuid: string;
  name: string;
  message: string;
};
