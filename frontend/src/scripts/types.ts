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

export type PipelineBlueprintDetail = {
  id: number;
  name: string;
  comment: PipelineStepStatus;
};

export type ExperimentUpload = {
  name: string;
  mail: string | undefined;
  comment: string | undefined;
  pipelineId: number;
};

export type GlobalDataDetails = {
  id: number;
  name: string;
  comment: string | null | undefined;
  creationTime: string;
};

export type GlobalDataFileDetails = {
  pathComponents: string[];
};

export type FileTreeNode = {
  id: string;
  label: string;
  children: FileTreeNode[];
  parents: string[];
  isFile: boolean;
  isUploaded: boolean;
};

export type ErrorResponse = {
  code: number;
  uuid: string;
  name: string;
  message: string;
};
