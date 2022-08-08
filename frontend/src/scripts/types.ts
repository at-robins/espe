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
  mail: string;
  comment: string;
  pipelineId: number;
};

export type ErrorResponse = {
  code: number;
  uuid: string;
  name: string;
  message: string;
};
