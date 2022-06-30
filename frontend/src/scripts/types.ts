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
  mail: string;
  comment: string;
  pipelineId: number;
};
