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
