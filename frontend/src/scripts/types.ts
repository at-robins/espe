import type { AxiosRequestConfig } from "axios";
import type { ComponentPublicInstance } from "vue";
import Poller from "../components/shared/Poller.vue";

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

export enum ExperimentExecutionStatus {
  Aborted = "Aborted",
  Failed = "Failed",
  Finished = "Finished",
  Running = "Running",
  Waiting = "Waiting",
  None = "None",
}

export type ExperimentStepLog = {
  stdout: string | null | undefined;
  stderr: string | null | undefined;
  exitCode: string | null | undefined;
};

export type ExperimentStepLogs = {
  build: ExperimentStepLog;
  run: ExperimentStepLog;
};

export type PollerPostData = {
  config: AxiosRequestConfig;
  data: any;
};

/**
 * An interface to interact with the Poller component.
 */
export interface PollerInterface
  extends ComponentPublicInstance<typeof Poller> {
  /**
   * Performs an immidiate polling action.
   */
  pollNow: () => void;
}
