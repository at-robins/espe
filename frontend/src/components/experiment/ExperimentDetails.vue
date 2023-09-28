<template>
  <div class="q-pa-md">
    <q-card>
      <div v-if="!loadingError">
        <div v-if="isLoadingDetails || !experiment">
          <q-card-section class="flex-center row">
            <q-spinner color="primary" size="xl" />
          </q-card-section>
        </div>
        <div v-else>
          <q-card-section>
            <entity-title
              :title="experiment.name"
              :entity-id="experiment.id"
              endpoint-type="experiments"
              @update:title="updateTitle"
            />
          </q-card-section>
          <q-card-section>
            <experiment-status-indicator :status="status" />
          </q-card-section>
          <q-card-section>
            <entity-comment
              :comment="experiment.comment"
              :entity-id="experiment.id"
              endpoint-type="experiments"
            />
          </q-card-section>
          <q-card-section>
            <entity-pipeline
              :pipeline-id="experiment.pipelineId"
              :entity-id="experiment.id"
              endpoint-type="experiments"
              @update:selected-pipeline="selectedPipeline = $event"
            />
          </q-card-section>
          <q-card-section v-if="selectedPipeline">
            <entity-pipeline-variables
              :pipeline="selectedPipeline"
              :entity-id="experiment.id"
              endpoint-type="experiments"
            />
          </q-card-section>
          <q-card-section>
            <div class="row">
              <div class="q-ma-md">
                <q-btn :icon="symOutlinedAccountTree" round>
                  <q-tooltip>Upload files as experiment input.</q-tooltip>
                </q-btn>
              </div>
              <q-separator vertical class="q-ml-md q-mr-md" />
              <file-tree
                v-model="fileNodes"
                :base-directory-label="experiment ? experiment.name : 'Root'"
                @added-file="uploadFile"
                @added-folder="uploadFolder"
                @deleted-path="deletePath"
                @deleted-all="deleteAll"
                class="col"
              />
            </div>
          </q-card-section>
          <q-card-section>
            <div class="row">
              <div class="q-ma-md">
                <q-btn
                  v-if="status == ExperimentExecutionStatus.None"
                  :icon="matPlayCircle"
                  label="Submit"
                  color="positive"
                  :loading="isSubmitting"
                  @click="submitExperiment"
                >
                  <q-tooltip>
                    Submit the experiment for execution with the specified
                    pipeline.
                  </q-tooltip>
                </q-btn>
                <q-btn
                  v-if="
                    status == ExperimentExecutionStatus.Running ||
                    status == ExperimentExecutionStatus.Waiting
                  "
                  :icon="matStopCircle"
                  label="Abort"
                  color="negative"
                  @click="aborExperiment"
                >
                  <q-tooltip>
                    Abort the current experiment execution.
                  </q-tooltip>
                </q-btn>
                <q-btn
                  v-if="
                    status == ExperimentExecutionStatus.Aborted ||
                    status == ExperimentExecutionStatus.Failed ||
                    status == ExperimentExecutionStatus.Finished
                  "
                  :icon="matRestartAlt"
                  label="Restart"
                  color="positive"
                >
                  <q-tooltip> Restart the experiment execution. </q-tooltip>
                </q-btn>
              </div>
            </div>
          </q-card-section>
        </div>
      </div>
      <div v-else>
        <error-popup :error-response="loadingError" />
      </div>
    </q-card>
    <q-dialog v-model="showServerError" v-if="serverError">
      <error-popup :error-response="serverError" />
    </q-dialog>
  </div>
</template>

<script setup lang="ts">
import {
  type ErrorResponse,
  type FileDetails,
  type FileTreeNode,
  type FilePath,
  type ExperimentDetails,
} from "@/scripts/types";
import axios from "axios";
import { ref, onMounted, type Ref } from "vue";
import ErrorPopup from "../ErrorPopup.vue";
import FileTree from "../FileTree.vue";
import EntityComment from "../shared/EntityComment.vue";
import EntityTitle from "../shared/EntityTitle.vue";
import { symOutlinedAccountTree } from "@quasar/extras/material-symbols-outlined";
import EntityPipeline from "../shared/EntityPipeline.vue";
import EntityPipelineVariables from "../shared/EntityPipelineVariables.vue";
import ExperimentStatusIndicator from "../shared/ExperimentStatusIndicator.vue";
import type { PipelineBlueprint } from "@/scripts/pipeline-blueprint";
import { ExperimentExecutionStatus } from "@/scripts/types";
import {
  matPlayCircle,
  matRestartAlt,
  matStopCircle,
} from "@quasar/extras/material-icons";
import { useRouter } from "vue-router";

// The intervall in which status updates are requested from the server.
const POLLING_INTERVALL_MILLISECONDS = 10000;

const files: Ref<Array<FileDetails>> = ref([]);
const experiment: Ref<ExperimentDetails | null> = ref(null);
const fileNodes: Ref<Array<FileTreeNode>> = ref([]);
const isLoadingDetails = ref(false);
const loadingError: Ref<ErrorResponse | null> = ref(null);
const serverError: Ref<ErrorResponse | null> = ref(null);
const showServerError = ref(false);
const selectedPipeline: Ref<PipelineBlueprint | null> = ref(null);
const isAborting = ref(false);
const isSubmitting = ref(false);
const status = ref(ExperimentExecutionStatus.None);
const isPolling = ref(false);
const pollingError: Ref<ErrorResponse | null> = ref(null);
const router = useRouter();
const this_route = router.currentRoute.value.fullPath;

const props = defineProps({
  id: { type: String, required: true },
});

function updateTitle(title: string) {
  if (experiment.value) {
    experiment.value.name = title;
  }
}

function getFileTreeNodes(files: FileDetails[]): FileTreeNode[] {
  const tree: FileTreeNode[] = [];
  for (const globalDataFile of files) {
    let currentNodes = tree;
    for (let i = 0; i < globalDataFile.pathComponents.length; i++) {
      const pathComponent = globalDataFile.pathComponents[i];
      const found = currentNodes.find((node) => node.label === pathComponent);
      if (found) {
        currentNodes = found.children;
      } else {
        const parents = globalDataFile.pathComponents.slice(0, i);
        const newNode: FileTreeNode = {
          id: parents.join("") + pathComponent,
          label: pathComponent,
          children: [],
          parents: parents,
          isFile:
            i === globalDataFile.pathComponents.length - 1 &&
            globalDataFile.isFile,
          isUploaded: false,
          error: null,
        };
        currentNodes.push(newNode);
        currentNodes = newNode.children;
      }
    }
  }
  return tree;
}

onMounted(() => {
  loadDetails();
  pollStatusChanges();
});

/**
 * Initial loading of details from the server.
 */
function loadDetails() {
  isLoadingDetails.value = true;
  loadingError.value = null;
  axios
    .get("/api/experiments/" + props.id)
    .then((response) => {
      experiment.value = response.data;
      return axios.get("/api/files/experiments/" + props.id);
    })
    .then((response) => {
      files.value = response.data;
      fileNodes.value = getFileTreeNodes(files.value);
    })
    .catch((error) => {
      experiment.value = null;
      files.value = [];
      fileNodes.value = [];
      loadingError.value = error.response.data;
    })
    .finally(() => {
      isLoadingDetails.value = false;
    });
}

function getFileNodeByNode(query: FileTreeNode): FileTreeNode | undefined {
  let foundFileNodes: FileTreeNode[] | undefined = fileNodes.value;
  for (const parentLabel of query.parents) {
    foundFileNodes = foundFileNodes?.find(
      (val) => val.label === parentLabel
    )?.children;
  }
  const foundFileNode: FileTreeNode | undefined = foundFileNodes?.find(
    (val) => val.id === query.id
  );
  return foundFileNode;
}

function uploadFolder(node: FileTreeNode) {
  if (!node.error) {
    const queryFolderNode = getFileNodeByNode(node);
    if (queryFolderNode) {
      queryFolderNode.isUploaded = true;
      const uploadInfo: FilePath = {
        pathComponents: [...node.parents],
      };
      uploadInfo.pathComponents.push(node.label);
      const config = {
        headers: {
          "content-type": "application/json",
        },
      };
      axios
        .post(
          "/api/folders/experiments/" + props.id,
          JSON.stringify(uploadInfo),
          config
        )
        .catch((error) => {
          queryFolderNode.error = error.response.data;
        })
        .finally(() => {
          if (queryFolderNode) {
            queryFolderNode.isUploaded = false;
          }
        });
    }
  }
}

function uploadFile(file: File, node: FileTreeNode) {
  if (!node.error) {
    const queryFileNode = getFileNodeByNode(node);
    if (queryFileNode) {
      queryFileNode.isUploaded = true;
      const formData = new FormData();
      formData.append("file", file);
      const uploadInfo: FilePath = {
        pathComponents: [...node.parents],
      };
      uploadInfo.pathComponents.push(file.name);
      formData.append("form", JSON.stringify(uploadInfo));
      const config = {
        headers: {
          "content-type": "multipart/form-data",
        },
      };
      axios
        .post("/api/files/experiments/" + props.id, formData, config)
        .catch((error) => {
          queryFileNode.error = error.response.data;
        })
        .finally(() => {
          if (queryFileNode) {
            queryFileNode.isUploaded = false;
          }
        });
    }
  }
}

function deletePath(node: FileTreeNode) {
  if (!node.error) {
    const pathComponents = [...node.parents];
    pathComponents.push(node.label);
    const pathUpload: FilePath = {
      pathComponents: pathComponents,
    };
    axios
      .delete("/api/files/experiments/" + props.id, {
        headers: {
          "content-type": "application/json",
        },
        data: JSON.stringify(pathUpload),
      })
      .catch((error) => {
        serverError.value = error.response.data;
        showServerError.value = true;
      });
  }
}

function deleteAll() {
  const pathUpload: FilePath = {
    pathComponents: [],
  };
  axios
    .delete("/api/files/experiments/" + props.id, {
      headers: {
        "content-type": "application/json",
      },
      data: JSON.stringify(pathUpload),
    })
    .catch((error) => {
      serverError.value = error.response.data;
      showServerError.value = true;
    });
}

function submitExperiment() {
  if (!isSubmitting.value) {
    isSubmitting.value = true;
    axios
      .post("/api/experiments/" + props.id)
      .then(() => {
        status.value = ExperimentExecutionStatus.Waiting;
      })
      .catch((error) => {
        serverError.value = error.response.data;
        showServerError.value = true;
      })
      .finally(() => {
        isSubmitting.value = false;
      });
  }
}

function aborExperiment() {
  if (!isAborting.value) {
    isAborting.value = true;
    axios
      .post("/api/experiments/" + props.id + "/abort")
      .then(() => {
        status.value = ExperimentExecutionStatus.Aborted;
      })
      .catch((error) => {
        serverError.value = error.response.data;
        showServerError.value = true;
      })
      .finally(() => {
        isAborting.value = false;
      });
  }
}

/**
 * Conitinuesly polls changes from the server.
 */
function pollStatusChanges() {
  if (
    !isPolling.value &&
    !loadingError.value &&
    !pollingError.value &&
    // Stop polling if the route changes.
    router.currentRoute.value.fullPath === this_route
  ) {
    pollingError.value = null;
    axios
      .get("/api/experiments/" + props.id + "/status")
      .then((response) => {
        status.value = response.data;
        setTimeout(pollStatusChanges, POLLING_INTERVALL_MILLISECONDS);
      })
      .catch((error) => {
        pollingError.value = error.response.data;
      })
      .finally(() => {
        isPolling.value = false;
      });
  }
}
</script>
