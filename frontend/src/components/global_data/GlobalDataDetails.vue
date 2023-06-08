<template>
  <div class="q-pa-md">
    <q-card>
      <div v-if="!loadingError">
        <div v-if="isLoadingGlobalDataDetails || !globalData">
          <q-card-section class="flex-center row">
            <q-spinner color="primary" size="xl" />
          </q-card-section>
        </div>
        <div v-else>
          <q-card-section>
            <global-data-title
              :title="globalData.name"
              :global-data-id="globalData.id"
              @update:title="updateTitle"
            />
          </q-card-section>
          <q-card-section>
            <global-data-comment
              :comment="globalData.comment"
              :global-data-id="globalData.id"
            />
          </q-card-section>
          <q-card-section>
            <div class="row">
              <div class="q-ma-md">
                <q-btn :icon="symOutlinedAccountTree" round>
                  <q-tooltip>Upload files to the global repository.</q-tooltip>
                </q-btn>
              </div>
              <q-separator vertical class="q-ml-md q-mr-md" />
              <file-tree
                v-model="fileNodes"
                :base-directory-label="globalData ? globalData.name : 'Root'"
                @added-file="uploadFile"
                @added-folder="uploadFolder"
                @deleted-path="deletePath"
                @deleted-all="deleteAll"
                class="col"
              />
            </div>
          </q-card-section>
        </div>
      </div>
      <div v-else>
        <error-popup :error-response="loadingError" />
      </div>
    </q-card>
    <q-dialog v-model="showDeletionError" v-if="deletionError">
      <error-popup :error-response="deletionError" />
    </q-dialog>
  </div>
</template>

<script setup lang="ts">
import {
  type ErrorResponse,
  type GlobalDataDetails,
  type GlobalDataFileDetails,
  type FileTreeNode,
  type GlobalDataFilePath,
} from "@/scripts/types";
import axios from "axios";
import { ref, onMounted, type Ref } from "vue";
import ErrorPopup from "../ErrorPopup.vue";
import FileTree from "../FileTree.vue";
import GlobalDataTitle from "./GlobalDataTitle.vue";
import GlobalDataComment from "./GlobalDataComment.vue";
import { symOutlinedAccountTree } from "@quasar/extras/material-symbols-outlined";

const files: Ref<Array<GlobalDataFileDetails>> = ref([]);
const globalData: Ref<GlobalDataDetails | null> = ref(null);
const fileNodes: Ref<Array<FileTreeNode>> = ref([]);
const isLoadingGlobalDataDetails = ref(false);
const loadingError: Ref<ErrorResponse | null> = ref(null);
const deletionError: Ref<ErrorResponse | null> = ref(null);
const showDeletionError = ref(false);

const props = defineProps({
  id: { type: String, required: true },
});

function updateTitle(title: string) {
  if (globalData.value) {
    globalData.value.name = title;
  }
}

function getFileTreeNodes(files: GlobalDataFileDetails[]): FileTreeNode[] {
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
  loadGlobalDataDetails();
});

/**
 * Initial loading of details from the server.
 */
function loadGlobalDataDetails() {
  isLoadingGlobalDataDetails.value = true;
  loadingError.value = null;
  axios
    .get("/api/globals/" + props.id)
    .then((response) => {
      globalData.value = response.data;
      return axios.get("/api/globals/" + props.id + "/files");
    })
    .then((response) => {
      files.value = response.data;
      fileNodes.value = getFileTreeNodes(files.value);
    })
    .catch((error) => {
      globalData.value = null;
      files.value = [];
      fileNodes.value = [];
      loadingError.value = error.response.data;
    })
    .finally(() => {
      isLoadingGlobalDataDetails.value = false;
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
      const uploadInfo: GlobalDataFilePath = {
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
          "/api/globals/" + props.id + "/folders",
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
      const uploadInfo: GlobalDataFilePath = {
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
        .post("/api/globals/" + props.id + "/files", formData, config)
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
    const pathUpload: GlobalDataFilePath = {
      pathComponents: pathComponents,
    };
    axios
      .delete("/api/globals/" + props.id + "/files", {
        headers: {
          "content-type": "application/json",
        },
        data: JSON.stringify(pathUpload),
      })
      .catch((error) => {
        deletionError.value = error.response.data;
        showDeletionError.value = true;
      });
  }
}

function deleteAll() {
  const pathUpload: GlobalDataFilePath = {
    pathComponents: [],
  };
  axios
    .delete("/api/globals/" + props.id + "/files", {
      headers: {
        "content-type": "application/json",
      },
      data: JSON.stringify(pathUpload),
    })
    .catch((error) => {
      deletionError.value = error.response.data;
      showDeletionError.value = true;
    });
}
</script>
