<template>
  <div class="q-pa-md">
    <q-card>
      <div v-if="!loadingError">
        <div v-if="isLoadingGlobalDataDetails">
          <q-card-section class="flex-center row">
            <q-spinner color="primary" size="xl" />
          </q-card-section>
        </div>
        <div v-else>
          <q-card-section>
            <div class="text-h6">{{ globalData?.name }}</div>
          </q-card-section>
          <q-card-section>
            <div>{{ globalData?.comment }}</div>
          </q-card-section>
          <q-card-section>
            <file-tree
              v-model="fileNodes"
              :base-directory-label="globalData ? globalData.name : 'Root'"
            />
          </q-card-section>
        </div>
      </div>
      <div v-else>
        <error-popup :error-response="loadingError" />
      </div>
    </q-card>
  </div>
</template>

<script setup lang="ts">
import {
  type ErrorResponse,
  type GlobalDataDetails,
  type GlobalDataFileDetails,
  type FileTreeNode,
} from "@/scripts/types";
import axios from "axios";
import { ref, onMounted, type Ref } from "vue";
import ErrorPopup from "./ErrorPopup.vue";
import FileTree from "./FileTree.vue";

const files: Ref<Array<GlobalDataFileDetails>> = ref([]);
const globalData: Ref<GlobalDataDetails | null> = ref(null);
const fileNodes: Ref<Array<FileTreeNode>> = ref([]);
const isLoadingGlobalDataDetails = ref(false);
const loadingError: Ref<ErrorResponse | null> = ref(null);

const props = defineProps({
  id: { type: String, required: true },
});

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
</script>
