<template>
  <div class="q-pa-md q-gutter-md">
    <q-table
      title="Global data repositories"
      :rows="loadedGlobalRepos"
      :columns="columns"
      row-key="id"
      :loading="isLoadingGlobals"
      no-data-label="No repositories have been declared yet."
      :rows-per-page-options="[10, 15, 20, 25, 50, 0]"
    >
      <template v-slot:body-cell="props">
        <q-td :props="props" @click="navigateToRepository(props.row.id)">
          {{ props.value }}
        </q-td>
        <q-menu touch-position context-menu>
          <q-list dense>
            <q-item
              clickable
              v-close-popup
              @click="navigateToRepository(props.row.id)"
            >
              <q-item-section
                ><div class="row flex-center no-wrap">
                  <q-icon
                    :name="matFileOpen"
                    color="primary"
                    size="xs"
                    class="q-mr-xs"
                  />
                  <div class="col-8">Open</div>
                </div></q-item-section
              >
            </q-item>
            <q-item
              clickable
              v-close-popup
              @click="deleteGlobalDataDetail(props.row.id)"
            >
              <q-item-section>
                <div class="row flex-center no-wrap">
                  <q-icon
                    :name="matDelete"
                    color="negative"
                    size="xs"
                    class="q-mr-xs"
                  />
                  <div class="col-8">Delete</div>
                </div>
              </q-item-section>
            </q-item>
          </q-list>
        </q-menu>
      </template>
    </q-table>
    <q-dialog v-model="showLoadingError" v-if="loadingError">
      <error-popup :error-response="loadingError" />
    </q-dialog>
  </div>
</template>

<script setup lang="ts">
import { type ErrorResponse, type GlobalDataDetails } from "@/scripts/types";
import axios from "axios";
import { ref, onMounted, type Ref } from "vue";
import ErrorPopup from "./ErrorPopup.vue";
import { DateTime } from "luxon";
import { type QTableProps } from "quasar";
import { matDelete, matFileOpen } from "@quasar/extras/material-icons";
import { useRouter } from "vue-router";

const loadedGlobalRepos: Ref<Array<GlobalDataDetails>> = ref([]);
const isLoadingGlobals = ref(false);
const loadingError: Ref<ErrorResponse | null> = ref(null);
const showLoadingError = ref(false);
const router = useRouter();

const columns: QTableProps["columns"] = [
  {
    name: "name",
    required: true,
    label: "Repository name",
    align: "left",
    field: (row: GlobalDataDetails) => row.name,
    sortable: true,
  },
  {
    name: "date",
    label: "Created",
    field: (row: GlobalDataDetails) => DateTime.fromISO(row.creationTime),
    sortable: true,
    format: (val: DateTime) => val.toLocaleString(DateTime.DATETIME_FULL),
  },
];

onMounted(() => {
  loadGlobalDataDetails();
});

/**
 * Initial loading of details from the server.
 */
function loadGlobalDataDetails() {
  isLoadingGlobals.value = true;
  loadingError.value = null;
  axios
    .get("/api/globals")
    .then((response) => {
      loadedGlobalRepos.value = response.data;
    })
    .catch((error) => {
      loadedGlobalRepos.value = [];
      loadingError.value = error.response.data;
      showLoadingError.value = true;
    })
    .finally(() => {
      isLoadingGlobals.value = false;
    });
}

/**
 * Navigates to the respective repository details page.
 */
function navigateToRepository(id: number) {
  router.push({ name: "globals_detail", params: { id: id } });
}

/**
 * Deletes the repository from the server.
 */
function deleteGlobalDataDetail(id: number) {
  axios
    .delete("/api/globals/" + id)
    .then(() => {
      loadedGlobalRepos.value = loadedGlobalRepos.value.filter(
        (repo) => repo.id !== id
      );
    })
    .catch((error) => {
      loadedGlobalRepos.value = [];
      loadingError.value = error.response.data;
      showLoadingError.value = true;
    })
    .finally(() => {
      isLoadingGlobals.value = false;
    });
}
</script>
