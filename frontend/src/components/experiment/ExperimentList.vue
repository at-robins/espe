<template>
  <div class="q-pa-md q-gutter-md">
    <q-table
      title="Existing experiments"
      :rows="loadedExperiments"
      :columns="columns"
      row-key="id"
      :loading="isLoading"
      no-data-label="No experiments have been declared yet."
      :rows-per-page-options="[10, 15, 20, 25, 50, 0]"
      flat
    >
      <template v-slot:body-cell="props">
        <q-td
          :props="props"
          @click="navigateToRepository(props.row.id)"
          class="cursor-pointer"
        >
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
              @click="deleteExperimentDetail(props.row.id)"
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
import { type ErrorResponse, type ExperimentDetails } from "@/scripts/types";
import axios from "axios";
import { ref, onMounted, type Ref } from "vue";
import ErrorPopup from "../ErrorPopup.vue";
import { DateTime } from "luxon";
import { type QTableProps } from "quasar";
import { matDelete, matFileOpen } from "@quasar/extras/material-icons";
import { useRouter } from "vue-router";

const loadedExperiments: Ref<Array<ExperimentDetails>> = ref([]);
const isLoading = ref(false);
const loadingError: Ref<ErrorResponse | null> = ref(null);
const showLoadingError = ref(false);
const router = useRouter();

const columns: QTableProps["columns"] = [
  {
    name: "name",
    required: true,
    label: "Experiment name",
    align: "left",
    field: (row: ExperimentDetails) => row.name,
    sortable: true,
  },
  {
    name: "pipeline",
    label: "Pipeline ID",
    align: "left",
    field: (row: ExperimentDetails) => row.pipelineId,
    sortable: true,
  },
  {
    name: "mail",
    label: "E-mail",
    align: "left",
    field: (row: ExperimentDetails) => row.mail,
    sortable: true,
  },
  {
    name: "date",
    label: "Created",
    field: (row: ExperimentDetails) => DateTime.fromISO(row.creationTime),
    sortable: true,
    format: (val: DateTime) => val.toLocaleString(DateTime.DATETIME_FULL),
  },
];

onMounted(() => {
  loadExperimentDetails();
});

/**
 * Initial loading of details from the server.
 */
function loadExperimentDetails() {
  isLoading.value = true;
  loadingError.value = null;
  axios
    .get("/api/experiments")
    .then((response) => {
      loadedExperiments.value = response.data;
    })
    .catch((error) => {
      loadedExperiments.value = [];
      loadingError.value = error.response.data;
      showLoadingError.value = true;
    })
    .finally(() => {
      isLoading.value = false;
    });
}

/**
 * Navigates to the respective experiment details page.
 */
function navigateToRepository(id: number) {
  router.push({ name: "experiments_detail", params: { id: id } });
}

/**
 * Deletes the experiment from the server.
 */
function deleteExperimentDetail(id: number) {
  axios
    .delete("/api/experiments/" + id)
    .then(() => {
      loadedExperiments.value = loadedExperiments.value.filter(
        (repo) => repo.id !== id
      );
    })
    .catch((error) => {
      loadedExperiments.value = [];
      loadingError.value = error.response.data;
      showLoadingError.value = true;
    })
    .finally(() => {
      isLoading.value = false;
    });
}
</script>
