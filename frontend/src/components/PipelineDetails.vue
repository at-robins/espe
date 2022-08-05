<template>
  <div class="q-pa-md q-gutter-md">
    <q-card>
      <q-card-section>
        <div class="text-h6">Pipeline {{ id }}</div>
      </q-card-section>
      <div class="q-pa-md gutter-md wrap row" v-if="!loadingError">
        <div
          v-for="(detail, index) in details"
          :key="detail.id"
          class="q-pb-md"
        >
          <q-btn rounded no-caps @click="selectStep(detail)">
            <q-spinner
              v-if="isStepRunning(detail)"
              :color="getDetailIconColour(detail)"
              class="on-left"
            ></q-spinner>
            <q-icon
              v-else
              :name="getDetailIcon(detail)"
              :color="getDetailIconColour(detail)"
              left
            ></q-icon>
            <div class="text-center">{{ detail.name }}</div>
          </q-btn>
          <q-icon
            v-if="index < details.length - 1"
            name="trending_flat"
          ></q-icon>
        </div>
      </div>
      <div v-else>
        <error-popup :error-response="loadingError" />
      </div>
    </q-card>
    <q-card>
      <q-card-section>
        <div v-if="selectedStep === null" class="text-h6">
          Select a step to display further information.
        </div>
        <div v-else class="text-h6">
          Step {{ selectedStep.id }} - {{ selectedStep.name }}
        </div>
      </q-card-section>
      <div v-if="selectedStep !== null" class="q-gutter-md q-pa-md col">
        <q-btn label="Display logs" class="row" />
        <q-btn label="Download output" class="row" />
      </div>
    </q-card>
    <q-dialog v-model="showPollingError" v-if="pollingError">
      <error-popup :error-response="pollingError" />
    </q-dialog>
  </div>
</template>

<script setup lang="ts">
import {
  PipelineStepStatus,
  type ErrorResponse,
  type PipelineStepDetail,
} from "@/scripts/types";
import axios from "axios";
import { ref, onMounted, type Ref } from "vue";
import ErrorPopup from "./ErrorPopup.vue";

// The intervall in which pipeline updates are requested from the server.
const POLLING_INTERVALL_MILLISECONDS = 30000;

const details: Ref<Array<PipelineStepDetail>> = ref([]);
const isLoadingPipelineDetails = ref(false);
const loadingError: Ref<ErrorResponse | null> = ref(null);
const isPollingPipelineDetails = ref(false);
const pollingError: Ref<ErrorResponse | null> = ref(null);
const selectedStep: Ref<PipelineStepDetail | null> = ref(null);
const showPollingError = ref(false);

const props = defineProps({
  id: { type: String, required: true },
});

onMounted(() => {
  loadPipelineDetails();
});

/**
 * Initial loading of details from the server.
 */
function loadPipelineDetails() {
  isLoadingPipelineDetails.value = true;
  loadingError.value = null;
  axios
    .get("/api/pipeline/instance/" + props.id)
    .then((response) => {
      setPipelineDetails(response.data);
      setTimeout(pollDetailsChanges, POLLING_INTERVALL_MILLISECONDS);
    })
    .catch((error) => {
      details.value = [];
      loadingError.value = error.response.data;
    })
    .finally(() => {
      isLoadingPipelineDetails.value = false;
    });
}

/**
 * Conitinuesly polls changes from the server.
 */
function pollDetailsChanges() {
  if (
    !isPollingPipelineDetails.value &&
    !loadingError.value &&
    !pollingError.value
  ) {
    pollingError.value = null;
    axios
      .get("/api/pipeline/instance/" + props.id)
      .then((response) => {
        setPipelineDetails(response.data);
        setTimeout(pollDetailsChanges, POLLING_INTERVALL_MILLISECONDS);
      })
      .catch((error) => {
        showPollingError.value = true;
        pollingError.value = error.response.data;
      })
      .finally(() => {
        isPollingPipelineDetails.value = false;
      });
  }
}

function setPipelineDetails(response: PipelineStepDetail[]) {
  details.value = response;
  if (selectedStep.value) {
    const update_selected = response.find(
      (detail) => detail.id === selectedStep.value?.id
    );
    selectedStep.value = !update_selected ? null : update_selected;
  }
}

/**
 * Returns the icon corresponding to the pipeline step status.
 */
function getDetailIcon(detail: PipelineStepDetail): string {
  switch (detail.status) {
    case PipelineStepStatus.Success:
      return "check_circle_outline";
    case PipelineStepStatus.Failed:
      return "highlight_off";
    case PipelineStepStatus.Pending:
      return "update";
    case PipelineStepStatus.Running:
      return "update";
    default:
      throw new Error("Unknown status: " + detail.status);
  }
}

/**
 * Returns the icon colour corresponding to the pipeline step status.
 */
function getDetailIconColour(detail: PipelineStepDetail): string {
  switch (detail.status) {
    case PipelineStepStatus.Success:
      return "positive";
    case PipelineStepStatus.Failed:
      return "negative";
    case PipelineStepStatus.Pending:
      return "primary";
    case PipelineStepStatus.Running:
      return "primary";
    default:
      throw new Error("Unknown status: " + detail.status);
  }
}

/**
 * Convinience method to check if the pipeline step is currently being run by the server.
 */
function isStepRunning(detail: PipelineStepDetail): boolean {
  return detail.status === PipelineStepStatus.Running;
}

/**
 * Selects the specified pipeline step to display related information.
 */
function selectStep(detail: PipelineStepDetail) {
  selectedStep.value = detail;
}
</script>
