<template>
  <div class="row no-wrap">
    <div class="q-ma-md">
      <q-btn :icon="symOutlinedRoute" :loading="isUpdatingPipeline" round>
        <q-tooltip>Select a pipeline.</q-tooltip>
      </q-btn>
    </div>
    <q-separator vertical class="q-ml-md q-mr-md" />

    <div class="col">
      <div class="q-ma-md row">
        <q-select
          v-model="pipelineModel"
          :options="pipelineOptions"
          :readonly="isUpdatingPipeline || locked"
          option-label="name"
          option-value="id"
          clearable
          @update:model-value="updatePipeline"
          :loading="isLoading || isUpdatingPipeline"
          label="Select a pipeline"
          :error="!!pipelineErrorMessage"
          :error-message="pipelineErrorMessage"
          bottom-slots
          class="col"
        >
          <q-tooltip v-if="locked">
            Altering the pipeline is not possible while the experiment is
            executed or output information is downloaded.
          </q-tooltip>
        </q-select>
        <div class="q-ml-md q-mt-md">
          <q-btn
            round
            outline
            :icon="symOutlinedSync"
            @click="updateAllPipelines"
            color="primary"
            :loading="isUpdatingAllPipelines"
          >
            <q-tooltip> Reload all pipelines. </q-tooltip>
          </q-btn>
        </div>
      </div>
      <div v-if="pipelineModel?.description" class="q-ma-md">
        <span v-html="pipelineModel.description" />
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import axios from "axios";
import { ref, type PropType, type Ref, onMounted } from "vue";
import {
  symOutlinedRoute,
  symOutlinedSync,
} from "@quasar/extras/material-symbols-outlined";
import { type ErrorResponse } from "@/scripts/types";
import { type PipelineBlueprint } from "@/scripts/pipeline-blueprint";
import { error_to_string } from "@/scripts/utilities";

const props = defineProps({
  pipelineId: {
    type: String as PropType<string | null>,
    default: null,
    required: false,
  },
  entityId: { type: Number, required: true },
  endpointType: { type: String, required: true },
  locked: { type: Boolean, required: false, default: false },
});

const isLoading = ref(false);
const isUpdatingAllPipelines = ref(false);
const isUpdatingPipeline = ref(false);
const pipelineErrorMessage = ref("");
const pipelineModel: Ref<PipelineBlueprint | null> = ref(null);
const pipelineOptions: Ref<Array<PipelineBlueprint>> = ref([]);

const emit = defineEmits<{
  (event: "update:selectedPipeline", value: PipelineBlueprint | null): void;
}>();

onMounted(() => {
  loadPipelineDetails();
});

/**
 * Initial loading of details from the server.
 */
function loadPipelineDetails() {
  isLoading.value = true;
  pipelineErrorMessage.value = "";
  axios
    .get("/api/" + props.endpointType + "/" + props.entityId + "/pipelines")
    .then((response) => {
      pipelineOptions.value = response.data;
      const selectedPipeline = pipelineOptions.value.find(
        (pipeline) => pipeline.id === props.pipelineId
      );
      pipelineModel.value = selectedPipeline ? selectedPipeline : null;
    })
    .then(() => {
      emit("update:selectedPipeline", pipelineModel.value);
    })
    .catch((error: ErrorResponse) => {
      pipelineOptions.value = [];
      pipelineErrorMessage.value = error.message;
    })
    .finally(() => {
      isLoading.value = false;
    });
}

function updatePipeline() {
  if (!isUpdatingPipeline.value) {
    isUpdatingPipeline.value = true;
    pipelineErrorMessage.value = "";
    const formData = JSON.stringify(
      pipelineModel.value ? pipelineModel.value.id : null
    );
    const config = {
      headers: {
        "Content-Type": "application/json",
      },
    };
    axios
      .patch(
        "/api/" + props.endpointType + "/" + props.entityId + "/pipeline",
        formData,
        config
      )
      .then(() => {
        emit("update:selectedPipeline", pipelineModel.value);
      })
      .catch((error) => {
        pipelineErrorMessage.value = error_to_string(error.response.data);
      })
      .finally(() => {
        isUpdatingPipeline.value = false;
      });
  }
}

function updateAllPipelines() {
  if (!isUpdatingAllPipelines.value) {
    isUpdatingAllPipelines.value = true;
    pipelineErrorMessage.value = "";
    axios
      .patch("/api/pipelines/blueprints")
      .then(() => {
        loadPipelineDetails();
      })
      .catch((error) => {
        pipelineErrorMessage.value = error_to_string(error.response.data);
      })
      .finally(() => {
        isUpdatingAllPipelines.value = false;
      });
  }
}
</script>
<style scoped lang="scss"></style>
