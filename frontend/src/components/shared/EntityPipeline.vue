<template>
  <div class="row no-wrap">
    <div class="q-ma-md">
      <q-btn :icon="symOutlinedRoute" :loading="isUpdatingPipeline" round>
        <q-tooltip>Select a pipeline.</q-tooltip>
      </q-btn>
    </div>
    <q-separator vertical class="q-ml-md q-mr-md" />

    <div class="col">
      <div class="q-ma-md">
        <q-select
          v-model="pipelineModel"
          :options="pipelineOptions"
          option-label="name"
          option-value="id"
          clearable
          @update:model-value="updatePipeline"
          :loading="isLoading || isUpdatingPipeline"
          label="Select a pipeline"
          :error="!!pipelineErrorMessage"
          :error-message="pipelineErrorMessage"
          bottom-slots
        />
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import axios from "axios";
import { ref, type PropType, type Ref, onMounted } from "vue";
import { symOutlinedRoute } from "@quasar/extras/material-symbols-outlined";
import {
  type ErrorResponse,
  type PipelineBlueprintDetail,
} from "@/scripts/types";
import { error_to_string } from "@/scripts/utilities";

const props = defineProps({
  pipelineId: {
    type: String as PropType<string | null>,
    default: null,
    required: false,
  },
  entityId: { type: Number, required: true },
  endpointType: { type: String, required: true },
});

const isLoading = ref(false);
const isUpdatingPipeline = ref(false);
const pipelineErrorMessage = ref("");
const pipelineModel: Ref<PipelineBlueprintDetail | null> = ref(null);
const pipelineOptions: Ref<Array<PipelineBlueprintDetail>> = ref([]);

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
    .get("/api/pipelines/blueprints")
    .then((response) => {
      pipelineOptions.value = response.data;
      const selectedPipeline = pipelineOptions.value.find(
        (pipeline) => pipeline.id === props.pipelineId
      );
      pipelineModel.value = selectedPipeline ? selectedPipeline : null;
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
      .catch((error) => {
        pipelineErrorMessage.value = error_to_string(error.response.data);
      })
      .finally(() => {
        isUpdatingPipeline.value = false;
      });
  }
}
</script>
<style scoped lang="scss"></style>
