<template>
  <div class="q-pa-md gutter-md">
    <q-input
      v-model="sampleName"
      label="Sample name"
      hint="Add the sample's name."
      :disable="isUploadingSample"
      counter
      :maxlength="128"
    />
    <q-input
      v-model="mail"
      label="E-mail (optional)"
      hint="E-mail address for notifications."
      :disable="isUploadingSample"
      counter
      :maxlength="128"
    />
    <q-input
      v-model="comment"
      label="Description (optional)"
      hint="Add a comment or description for this sample."
      :disable="isUploadingSample"
      counter
      :maxlength="256"
    />
    <q-select
      v-model="pipeline"
      :options="pipelineOptions"
      :disable="isUploadingSample"
      :loading="isLoadingPipelineBlueprints"
      option-value="id"
      option-label="name"
      label="Pipeline"
      hint="Select the pipeline to process the sample with."
      :error="!!loadPipelineBlueprintsError"
      :error-message="getLoadPipelineBlueprintsErrorMessage()"
    ></q-select>
    <q-file
      color="teal"
      v-model="sample"
      label="Upload sample"
      hint="Upload the sample data in the FASTQ format."
      :loading="isUploadingSample"
      :readonly="isUploadingSample"
      max-files="1"
    >
      <template v-slot:prepend>
        <q-icon name="cloud_upload" />
      </template>
    </q-file>
    <q-btn
      color="primary"
      label="Submit"
      :disable="!isUploadAllowed()"
      @click="uploadSample"
      class="q-mt-md"
    />
    <q-dialog v-model="showUploadError" v-if="uploadError">
      <error-popup :error-response="uploadError" />
    </q-dialog>
  </div>
</template>

<script setup lang="ts">
import ErrorPopup from "@/components/ErrorPopup.vue";
import type {
  ErrorResponse,
  ExperimentUpload,
  PipelineBlueprintDetail,
} from "@/scripts/types";
import axios from "axios";
import { type Ref, ref, onMounted } from "vue";

const isLoadingPipelineBlueprints = ref(false);
const sampleName = ref("");
const mail = ref("");
const comment = ref("");
const sample: Ref<File | null> = ref(null);
const isUploadingSample = ref(false);
const uploadError: Ref<ErrorResponse | null> = ref(null);
const loadPipelineBlueprintsError: Ref<ErrorResponse | null> = ref(null);
const pipeline: Ref<PipelineBlueprintDetail | null> = ref(null);
const pipelineOptions: Ref<Array<PipelineBlueprintDetail>> = ref([]);
const showUploadError = ref(false);

onMounted(() => {
  loadPipelineBlueprints();
});

function isUploadAllowed(): boolean {
  return (
    !!sample.value &&
    !!sampleName.value &&
    !!pipeline.value &&
    !isUploadingSample.value &&
    !isLoadingPipelineBlueprints.value
  );
}

function getLoadPipelineBlueprintsErrorMessage(): string {
  return loadPipelineBlueprintsError.value
    ? "The pipelines could not be loaded due to an error with ID: " +
        loadPipelineBlueprintsError.value.uuid
    : "";
}

function uploadSample() {
  if (isUploadAllowed() && sample.value && pipeline.value) {
    isUploadingSample.value = true;
    uploadError.value = null;
    const formData = new FormData();
    formData.append("file", sample.value);
    const uploadInfo: ExperimentUpload = {
      name: sampleName.value,
      mail: mail.value,
      comment: comment.value,
      pipelineId: pipeline.value.id,
    };
    formData.append("form", JSON.stringify(uploadInfo));
    const config = {
      headers: {
        "content-type": "multipart/form-data",
      },
    };
    axios
      .post("/api/experiment", formData, config)
      .catch((error) => {
        showUploadError.value = true;
        uploadError.value = error.response.data;
      })
      .finally(() => {
        sample.value = null;
        isUploadingSample.value = false;
      });
  }
}

/**
 * Initial loading of all available pipelines.
 */
function loadPipelineBlueprints() {
  isLoadingPipelineBlueprints.value = true;
  loadPipelineBlueprintsError.value = null;
  axios
    .get("/api/pipeline/blueprint")
    .then((response) => {
      pipelineOptions.value = response.data;
    })
    .catch((error) => {
      pipelineOptions.value = [];
      loadPipelineBlueprintsError.value = error.response.data;
    })
    .finally(() => {
      isLoadingPipelineBlueprints.value = false;
    });
}
</script>
