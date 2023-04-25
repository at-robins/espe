<template>
  <div class="q-pa-md gutter-md">
    <q-form @submit="uploadSample" :greedy="true">
      <q-input
        v-model="sampleName"
        for="sample-submission-input-name"
        label="Sample name"
        hint="Add the sample's name."
        :disable="isUploadingSample"
        counter
        :maxlength="MAX_SAMPLE_NAME_LENGTH"
        lazy-rules
        :rules="[
          (val: string | null) => (val && val.length > 0) || 'A sample name is required.',
        ]"
      />
      <q-input
        v-model="mail"
        for="sample-submission-input-mail"
        label="E-mail (optional)"
        hint="E-mail address for notifications."
        :disable="isUploadingSample"
        counter
        lazy-rules
        :rules="[validateEMail]"
        :maxlength="MAX_SAMPLE_NAME_LENGTH"
      />
      <q-input
        v-model="comment"
        for="sample-submission-input-comment"
        label="Description (optional)"
        hint="Add a comment or description for this sample."
        :disable="isUploadingSample"
        counter
        :maxlength="MAX_SAMPLE_NAME_LENGTH * 2"
      />
      <q-select
        v-model="pipeline"
        for="sample-submission-input-pipeline"
        :options="pipelineOptions"
        :disable="isUploadingSample"
        :loading="isLoadingPipelineBlueprints"
        option-value="id"
        option-label="name"
        label="Pipeline"
        hint="Select the pipeline to process the sample with."
        :error="!!loadPipelineBlueprintsError"
        :error-message="getLoadPipelineBlueprintsErrorMessage()"
        lazy-rules
        :rules="[(val: PipelineBlueprintDetail | null) => !!val || 'A pipeline needs to be selected.']"
      ></q-select>
      <q-file
        for="sample-submission-input-file"
        color="teal"
        v-model="sample"
        label="Upload sample"
        hint="Upload the sample data in the FASTQ format."
        :loading="isUploadingSample"
        :readonly="isUploadingSample"
        max-files="1"
        accept=".fastq.gz"
        lazy-rules
        :rules="[
          (val: File | null) =>
            (val && val.name.endsWith('.fastq.gz')) ||
            'A compressed FASTQ file is required (.fastq.gz).',
        ]"
        @rejected="sample = null"
        @update:model-value="setSampleNameIfNone"
      >
        <template v-slot:prepend>
          <q-icon name="cloud_upload" />
        </template>
      </q-file>
      <q-btn
        color="primary"
        label="Submit"
        type="submit"
        class="q-mt-md"
        :loading="isUploadingSample"
      />
    </q-form>
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

const MAX_SAMPLE_NAME_LENGTH = 128;

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

/**
 * Validates the user specified e-mail address.
 */
function validateEMail(value: string | null | undefined): boolean | string {
  if (!value) {
    return true;
  } else {
    const validationPattern =
      // eslint-disable-next-line no-control-regex
      /(?:[a-z0-9!#$%&'*+/=?^_`{|}~-]+(?:\.[a-z0-9!#$%&'*+/=?^_`{|}~-]+)*|"(?:[\x01-\x08\x0b\x0c\x0e-\x1f\x21\x23-\x5b\x5d-\x7f]|\\[\x01-\x09\x0b\x0c\x0e-\x7f])*")@(?:(?:[a-z0-9](?:[a-z0-9-]*[a-z0-9])?\.)+[a-z0-9](?:[a-z0-9-]*[a-z0-9])?|\[(?:(?:(2(5[0-5]|[0-4][0-9])|1[0-9][0-9]|[1-9]?[0-9]))\.){3}(?:(2(5[0-5]|[0-4][0-9])|1[0-9][0-9]|[1-9]?[0-9])|[a-z0-9-]*[a-z0-9]:(?:[\x01-\x08\x0b\x0c\x0e-\x1f\x21-\x5a\x53-\x7f]|\\[\x01-\x09\x0b\x0c\x0e-\x7f])+)\])/;
    return validationPattern.test(value) || "The specified e-mail is invalid.";
  }
}

/**
 * Checks if the upload is currently allowed.
 */
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
      mail: mail.value ? mail.value : undefined,
      comment: comment.value ? comment.value : undefined,
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
      loadPipelineBlueprintsError.value = error?.response?.data;
    })
    .finally(() => {
      isLoadingPipelineBlueprints.value = false;
    });
}

/**
 * Sets the sample name to the file name without extension
 * if no sample name was specified previously.
 */
function setSampleNameIfNone(file: File | undefined | null) {
  if (file && !sampleName.value) {
    let fileNameWithoutExtension = file.name.split(".")[0];
    if (fileNameWithoutExtension.length > MAX_SAMPLE_NAME_LENGTH) {
      fileNameWithoutExtension = fileNameWithoutExtension.substring(
        0,
        MAX_SAMPLE_NAME_LENGTH
      );
    }
    sampleName.value = fileNameWithoutExtension;
  }
}
</script>
