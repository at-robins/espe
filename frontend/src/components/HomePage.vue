<template>
  <div class="q-pa-md gutter-md">
    <q-input
      v-model="sampleName"
      label="Sample name"
      hint="Add the sample's name."
      counter
      :maxlength="128"
    />
    <q-input
      v-model="mail"
      label="E-mail (optional)"
      hint="E-mail address for notifications."
      counter
      :maxlength="128"
    />
    <q-input
      v-model="comment"
      label="Description (optional)"
      hint="Add a comment or description for this sample."
      counter
      :maxlength="256"
    />
    <q-select
      v-model="pipeline"
      :options="pipelineOptions"
      label="Pipeline"
      hint="Select the pipeline to process the sample with."
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
      @click="uploadSample"
      class="q-mt-md"
    />
    <q-dialog v-model="showError" v-if="uploadError">
      <error-popup :error-response="uploadError" />
    </q-dialog>
  </div>
</template>

<script setup lang="ts">
import ErrorPopup from "@/components/ErrorPopup.vue";
import type { ErrorResponse, ExperimentUpload } from "@/scripts/types";
import axios from "axios";
import { type Ref, ref } from "vue";

const sampleName = ref("");
const mail = ref("");
const comment = ref("");
const sample: Ref<File | null> = ref(null);
const isUploadingSample = ref(false);
const uploadError: Ref<ErrorResponse | null> = ref(null);
const pipeline = ref(0);
const pipelineOptions = ref(["ATACseq", "ChIPseq", "RNAseq"]);
const showError = ref(false);

function uploadSample() {
  if (sample.value && sampleName.value) {
    isUploadingSample.value = true;
    uploadError.value = null;
    const formData = new FormData();
    formData.append("file", sample.value);
    const uploadInfo: ExperimentUpload = {
      name: sampleName.value,
      mail: mail.value,
      comment: comment.value,
      pipelineId: pipeline.value,
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
        showError.value = true;
        uploadError.value = error.response.data;
      })
      .finally(() => {
        sample.value = null;
        isUploadingSample.value = false;
      });
  }
}
</script>
