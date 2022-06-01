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
      :error="!!uploadErrorMessage"
      :error-message="uploadErrorMessage"
      max-files="1"
    >
      <template v-slot:prepend>
        <q-icon name="cloud_upload" />
      </template>
    </q-file>
  </div>
</template>

<script setup lang="ts">
import { type Ref, ref } from "vue";

const sampleName = ref("");
const mail = ref("");
const comment = ref("");
const sample: Ref<File | null> = ref(null);
const isUploadingSample = ref(false);
const uploadErrorMessage = ref("");
const pipeline = ref("");
const pipelineOptions = ref(["ATACseq", "ChIPseq", "RNAseq"]);
</script>
