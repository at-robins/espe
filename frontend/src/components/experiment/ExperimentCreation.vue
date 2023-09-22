<template>
  <q-input
    outlined
    bottom-slots
    v-model="name"
    label="Experiment name"
    counter
    maxlength="512"
    :readonly="isUploadingExperiment"
    :error="!!uploadErrorMessage"
    :error-message="uploadErrorMessage"
    @keydown.enter="uploadExperiment"
  >
    <template v-slot:before>
      <q-icon :name="symOutlinedScience" color="primary" />
    </template>

    <template v-slot:append>
      <q-icon
        v-if="!!name"
        name="close"
        @click="name = ''"
        class="cursor-pointer"
      />
    </template>

    <template v-slot:hint>
      Specifiy a name to create a new experiment.
    </template>

    <template v-slot:after>
      <q-btn
        round
        color="primary"
        icon="add"
        :disable="isUploadingExperiment || !name"
        @click="uploadExperiment"
      />
    </template>
  </q-input>
</template>

<script setup lang="ts">
import axios from "axios";
import { ref, type Ref } from "vue";
import { useRouter } from "vue-router";
import { symOutlinedScience } from "@quasar/extras/material-symbols-outlined";

const router = useRouter();
const name: Ref<string> = ref("");
const isUploadingExperiment = ref(false);
const uploadErrorMessage = ref("");

function uploadExperiment() {
  if (name.value) {
    isUploadingExperiment.value = true;
    uploadErrorMessage.value = "";
    const formData = JSON.stringify(name.value);
    const config = {
      headers: {
        "Content-Type": "application/json",
      },
    };
    axios
      .post("/api/experiments", formData, config)
      .then((response) => {
        return router.push({
          name: "experiments_detail",
          params: { id: response.data },
        });
      })
      .catch((error) => {
        uploadErrorMessage.value = error;
      })
      .finally(() => {
        name.value = "";
        isUploadingExperiment.value = false;
      });
  } else {
    uploadErrorMessage.value = "A valid name is required.";
  }
}
</script>
<style scoped lang="scss"></style>
