<template>
  <q-input
    outlined
    bottom-slots
    v-model="name"
    label="Global data repository name"
    counter
    maxlength="512"
    :readonly="isUploadingGlobalData"
    :error="!!uploadErrorMessage"
    :error-message="uploadErrorMessage"
    @keydown.enter="uploadGlobalData"
  >
    <template v-slot:before>
      <q-icon :name="matPublic" color="primary" />
    </template>

    <template v-slot:append>
      <q-icon
        v-if="!!name"
        :name="matClose"
        @click="name = ''"
        class="cursor-pointer"
      />
    </template>

    <template v-slot:hint>
      Specifiy a name to create a new global data repository.
    </template>

    <template v-slot:after>
      <q-btn
        round
        color="primary"
        :icon="matAdd"
        :disable="isUploadingGlobalData || !name"
        @click="uploadGlobalData"
      />
    </template>
  </q-input>
</template>

<script setup lang="ts">
import axios from "axios";
import { ref, type Ref } from "vue";
import { useRouter } from "vue-router";
import { matAdd, matClose, matPublic } from "@quasar/extras/material-icons";

const router = useRouter();
const name: Ref<string> = ref("");
const isUploadingGlobalData = ref(false);
const uploadErrorMessage = ref("");

function uploadGlobalData() {
  if (name.value) {
    isUploadingGlobalData.value = true;
    uploadErrorMessage.value = "";
    const formData = JSON.stringify(name.value);
    const config = {
      headers: {
        "Content-Type": "application/json",
      },
    };
    axios
      .post("/api/globals", formData, config)
      .then((response) => {
        return router.push({
          name: "globals_detail",
          params: { id: response.data },
        });
      })
      .catch((error) => {
        uploadErrorMessage.value = error;
      })
      .finally(() => {
        name.value = "";
        isUploadingGlobalData.value = false;
      });
  } else {
    uploadErrorMessage.value = "A valid name is required.";
  }
}
</script>
<style scoped lang="scss"></style>
