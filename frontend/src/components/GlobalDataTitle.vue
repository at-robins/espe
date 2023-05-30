<template>
  <div class="row no-wrap">
    <q-item-label v-if="!editTitleMode" class="text-h6 shrink q-ma-sm"
      >{{ titleModel }}
      <q-menu touch-position context-menu>
        <q-list dense>
          <q-item clickable v-close-popup @click="editTitleField">
            <q-item-section
              ><div class="row flex-center no-wrap">
                <q-icon
                  :name="matEdit"
                  color="primary"
                  size="xs"
                  class="q-mr-xs"
                />
                <div class="col">Edit...</div>
              </div></q-item-section
            >
          </q-item>
        </q-list>
      </q-menu>
    </q-item-label>
    <q-input
      v-else
      ref="titleInputRef"
      v-model="titleModel"
      @keydown.enter="updateTitle"
      @blur="updateTitle"
      :readonly="isUpdatingTitle"
      class="text-h6 col-grow q-ma-sm"
    />

    <q-tooltip
      v-if="updateTitleErrorMessage"
      class="bg-red"
      anchor="bottom start"
      self="top start"
    >
      {{ updateTitleErrorMessage }}
    </q-tooltip>
  </div>
</template>

<script setup lang="ts">
import axios from "axios";
import { nextTick, ref, type Ref } from "vue";
import { matEdit } from "@quasar/extras/material-icons";

const props = defineProps({
  title: { type: String, required: true },
  globalDataId: { type: Number, required: true },
});

const emit = defineEmits<{
  (event: "update:title", title: string): void;
}>();

const isUpdatingTitle = ref(false);
const updateTitleErrorMessage = ref("");
const editTitleMode = ref(false);
const titleModel = ref(props.title);
const titleInputRef: Ref<HTMLInputElement | null> = ref(null);

function updateTitle() {
  isUpdatingTitle.value = true;
  updateTitleErrorMessage.value = "";
  const formData = JSON.stringify(titleModel.value);
  const config = {
    headers: {
      "Content-Type": "application/json",
    },
  };
  axios
    .patch("/api/globals/" + props.globalDataId + "/name", formData, config)
    .then(() => {
      emit("update:title", titleModel.value);
    })
    .catch((error) => {
      updateTitleErrorMessage.value = error;
    })
    .finally(() => {
      isUpdatingTitle.value = false;
      editTitleMode.value = false;
    });
}

function editTitleField() {
  editTitleMode.value = true;
  nextTick(() => {
    if (titleInputRef.value) {
      titleInputRef.value.focus();
    }
  });
}
</script>
<style scoped lang="scss"></style>
