<template>
  <div class="row no-wrap">
    <div class="q-ma-md">
      <q-btn
        @click="saveOrChangeToEditMode"
        :icon="symOutlinedDescription"
        :loading="isUpdatingComment"
        round
      >
        <q-tooltip>Specify a description.</q-tooltip>
      </q-btn>
    </div>
    <q-separator vertical class="q-ml-md q-mr-md" />
    <q-menu v-if="!editCommentMode" touch-position context-menu>
      <q-list dense>
        <q-item clickable v-close-popup @click="editCommentMode = true">
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

    <div class="col">
      <div v-show="!editCommentMode" v-html="commentModel" class="q-ma-md" />
      <div v-show="editCommentMode" style="min-width: 100%" class="q-ma-md">
        <q-editor
          @update:model-value="onEditorUpdate"
          :readonly="!editCommentMode"
          style="min-width: 100%"
          v-model="commentModel"
          :dense="$q.screen.lt.md"
          :definitions="{
            save: {
              tip: 'Save changes',
              icon: 'save',
              label: 'Save',
              handler: saveOrChangeToEditMode,
            },
          }"
          :toolbar="[
            [
              {
                label: $q.lang.editor.align,
                icon: $q.iconSet.editor.align,
                fixedLabel: true,
                options: ['left', 'center', 'right', 'justify'],
              },
            ],
            [
              'bold',
              'italic',
              'strike',
              'underline',
              'subscript',
              'superscript',
            ],
            ['link'],
            [
              {
                label: $q.lang.editor.formatting,
                icon: $q.iconSet.editor.formatting,
                list: 'no-icons',
                options: ['p', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'code'],
              },
              {
                label: $q.lang.editor.fontSize,
                icon: $q.iconSet.editor.fontSize,
                fixedLabel: true,
                fixedIcon: true,
                list: 'no-icons',
                options: [
                  'size-1',
                  'size-2',
                  'size-3',
                  'size-4',
                  'size-5',
                  'size-6',
                  'size-7',
                ],
              },
              {
                label: $q.lang.editor.defaultFont,
                icon: $q.iconSet.editor.font,
                fixedIcon: true,
                list: 'no-icons',
                options: [
                  'default_font',
                  'arial',
                  'comic_sans',
                  'courier_new',
                  'impact',
                  'lucida_grande',
                  'times_new_roman',
                  'verdana',
                ],
              },
              'removeFormat',
            ],
            ['unordered', 'ordered'],

            ['undo', 'redo'],
            ['viewsource', 'save'],
          ]"
          :fonts="{
            arial: 'Arial',
            arial_black: 'Arial Black',
            comic_sans: 'Comic Sans MS',
            courier_new: 'Courier New',
            impact: 'Impact',
            lucida_grande: 'Lucida Grande',
            times_new_roman: 'Times New Roman',
            verdana: 'Verdana',
          }"
        />
      </div>
    </div>

    <q-tooltip
      v-if="updateCommentErrorMessage"
      class="bg-red"
      anchor="bottom start"
      self="top start"
    >
      {{ updateCommentErrorMessage }}
    </q-tooltip>
  </div>
</template>

<script setup lang="ts">
import axios from "axios";
import { ref, type PropType, type Ref } from "vue";
import { matEdit } from "@quasar/extras/material-icons";
import { symOutlinedDescription } from "@quasar/extras/material-symbols-outlined";

type Timer = ReturnType<typeof setTimeout>;

const props = defineProps({
  comment: {
    type: String as PropType<string | null>,
    default: null,
    required: false,
  },
  entityId: { type: Number, required: true },
  endpointType: { type: String, required: true },
});

const isUpdatingComment = ref(false);
const updateCommentErrorMessage = ref("");
const editCommentMode = ref(false);
const commentModel = ref(props.comment ? props.comment : "");

const updateTimer: Ref<Timer | null> = ref(null);
const UPDATE_DELAY = 3000.0;

function onEditorUpdate() {
  if (updateTimer.value) {
    clearTimeout(updateTimer.value);
  }
  updateTimer.value = setTimeout(clearTimerAndUpdateInstructions, UPDATE_DELAY);
}

function clearTimerAndUpdateInstructions() {
  updateTimer.value = null;
  updateComment();
}

function saveOrChangeToEditMode() {
  if (editCommentMode.value) {
    editCommentMode.value = false;
    clearTimerAndUpdateInstructions();
  } else {
    editCommentMode.value = true;
  }
}

function updateComment() {
  if (
    !isUpdatingComment.value &&
    (props.comment !== commentModel.value || !!updateCommentErrorMessage.value)
  ) {
    isUpdatingComment.value = true;
    updateCommentErrorMessage.value = "";
    const formData = JSON.stringify(commentModel.value);
    const config = {
      headers: {
        "Content-Type": "application/json",
      },
    };
    axios
      .patch(
        "/api/" + props.endpointType + "/" + props.entityId + "/comment",
        formData,
        config
      )
      .catch((error) => {
        updateCommentErrorMessage.value = error;
      })
      .finally(() => {
        isUpdatingComment.value = false;
      });
  }
}
</script>
<style scoped lang="scss"></style>
