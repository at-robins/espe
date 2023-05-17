<template>
  <div>
    <q-tree
      ref="treeReference"
      :nodes="fileTree"
      node-key="id"
      no-nodes-label="No files present."
      text-color="primary"
      :selected-color="showContextMenu ? 'accent' : 'primary'"
      v-model:selected="selectedNode"
      default-expand-all
    >
      <template v-slot:default-header="prop">
        <div class="col items-center">
          <div
            @click="
              treeReference?.setExpanded(
                prop.node.id,
                !treeReference?.isExpanded(prop.node.id)
              )
            "
            @contextmenu="selectedNode = prop.node.id"
            class="col"
          >
            <q-icon
              v-if="!prop.node.isUploaded"
              :name="prop.node.isFile ? matFilePresent : matFolder"
            />
            <q-spinner v-else />
            {{ prop.node.label }}
          </div>
          <q-menu
            @update:model-value="(value) => (showContextMenu = value)"
            touch-position
            context-menu
          >
            <q-list dense>
              <q-item
                v-if="!prop.node.isFile"
                clickable
                v-close-popup
                @click="openNewFolderPopup"
              >
                <q-item-section
                  ><div class="row flex-center no-wrap">
                    <q-icon
                      :name="matCreateNewFolder"
                      color="primary"
                      size="xs"
                      class="q-mr-xs"
                    />
                    <div class="col">Add folder...</div>
                  </div></q-item-section
                >
              </q-item>
              <q-item
                v-if="!prop.node.isFile"
                clickable
                v-close-popup
                @click="pickFile"
              >
                <q-item-section
                  ><div class="row flex-center no-wrap">
                    <q-icon
                      :name="matUploadFile"
                      color="primary"
                      size="xs"
                      class="q-mr-xs"
                    />
                    <div class="col">Add file...</div>
                  </div></q-item-section
                >
              </q-item>
              <q-item v-if="prop.node.id !== ROOT_ID" clickable v-close-popup>
                <q-item-section>
                  <div class="row flex-center no-wrap">
                    <q-icon
                      :name="matDelete"
                      color="negative"
                      size="xs"
                      class="q-mr-xs"
                    />
                    <div class="col">Delete</div>
                  </div>
                </q-item-section>
              </q-item>
            </q-list>
          </q-menu>
        </div>
      </template>
    </q-tree>
    <q-dialog v-model="showFolderNamePopup" persistent>
      <q-card style="min-width: 350px">
        <q-card-section>
          <div class="text-h6">Folder name</div>
        </q-card-section>
        <q-form @submit="createNewFolder" :greedy="true">
          <q-card-section class="q-pt-none">
            <q-input
              dense
              v-model="folderName"
              autofocus
              :rules="componentValidationRules"
            />
          </q-card-section>

          <q-card-actions align="right" class="text-primary">
            <q-btn
              flat
              label="Cancel"
              @click="folderName = null"
              v-close-popup
            />
            <q-btn flat label="Add folder" type="submit" />
          </q-card-actions>
        </q-form>
      </q-card>
    </q-dialog>
    <q-file
      v-model="fileModel"
      ref="fileReference"
      @update:model-value="createNewFile"
      style="display: none"
    />
  </div>
</template>

<script setup lang="ts">
import { type FileTreeNode } from "@/scripts/types";
import { computed, ref, type PropType, type Ref, nextTick } from "vue";
import {
  matDelete,
  matCreateNewFolder,
  matFolder,
  matFilePresent,
  matUploadFile,
} from "@quasar/extras/material-icons";
import type { QFile, QTree } from "quasar";

const ROOT_ID = ".";
const ILLEGAL_COMPONENTS = [".", "..", "~"];
const ILLEGAL_COMPONENT_CHARACTERS = [
  "/",
  "\\",
  "<",
  ">",
  ":",
  "*",
  "?",
  "|",
  '"',
  "\x00",
  "\x01",
  "\x02",
  "\x03",
  "\x04",
  "\x05",
  "\x06",
  "\x07",
  "\x08",
  "\x09",
  "\x0A",
  "\x0B",
  "\x0C",
  "\x0D",
  "\x0E",
  "\x0F",
  "\x10",
  "\x11",
  "\x12",
  "\x13",
  "\x14",
  "\x15",
  "\x16",
  "\x17",
  "\x18",
  "\x19",
  "\x1A",
  "\x1B",
  "\x1C",
  "\x1D",
  "\x1E",
  "\x1F",
  "\x7F",
];

const props = defineProps({
  modelValue: { type: Object as PropType<FileTreeNode[]>, required: true },
  baseDirectoryLabel: { type: String, default: "Root", required: false },
});

const selectedNode: Ref<string | null> = ref(null);
const showFolderNamePopup = ref(false);
const folderName: Ref<string | null> = ref(null);
const treeReference: Ref<QTree | null> = ref(null);
const fileModel: Ref<File | null> = ref(null);
const fileReference: Ref<QFile | null> = ref(null);
const showContextMenu = ref(false);

const componentValidationRules = [
  (val: string) => !!val || "The name cannot be empty.",
  (val: string) =>
    !treeReference.value
      ?.getNodeByKey(selectedNode.value)
      .children.some((node: FileTreeNode) => node.label === val) ||
    "The name already exists.",
  (val: string) =>
    !ILLEGAL_COMPONENTS.some((illegal_comp) => illegal_comp === val) ||
    'The name cannot be any of the following ""' +
      ILLEGAL_COMPONENTS.join(", ") +
      '"".',
  (val: string) =>
    !ILLEGAL_COMPONENT_CHARACTERS.some((illegalChar) =>
      val.includes(illegalChar)
    ) ||
    "The name contains one of the following illegal characters: " +
      ILLEGAL_COMPONENT_CHARACTERS.join(", "),
];

const emit = defineEmits<{
  (event: "update:modelValue", nodes: FileTreeNode[]): void;
  (event: "addedFile", file: File, node: FileTreeNode): void;
}>();

const fileTree = computed(() => {
  const computedValue: FileTreeNode[] = [
    {
      id: ROOT_ID,
      label: props.baseDirectoryLabel,
      children: [...props.modelValue],
      parents: [],
      isFile: false,
      isUploaded: false,
    },
  ];
  return computedValue;
});

function openNewFolderPopup() {
  showFolderNamePopup.value = true;
}

function pickFile() {
  fileReference.value?.pickFiles();
}

function createNewFolder() {
  showFolderNamePopup.value = false;
  const parent: FileTreeNode | null = treeReference.value?.getNodeByKey(
    selectedNode.value
  );
  const label = folderName.value;
  if (parent && label) {
    const createdNode = createNewNode(label, parent, false);
    // Expands the parent folder.
    if (createdNode) {
      nextTick(() => treeReference.value?.setExpanded(parent.id, true));
    }
  }
  folderName.value = null;
}

function createNewFile(value: File | null) {
  const parent: FileTreeNode | null = treeReference.value?.getNodeByKey(
    selectedNode.value
  );
  if (value && parent) {
    const createdNode = createNewNode(value.name, parent, true);
    if (createdNode) {
      emit("addedFile", value, createdNode);
      nextTick(() => treeReference.value?.setExpanded(parent.id, true));
    }
  }
}

function createNewNode(
  label: string,
  parent: FileTreeNode,
  isFile: boolean
): FileTreeNode | null {
  let parents: string[] = [];
  if (parent.id !== ROOT_ID) {
    parents = [...parent.parents];
    parents.push(parent.id);
  }
  const newValue = [...props.modelValue];
  let currentNodes = newValue;
  for (const parent of parents) {
    const found = currentNodes.find((node) => node.id === parent)?.children;
    if (found) {
      currentNodes = found;
    }
  }
  if (!parent.children.some((node) => node.label === label)) {
    const id = parents.join("") + label;
    const newNode: FileTreeNode = {
      id: id,
      label: label,
      children: [],
      parents: parents,
      isFile: isFile,
      isUploaded: false,
    };
    currentNodes.push(newNode);
    emit("update:modelValue", newValue);
    return newNode;
  }
  return null;
}
</script>
