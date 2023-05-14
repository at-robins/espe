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
            <q-icon :name="prop.node.isFile ? matFilePresent : matFolder" />
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

        <q-card-section class="q-pt-none">
          <q-input
            dense
            v-model="folderName"
            autofocus
            @keyup.enter="createNewFolder"
          />
        </q-card-section>

        <q-card-actions align="right" class="text-primary">
          <q-btn flat label="Cancel" v-close-popup />
          <q-btn
            flat
            label="Add folder"
            @click="createNewFolder"
            v-close-popup
          />
        </q-card-actions>
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

const emit = defineEmits<{
  (event: "update:modelValue", nodes: FileTreeNode[]): void;
}>();

const fileTree = computed(() => {
  const computedValue: FileTreeNode[] = [
    {
      id: ROOT_ID,
      label: props.baseDirectoryLabel,
      children: [...props.modelValue],
      parents: [],
      isFile: false,
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
    const createdId = createNewNode(label, parent, false);
    // Expands the parent folder.
    if (createdId) {
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
    createNewNode(value.name, parent, true);
  }
}

function createNewNode(
  label: string,
  parent: FileTreeNode,
  isFile: boolean
): string | null {
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
  if (parent.children.some((node) => node.label === label)) {
    // TODO: show error
  } else {
    const id = parents.join("") + label;
    const newNode: FileTreeNode = {
      id: id,
      label: label,
      children: [],
      parents: parents,
      isFile: isFile,
    };
    currentNodes.push(newNode);
    emit("update:modelValue", newValue);
    return id;
  }
  return null;
}
</script>
