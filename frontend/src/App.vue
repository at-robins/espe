<script setup lang="ts">
import { matPublic } from "@quasar/extras/material-icons";
import {
  symOutlinedInfo,
  symOutlinedLibraryBooks,
  symOutlinedScience,
} from "@quasar/extras/material-symbols-outlined";
import { ref } from "vue";
import { RouterView, useRouter } from "vue-router";

const leftMenuOpen = ref(false);
const router = useRouter();

function navigateToGlobals() {
  router.push({ name: "globals" });
}

function navigateToExperiments() {
  router.push({ name: "experiments" });
}

function navigateToInfo() {
  router.push({ name: "info" });
}

function navigateToDocumentation() {
  router.push({ name: "documentation" });
}

function triggerSidebar() {
  leftMenuOpen.value = !leftMenuOpen.value;
}
</script>

<template>
  <q-layout view="hHh lpr fFf">
    <q-header elevated class="bg-primary text-white">
      <q-toolbar>
        <q-btn dense flat round icon="menu" @click="triggerSidebar" />

        <q-toolbar-title>
          <q-avatar class="q-mr-md">
            <img src="/icon_main.svg" />
          </q-avatar>
          <b>E</b>ncapsulated <b>S</b>equencing <b>P</b>ipeline <b>E</b>ngine
        </q-toolbar-title>
      </q-toolbar>
    </q-header>

    <q-drawer
      v-model="leftMenuOpen"
      side="left"
      overlay
      elevated
      @mouseleave="leftMenuOpen = false"
    >
      <q-list separator>
        <q-item clickable v-ripple @click="navigateToExperiments">
          <q-item-section avatar>
            <q-icon :name="symOutlinedScience" />
          </q-item-section>
          <q-item-section>Experiments</q-item-section>
        </q-item>
        <q-item clickable v-ripple @click="navigateToGlobals">
          <q-item-section avatar>
            <q-icon :name="matPublic" />
          </q-item-section>
          <q-item-section>Global data repositories</q-item-section>
        </q-item>
        <q-item clickable v-ripple @click="navigateToDocumentation">
          <q-item-section avatar>
            <q-icon :name="symOutlinedLibraryBooks" />
          </q-item-section>
          <q-item-section>Documentation</q-item-section>
        </q-item>
        <q-item clickable v-ripple @click="navigateToInfo">
          <q-item-section avatar>
            <q-icon :name="symOutlinedInfo" />
          </q-item-section>
          <q-item-section>Information</q-item-section>
        </q-item>
      </q-list>
    </q-drawer>

    <q-page-container>
      <router-view
        :key="router.currentRoute.value.path"
        @trigger-sidebar="triggerSidebar"
      />
    </q-page-container>
  </q-layout>
</template>

<style>
@import "@/assets/base.css";
</style>
