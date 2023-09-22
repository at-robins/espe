<script setup lang="ts">
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
</script>

<template>
  <q-layout view="hHh lpr fFf">
    <q-header elevated class="bg-primary text-white">
      <q-toolbar>
        <q-btn
          dense
          flat
          round
          icon="menu"
          @click="leftMenuOpen = !leftMenuOpen"
        />

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
          <q-item-section>Experiments</q-item-section>
        </q-item>
        <q-item clickable v-ripple @click="navigateToGlobals">
          <q-item-section>Global data repositories</q-item-section>
        </q-item>
      </q-list>
    </q-drawer>

    <q-page-container>
      <router-view :key="router.currentRoute.value.path" />
    </q-page-container>
  </q-layout>
</template>

<style>
@import "@/assets/base.css";
</style>
