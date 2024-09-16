<script setup lang="ts">
import { matMenu } from "@quasar/extras/material-icons";
import {
  symOutlinedLibraryBooks,
  symOutlinedScience,
} from "@quasar/extras/material-symbols-outlined";
import { useRouter } from "vue-router";

const router = useRouter();
const emit = defineEmits<{
  (event: "triggerSidebar"): void;
}>();

function navigateToDocumentation() {
  router.push({ name: "documentation" });
}

function navigateToExperiments() {
  router.push({ name: "experiments" });
}

function clickSidebarButton() {
  emit("triggerSidebar");
}
</script>

<template>
  <main>
    <div class="q-pa-xl background-image-container row flex-center">
      <img class="background-image" src="/icon_main.svg" />
    </div>
    <div class="q-pa-lg row flex-center text-description-container">
      <div class="text-description">
        ESPE is a sequencing pipeline engine to automate your biological
        workflows.<br />
        Get started by creating a new
        <span class="q-pa-sm">
          <q-btn
            @click="navigateToExperiments"
            :icon="symOutlinedScience"
            label="experiment"
            no-caps
          />
        </span>
        or reading the
        <span class="q-pa-sm">
          <q-btn
            @click="navigateToDocumentation"
            :icon="symOutlinedLibraryBooks"
            label="documentation"
            no-caps
          />
        </span>
        .
        <br />
        <span class="q-pa-sm">
          Everything can also be accessed conveniently via the
          <q-btn
            color="primary"
            :icon="matMenu"
            label="sidebar"
            no-caps
            @click="clickSidebarButton"
          />
        </span>
        .
      </div>
    </div>
  </main>
</template>

<style scoped lang="scss">
.background-image {
  height: 40vh;
  opacity: 0.8;
}

.background-image-container {
  width: 100%;
  height: 100%;
  animation: fadeIn-image 5s;
  -webkit-animation: fadeIn-image 5s;
  -moz-animation: fadeIn-image 5s;
  -o-animation: fadeIn-image 5s;
  -ms-animation: fadeIn-image 5s;
}

.text-description-container {
  animation: fadeIn 2s;
  -webkit-animation: fadeIn 2s;
  -moz-animation: fadeIn 2s;
  -o-animation: fadeIn 2s;
  -ms-animation: fadeIn 2s;
}

.text-description {
  line-height: 3em;
  font-size: 1.75em;
  text-align: center;
}

@keyframes fadeIn-image {
  0% {
    opacity: 0;
  }
  100% {
    opacity: 1;
  }
}

@-moz-keyframes fadeIn-image {
  0% {
    opacity: 0;
  }
  100% {
    opacity: 1;
  }
}

@-webkit-keyframes fadeIn-image {
  0% {
    opacity: 0;
  }
  100% {
    opacity: 1;
  }
}

@-o-keyframes fadeIn-image {
  0% {
    opacity: 0;
  }
  100% {
    opacity: 1;
  }
}

@-ms-keyframes fadeIn-image {
  0% {
    opacity: 0;
  }
  100% {
    opacity: 1;
  }
}

@keyframes fadeIn {
  0% {
    transform: translateY(100%);
    opacity: 0;
  }
  100% {
    transform: translateY(0%);
    opacity: 1;
  }
}

@-moz-keyframes fadeIn {
  0% {
    transform: translateY(100%);
    opacity: 0;
  }
  100% {
    transform: translateY(0%);
    opacity: 1;
  }
}

@-webkit-keyframes fadeIn {
  0% {
    transform: translateY(100%);
    opacity: 0;
  }
  100% {
    transform: translateY(0%);
    opacity: 1;
  }
}

@-o-keyframes fadeIn {
  0% {
    transform: translateY(100%);
    opacity: 0;
  }
  100% {
    transform: translateY(0%);
    opacity: 1;
  }
}

@-ms-keyframes fadeIn {
  0% {
    transform: translateY(100%);
    opacity: 0;
  }
  100% {
    transform: translateY(0%);
    opacity: 1;
  }
}
</style>
