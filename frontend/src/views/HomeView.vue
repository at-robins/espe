<script setup lang="ts">
import { useRouter } from "vue-router";

const router = useRouter();

function navigateToDocumentation() {
  router.push({ name: "documentation" });
}

function navigateToExperiments() {
  router.push({ name: "experiments" });
}
</script>

<template>
  <main>
    <!-- <div class="q-ma-lg row flex-center">
      <img class="fade-in-text" src="/icon_main.svg" style="width: 80%; position: fixed;" />
    </div> -->
    <div class="q-ma-xl row flex-center fade-in-text">
      <div class="text-description q-pa-lg">
        ESPE is a sequencing pipeline engine to automate your biological
        workflows.<br />
        Get started by creating a new
        <q-btn @click="navigateToExperiments" label="experiment" no-caps />
        or reading the
        <q-btn
          @click="navigateToDocumentation"
          label="documentation"
          no-caps
        />.
        <br />
        Everything can also be accessed via the sidebar.
      </div>
    </div>
  </main>
</template>

<style scoped lang="scss">
.fade-in-text {
  animation: fadeIn 2s;
  -webkit-animation: fadeIn 2s;
  -moz-animation: fadeIn 2s;
  -o-animation: fadeIn 2s;
  -ms-animation: fadeIn 2s;
  background-image: url("/icon_main.svg");
  background-position: center;
  background-repeat: no-repeat;
  background-size: contain;
  height: 80vh;
}

.text-description {
  font-size: 5vh;
  text-align: center;
  background-color: rgb(255, 255, 255, 0.9);
  border: black;
  border-style: solid;
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
