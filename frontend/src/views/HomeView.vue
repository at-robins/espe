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
    <div class="background-image-container row flex-center">
      <img class="background-image" src="/icon_main.svg" />
    </div>
    <div class="row flex-center text-description-container">
      <div class="text-description q-pa-xl">
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
.background-image {
  height: 90%;
  opacity: 0.3;
}

.background-image-container {
  position: absolute;
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
  height: 80vh;
}

.text-description {
  line-height: 3em;
  font-size: 1.5em;
  text-align: center;
  background-color: rgb(255, 255, 255, 0.9);
  border: black;
  border-style: solid;
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
