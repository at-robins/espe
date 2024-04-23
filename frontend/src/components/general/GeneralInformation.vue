<template>
  <div class="q-pa-lg">
    <div class="q-pb-xl">
      <span class="row" style="align-items: center">
        <q-avatar class="row q-pr-md">
          <q-icon :name="symOutlinedGroups" />
        </q-avatar>
        <h3 class="">Team</h3>
      </span>
      <q-separator inset class="q-mt-md q-mb-md" />
      <author-item
        name="Robin Patrick Schenk"
        image-url="https://www.immunology-med2.de/-_-/res/841d1fe7-1eed-4014-9ded-23b8a68a7178/images/files/841d1fe7-1eed-4014-9ded-23b8a68a7178/f5b48615-6425-4370-90c8-3f0cd5473dcf/480-480/99dfa70b33654583f760ce4f390ff9dbf7524ab6"
        :infos="[
          {
            name: 'Role',
            value: 'Lead developer',
            url: null,
          },
          {
            name: 'ORCID',
            value: '0000-0001-6823-8128',
            url: 'https://orcid.org/0000-0001-6823-8128',
          },
          {
            name: 'Affiliation',
            value: 'Wiedemann lab',
            url: 'https://www.immunology-med2.de/wiedemann-lab',
          },
        ]"
      />
      <author-item
        name="PD Dr. Gabriela M. Wiedemann"
        image-url="https://www.immunology-med2.de/-_-/res/841d1fe7-1eed-4014-9ded-23b8a68a7178/images/files/841d1fe7-1eed-4014-9ded-23b8a68a7178/1be41883-d2fb-4918-9f95-b5d22c3c196c/768-768/5480952c7d0b69c80c0ed82e5547b945e3b614ab"
        :infos="[
          {
            name: 'Role',
            value: 'Principle investigator',
            url: null,
          },
          {
            name: 'ORCID',
            value: '0000-0001-6138-9791',
            url: 'https://orcid.org/0000-0001-6138-9791',
          },
          {
            name: 'Affiliation',
            value: 'Wiedemann lab',
            url: 'https://www.immunology-med2.de/wiedemann-lab',
          },
        ]"
      />
      <author-item
        name="Kevin Schmid"
        image-url="https://www.immunology-med2.de/-_-/res/841d1fe7-1eed-4014-9ded-23b8a68a7178/images/files/841d1fe7-1eed-4014-9ded-23b8a68a7178/98211c54-8678-485c-953a-1b7e7aa4aa21/768-768/b465f8273383b30fe74a4bbbd8587ba8240e37b3"
        :infos="[
          {
            name: 'Role',
            value: 'Tester',
            url: null,
          },
          {
            name: 'Affiliation',
            value: 'Wiedemann lab',
            url: 'https://www.immunology-med2.de/wiedemann-lab',
          },
        ]"
      />
      <author-item
        name="Sebastian Deschler"
        image-url="https://www.immunology-med2.de/-_-/res/841d1fe7-1eed-4014-9ded-23b8a68a7178/images/files/841d1fe7-1eed-4014-9ded-23b8a68a7178/6a2b099b-bc54-47aa-8559-c2b3a9430ac9/768-751/acd241b5f590b446a772d0ecfff48701df3cb0a4"
        :infos="[
          {
            name: 'Role',
            value: 'Tester',
            url: null,
          },
          {
            name: 'ORCID',
            value: '0000-0003-4680-9609',
            url: 'https://orcid.org/0000-0003-4680-9609',
          },
          {
            name: 'Affiliation',
            value: 'Böttcher lab',
            url: 'https://www.immunology-med2.de/bottcher-lab',
          },
        ]"
      />
    </div>
    <div class="q-pb-xl">
      <span class="row" style="align-items: center">
        <q-avatar class="row q-pr-md">
          <q-icon :name="symOutlinedCode" />
        </q-avatar>
        <h3 class="">Source code</h3>
      </span>
      <q-separator inset class="q-mt-md q-mb-md" />
      <div class="q-pl-md">
        <p>
          The ESPE source code is freely available via
          <a href="https://github.com/at-robins/espe">GitHub</a>.
        </p>
      </div>
    </div>
    <div class="q-pb-xl">
      <span class="row" style="align-items: center">
        <q-avatar class="row q-pr-md">
          <q-icon :name="symOutlinedBugReport" />
        </q-avatar>
        <h3 class="">Reporting issues</h3>
      </span>
      <q-separator inset class="q-mt-md q-mb-md" />
      <div class="q-pl-md">
        <p>
          If you encounter any issues or bugs, you can create a ticket using the
          <a href="https://github.com/at-robins/espe/issues">issue tracker</a>.
        </p>
      </div>
    </div>
    <div class="q-pb-xl">
      <span class="row" style="align-items: center">
        <q-avatar class="row q-pr-md">
          <q-icon :name="symOutlinedBookmark" />
        </q-avatar>
        <h3 class="">Cite</h3>
      </span>
      <q-separator inset class="q-mt-md q-mb-md" />
      <div class="q-pl-md">
        <p>
          There is currently no publication for our pipeline engine. In the
          meantime you can cite it as follows:
        </p>
        <p>
          Schenk, R. P., & Wiedemann, G. M. (2024). ESPE: Encapsulated
          Sequencing Pipeline Engine (0.1).
          <a href="https://github.com/at-robins/espe"
            >https://github.com/at-robins/espe</a
          >
        </p>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
import {
  symOutlinedBookmark,
  symOutlinedBugReport,
  symOutlinedCode,
  symOutlinedGroups,
} from "@quasar/extras/material-symbols-outlined";
import AuthorItem from "./AuthorItem.vue";
</script>
<style scoped lang="scss"></style>
