import { createApp } from "vue";
import App from "./App.vue";
import router from "./router";
import "@quasar/extras/material-icons/material-icons.css";
import "quasar/src/css/index.sass";

import { Quasar } from "quasar";
import quasarLang from "quasar/lang/en-GB";
import quasarIconSet from "quasar/icon-set/svg-material-icons";

const app = createApp(App);

app.use(router);

app.use(Quasar, {
  plugins: {},
  lang: quasarLang,
  iconSet: quasarIconSet,
});

app.mount("#app");
