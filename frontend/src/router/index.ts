import { createRouter, createWebHistory } from "vue-router";

const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  routes: [
    {
      path: "/",
      redirect: "/ui/",
    },
    {
      path: "/ui/",
      name: "home",
      component: () => import("../views/HomeView.vue"),
    },
    {
      path: "/ui/info",
      name: "info",
      component: () => import("../views/GeneralInformationView.vue"),
    },
    {
      path: "/ui/globals",
      name: "globals",
      component: () => import("@/views/GlobalDataView.vue"),
    },
    {
      path: "/ui/globals/:id",
      name: "globals_detail",
      component: () => import("@/views/GlobalDataDetailsView.vue"),
      props: true,
    },
    {
      path: "/ui/experiments",
      name: "experiments",
      component: () => import("@/views/ExperimentView.vue"),
    },
    {
      path: "/ui/experiments/:id",
      name: "experiments_detail",
      component: () => import("@/views/ExperimentDetailsView.vue"),
      props: true,
    },
    {
      path: "/ui/experiments/:id/run",
      name: "experiments_run_detail",
      component: () => import("@/views/ExperimentRunDetailsView.vue"),
      props: true,
    },
  ],
});

export default router;
