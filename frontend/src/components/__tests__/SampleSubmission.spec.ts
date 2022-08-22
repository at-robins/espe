import SampleSubmission from "../SampleSubmission.vue";
import { describe, it, expect } from "vitest";

import { mount } from "@vue/test-utils";
import { Quasar } from "quasar";

describe("error if name is empty", () => {
  it("renders properly", async () => {
    const wrapper = mount(SampleSubmission, {
      props: {},
      global: {
        plugins: [Quasar],
      },
    });
    expect(wrapper.findAll(".q-field--error")).toHaveLength(0);
    const submit = wrapper.find(
      "[data-testid=sample-submission-button-submit]"
    );
    await submit.trigger("submit");
    expect(
      wrapper.findAll("label").at(0)?.classes().includes("q-field--error")
    ).toBe(true);
  });
});

describe("error if pipeline is not set", () => {
  it("renders properly", async () => {
    const wrapper = mount(SampleSubmission, {
      props: {},
      global: {
        plugins: [Quasar],
      },
    });
    expect(wrapper.findAll(".q-field--error")).toHaveLength(0);
    const inputName = wrapper.find(
      "[data-testid=sample-submission-input-name]"
    );
    await inputName.setValue("Dummy name");
    const submit = wrapper.find(
      "[data-testid=sample-submission-button-submit]"
    );
    await submit.trigger("submit");
    expect(
      wrapper.findAll("label").at(3)?.classes().includes("q-field--error")
    ).toBe(true);
  });
});

describe("error if e-mail incorrect", () => {
  it("renders properly", async () => {
    const wrapper = mount(SampleSubmission, {
      props: {},
      global: {
        plugins: [Quasar],
      },
    });
    expect(wrapper.findAll(".q-field--error")).toHaveLength(0);
    const inputName = wrapper.find(
      "[data-testid=sample-submission-input-name]"
    );
    await inputName.setValue("not an e-mail address");
    const inputMail = wrapper.find(
      "[data-testid=sample-submission-input-mail]"
    );
    await inputMail.setValue("Dummy name");
    const submit = wrapper.find(
      "[data-testid=sample-submission-button-submit]"
    );
    await submit.trigger("submit");
    expect(
      wrapper.findAll("label").at(1)?.classes().includes("q-field--error")
    ).toBe(true);
  });
});
