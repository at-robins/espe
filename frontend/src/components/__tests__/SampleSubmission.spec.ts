import SampleSubmission from "../SampleSubmission.vue";
import { describe, it, expect } from "vitest";

import { VueWrapper, mount } from "@vue/test-utils";
import { Quasar } from "quasar";
import { type ComponentPublicInstance } from "vue";

function mountWrapper(): VueWrapper<ComponentPublicInstance> {
  return mount(SampleSubmission, {
    props: {},
    global: {
      plugins: [Quasar],
    },
  });
}

function findSubmitButton(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.find("button");
}

function findFileInput(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.find("#sample-submission-input-file");
}

function findMailInput(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.find("#sample-submission-input-mail");
}

function findNameInput(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.find("#sample-submission-input-name");
}

function findPipelineInput(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.find("#sample-submission-input-pipeline");
}

function findAllErrors(wrapper: VueWrapper<ComponentPublicInstance>) {
  return wrapper.findAll(".q-field--error");
}

function triggerSubmit(
  wrapper: VueWrapper<ComponentPublicInstance>
): Promise<void> {
  return findSubmitButton(wrapper).trigger("submit");
}

describe("no error if everything is correct", () => {
  it("renders properly", async () => {
    const wrapper = mountWrapper();

    expect(findAllErrors(wrapper)).toHaveLength(0);

    const inputName = findNameInput(wrapper);
    await inputName.setValue("Dummy name");

    const inputMail = findMailInput(wrapper);
    await inputMail.setValue("a.valid@mail.com");

    const inputPipeline = findPipelineInput(wrapper);
    await inputPipeline.setValue({
      id: 0,
      name: "test",
      comment: "",
    });

    const inputFile = findFileInput(wrapper);
    const fileElement = inputFile.element as HTMLInputElement;
    const testFile = new File(["foo"], "programmatically_created.fastq.gz");
    const list = new FileList();
    list[0] = testFile;
    fileElement.files = list;
    await inputFile.trigger("change");

    await triggerSubmit(wrapper);
    expect(findAllErrors(wrapper)).toHaveLength(0);
  });
});

describe("error if name is empty", () => {
  it("renders properly", async () => {
    const wrapper = mountWrapper();
    expect(findAllErrors(wrapper)).toHaveLength(0);
    await triggerSubmit(wrapper);
    expect(
      wrapper.findAll("label").at(0)?.classes().includes("q-field--error")
    ).toBe(true);
  });
});

describe("error if pipeline is not set", () => {
  it("renders properly", async () => {
    const wrapper = mountWrapper();
    expect(findAllErrors(wrapper)).toHaveLength(0);
    const inputName = findNameInput(wrapper);
    await inputName.setValue("Dummy name");
    await triggerSubmit(wrapper);
    expect(
      wrapper.findAll("label").at(3)?.classes().includes("q-field--error")
    ).toBe(true);
  });
});

describe("error if e-mail incorrect", () => {
  it("renders properly", async () => {
    const wrapper = mountWrapper();
    expect(findAllErrors(wrapper)).toHaveLength(0);
    const inputName = findNameInput(wrapper);
    await inputName.setValue("Dummy name");
    const inputMail = findMailInput(wrapper);
    await inputMail.setValue("not an e-mail address");
    await triggerSubmit(wrapper);
    expect(
      wrapper.findAll("label").at(1)?.classes().includes("q-field--error")
    ).toBe(true);
  });
});
