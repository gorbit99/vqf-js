import { defaultParams, VQF, VQFParams } from "./vqf";

import "./style.scss";

const vqfParamsElement = document.querySelector(".vqf-params");
const fileInputElement = document.querySelector(
  "#fileInput",
) as HTMLInputElement;
const resultsElement = document.querySelector(
  ".result-download",
) as HTMLButtonElement;

const sampleRateElement = document.querySelector(
  "#sampleRate",
) as HTMLInputElement;
const gyrHzElement = document.querySelector("#gyrHz") as HTMLInputElement;
const accHzElement = document.querySelector("#accHz") as HTMLInputElement;
const magHzElement = document.querySelector("#magHz") as HTMLInputElement;
const useMagElement = document.querySelector("#useMag") as HTMLInputElement;
const magHzWrapper = document.querySelector(
  ".mag-hz-wrapper",
) as HTMLDivElement;

const hideResult = () => resultsElement.classList.add("hidden");
sampleRateElement.addEventListener("change", hideResult);
gyrHzElement.addEventListener("change", hideResult);
accHzElement.addEventListener("change", hideResult);
magHzWrapper.addEventListener("change", hideResult);
useMagElement.addEventListener("change", hideResult);

useMagElement.addEventListener("change", () => {
  magHzWrapper.classList.toggle("hidden", !useMagElement.checked);
});
magHzWrapper.classList.toggle("hidden", !useMagElement.checked);

const vqfParams = defaultParams();

vqfParamsElement?.append(
  ...Object.keys(vqfParams).map((key) =>
    typeof vqfParams[key as keyof VQFParams] === "boolean"
      ? createBooleanInput(key as keyof VQFParams)
      : createNumberInput(key as keyof VQFParams),
  ),
);

function createBooleanInput(paramName: keyof VQFParams) {
  const wrapper = document.createElement("div");
  wrapper.classList.add("input__wrapper", "input__wrapper--boolean");
  const label = document.createElement("label");
  label.classList.add("input__label");
  label.htmlFor = "input-" + paramName;
  label.innerText = paramName;
  const input = document.createElement("input");
  input.type = "checkbox";
  input.name = paramName;
  input.id = "input-" + paramName;
  input.checked = vqfParams[paramName] as boolean;
  input.addEventListener("change", () => {
    (vqfParams as any)[paramName] = input.checked;
    hideResult();
  });

  wrapper.append(label, input);
  return wrapper;
}

function createNumberInput(paramName: keyof VQFParams) {
  const wrapper = document.createElement("div");
  wrapper.classList.add("input__wrapper", "input__wrapper--number");
  const label = document.createElement("label");
  label.classList.add("input__label");
  label.htmlFor = "input-" + paramName;
  label.innerText = paramName;
  const input = document.createElement("input");
  input.type = "number";
  input.name = paramName;
  input.id = "input-" + paramName;
  input.value = (vqfParams[paramName] as number).toString();
  input.addEventListener("change", () => {
    (vqfParams as any)[paramName] = parseFloat(input.value);
    hideResult();
  });

  wrapper.append(label, input);
  return wrapper;
}

fileInputElement.addEventListener("change", async () => {
  const file = fileInputElement.files?.[0];

  if (!file) {
    resultsElement.classList.add("hidden");
    return;
  }

  const gyroRate = parseFloat(gyrHzElement.value);
  const accRate = parseFloat(accHzElement.value);
  const magRate = parseFloat(magHzElement.value);

  const vqf = new VQF(
    1 / gyroRate,
    1 / accRate,
    useMagElement.checked ? 1 / magRate : -1,
    vqfParams,
  );

  const text = await file.text();

  const output: [number, number, number, number][] = [];

  let gyroSampleCount = 0;

  text.split("\n").map((line) => {
    const [type, x, y, z] = line.split(";").map((value) => parseFloat(value));

    if (type === 0) {
      vqf.updateGyr([x, y, z]);
      gyroSampleCount++;
      if (
        gyroSampleCount / gyroRate >=
        1 / parseFloat(sampleRateElement.value)
      ) {
        output.push(vqf.getQuat6D());
        gyroSampleCount = 0;
      }
    } else if (type === 1) {
      vqf.updateAcc([x, y, z]);
    } else {
      vqf.updateMag([x, y, z]);
    }
  });

  const resultText = [
    "Time;W;X;Y;Z",
    ...output.map((line, index) => {
      const time = (1 / parseFloat(sampleRateElement.value)) * index;
      return [time, ...line].join(";");
    }),
  ].join("\n");

  resultsElement.classList.remove("hidden");

  resultsElement.addEventListener("click", () => {
    const blob = new Blob([resultText], { type: "text/csv" });
    const elem = window.document.createElement("a");
    elem.href = window.URL.createObjectURL(blob);
    elem.download = "result.csv";
    document.body.appendChild(elem);
    elem.click();
    document.body.removeChild(elem);
  });
});
