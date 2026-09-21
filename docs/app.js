import { computeMetrics } from './calculations.mjs';

const form = document.querySelector('#cmr-form');
const errorBox = document.querySelector('#error');
const themeButton = document.querySelector('#theme-toggle');
const resultFields = ['lvef', 'sv', 'edvi', 'esvi', 'svi', 'co', 'ci', 'lvmi', 'gls', 'gcs', 'grs', 'ecv'];

function formValues() {
  return Object.fromEntries(new FormData(form).entries());
}

function fmt(value, digits = 1) {
  return value == null ? '—' : Number(value).toFixed(digits);
}

function render(result) {
  const values = {
    lvef: `${fmt(result.lvef)}%`,
    sv: `${fmt(result.sv)} mL`,
    edvi: `${fmt(result.edvi)} mL/m²`,
    esvi: `${fmt(result.esvi)} mL/m²`,
    svi: `${fmt(result.svi)} mL/m²`,
    co: `${fmt(result.co, 2)} L/min`,
    ci: `${fmt(result.ci, 2)} L/min/m²`,
    lvmi: `${fmt(result.lvmi)} g/m²`,
    gls: `${fmt(result.gls)}%`,
    gcs: `${fmt(result.gcs)}%`,
    grs: `${fmt(result.grs)}%`,
    ecv: result.ecv == null ? 'Not calculated' : `${fmt(result.ecv)}%`,
  };
  resultFields.forEach((key) => {
    document.querySelector(`[data-result="${key}"]`).textContent = values[key];
  });
}

function analyze() {
  errorBox.hidden = true;
  try {
    render(computeMetrics(formValues()));
  } catch (error) {
    errorBox.textContent = error.message;
    errorBox.hidden = false;
  }
}

form.addEventListener('submit', (event) => {
  event.preventDefault();
  analyze();
});

document.querySelector('#reset').addEventListener('click', () => {
  form.reset();
  errorBox.hidden = true;
  analyze();
});

let savedTheme = null;
try { savedTheme = localStorage.getItem('cmr-theme'); } catch (_) { /* Storage may be unavailable. */ }
if (savedTheme === 'dark') document.documentElement.dataset.theme = 'dark';
themeButton.setAttribute('aria-pressed', String(savedTheme === 'dark'));
themeButton.addEventListener('click', () => {
  const dark = document.documentElement.dataset.theme !== 'dark';
  document.documentElement.dataset.theme = dark ? 'dark' : 'light';
  try { localStorage.setItem('cmr-theme', dark ? 'dark' : 'light'); } catch (_) { /* Theme still works for this page view. */ }
  themeButton.setAttribute('aria-pressed', String(dark));
});

analyze();
