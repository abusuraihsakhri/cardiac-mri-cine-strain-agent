# Cardiac MRI Cine & Strain Calculator

Quantitative cardiac MRI calculator for already-measured cine volumes, global strain values, and T1/ECV inputs. It derives volumetric and hemodynamic indices plus ECV without external runtime dependencies.

## What it calculates

- LVEF and stroke volume from LVEDV/LVESV
- Indexed LV volumes, stroke volume, and LV mass using BSA
- Cardiac output and cardiac index
- GLS/GCS/GRS summary values supplied by the user
- Extracellular volume fraction from paired myocardial/blood T1 values and hematocrit
- CSV batch processing through the Python CLI

The project does **not** segment cine images, perform feature tracking from DICOM data, infer a disease diagnosis, or recommend treatment. CMR strain and tissue-mapping reference intervals are method dependent; interpretation should use an appropriate local or method-specific reference interval.

## Use

Python 3.9+:

```bash
python -m pip install -e .
cardio-mri-strain evaluate --edv 145 --esv 58 --hr 72 --bsa 1.85 --gls -20.5 --gcs -22 --grs 42
cardio-mri-strain batch -i sample.csv -o results.csv
```

Run tests:

```bash
python -m pip install pytest
python -m pytest -q
```

## Browser application

The GitHub Pages application performs the core arithmetic locally in the browser. It has no server component and does not upload entered values. Do not enter direct patient identifiers into a public or shared browser session.

The browser version uses a small JavaScript implementation of the same formulas rather than Pyodide. For this calculator, loading a full Python/WASM runtime would add substantial startup and network overhead without adding functionality.

## Technology

- Python standard library for the CLI and batch processor
- Static HTML/CSS/JavaScript for GitHub Pages
- GitHub Actions for CI and Pages deployment

Current evergreen versions of Chrome, Edge, Firefox, and Safari are targeted. The responsive layout supports desktop and mobile screens, with light mode by default and an optional dark theme.

## License

MIT. See [LICENSE](LICENSE).
