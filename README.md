# Kidney function and albuminuria assessment in type 2 diabetes

Analysis code for kidney-marker availability and paired eGFR–UACR categories
among adults hospitalized in Endocrinology and Diabetology, 2017–September 2026.
The code follows the v27 study analysis.

## Data

**No real patient data or study results are included.** `mock_data.py` generates
four linked tables in memory from wholly artificial values. Records are not
sampled or reconstructed from patients; identifiers begin with `SYNTHETIC-`.
The scripts read no clinical files. Example outputs demonstrate execution and
do not reproduce the paper's results.

## Files and use

- `derive.py`: eligibility, laboratory selection, renal categories and medications.
- `run_analysis.py`: descriptive tables and sensitivity analyses.
- `mock_data.py`: artificial input tables.
- `validate.py`: independent selection, calculation and denominator checks.
- `variable_dictionary.csv`: input fields and principal derived variables.
- `code_lists.csv`: clinical concepts, units and classification rules.

Use Python 3.14 with the pinned dependencies:

```sh
python3.14 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python run_analysis.py
python validate.py
```

Aggregate synthetic outputs are written to `outputs/`, which is excluded from
Git. Patient-level examples remain in memory. The study flowchart is maintained
in BioRender; this repository does not generate figures.

## Study definitions

Inputs represent the reviewed T2D cohort. Source-level diabetes adjudication and
clinical extraction precede this code. The analysis selects discharged adults,
excludes current pregnancy and uses each patient's first eligible admission.
Missing markers do not exclude patients from the availability denominator.

Paired categories use the earliest preferred reported UACR with a creatinine
within seven calendar days. eGFR uses CKD-EPI 2021. Categories split at eGFR
60 mL/min/1.73 m² and UACR 30 mg/g; maintenance dialysis and kidney transplant
history exclude patients from paired classification. Single-admission categories
do not establish chronic kidney disease. Detailed rules are in the code and
accompanying dictionaries.

## Analyses

Outputs cover marker availability, renal categories, sample selection and
contemporary discharge recommendations. Proportions use Wilson 95% confidence
intervals; medication estimates use 10,000 patient-clustered bootstrap resamples.

Sensitivity analyses vary pairing time, creatinine selection, age bounds,
reported eGFR, qualifying admission and the observation window. The laboratory's
reported eGFR equation and indexing convention remain unspecified. There are no
regression models, hypothesis tests or imputations.

Code and artificial examples use the MIT License.
