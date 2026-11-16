# Kidney function and albuminuria assessment in type 2 diabetes

Analysis code for the E11 clinic study of kidney-marker availability and paired
eGFR–UACR categories at the first eligible hospitalization.

**No real patient data or study results are included.** `mock_data.py` creates
wholly artificial linked examples in memory from fixed values and arithmetic
sequences. They are not sampled, altered or reconstructed from patients. All
identifiers begin with `SYNTHETIC-`. No clinical files or excerpts are read.
Example outputs demonstrate the code and do not reproduce the paper's results.

## Files and use

- `derive.py`: eligibility, laboratory selection, renal categories and medications.
- `run_analysis.py`: descriptive tables and all planned sensitivity analyses.
- `mock_data.py`: four artificial input tables.
- `validate.py`: independent pairing, eGFR, interval and denominator checks.
- `variable_dictionary.csv`: input fields and principal derived variables.
- `code_lists.csv`: clinical concepts, units and classification rules.

Use Python 3.14 and the pinned packages:

```sh
python3.14 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python run_analysis.py
python validate.py
```

Only aggregate synthetic outputs are written, to `outputs/`; patient-level
examples remain in memory. The `.gitignore` allows only the ten reviewed source
and documentation files. No dataset, patient audit, manuscript or figure is
included. The numerical code is adapted directly from the study scripts;
clinical file loading and patient-level exports have been removed.

## Selection and variables

The input represents the released cohort after source-level review of diabetes
type. E11-coded admissions were reviewed using diagnosis and clinical
information. Documented non-T2D and unresolved type conflicts were excluded
upstream; E10 wording alone did not override explicit adjudicated T2D. This
review and the extraction of structured clinical assertions are not reconstructed
from raw medical records by the analysis scripts.

Select discharged adults with admissions beginning 1 January 2017–3 September
2026. Exclude current pregnancy; historical pregnancy alone does not exclude.
Use each patient's first eligible admission, ordered by admission date,
discharge date and identifier. Missing markers or covariates do not exclude
patients from the availability analysis.

Retain dialysis/transplant patients in the availability denominator. Exclude
current maintenance dialysis and kidney transplant history from paired categories.
Historical or acute dialysis alone does not exclude. Status is admission-specific;
unresolved renal status stops classification.

Laboratory results must be linked, numeric, exact, in the expected units and
within the admission dates. Use collection time, otherwise execution time;
exclude results executed after discharge. eGFR uses CKD-EPI 2021, released sex
and calendar-year age. Select the earliest preferred reported UACR with a
creatinine within seven calendar days. Rank creatinine matches by calendar-day
gap, elapsed hours, timestamp and result identifier. UACR timestamp ties use
its result identifier first. Reported zero UACR is retained. Derived alternatives
are not counted twice or substituted for reported values. The study release
had no preferred derived-only UACR. Missing UACR is never classified as normal.

The four groups split at eGFR 60 mL/min/1.73 m² and UACR 30 mg/g. Single-admission
categories do not establish chronic kidney disease. BMI uses the released
height/weight calculation, otherwise documented BMI. HbA1c is the first usable
inpatient result. No covariates are imputed.

## Analyses

Primary analysis: joint-marker availability during the first eligible admission,
regardless of measurement interval, with clinic and clinic/year descriptions.
Secondary analysis: paired renal categories within seven days, conditional
marker proportions, detailed G-by-A categories and selection into the paired
sample. Report n/N (%), medians [IQR] and Wilson 95% confidence intervals.

Supplementary medication descriptions use all qualifying Diabetology admissions
beginning 14 July 2023, with an available discharge regimen and the same renal
exclusions. Identify active recommendations using eligible discharge mentions,
ingredients and classes, including combination products. Sacubitril/valsartan
counts as ARB exposure and is also counted separately. Finerenone is a count.
SGLT2 and ACEi/ARB intervals use 10,000 patient-clustered bootstrap resamples
(seed 20261115), retaining each sampled patient's admissions and categories.
These are recommendations, not verified medication use or treatment eligibility.

Sensitivities cover same-day pairing; first and last creatinine with UACR fixed;
reported eGFR matched to the selected creatinine execution event; minimum and
maximum age; first qualifying paired admission, including later admissions;
and admissions through 31 July 2026. Report sample changes separately from
reclassification. First/last-creatinine checks may exceed seven days. The
reported laboratory eGFR equation and body-surface-area convention were
unspecified. There are no regression models, hypothesis tests or imputations.

The repository follows the code, variable dictionary, code-list and synthetic
example structure used for the companion adoption and dyslipidemia papers.
Code and artificial examples use the MIT License.
