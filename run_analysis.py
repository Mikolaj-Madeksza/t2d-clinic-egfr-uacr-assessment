"""Descriptive estimates, planned sensitivities and synthetic example tables."""
import json
import numpy as np
import pandas as pd
from derive import OUT, ID, PID, CATEGORIES, build, category, egfr, first_admission

SEED = 20261115
BOOTSTRAPS = 10000


def wilson(k, n):
    if n == 0:
        return np.nan, np.nan
    z = 1.959963984540054
    p = k / n
    den = 1 + z*z/n
    centre = (p + z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return 100*(centre-half), 100*(centre+half)


def proportion(k, n):
    return f'{k}/{n} ({100*k/n:.1f}%)' if n else '—'


def interval(k, n):
    lo, hi = wilson(k, n)
    return f'{lo:.1f}–{hi:.1f}' if n else '—'


def continuous(s):
    x = s.dropna()
    if x.empty:
        return '—; n=0'
    q = x.quantile([.25, .5, .75])
    return f'{q.loc[.5]:.1f} [{q.loc[.25]:.1f}–{q.loc[.75]:.1f}]; n={len(x)}'


def characteristics(groups):
    rows = []
    specs = [('Patients', lambda d: str(len(d))), ('Age, years', lambda d: continuous(d.age)),
             ('Female', lambda d: proportion(int(d.sex.eq('female').sum()), len(d))),
             ('Diabetology', lambda d: proportion(int(d.clinic.eq('Diabetology').sum()), len(d))),
             ('Endocrinology', lambda d: proportion(int(d.clinic.eq('Endocrinology').sum()), len(d))),
             ('Admission year', lambda d: f'{d.year.median():.0f} [{d.year.quantile(.25):.0f}–{d.year.quantile(.75):.0f}]'),
             ('Hospital stay, days', lambda d: continuous(d.stay_days)),
             ('BMI, kg/m²', lambda d: continuous(d.bmi)),
             ('BMI missing', lambda d: proportion(int(d.bmi.isna().sum()), len(d))),
             ('HbA1c, %', lambda d: continuous(d.hba1c)),
             ('HbA1c missing', lambda d: proportion(int(d.hba1c.isna().sum()), len(d)))]
    for label, fn in specs:
        rows.append({'Characteristic': label, **{name: fn(d) for name, d in groups.items()}})
    return pd.DataFrame(rows)


def phenotype_rows(d, label):
    rows = []
    for cat in CATEGORIES:
        g = d[d.category.eq(cat)]
        rows.append({'analysis': label, 'category': cat, 'n': len(g), 'N': len(d), 'percent': 100*len(g)/len(d) if len(d) else np.nan,
                     'ci_low': wilson(len(g), len(d))[0], 'ci_high': wilson(len(g), len(d))[1]})
    return rows


def cluster_intervals(d):
    patients = sorted(d[PID].unique())
    # Resample patients from the whole subgroup, preserving all their admissions and categories.
    rng = np.random.default_rng(SEED)
    weights = rng.multinomial(len(patients), np.repeat(1/len(patients), len(patients)), size=BOOTSTRAPS)
    rows = []
    for cat in CATEGORIES:
        g = d[d.category.eq(cat)]
        den = g.groupby(PID).size().reindex(patients, fill_value=0).to_numpy()
        sample_den = weights @ den
        for drug in ['sglt2', 'ace_arb']:
            num = g.groupby(PID)[drug].sum().reindex(patients, fill_value=0).to_numpy()
            sample_num = weights @ num
            valid = sample_den > 0
            lo, hi = np.quantile(100*sample_num[valid]/sample_den[valid], [.025, .975])
            rows.append({'category': cat, 'drug': drug, 'n': int(g[drug].sum()), 'N': len(g), 'patients': g[PID].nunique(),
                         'percent': 100*g[drug].mean(), 'ci_low': lo, 'ci_high': hi, 'valid_bootstraps': int(valid.sum())})
    return pd.DataFrame(rows)


def sensitivity_analyses(e, index, main, p7, p0, cr, good):
    variants = {'Primary': main.copy()}
    variants['Same-day pair'] = index[~index.renal_excluded].merge(p0, on=[ID, PID], validate='one_to_one')
    for which in ['first', 'last']:
        selected = cr.sort_values(['time', 'lab_result_id']).drop_duplicates(ID, keep=which)
        d = main.drop(columns=['egfr', 'creatinine', 'creatinine_time', 'creatinine_id', 'gap_days']).merge(
            selected[[ID, 'egfr', 'value', 'time', 'lab_result_id']], on=ID, validate='one_to_one').rename(columns={'value':'creatinine','time':'creatinine_time','lab_result_id':'creatinine_id'})
        d['gap_days'] = (d.uacr_time.dt.normalize() - d.creatinine_time.dt.normalize()).dt.days.abs()
        d['category'] = category(d.egfr, d.uacr)
        variants[which.title() + ' creatinine'] = d
    # Match the reported estimate to the selected creatinine's execution event.
    reported = good[good.analyte.eq('egfr_reported_ml_min')].copy()
    keys = [ID, PID, 'executed_at']
    conflicts = reported.groupby(keys).value.nunique()
    assert not conflicts.gt(1).any(), 'Conflicting reported eGFR at one execution event.'
    reported = reported.sort_values('lab_result_id').drop_duplicates(keys)
    d = main.merge(reported[keys + ['value', 'lab_result_id']], left_on=[ID, PID, 'creatinine_executed_at'], right_on=keys, validate='one_to_one')
    d['egfr'] = d.value
    d['category'] = category(d.egfr, d.uacr)
    variants['Reported eGFR'] = d
    for field, label in [('age_at_admission_min_years', 'Minimum age'), ('age_at_admission_max_years', 'Maximum age')]:
        d = main.copy()
        d['egfr'] = egfr(d.creatinine, d[field], d.sex)
        d['category'] = category(d.egfr, d.uacr)
        variants[label] = d
    all_pairs = e[~e.renal_excluded].merge(p7, on=[ID, PID], validate='one_to_one')
    variants['First qualifying admission'] = first_admission(all_pairs)
    distribution, changes, medians = [], [], []
    for label, d in variants.items():
        distribution.extend(phenotype_rows(d, label))
        common = main[[PID, ID, 'category', 'egfr', 'uacr']].merge(d[[PID, ID, 'category', 'egfr', 'uacr']], on=PID, suffixes=('_primary', '_alternative'), validate='one_to_one')
        changes.append({'Analysis': label, 'Patients': len(d), 'Common patients': len(common),
                        'Changed category': int(common.category_primary.ne(common.category_alternative).sum()),
                        'Added patients': int((~d[PID].isin(main[PID])).sum()),
                        'Removed patients': int((~main[PID].isin(d[PID])).sum()),
                        'Different admission': int(common[ID+'_primary'].ne(common[ID+'_alternative']).sum()),
                        'Pairs over 7 days': int(d.gap_days.gt(7).sum())})
        medians.append({'Analysis': label, 'eGFR': continuous(d.egfr), 'UACR': continuous(d.uacr)})
    pd.DataFrame(distribution).to_csv(OUT / 'sensitivity_estimates.csv', index=False)
    pd.DataFrame(changes).to_csv(OUT / 'sensitivity_changes.csv', index=False)
    pd.DataFrame(medians).to_csv(OUT / 'sensitivity_continuous.csv', index=False)
    return variants, pd.DataFrame(distribution), pd.DataFrame(changes)


def run():
    e, index, main, p7, p0, cr, good = build()
    tables = []
    def add(key, title, d, note, abbreviations=''):
        d.to_csv(OUT / (key + '.csv'), index=False)
        tables.append({'key': key, 'title': 'SYNTHETIC EXAMPLE — ' + title, 'columns': list(d.columns), 'rows': d.fillna('—').astype(str).values.tolist(), 'note': note, 'abbreviations': abbreviations})

    t1 = characteristics({'All patients': index, 'Both markers': index[index.both], 'Without both': index[~index.both]})
    availability = []
    for label, mask in [('Creatinine sufficient for eGFR', index.creatinine_available), ('UACR', index.uacr_available), ('Both markers', index.both), ('eGFR only', index.availability.eq('eGFR only')), ('UACR only', index.availability.eq('UACR only')), ('Neither', index.availability.eq('Neither')), ('Pair within 7 days', index.pair7), ('Same-day pair', index.pair0)]:
        k = int(mask.sum()); availability.append({'Measure': label, 'n': k, 'N': len(index), 'percent': 100*k/len(index), 'ci_low':wilson(k,len(index))[0], 'ci_high':wilson(k,len(index))[1]})
    pd.DataFrame(availability).to_csv(OUT / 'availability_estimates.csv', index=False)
    extra = []
    for a in availability:
        extra.append({'Characteristic': a['Measure'], 'All patients': proportion(a['n'], a['N']), 'Both markers':'—','Without both':'—'})
    for clinic, g in index.groupby('clinic'):
        extra.append({'Characteristic': 'Both markers: '+clinic, 'All patients': proportion(int(g.both.sum()),len(g)), 'Both markers':'—','Without both':'—'})
    t1 = pd.concat([t1, pd.DataFrame(extra)], ignore_index=True)
    add('table_1', 'Table 1. Patient characteristics and kidney-marker availability at the index admission', t1,
        f'Values are n/N (%) or median [IQR], with the number of available observations. All eligible patients are included, including those receiving maintenance dialysis or with a kidney transplant. Both markers were available in {proportion(int(index.both.sum()), len(index))}; 95% CI {interval(int(index.both.sum()), len(index))}%. Clinical characteristics are not required for cohort inclusion. BMI uses released height/weight calculations where available, otherwise the documented value. HbA1c is the first usable inpatient result. Counts refer to patients.',
        'BMI, body mass index; CI, confidence interval; eGFR, estimated glomerular filtration rate; HbA1c, glycated hemoglobin; IQR, interquartile range; UACR, urine albumin-to-creatinine ratio')

    rows = []
    for cat in CATEGORIES:
        g = main[main.category.eq(cat)]
        rows.append({'Category or measure': cat, 'Patients, n/N (%)':proportion(len(g),len(main)), '95% CI, %':interval(len(g),len(main)), 'eGFR, median [IQR]':continuous(g.egfr).split(';')[0], 'UACR, median [IQR]':continuous(g.uacr).split(';')[0]})
    rows.append({'Category or measure':'All paired patients','Patients, n/N (%)':proportion(len(main),len(main)),'95% CI, %':'—','eGFR, median [IQR]':continuous(main.egfr).split(';')[0],'UACR, median [IQR]':continuous(main.uacr).split(';')[0]})
    outcomes = [('At least one abnormal marker', (main.egfr.lt(60)|main.uacr.ge(30)), pd.Series(True,index=main.index)),
                ('UACR ≥30 among eGFR ≥60', main.uacr.ge(30)&main.egfr.ge(60), main.egfr.ge(60)),
                ('UACR <30 among eGFR <60', main.uacr.lt(30)&main.egfr.lt(60), main.egfr.lt(60))]
    estimates = []
    for label, numerator, denominator in outcomes:
        k,n=int(numerator.sum()),int(denominator.sum())
        rows.append({'Category or measure':label,'Patients, n/N (%)':proportion(k,n),'95% CI, %':interval(k,n),'eGFR, median [IQR]':'—','UACR, median [IQR]':'—'})
        estimates.append({'measure':label,'n':k,'N':n,'percent':100*k/n,'ci_low':wilson(k,n)[0],'ci_high':wilson(k,n)[1]})
    pd.DataFrame(estimates).to_csv(OUT / 'complementarity_estimates.csv',index=False)
    add('table_2', 'Table 2. Renal-marker categories among patients with paired measurements', pd.DataFrame(rows),
        'The earliest UACR with a creatinine result within seven calendar days was selected; the nearest creatinine was used. eGFR was calculated with CKD-EPI 2021. Units are mL/min/1.73 m² for eGFR and mg/g for UACR. G1–2 denotes eGFR ≥60; G3–5 denotes eGFR <60; A1 denotes UACR <30; A2–3 denotes UACR ≥30. Current maintenance dialysis and kidney transplantation were excluded. Confidence intervals use the Wilson method. These are single-admission marker categories, not confirmed CKD.',
        'CI, confidence interval; CKD, chronic kidney disease; CKD-EPI, Chronic Kidney Disease Epidemiology Collaboration; eGFR, estimated glomerular filtration rate; IQR, interquartile range; UACR, urine albumin-to-creatinine ratio')

    definitions = [
        ['Eligibility','Adults with adjudicated T2D, discharged admissions beginning 1 January 2017–3 September 2026; current pregnancy excluded. Historical pregnancy alone does not exclude.'],
        ['Index admission','First eligible admission by admission date, discharge date and encounter identifier.'],
        ['Availability','Usable linked serum creatinine and/or quantitative UACR during the index admission. Missing covariates do not exclude patients.'],
        ['Laboratory timing','Collection time first; execution time if collection time is missing. Dates must be within the admission; no result executed after discharge is used.'],
        ['UACR','Preferred numeric reported result in mg/g. Reported and derived alternatives are not counted twice. No derived-only UACR was selected. Reported zero is retained; censored results are not treated as exact.'],
        ['Pair selection','Earliest UACR with a creatinine within seven calendar days. Ties: UACR identifier; then creatinine day gap, absolute time gap, timestamp and identifier.'],
        ['eGFR','CKD-EPI 2021 creatinine equation; age is the released calendar-year age at admission. Minimum/maximum age bounds are assessed separately.'],
        ['Renal exclusions','Current maintenance dialysis or kidney transplant history documented during the admission. Historical or acute dialysis alone does not exclude.'],
        ['BMI and HbA1c','BMI calculated from released index height/weight, otherwise documented BMI. First usable inpatient HbA1c. No imputation.'],
        ['Drug recommendations','Active eligible discharge mentions. Ingredients and combination products identify SGLT2 inhibitors, ACE inhibitors/ARBs and finerenone. ARB-containing sacubitril/valsartan is included and counted separately in the medication table note.'],
        ['Intervals','Wilson 95% confidence intervals for patient-level proportions. Medication intervals: 10,000 patient-clustered bootstrap samples, seed 20261115, percentile limits.'],
        ['Reported eGFR','Matched to the selected creatinine by patient, admission and execution timestamp. The laboratory equation and BSA-indexing convention are not specified in the release; this is a sensitivity comparison.'],
    ]
    add('table_s1','Table S1. Analysis definitions',pd.DataFrame(definitions,columns=['Measure','Definition']),
        'Synthetic demonstration of the v27 analysis rules. Prior adjudicated diabetes classifications were retained; structured E10 code assertions did not automatically override explicit T2D documentation.',
        'ACE, angiotensin-converting enzyme; ARB, angiotensin receptor blocker; BMI, body mass index; BSA, body surface area; CKD-EPI, Chronic Kidney Disease Epidemiology Collaboration; eGFR, estimated glomerular filtration rate; HbA1c, glycated hemoglobin; SGLT2, sodium-glucose cotransporter 2; T2D, type 2 diabetes; UACR, urine albumin-to-creatinine ratio')

    rows = []
    for (clinic, year),g in index.groupby(['clinic','year']):
        rows.append({'Clinic':clinic,'Year':str(year),'Patients':len(g),'eGFR available':proportion(int(g.creatinine_available.sum()),len(g)),'UACR available':proportion(int(g.uacr_available.sum()),len(g)),'Both markers':proportion(int(g.both.sum()),len(g))})
    add('table_s2','Table S2. Kidney-marker availability by clinic and year',pd.DataFrame(rows),
        'Values are n/N (%). Each patient contributes only their first eligible admission in the entire cohort. Diabetology begins on 14 July 2023; 2026 includes admissions beginning through 3 September. Clinic comparisons also reflect different periods and patient populations.',
        'eGFR, estimated glomerular filtration rate; UACR, urine albumin-to-creatinine ratio')

    main['G'] = pd.cut(main.egfr,[-np.inf,15,30,45,60,90,np.inf],labels=['G5','G4','G3b','G3a','G2','G1'],right=False)
    main['A'] = pd.cut(main.uacr,[-np.inf,30,300,np.inf],labels=['A1','A2','A3'],right=False)
    ga = pd.crosstab(main.G,main.A,dropna=False).reindex(['G1','G2','G3a','G3b','G4','G5'])
    rows = []
    for g in ga.index:
        rows.append({'eGFR category':g,**{a:proportion(int(ga.loc[g,a]),len(main)) for a in ['A1','A2','A3']},'Total':proportion(int(ga.loc[g].sum()),len(main))})
    add('table_s3','Table S3. Detailed eGFR and UACR categories',pd.DataFrame(rows),
        f'Values are n/N (%), with N={len(main)} for every cell. eGFR thresholds: G1 ≥90, G2 60–<90, G3a 45–<60, G3b 30–<45, G4 15–<30, G5 <15 mL/min/1.73 m². UACR thresholds: A1 <30, A2 30–<300, A3 ≥300 mg/g. Current maintenance dialysis and kidney transplantation are excluded. Single-admission categories do not establish CKD.',
        'CKD, chronic kidney disease; eGFR, estimated glomerular filtration rate; UACR, urine albumin-to-creatinine ratio')

    excluded=index[index.pair7&index.renal_excluded]

    rows=[]
    for field,label in [('maintenance_dialysis','Current maintenance dialysis'),('transplant','Kidney transplant history'),('renal_excluded','Either exclusion'),('historical_dialysis','Historical dialysis'),('acute_dialysis','Acute dialysis'),('renal_status_unresolved','Unresolved status')]:
        g=index[index[field]]
        rows.append({'Documented status':label,'Index patients':len(g),'Both markers':int(g.both.sum()),'Pair within 7 days':int(g.pair7.sum())})
    add('table_s4','Table S4. Dialysis and kidney transplant status at the index admission',pd.DataFrame(rows),
        'Counts refer to distinct patients. Rows can overlap. Status is based on affirmative, encounter-specific source assertions. Historical and acute dialysis do not by themselves establish current maintenance dialysis. Transplant history includes failed grafts.', '')

    recent=e[e.clinic.isin(['Diabetology','Endocrinology'])&e.admission_date.ge('2023-01-01')&e.regimen_available&~e.renal_excluded].merge(p7,on=[ID,PID],validate='one_to_one')
    boot=cluster_intervals(recent);boot.to_csv(OUT/'medication_estimates.csv',index=False)
    rows=[]
    for cat in CATEGORIES:
        g=recent[recent.category.eq(cat)];row={'Category':cat,'Admissions / patients':f'{len(g)} / {g[PID].nunique()}'}
        for drug,label in [('sglt2','SGLT2 inhibitor'),('ace_arb','ACEi/ARB')]:
            b=boot[boot.category.eq(cat)&boot.drug.eq(drug)].iloc[0]
            row[label]=proportion(int(b.n),int(b.N))+f'; CI {b.ci_low:.1f}–{b.ci_high:.1f}'
        row['Finerenone, n']=int(g.finerenone.sum());rows.append(row)
    add('table_s5','Table S5. Discharge medication recommendations by renal-marker category',pd.DataFrame(rows),
        f'Endocrinology and Diabetology admissions from 1 January 2023 to 3 September 2026 with paired renal markers and an available discharge regimen; maintenance dialysis and kidney transplantation excluded. Total: {len(recent)} admissions from {recent[PID].nunique()} patients. Values are n/N (%) with patient-clustered 95% bootstrap confidence intervals; finerenone is a count. Patients may contribute to more than one category. ACEi/ARB includes sacubitril/valsartan ({int(recent.sacubitril_valsartan.sum())} admissions).',
        'ACEi, angiotensin-converting enzyme inhibitor; ARB, angiotensin receptor blocker; CI, confidence interval; SGLT2, sodium-glucose cotransporter 2')

    variants,dist,changes=sensitivity_analyses(e,index,main,p7,p0,cr,good)
    rows=[]
    for label,d in variants.items():
        rows.append({'Analysis':label,'N':len(d),**{cat:f'{int(d.category.eq(cat).sum())} ({100*d.category.eq(cat).mean():.1f}%)' for cat in CATEGORIES}})
    add('table_s6','Table S6. Renal-marker categories in sensitivity analyses',pd.DataFrame(rows),
        'Values are n (% within each analysis). The same renal exclusions apply throughout. First/last creatinine analyses keep the primary UACR and admission fixed. The first qualifying admission may be later than the index admission. Reported eGFR uses the common set with a matched reported result; its original equation is unspecified. Age bounds use the released minimum and maximum ages.',
        'eGFR, estimated glomerular filtration rate; UACR, urine albumin-to-creatinine ratio')
    add('table_s7','Table S7. Sample changes and reclassification in sensitivity analyses',changes.drop(columns=['Patients','Different admission']),
        'Changed category is assessed only among patients common to the primary and alternative analysis. Added/removed patients are relative to the primary analysis. Pairs over seven days occur only in the alternative first/last-creatinine checks. Changes in sample composition are separate from reclassification.', '')
    (OUT/'tables.json').write_text(json.dumps(tables,ensure_ascii=False,indent=2))
    summary={'data':'SYNTHETIC ONLY','admissions':len(e),'patients':len(index),'both':int(index.both.sum()),'both_ci':wilson(int(index.both.sum()),len(index)),
             'paired_before_exclusions':int(index.pair7.sum()),'paired_exclusions':len(excluded),'main_pairs':len(main),'categories':main.category.value_counts().to_dict(),
             'medication_admissions':len(recent),'medication_patients':recent[PID].nunique(),'finerenone':int(recent.finerenone.sum()),
             'execution_fallback_uacr':int(main.uacr_execution_fallback.sum()),'execution_fallback_creatinine':int(main.creatinine_execution_fallback.sum()),
             'reported_zero_uacr':int(main.uacr.eq(0).sum()),'seed':SEED,'bootstrap_samples':BOOTSTRAPS}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,default=lambda x:int(x) if isinstance(x,np.integer) else str(x)))
    print(json.dumps(summary,indent=2,default=lambda x:int(x) if isinstance(x,np.integer) else str(x)))
    return summary


if __name__ == '__main__':
    run()
