"""Study derivations demonstrated on wholly artificial linked tables."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'outputs'
ID = 'encounter_study_id'
PID = 'patient_study_id'
CATEGORIES = ['G1–2/A1', 'G1–2/A2–3', 'G3–5/A1', 'G3–5/A2–3']


def first_admission(d):
    return d.sort_values(['admission_date', 'discharge_date', ID]).drop_duplicates(PID).copy()


def egfr(creatinine, age, sex):
    female = sex.eq('female')
    ratio = creatinine / np.where(female, .7, .9)
    return 142 * np.minimum(ratio, 1)**np.where(female, -.241, -.302) * np.maximum(ratio, 1)**-1.2 * .9938**age * np.where(female, 1.012, 1)


def category(gfr, uacr):
    return pd.Series(np.select([gfr.ge(60) & uacr.lt(30), gfr.ge(60) & uacr.ge(30), gfr.lt(60) & uacr.lt(30)], CATEGORIES[:3], default=CATEGORIES[3]), index=gfr.index).where(gfr.notna() & uacr.notna())


def read_release():
    from mock_data import make_tables
    return make_tables()


def renal_status(e, h):
    positive = h[h.assertion.eq('present')]
    flags = {
        'maintenance_dialysis': positive.concept.eq('renal_replacement_therapy') & positive.dialysis_context.eq('current_maintenance'),
        'transplant': positive.concept.eq('kidney_transplant_history'),
        'historical_dialysis': positive.concept.eq('renal_replacement_therapy') & positive.dialysis_context.eq('historical'),
        'acute_dialysis': positive.concept.eq('renal_replacement_therapy') & positive.dialysis_context.eq('acute'),
        'documented_aki': positive.concept.eq('acute_kidney_injury'),
    }
    for name, mask in flags.items():
        e[name] = e[ID].isin(positive.loc[mask, ID])
    unresolved = h[h.concept.isin(['renal_replacement_therapy', 'kidney_transplant_history']) & h.assertion.eq('unclear')]
    unknown_context = positive[positive.concept.eq('renal_replacement_therapy') & ~positive.dialysis_context.isin(['current_maintenance', 'historical', 'acute'])]
    e['renal_status_unresolved'] = e[ID].isin(pd.concat([unresolved, unknown_context])[ID])
    assert not e.renal_status_unresolved.any(), 'Resolve renal status before classification.'
    e['renal_excluded'] = e.maintenance_dialysis | e.transplant
    return e


def laboratory_results(e, labs):
    names = ['serum_creatinine_mg_dl', 'uacr_mg_g', 'egfr_reported_ml_min', 'hba1c_percent', 'urine_albumin_mg_l', 'urine_creatinine_mg_dl']
    q = labs[labs.analyte.isin(names)].merge(e[[ID, PID, 'admission_date', 'discharge_date', 'age', 'sex']], on=[ID, PID], how='left', validate='many_to_one')
    q['time'] = q.collected_at.fillna(q.executed_at)
    q['day'] = q.time.dt.normalize()
    q['execution_fallback'] = q.collected_at.isna()
    units = dict(zip(names, ['mg/dL', 'mg/g', 'mL/min', '%', 'mg/L', 'mg/dL']))
    q['reason'] = 'included'
    checks = [
        ('outside_eligible_admissions', q.admission_date.isna()),
        ('unlinked', ~q.encounter_link_status.eq('linked')),
        ('invalid_result', ~q.result_status.isin(['parsed', 'derived']) | q.value.isna()),
        ('censored', ~q.comparison_operator.fillna('=').eq('=')),
        ('incorrect_unit', ~q.unit.eq(q.analyte.map(units))),
        ('invalid_value', q.value.lt(0) | (q.value.eq(0) & ~q.analyte.eq('uacr_mg_g'))),
        ('missing_time', q.time.isna()),
        ('outside_admission_dates', ~q.day.between(q.admission_date, q.discharge_date)),
        ('after_discharge', q.executed_at.dt.normalize().gt(q.discharge_date)),
        ('uacr_alternative', q.analyte.eq('uacr_mg_g') & ~q.is_preferred_for_analysis.eq(True)),
    ]
    for reason, mask in checks:
        q.loc[q.reason.eq('included') & mask, 'reason'] = reason
    q.groupby(['analyte', 'reason', 'execution_fallback'], dropna=False).size().rename('records').reset_index().to_csv(OUT / 'laboratory_checks.csv', index=False)
    good = q[q.reason.eq('included')].copy()
    cr = good[good.analyte.eq('serum_creatinine_mg_dl')].copy()
    cr['egfr'] = egfr(cr.value, cr.age, cr.sex)
    uacr = good[good.analyte.eq('uacr_mg_g')].copy()
    # v27 has no preferred derived-only UACR; reported and derived alternatives stay separate.
    assert uacr.result_origin.eq('reported').all()
    assert not labs[labs.analyte.eq('uacr_mg_g') & labs.result_status.eq('parsed') & labs.value.notna()].is_preferred_for_analysis.eq(False).any()
    for field, d in [('creatinine_available', cr), ('uacr_available', uacr)]:
        e[field] = e[ID].isin(d[ID])
    e['both'] = e.creatinine_available & e.uacr_available
    e['availability'] = np.select([e.both, e.creatinine_available, e.uacr_available], ['Both', 'eGFR only', 'UACR only'], default='Neither')
    hb = good[good.analyte.eq('hba1c_percent')].sort_values(['time', 'executed_at', 'lab_result_id']).drop_duplicates(ID)
    e['hba1c'] = e[ID].map(hb.set_index(ID).value)
    return e, cr, uacr, good


def pairs(cr, uacr, maximum_days=7):
    u = uacr[[ID, PID, 'lab_result_id', 'time', 'day', 'value', 'execution_fallback']].rename(columns={'lab_result_id': 'uacr_id', 'time': 'uacr_time', 'day': 'uacr_day', 'value': 'uacr', 'execution_fallback': 'uacr_execution_fallback'})
    c = cr[[ID, PID, 'lab_result_id', 'time', 'day', 'value', 'egfr', 'executed_at', 'execution_fallback']].rename(columns={'lab_result_id': 'creatinine_id', 'time': 'creatinine_time', 'day': 'creatinine_day', 'value': 'creatinine', 'execution_fallback': 'creatinine_execution_fallback', 'executed_at': 'creatinine_executed_at'})
    p = u.merge(c, on=[ID, PID], validate='many_to_many')
    p['gap_days'] = (p.uacr_day - p.creatinine_day).dt.days.abs()
    p['gap_hours'] = (p.uacr_time - p.creatinine_time).dt.total_seconds().abs() / 3600
    # UACR ties are resolved before ranking its possible creatinine matches.
    p = p[p.gap_days.le(maximum_days)].sort_values(['uacr_time', 'uacr_id', 'gap_days', 'gap_hours', 'creatinine_time', 'creatinine_id']).drop_duplicates(ID)
    p['category'] = category(p.egfr, p.uacr)
    return p


def medications(e, meds):
    q = meds[meds[ID].isin(e[ID])].copy()
    q['active'] = q.discharge_exposure_eligible & q.medication_context.eq('discharge_recommendation')
    d = q[q.active].copy()
    assert not d.medication_use_status.isin(['stopped', 'temporarily_held', 'replaced', 'reference_only', 'conditional_replacement']).any()
    ingredient = d.ingredient.fillna('').str.lower()
    classes = d.medication_class.fillna('')
    ace = r'\b(?:captop|enalap|lisinop|perindop|ramip|quinap|trandolap|zofenop|cilazap|benazep|fosinop)'
    arb = r'\b(?:losartan|valsartan|walsartan|telmisartan|candesartan|kandesartan|irbesartan|olmesartan|eprosartan|azilsartan)'
    d['sglt2'] = ingredient.str.contains('gliflozin|gliflozyn') | classes.str.contains('sglt2')
    d['ace_arb'] = ingredient.str.contains(ace + '|' + arb, regex=True) | classes.str.contains('^ace_inhibitor|^arb|angiotensin_receptor_blocker')
    d['finerenone'] = ingredient.str.contains('fineren')
    d['sacubitril_valsartan'] = ingredient.str.contains('sacubitril') & ingredient.str.contains('valsartan')
    for field in ['sglt2', 'ace_arb', 'finerenone', 'sacubitril_valsartan']:
        e[field] = e[ID].isin(d.loc[d[field], ID])
    e['regimen_available'] = e.discharge_medication_context.eq('discharge_recommendation')
    q.groupby(['active', 'medication_use_status'], dropna=False).size().rename('records').reset_index().to_csv(OUT / 'medication_checks.csv', index=False)
    return e


def build():
    OUT.mkdir(exist_ok=True)
    t = read_release()
    e = t['encounters'].copy()
    for field in ['admission_date', 'discharge_date']:
        e[field] = pd.to_datetime(e[field])
    assert e.discharge_date.ge(e.admission_date).all()
    assert e.groupby(PID).sex.nunique().eq(1).all()
    pregnancy = set(t['clinical_history'].loc[t['clinical_history'].concept.eq('pregnancy') & t['clinical_history'].assertion.eq('present') & ~t['clinical_history'].temporality.eq('historical'), ID])
    e['eligible'] = e.hospitalization_status.eq('discharged') & e.age_at_admission_min_years.ge(18) & ~e[ID].isin(pregnancy) & e.admission_date.between('2017-01-01', '2026-09-03')
    e = e[e.eligible].copy()
    assert e.sex.isin(['female', 'male']).all()
    e['age'] = e.age_at_admission_calendar_years
    e['year'] = e.admission_date.dt.year
    e['stay_days'] = (e.discharge_date - e.admission_date).dt.days + 1
    e['bmi'] = e.calculated_bmi_kg_m2.fillna(e.documented_bmi_kg_m2)
    assert e.bmi.dropna().gt(0).all()
    e = renal_status(e, t['clinical_history'])
    e, cr, uacr, good = laboratory_results(e, t['labs'])
    e = medications(e, t['medications'])
    p7, p0 = pairs(cr, uacr), pairs(cr, uacr, 0)
    e['pair7'] = e[ID].isin(p7[ID])
    e['pair0'] = e[ID].isin(p0[ID])
    index = first_admission(e)
    main = index[~index.renal_excluded].merge(p7, on=[ID, PID], validate='one_to_one')
    flow = [
        ('Released cohort', len(t['encounters']), t['encounters'][PID].nunique()),
        ('Documented pregnancy excluded', len(pregnancy), t['encounters'].loc[t['encounters'][ID].isin(pregnancy), PID].nunique()),
        ('Eligible cohort', len(e), e[PID].nunique()),
        ('First eligible admission', len(index), len(index)),
        ('Both markers', int(index.both.sum()), int(index.both.sum())),
        ('Pair within seven days', int(index.pair7.sum()), int(index.pair7.sum())),
        ('Dialysis/transplant excluded from paired categories', int((index.pair7 & index.renal_excluded).sum()), int((index.pair7 & index.renal_excluded).sum())),
        ('Main paired categories', len(main), len(main)),
    ]
    pd.DataFrame(flow, columns=['step', 'admissions', 'patients']).to_csv(OUT / 'cohort_flow.csv', index=False)
    return e, index, main, p7, p0, cr, good
