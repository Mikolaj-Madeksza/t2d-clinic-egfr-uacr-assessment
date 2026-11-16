"""Independent checks of cohort membership, pairing, estimates and exclusions."""
import json
import math
import pandas as pd
from derive import OUT, ID, PID, build
from mock_data import make_tables


def run():
    source = make_tables()
    e, index, main, p7, p0, cr, good = build()
    raw = source['encounters'].copy()
    for field in ['admission_date','discharge_date']:
        raw[field] = pd.to_datetime(raw[field])
    h = source['clinical_history']
    pregnant = set(h.loc[h.concept.eq('pregnancy') & h.assertion.eq('present') & h.temporality.ne('historical'), ID])
    eligible = raw[raw.age_at_admission_min_years.ge(18) & raw.hospitalization_status.eq('discharged') & ~raw[ID].isin(pregnant) & raw.admission_date.between('2017-01-01','2026-09-03')]
    expected_index = eligible.sort_values(['admission_date','discharge_date',ID]).groupby(PID,sort=False).head(1)
    assert len(e)==157 and len(index)==80 and len(main)==65
    assert 'SYNTHETIC-ADMISSION-002' not in set(e[ID])
    assert 'SYNTHETIC-ADMISSION-003' in set(e[ID])
    assert set(expected_index[ID]) == set(index[ID])
    assert index[PID].is_unique and main[PID].is_unique
    l = source['labs'].set_index('lab_result_id',drop=False)
    independent_both = 0
    expected_pairs = {}
    for _,enc in expected_index.iterrows():
        q = l[l[ID].eq(enc[ID]) & l[PID].eq(enc[PID])].copy()
        q['t'] = q.collected_at.combine_first(q.executed_at)
        q = q[q.encounter_link_status.eq('linked') & q.result_status.eq('parsed') & q.value.notna() & q.comparison_operator.fillna('=').eq('=') & q.t.dt.normalize().between(enc.admission_date,enc.discharge_date) & (q.executed_at.isna() | q.executed_at.dt.normalize().le(enc.discharge_date))]
        c = q[q.analyte.eq('serum_creatinine_mg_dl') & q.unit.eq('mg/dL') & q.value.gt(0)]
        u = q[q.analyte.eq('uacr_mg_g') & q.unit.eq('mg/g') & q.value.ge(0) & q.is_preferred_for_analysis.eq(True)]
        independent_both += bool(len(c) and len(u))
        candidates = []
        for ur in u.itertuples(index=False):
            for cr in c.itertuples(index=False):
                gap = abs((ur.t.normalize()-cr.t.normalize()).days)
                if gap <= 7:
                    candidates.append((ur.t,ur.lab_result_id,gap,abs((ur.t-cr.t).total_seconds()),cr.t,cr.lab_result_id))
        if candidates:
            chosen = min(candidates)
            expected_pairs[enc[ID]] = (chosen[1],chosen[5])
    assert independent_both == int(index.both.sum())
    assert set(expected_pairs) == set(index.loc[index.pair7,ID])
    pos = h[h.assertion.eq('present')]
    renal = set(pos.loc[pos.concept.eq('kidney_transplant_history') | (pos.concept.eq('renal_replacement_therapy') & pos.dialysis_context.eq('current_maintenance')),ID])
    assert set(main[ID]) == set(expected_pairs)-renal
    for row in main.itertuples(index=False):
        assert expected_pairs[getattr(row,ID)] == (row.uacr_id,row.creatinine_id)
        cr = l.loc[row.creatinine_id,'value']
        # Independent piecewise form of CKD-EPI 2021.
        if row.sex == 'female':
            k, alpha, multiplier = .7,-.241,1.012
        else:
            k, alpha, multiplier = .9,-.302,1
        exponent = alpha if cr <= k else -1.2
        expected = 142*(cr/k)**exponent * .9938**row.age * multiplier
        assert math.isclose(expected,row.egfr,rel_tol=1e-12)
    estimates = pd.read_csv(OUT/'availability_estimates.csv')
    for r in estimates.itertuples():
        # Wilson intervals cross-checked by inversion of the binomial score statistic.
        for bound in [r.ci_low/100,r.ci_high/100]:
            if 0 < bound < 1:
                score=(r.n-r.N*bound)**2/(r.N*bound*(1-bound))
                assert math.isclose(score,1.959963984540054**2,rel_tol=1e-9)
    meds = e[e.clinic.eq('Diabetology') & e.admission_date.ge('2023-07-14') & e.regimen_available & ~e.renal_excluded].merge(p7,on=[ID,PID])
    assert meds[ID].is_unique and not meds[ID].isin(renal).any()
    assert meds.clinic.eq('Diabetology').all() and meds.admission_date.ge('2023-07-14').all()
    active = source['medications']
    active = active[active.discharge_exposure_eligible & active.medication_context.eq('discharge_recommendation')]
    for drug, ingredients in [('sglt2',['empagliflozin / metformin']),('ace_arb',['ramipril','sacubitril / valsartan']),('finerenone',['finerenone'])]:
        expected = set(active.loc[active.ingredient.isin(ingredients),ID])
        assert set(meds.loc[meds[drug],ID]) == set(meds[ID]) & expected
    boot = pd.read_csv(OUT/'medication_estimates.csv')
    assert boot.valid_bootstraps.eq(10000).all()
    assert boot.ci_low.between(0,100).all() and boot.ci_high.between(0,100).all() and boot.ci_low.le(boot.ci_high).all()
    changes = pd.read_csv(OUT/'sensitivity_changes.csv')
    assert changes['Changed category'].le(changes['Common patients']).all()
    assert (changes.Patients == len(main)+changes['Added patients']-changes['Removed patients']).all()
    result = {'status':'PASS','independent_index_patients':len(index),'independent_joint_availability':independent_both,
              'independently_checked_pairs':len(main),'medication_denominators':'PASS','sensitivity_denominators':'PASS'}
    (OUT/'validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    run()
