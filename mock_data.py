"""Wholly artificial examples; no files, clinical records or distributions are read."""
import numpy as np
import pandas as pd

ID = 'encounter_study_id'
PID = 'patient_study_id'


def make_tables():
    tables = {name: [] for name in ['encounters', 'clinical_history', 'labs', 'medications']}
    for i in range(160):
        patient = i % 80
        date = pd.Timestamp(2023 if i < 80 else 2025, 8, 1) + pd.Timedelta(days=i % 20)
        if i == 159:
            date = pd.Timestamp('2026-08-01')
        keys = {ID: f'SYNTHETIC-ADMISSION-{i:03}', PID: f'SYNTHETIC-PATIENT-{patient:03}'}
        age = 40 + patient % 40
        tables['encounters'].append({**keys, 'admission_date': date, 'discharge_date': date + pd.Timedelta(days=12),
            'clinic': 'Endocrinology' if i % 7 == 0 else 'Diabetology',
            'hospitalization_status': 'ongoing' if i == 0 else 'discharged',
            'sex': 'female' if patient % 2 else 'male',
            'age_at_admission_min_years': 16 if i == 1 else age-1,
            'age_at_admission_calendar_years': 17 if i == 1 else age,
            'age_at_admission_max_years': 18 if i == 1 else age,
            'calculated_bmi_kg_m2': np.nan if i % 5 == 0 else 24 + i % 12,
            'documented_bmi_kg_m2': 28.0 if i % 7 == 0 else np.nan,
            'discharge_medication_context': 'unavailable' if i % 19 == 0 else 'discharge_recommendation'})
        concept, assertion, temporality, context = 'pregnancy', 'explicitly_absent', 'current', None
        if i == 2: assertion = 'present'
        if i == 3: assertion, temporality = 'present', 'historical'
        if i in [4,5,6]:
            concept, assertion = 'renal_replacement_therapy', 'present'
            context = {4:'current_maintenance',5:'historical',6:'acute'}[i]
        if i == 7: concept, assertion = 'kidney_transplant_history', 'present'
        tables['clinical_history'].append({**keys,'history_record_id':f'SYNTHETIC-HISTORY-{i:03}',
            'concept':concept,'assertion':assertion,'temporality':temporality,'dialysis_context':context})

        def lab(analyte, value, unit, day=1, **extra):
            time = date + pd.Timedelta(days=day, hours=8)
            row = {**keys,'lab_result_id':f'SYNTHETIC-LAB-{len(tables["labs"]):04}',
                'analyte':analyte,'value':float(value),'unit':unit,'collected_at':time,'executed_at':time,
                'encounter_link_status':'linked','result_status':'parsed','comparison_operator':'=',
                'is_preferred_for_analysis':True,'result_origin':'reported'}
            row.update(extra); tables['labs'].append(row)
        # Independent fixed values span all four renal categories.
        cr = .8 if patient % 4 < 2 else 2.0
        uacr = 12 if patient % 2 == 0 else 120
        if i % 11 != 0:
            lab('serum_creatinine_mg_dl',cr,'mg/dL',collected_at=pd.NaT if i % 17 == 0 else date+pd.Timedelta(days=1,hours=8))
            lab('egfr_reported_ml_min',90 if cr == .8 else 35,'mL/min')
            lab('serum_creatinine_mg_dl',cr*1.2,'mg/dL',10)
        if i % 13 != 0:
            lab('uacr_mg_g',0 if i == 8 else uacr,'mg/g',2 if i % 3 == 0 else 1)
            lab('uacr_mg_g',uacr,'mg/g',result_origin='derived',result_status='derived',is_preferred_for_analysis=False)
        if i % 9 != 0: lab('hba1c_percent',6+i%5,'%',2)
        lab('serum_creatinine_mg_dl',9,'mg/dL',20)  # Outside admission.
        lab('uacr_mg_g',3,'mg/g',comparison_operator='<')  # Censored, not exact.
        for ingredient, drugclass, active in [('empagliflozin / metformin','combination',i%3!=0),
                ('ramipril','ace_inhibitor',i%5!=0),('finerenone','other',i%23==0),
                ('sacubitril / valsartan','combination',i%29==0)]:
            tables['medications'].append({**keys,'medication_record_id':f'SYNTHETIC-MEDICATION-{len(tables["medications"]):04}',
                'ingredient':ingredient,'medication_class':drugclass,'medication_context':'discharge_recommendation',
                'discharge_exposure_eligible':active,'medication_use_status':'active' if active else 'stopped'})
    return {name:pd.DataFrame(rows) for name,rows in tables.items()}
