
import sys
import hail as hl

codes = {
    'albumin': '30600',
    'alkaline_phosphatase': '30610',
    'alanine_aminotransferase': '30620',
    'apoliprotein_A': '30630',
    'apoliprotein_B': '30640',
    'aspartate_aminotransferase': '30650',
    'direct_bilirubin': '30660',
    'urea': '30670',
    'calcium': '30680',
    'cholesterol': '30690',
    'creatinine': '30700',
    'C_reactive_protein': '30710',
    'cystatin_C': '30720',
    'gamma_glutamyltransferase': '30730',
    'glucose': '30740',
    'glycated_haemoglobin': '30750',
    'hdl_cholesterol': '30760',
    'igf_1': '30770',
    'ldl': '30780',
    'lipoprotein_A': '30790',
    'oestradiol': '30800',
    'phosphate': '30810',
    'rheumatoid_factor': '30820',
    'shbg': '30830',
    'total_bilirubin': '30840',
    'testosterone': '30850',
    'total_protein': '30860',
    'triglycerides': '30870',
    'urate': '30880',
    'vitamin_D': '30890',
    'estimated_sample_dilution_factor': '30897'}

descriptions = {
    'albumin': 'Albumin',
    'alkaline_phosphatase': 'Alkaline phosphatase',
    'alanine_aminotransferase': 'Alanine aminotransferase',
    'apoliprotein_A': 'Apoliprotein A',
    'apoliprotein_B': 'Apoliprotein B',
    'aspartate_aminotransferase': 'Aspartate aminotransferase',
    'direct_bilirubin': 'Direct bilirubin',
    'urea': 'Urea',
    'calcium': 'Calcium',
    'cholesterol': 'Cholesterol',
    'creatinine': 'Creatinine',
    'C_reactive_protein': 'C-reactive protein',
    'cystatin_C': 'Cystatin C',
    'gamma_glutamyltransferase': 'Gamma glutamyltransferase',
    'glucose': 'Glucose',
    'glycated_haemoglobin': 'Glycated haemoglobin',
    'hdl_cholesterol': 'HDL cholesterol',
    'igf_1': 'IGF-1',
    'ldl': 'LDL direct',
    'lipoprotein_A': 'Lipoprotein A',
    'oestradiol': 'Oestradiol',
    'phosphate': 'Phosphate',
    'rheumatoid_factor': 'Rheumatoid factor',
    'shbg': 'SHBG',
    'total_bilirubin': 'Total bilirubin',
    'testosterone': 'Testosterone',
    'total_protein': 'Total protein',
    'triglycerides': 'Triglycerides',
    'urate': 'Urate',
    'vitamin_D': 'Vitamin D',
    'estimated_sample_dilution_factor': 'Estimated sample dilution factor'}

units = {
    'albumin': 'g/L',
    'alkaline_phosphatase': 'U/L',
    'alanine_aminotransferase': 'U/L',
    'apoliprotein_A': 'g/L',
    'apoliprotein_B': 'g/L',
    'aspartate_aminotransferase': 'U/L',
    'direct_bilirubin': 'umol/L',
    'urea': 'mmol/L',
    'calcium': 'mmol/L',
    'cholesterol': 'mmol/L',
    'creatinine': 'umol/L',
    'C_reactive_protein': 'mg/L',
    'cystatin_C': 'mg/L',
    'gamma_glutamyltransferase': 'U/L',
    'glucose': 'mmol/L',
    'glycated_haemoglobin': 'mmol/mol',
    'hdl_cholesterol': 'mmol/L',
    'igf_1': 'nmol/L',
    'ldl': 'mmol/L',
    'lipoprotein_A': 'nmol/L',
    'oestradiol': 'pmol/L',
    'phosphate': 'mmol/L',
    'rheumatoid_factor': 'IU/ml',
    'shbg': 'nmol/L',
    'total_bilirubin': 'umol/L',
    'testosterone': 'nmol/L',
    'total_protein': 'g/L',
    'triglycerides': 'mmol/L',
    'urate': 'umol/L',
    'vitamin_D': 'nmol/L',
    'estimated_sample_dilution_factor': 'factor'}

sex = sys.argv[1]

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError(f'Invalid sex argument "{sex}" - must be one of {{"both_sexes", "female", "male"}}.')

ht = hl.read_table(f'gs://ukb31063/hail/ukb31063.biomarkers_gwas.{sex}.ht')

n_missing = ht.aggregate(hl.struct(**{x: hl.agg.count_where(~hl.is_defined(ht[x])) for x in ht.row_value}))
n_non_missing = ht.aggregate(hl.struct(**{x: hl.agg.count_where(hl.is_defined(ht[x])) for x in ht.row_value}))

with hl.hadoop_open(f'gs://ukb31063-mega-gwas/biomarkers/phenotype_summary.biomarkers.{sex}.tsv.bgz', 'w') as f:
    f.write('\t'.join([
        'phenotype',
        'description',
        'variable_type',
        'source',
        'n_non_missing',
        'n_missing',
        'n_controls',
        'n_cases',
        'PHESANT_transformation',
        'notes']) + '\n')
    for biomarker, code in codes.items():
        f.write('\t'.join([
            code + '_raw',
            f'{descriptions[biomarker]} ({units[biomarker]})',
            'continuous_raw',
            'biomarkers',
            str(n_non_missing[biomarker + '_raw']),
            str(n_missing[biomarker + '_raw']),
            'NA',
            'NA',
            'NA',
            'NA']) + '\n')
        f.write('\t'.join([
            code + '_irnt',
            f'{descriptions[biomarker]} (quantile)',
            'continuous_irnt',
            'biomarkers',
            str(n_non_missing[biomarker + '_irnt']),
            str(n_missing[biomarker + '_irnt']),
            'NA',
            'NA',
            'NA',
            'NA']) + '\n')
