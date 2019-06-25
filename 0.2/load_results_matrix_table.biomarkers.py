
import sys
import hail as hl

sex = sys.argv[1]

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError(f'Invalid sex argument "{sex}" - must be one of {{"both_sexes", "female", "male"}}.')

try:
    dilution = sys.argv[2]
except:
    dilution = False
else:
    dilution = True

def load_results(sex, pipeline, dilution):
    if dilution:
        ht_autosomes = hl.read_table(f'gs://ukb31063-mega-gwas/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.autosomes.pipeline_{pipeline}.dilution_factor.ht')
        ht_chrX = hl.read_table(f'gs://ukb31063-mega-gwas/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.chrX.pipeline_{pipeline}.dilution_factor.ht')
    else:
        ht_autosomes = hl.read_table(f'gs://ukb31063-mega-gwas/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.autosomes.pipeline_{pipeline}.ht')
        ht_chrX = hl.read_table(f'gs://ukb31063-mega-gwas/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.chrX.pipeline_{pipeline}.ht')

    return hl.Table.union(ht_autosomes, ht_chrX)

ht_results = load_results(sex, 0, dilution)

ht_results = ht_results.annotate_globals(
    columns=hl.map(
        lambda i: hl.struct(phenotype=ht_results['phenotypes'][i]),
        hl.range(0, hl.len(ht_results['phenotypes']))))

ht_results = ht_results.annotate(entries=hl.map(
        lambda i: hl.struct(
            n=ht_results['n'][i],
            sum_x=ht_results['sum_x'][i],
            y_transpose_x=ht_results['y_transpose_x'][i][0],
            beta=ht_results['beta'][i][0],
            standard_error=ht_results['standard_error'][i][0],
            t_stat=ht_results['t_stat'][i][0],
            p_value=ht_results['p_value'][i][0]),
        hl.range(0, hl.len(ht_results['phenotypes']))))

ht_results = ht_results.select('varid', 'rsid', 'entries')
ht_results = ht_results.select_globals('columns')

for i in range(1, 7):
    ht = load_results(sex, i, dilution)
    new_phenotypes = list(ht['phenotypes'].collect()[0])

    ht_results = ht_results.annotate_globals(
        columns=ht_results['columns'].extend(hl.array([
            hl.struct(phenotype=x) for x in new_phenotypes])))

    ht_results = ht_results.annotate(
        entries=ht_results['entries'].extend(hl.rbind(
            ht[ht_results.key],
            lambda new_results: hl.array([
                hl.struct(
                    n=new_results['n'][i],
                    sum_x=new_results['sum_x'][i],
                    y_transpose_x=new_results['y_transpose_x'][i][0],
                    beta=new_results['beta'][i][0],
                    standard_error=new_results['standard_error'][i][0],
                    t_stat=new_results['t_stat'][i][0],
                    p_value=new_results['p_value'][i][0])
                for i in range(len(new_phenotypes))]))))

mt = ht_results._unlocalize_entries('entries', 'columns', ['phenotype'])

codes = hl.literal({
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
    'estimated_sample_dilution_factor': '30897'})

descriptions = hl.literal({
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
    'estimated_sample_dilution_factor': 'Estimated sample dilution factor'})

units = hl.literal({
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
    'estimated_sample_dilution_factor': 'factor'})

ht_summary = hl.import_table(f'gs://ukb31063-mega-gwas/biomarkers/phenotype_summary.biomarkers.{sex}.tsv.bgz',
    key='phenotype', impute=True)

mt = mt.annotate_cols(
    biomarker=hl.delimit(mt.phenotype.split('_')[:-1], '_'))

mt = mt.annotate_cols(
    code=codes[mt.biomarker] + '_' + mt.phenotype.split('_')[-1],
    description=hl.rbind(
        mt.phenotype.split('_')[-1],
        lambda variable_type: (
            hl.case()
              .when(variable_type == 'raw', descriptions[mt.biomarker] + ' (' + units[mt.biomarker] + ')')
              .when(variable_type == 'irnt', descriptions[mt.biomarker] + ' (quantile)')
              .default(''))))

mt = mt.annotate_cols(
    n_missing=ht_summary[mt.code]['n_missing'],
    n_non_missing=ht_summary[mt.code]['n_non_missing'])

mt = mt.repartition(500)

if dilution:
    mt.write(f'gs://ukb31063-mega-gwas/biomarkers/results-matrix-tables/ukb31063.biomarker_gwas_results.{sex}.dilution_factor.mt', overwrite=True)
else:
    mt.write(f'gs://ukb31063-mega-gwas/biomarkers/results-matrix-tables/ukb31063.biomarker_gwas_results.{sex}.mt', overwrite=True)

