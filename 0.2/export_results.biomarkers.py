
import sys
import hail as hl
from hail.utils import new_temp_file

sex = sys.argv[1]
pipeline = sys.argv[2]

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError(f'Invalid sex argument "{sex}" - must be one of {{"both_sexes", "female", "male"}}.')

try:
    dilution = sys.argv[3]
except:
    dilution = False
else:
    dilution = True

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

if dilution:
    ht_autosomes = hl.read_table(f'gs://ukb31063-mega-gwas/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.autosomes.pipeline_{pipeline}.dilution_factor.ht')
    ht_chrX = hl.read_table(f'gs://ukb31063-mega-gwas/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.chrX.pipeline_{pipeline}.dilution_factor.ht')
else:
    ht_autosomes = hl.read_table(f'gs://ukb31063-mega-gwas/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.autosomes.pipeline_{pipeline}.ht')
    ht_chrX = hl.read_table(f'gs://ukb31063-mega-gwas/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.chrX.pipeline_{pipeline}.ht')

ht_results = hl.Table.union(ht_autosomes, ht_chrX)
ht_results = ht_results.annotate(variant=hl.delimit(hl.array([
    ht_results['locus'].contig,
    hl.str(ht_results['locus'].position),
    ht_results['alleles'][0],
    ht_results['alleles'][1]]), delimiter=':'))
ht_results = ht_results.key_by('variant')
ht_results = ht_results.repartition(116)
ht_results = ht_results.cache()

phenotypes = ht_results['phenotypes'].collect()[0]
for i, phenotype in enumerate(phenotypes):
    variable_type = phenotype.split('_')[-1]
    code = codes[phenotype.replace('_raw', '').replace('_irnt', '')]
    ht_export = ht_results.annotate(
        n_complete_samples=ht_results['n'][i],
        AC=ht_results['sum_x'][i],
        ytx=ht_results['y_transpose_x'][i][0],
        beta=ht_results['beta'][i][0],
        se=ht_results['standard_error'][i][0],
        tstat=ht_results['t_stat'][i][0],
        pval=ht_results['p_value'][i][0])
    ht_export = ht_export.annotate(
        AF=ht_export['AC'] / (2 * ht_export['n_complete_samples']))
    ht_export = ht_export.annotate(
        minor_AF=hl.cond(ht_export['AF'] <= 0.5, ht_export['AF'], 1.0 - ht_export['AF']),
        minor_allele=hl.cond(ht_export['AF'] <= 0.5, ht_export['alleles'][1], ht_export['alleles'][0]))
    ht_export = ht_export.annotate(
        low_confidence_variant=ht_export['minor_AF'] < 0.001)
    ht_export = ht_export.select(
        'minor_allele',
        'minor_AF',
        'low_confidence_variant',
        'n_complete_samples',
        'AC',
        'ytx',
        'beta',
        'se',
        'tstat',
        'pval')
    if dilution:
        ht_export.export(f'gs://ukb31063-mega-gwas/biomarkers/results-tsvs/{code}_{variable_type}.gwas.imputed_v3.dilution_factor.{sex}.tsv.bgz')   
    else:
        ht_export.export(f'gs://ukb31063-mega-gwas/biomarkers/results-tsvs/{code}_{variable_type}.gwas.imputed_v3.{sex}.tsv.bgz')

