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



ht_hm3 = hl.read_table('gs://ukbb-ldsc-dev/ukb_hm3_snplist/hm3.r3.b37.auto_bi_af.ukbb_gwas_qcpos.no_mhc.ht')
ht_hm3 = ht_hm3.select('variant').key_by('variant').drop('rsid','locus','alleles')
ht_join = ht_results.join(ht_hm3, how='inner')
ht_join = ht_join.annotate(SNP=ht_join['rsid']).key_by('SNP').cache()



count = 1
phenotypes = ht_join['phenotypes'].collect()[0]
for i, phenotype in enumerate(phenotypes):
    variable_type = phenotype.split('_')[-1]
    code = codes[phenotype.replace('_raw', '').replace('_irnt', '')]
    print(f'Exporting LDSC sumstats for trait {code} ({count})...')
    ht_export = ht_join.annotate(
        A1 = ht_join.alleles[1],
        A2 = ht_join.alleles[0],
        N = ht_join['n'][i],
        Z = ht_join['t_stat'][i][0])
    ht_export = ht_export.select('A1','A2','N','Z')
    if dilution:
        ht_export.export(f'gs://ukb-mega-gwas-results/round2/additive-ldsc-sumstats/biomarkers/{code}_{variable_type}.imputed_v3.ldsc.dilution_factor.{sex}.tsv.gz')
    else:
        ht_export.export(f'gs://ukb-mega-gwas-results/round2/additive-ldsc-sumstats/biomarkers/{code}_{variable_type}.imputed_v3.ldsc.{sex}.tsv.gz')
    count += 1



print('#######################')
print('## COMPLETED     ######')
print('## Sex: {}'.format(sex))
print('## Pipeline number: {}'.format(pipeline))
print('## Dilution: {}'.format(dilution))
print('#######################')
