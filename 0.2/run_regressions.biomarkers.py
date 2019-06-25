
import sys
import hail as hl

sex = sys.argv[1]
contig = sys.argv[2]
pipeline = sys.argv[3]

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError(f'Invalid sex argument "{sex}" - must be one of {{"both_sexes", "female", "male"}}.')
if contig not in set(['autosomes', 'chrX', 'chrXY']):
    raise ValueError(f'Invalid contig argument "{contig}" - must be one of {{"autosomes", "chrX", "chrXY"}}.')

try:
    dilution = sys.argv[4]
except:
    dilution = False
else:
    dilution = True

ht_phenotypes = hl.read_table(f'gs://ukb31063-mega-gwas/biomarkers/pipelines/ukb31063.biomarkers_gwas.{sex}.pipeline_{pipeline}.ht')
ht_covariates = hl.read_table(f'gs://ukb31063/hail/ukb31063.neale_gwas_covariates.{sex}.ht')
ht_variants = hl.read_table('gs://ukb31063/hail/ukb31063.neale_gwas_variants.ht')

if dilution:
    ht = hl.read_table(f'gs://ukb31063/hail/ukb31063.biomarkers_gwas.{sex}.ht')
    ht = ht.select('estimated_sample_dilution_factor_raw')
    ht_covariates = ht_covariates.annotate(
        estimated_sample_dilution_factor=ht[ht_covariates.s]['estimated_sample_dilution_factor_raw'])

if contig == 'autosomes':
    contig_expr = 'chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'
else:
    contig_expr = contig

mt = hl.import_bgen(
    path=f'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_{contig_expr}_v3.bgen',
    sample_file=f'gs://ukb31063/ukb31063.{contig}.sample',
    entry_fields=['dosage'],
    variants=ht_variants)

mt = mt.annotate_cols(
    phenotypes=ht_phenotypes[mt.s],
    covariates=ht_covariates[mt.s])

phenotypes = list(mt['phenotypes'].keys())
 
ht = hl.linear_regression_rows(
    y=[[mt['phenotypes'][y]] for y in phenotypes],
    x=mt.dosage,
    covariates=[1, *[mt['covariates'][x] for x in list(mt['covariates'].keys())]],
    pass_through=['varid', 'rsid'])

ht = ht.annotate_globals(phenotypes=phenotypes)

if dilution:
    ht.write(f'gs://ukb31063-mega-gwas/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.{contig}.pipeline_{pipeline}.dilution_factor.ht',
             overwrite=True)
else:
    ht.write(f'gs://ukb31063-mega-gwas/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.{contig}.pipeline_{pipeline}.ht',
             overwrite=True)
