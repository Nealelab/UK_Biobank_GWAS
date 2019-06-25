
import sys
from hail import *
hc = HailContext()

sex = sys.argv[1]

kt_results_auto = hc.read_table('gs://ukb31063-mega-gwas/results-tables/ukb31063.isFemale_age.{0}.autosomes.kt'.format(sex))
kt_results_chrX = hc.read_table('gs://ukb31063-mega-gwas/results-tables/ukb31063.isFemale_age.{0}.chrX.kt'.format(sex))

kt_results = KeyTable.union(kt_results_auto, kt_results_chrX)
kt_results = kt_results.annotate('variant = Variant(v.contig.replace("^0", ""), v.start, v.ref, v.alt())')
kt_results = (kt_results.key_by('variant')
                        .order_by('variant')
                        .filter('isDefined(va.AF)')
                        .drop('v'))

kt_hm3 = hc.import_table('gs://ukb31063-mega-gwas/hm3_variants.tsv.bgz', key='v_ukb')
kt_hm3 = kt_hm3.select('v_ukb')
kt_hm3 = kt_hm3.annotate('v_ukb = Variant(v_ukb.replace("^0", ""))')
kt_results = kt_results.join(kt_hm3, how='inner')
print kt_results.count()

if sex == 'both_sexes':
    kt_export_sex = kt_results.annotate(['SNP = va.rsid',
                                         'A1 = variant.alt()',
                                         'A2 = variant.ref',
                                         'N = va.isFemale_results.nCompleteSamples',
                                         'Z = va.isFemale_results.tstat[0]'])
    kt_export_sex = kt_export_sex.select(['SNP', 'A1', 'A2', 'N', 'Z'])
    kt_export_sex.export('gs://ukb31063-mega-gwas/ldsc-export/sumstats-files/is_female.ldsc.imputed_v3.tsv.bgz')

kt_export_age = kt_results.annotate(['SNP = va.rsid',
                                     'A1 = variant.alt()',
                                     'A2 = variant.ref',
                                     'N = va.age_results.nCompleteSamples',
                                     'Z = va.age_results.tstat[0]'])
kt_export_age = kt_export_age.select(['SNP', 'A1', 'A2', 'N', 'Z'])
kt_export_age.export('gs://ukb31063-mega-gwas/ldsc-export/sumstats-files/age.ldsc.imputed_v3.{0}.tsv.bgz'.format(sex))
