
import sys
from hail import *
hc = HailContext()

sex = sys.argv[1]    

print '#####################'
print '## Starting... ######'
print '## Sex: {}'.format(sex)
print '#####################'

kt_results_auto = hc.read_table('gs://ukb31063-mega-gwas/results-tables/ukb31063.isFemale_age.{0}.autosomes.kt'.format(sex))
kt_results_chrX = hc.read_table('gs://ukb31063-mega-gwas/results-tables/ukb31063.isFemale_age.{0}.chrX.kt'.format(sex))

kt_results = KeyTable.union(kt_results_auto, kt_results_chrX)
kt_results = kt_results.annotate('variant = Variant(v.contig.replace("^0", ""), v.start, v.ref, v.alt())')
kt_results = (kt_results.key_by('variant')
                        .order_by('variant')
                        .filter('isDefined(va.AF)')
                        .drop('v'))

from pprint import pprint
pprint(kt_results.schema)

if sex == 'both_sexes':
    kt_export_sex = kt_results.annotate(['n_complete_samples = va.isFemale_results.nCompleteSamples',
                                         'AC = va.isFemale_results.AC',
                                         'ytx = va.isFemale_results.ytx[0]',
                                         'beta = va.isFemale_results.beta[0]',
                                         'se = va.isFemale_results.se[0]',
                                         'tstat = va.isFemale_results.tstat[0]',
                                         'pval = va.isFemale_results.pval[0]'])
    kt_export_sex = kt_export_sex.annotate('AF = AC/(2.0 * n_complete_samples)')
    kt_export_sex = kt_export_sex.annotate(['minor_allele = if (AF <= 0.5) variant.alt() else variant.ref',
                                            'minor_AF = if (AF <= 0.5) AF else 1.0 - AF'])
    kt_export_sex = kt_export_sex.annotate('expected_case_minor_AC = 2.0 * minor_AF * 194174')
    kt_export_sex = kt_export_sex.annotate('low_confidence_variant = (expected_case_minor_AC < 25) || (minor_AF < 0.001)')
    kt_export_sex = kt_export_sex.select(['variant',
                                          'minor_allele',
                                          'minor_AF',
                                          'expected_case_minor_AC',
                                          'low_confidence_variant',
                                          'n_complete_samples',
                                          'AC',
                                          'ytx',
                                          'beta',
                                          'se',
                                          'tstat',
                                          'pval'])
    kt_export_sex.export('gs://ukb31063-mega-gwas/export2-results-tsvs/is_female.gwas.imputed_v3.tsv.bgz')

kt_export_age = kt_results.annotate(['n_complete_samples = va.age_results.nCompleteSamples',
                                     'AC = va.age_results.AC',
                                     'ytx = va.age_results.ytx[0]',
                                     'beta = va.age_results.beta[0]',
                                     'se = va.age_results.se[0]',
                                     'tstat = va.age_results.tstat[0]',
                                     'pval = va.age_results.pval[0]'])
kt_export_age = kt_export_age.annotate('AF = AC/(2.0 * n_complete_samples)')
kt_export_age = kt_export_age.annotate(['minor_allele = if (AF <= 0.5) variant.alt() else variant.ref',
                                        'minor_AF = if (AF <= 0.5) AF else 1.0 - AF'])
kt_export_age = kt_export_age.annotate('low_confidence_variant = (minor_AF < 0.001)')
kt_export_age = kt_export_age.select(['variant',
                                      'minor_allele',
                                      'minor_AF',
                                      'low_confidence_variant',
                                      'n_complete_samples',
                                      'AC',
                                      'ytx',
                                      'beta',
                                      'se',
                                      'tstat',
                                      'pval'])
kt_export_age.export('gs://ukb31063-mega-gwas/export2-results-tsvs/age.gwas.imputed_v3.{0}.tsv.bgz'.format(sex))

