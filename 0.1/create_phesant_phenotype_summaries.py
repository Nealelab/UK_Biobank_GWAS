

from hail import *
hc = HailContext()

# both_sexes
kt_cat = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/still-more-phesant/neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_both_sexes.*_phenosummary.tsv',
                         comment='Field', no_header=True, missing='NA')
kt_cat = kt_cat.rename({'f0': 'phenotype',
                        'f1': 'description',
                        'f2': 'n_non_missing',
                        'f3': 'n_missing',
                        'f4': 'n_cases',
                        'f5': 'n_controls',
                        'f6': 'notes',
                        'f7': 'PHESANT_transformation'})
kt_cat = kt_cat.annotate('type = "cat"')
kt_cat = kt_cat.select(['phenotype',
                        'description',
                        'notes',
                        'PHESANT_transformation',
                        'n_non_missing',
                        'n_missing',
                        'n_controls',
                        'n_cases',
                        'type'])
cat_phenotypes = kt_cat.query('phenotype.collect()')

kt_cts = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb11214_final_july_reference_QC_more_phenos_and_corrected_restricted*_phenosummary.tsv',
                         comment='Field', no_header=True, missing='')
kt_cts = kt_cts.rename({'f0': 'phenotype',
                        'f1': 'description',
                        'f2': 'n_non_missing',
                        'f3': 'n_missing',
                        'f4': 'n_cases',
                        'f5': 'n_controls',
                        'f6': 'notes',
                        'f7': 'PHESANT_transformation'})
kt_cts = kt_cts.annotate('type = "cts"')
kt_cts = kt_cts.select(['phenotype',
                        'description',
                        'notes',
                        'PHESANT_transformation',
                        'n_non_missing',
                        'n_missing',
                        'n_controls',
                        'n_cases',
                        'type'])
kt_cts = kt_cts.filter('let c = ["{}"].toSet() in c.contains(phenotype)'.format('","'.join(cat_phenotypes)), keep=False)

kt = kt_cat.union(kt_cts)
kt = kt.key_by('phenotype')
print 'n_both_sexes: {:,}'.format(kt.count())

kt.write('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb31063.phesant_phenosummaries.both_sexes.kt', overwrite=True)
kt.export('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb31063.phesant_phenosummaries.both_sexes.tsv')

## females
kt_cat = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/still-more-phesant/neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_females.*_phenosummary.tsv',
                         comment='Field', no_header=True, missing='NA')
kt_cat = kt_cat.rename({'f0': 'phenotype',
                        'f1': 'description',
                        'f2': 'n_non_missing',
                        'f3': 'n_missing',
                        'f4': 'n_cases',
                        'f5': 'n_controls',
                        'f6': 'notes',
                        'f7': 'PHESANT_transformation'})
kt_cat = kt_cat.annotate('type = "cat"')
kt_cat = kt_cat.select(['phenotype',
                        'description',
                        'notes',
                        'PHESANT_transformation',
                        'n_non_missing',
                        'n_missing',
                        'n_controls',
                        'n_cases',
                        'type'])

kt_cts = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/still-more-phesant/neale_lab_parsed_and_restricted_to_QCed_samples_cts_irnt_females_phenosummary.tsv',
                         comment='Field', no_header=True, missing='NA')
kt_cts = kt_cts.rename({'f0': 'phenotype',
                        'f1': 'description',
                        'f2': 'n_non_missing',
                        'f3': 'n_missing',
                        'f4': 'n_cases',
                        'f5': 'n_controls',
                        'f6': 'notes',
                        'f7': 'PHESANT_transformation'})
kt_cts = kt_cts.annotate('type = "cts"')
kt_cts = kt_cts.select(['phenotype',
                        'description',
                        'notes',
                        'PHESANT_transformation',
                        'n_non_missing',
                        'n_missing',
                        'n_controls',
                        'n_cases',
                        'type'])

kt = kt_cat.union(kt_cts)
kt = kt.key_by('phenotype')
print 'n_females: {:,}'.format(kt.count())

kt.write('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb31063.phesant_phenosummaries.female.kt', overwrite=True)
kt.export('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb31063.phesant_phenosummaries.female.tsv')

# males
kt_cat = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/still-more-phesant/neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_males.*_phenosummary.tsv',
                         comment='Field', no_header=True, missing='NA')
kt_cat = kt_cat.rename({'f0': 'phenotype',
                        'f1': 'description',
                        'f2': 'n_non_missing',
                        'f3': 'n_missing',
                        'f4': 'n_cases',
                        'f5': 'n_controls',
                        'f6': 'notes',
                        'f7': 'PHESANT_transformation'})
kt_cat = kt_cat.annotate('type = "cat"')
kt_cat = kt_cat.select(['phenotype',
                        'description',
                        'notes',
                        'PHESANT_transformation',
                        'n_non_missing',
                        'n_missing',
                        'n_controls',
                        'n_cases',
                        'type'])

kt_cts = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/still-more-phesant/neale_lab_parsed_and_restricted_to_QCed_samples_cts_irnt_males_phenosummary.tsv',
                         comment='Field', no_header=True, missing='NA')
kt_cts = kt_cts.rename({'f0': 'phenotype',
                        'f1': 'description',
                        'f2': 'n_non_missing',
                        'f3': 'n_missing',
                        'f4': 'n_cases',
                        'f5': 'n_controls',
                        'f6': 'notes',
                        'f7': 'PHESANT_transformation'})
kt_cts = kt_cts.annotate('type = "cts"')
kt_cts = kt_cts.select(['phenotype',
                        'description',
                        'notes',
                        'PHESANT_transformation',
                        'n_non_missing',
                        'n_missing',
                        'n_controls',
                        'n_cases',
                        'type'])

kt = kt_cat.union(kt_cts)
kt = kt.key_by('phenotype')
print 'n_males: {:,}'.format(kt.count())

kt.write('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb31063.phesant_phenosummaries.male.kt', overwrite=True)
kt.export('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb31063.phesant_phenosummaries.male.tsv')
