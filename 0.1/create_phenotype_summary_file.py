
from pprint import pprint
from hail import *
hc = HailContext()

n_samples = {'both_sexes': 361194, 'female': 194174, 'male': 167020}

icd10_descriptions = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/icd10/coding19.tsv', key='coding')
icd10_descriptions = icd10_descriptions.rename({'coding': 'phenotype', 'meaning': 'description'})
icd10_descriptions = icd10_descriptions.select(['phenotype', 'description'])
icd10_descriptions = icd10_descriptions.annotate(['description = "Diagnoses - main ICD10: " + description'])
icd10_descriptions = icd10_descriptions.filter('phenotype.length() == 3')

finngen_descriptions = hc.import_table('gs://ukb31063-mega-gwas/phenotype-files/curated-phenotypes/2018-04-07_ukb-finngen-all-pheno-counts.tsv', key='NAME', missing='')
finngen_descriptions = finngen_descriptions.rename({'NAME': 'phenotype', 'LONGNAME': 'description'})
finngen_descriptions = finngen_descriptions.filter('!isDefined(HD_ICD_10)')
finngen_descriptions = finngen_descriptions.select(['phenotype', 'description'])

to_remove = {}
with hadoop_read('gs://ukb31063-mega-gwas/phenotype-files/still-more-phesant/to_remove_both_sexes.tsv') as f:
    to_remove['both_sexes'] = [x.strip() for x in f.readlines()]
with hadoop_read('gs://ukb31063-mega-gwas/phenotype-files/still-more-phesant/to_remove_females.tsv') as f:
    to_remove['female'] = [x.strip() for x in f.readlines()]
with hadoop_read('gs://ukb31063-mega-gwas/phenotype-files/still-more-phesant/to_remove_males.tsv') as f:
    to_remove['male'] = [x.strip() for x in f.readlines()]

phesant_to_copy = []
icd10_finngen_to_remove = []

for sex in ['both_sexes', 'female', 'male']:

    kt_icd10 = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/icd10/ukb31063.{}.icd10.phenosummary.pipeline.*.tsv'.format(sex), 
                               key='code', missing='', types={'n_controls': TInt(), 'n_cases': TInt()})
    kt_icd10 = kt_icd10.rename({'code': 'phenotype'})
    kt_icd10 = kt_icd10.join(icd10_descriptions, how='inner')
    kt_icd10 = kt_icd10.annotate(['n_missing = 0',
                                  'n_non_missing = {:}'.format(n_samples[sex]), 
                                  'source = "icd10"', 
                                  'notes = "NA"', 
                                  'PHESANT_transformation = "NA"'])
    kt_icd10 = kt_icd10.select(['phenotype',
                                'description',
                                'source',
                                'n_non_missing',
                                'n_missing',
                                'n_controls',
                                'n_cases',
                                'PHESANT_transformation',
                                'notes'])
    print 'n_phenotypes, icd10, {0}: {1:}'.format(sex, kt_icd10.count())

    kt_finngen = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/finngen/ukb31063.{}.finngen.phenosummary.pipeline.*.tsv'.format(sex), 
                                 key='code', missing='', types={'n_controls': TInt(), 'n_cases': TInt()})
    kt_finngen = kt_finngen.rename({'code': 'phenotype'})
    kt_finngen = kt_finngen.join(finngen_descriptions, how='inner')
    kt_finngen = kt_finngen.annotate(['n_missing = 0', 
                                      'n_non_missing = {:}'.format(n_samples[sex]), 
                                      'source = "finngen"', 
                                      'notes = "NA"', 
                                      'PHESANT_transformation = "NA"'])
    kt_finngen = kt_finngen.select(['phenotype',
                                    'description',
                                    'source',
                                    'n_non_missing',
                                    'n_missing',
                                    'n_controls',
                                    'n_cases',
                                    'PHESANT_transformation',
                                    'notes'])
    print 'n_phenotypes, finngen, {0}: {1:}'.format(sex, kt_finngen.count())

    kt_phesant = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb31063.phesant_phenosummaries.{}.tsv'.format(sex), 
                                 key='phenotype', missing='NA', types={'n_non_missing': TInt(), 'n_missing': TInt(), 'n_controls': TInt(), 'n_cases': TInt()})
    kt_phesant = kt_phesant.annotate('source = "phesant"')
    kt_phesant = kt_phesant.select(['phenotype',
                                    'description',
                                    'source',
                                    'n_non_missing',
                                    'n_missing',
                                    'n_controls',
                                    'n_cases',
                                    'PHESANT_transformation',
                                    'notes'])

    test = kt_phesant.filter('phenotype == "100001"')
    print test.take(1)

    if sex == 'both_sexes':
        plural_sex = 'both_sexes'
    if sex == 'female':
        plural_sex = 'females'
    if sex == 'male':
        plural_sex = 'males'

    kt_cts_irnt = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/still-more-phesant-cts-irnt/ukb31063.{}.still_more_phesant_cts_irnt.pipeline.0.kt'.format(plural_sex))
    kt_cts_raw = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/still-more-phesant-cts-raw/ukb31063.{}.still_more_phesant_cts_raw.pipeline.0.kt'.format(plural_sex))
    
    irnt_phenotypes = [x for x in kt_cts_irnt.columns if x != 's']
    raw_phenotypes = [x for x in kt_cts_raw.columns if x != 's']

    kt_phesant = kt_phesant.annotate("""types = let ip = ["{0}"].toSet() and rp = ["{1}"].toSet() in
                                        if (ip.contains(phenotype) && rp.contains(phenotype)) ["irnt", "raw"]
                                        else if (ip.contains(phenotype)) ["irnt"]
                                        else if (rp.contains(phenotype)) ["raw"]
                                        else ["cat"]""".format('","'.join(irnt_phenotypes), '","'.join(raw_phenotypes)))
    kt_phesant = kt_phesant.explode('types')
    kt_phesant = kt_phesant.annotate("""phenotype = if (types == "irnt") phenotype + "_irnt"
                                                    else if (types == "raw") phenotype + "_raw"
                                                    else phenotype""")
    kt_phesant = kt_phesant.drop('types') 
    print 'n_phenotypes, phesant, {0}: {1:}'.format(sex, kt_phesant.count())                                  

    kt_union = KeyTable.union(kt_icd10, kt_finngen, kt_phesant)
    
    kt_keep = kt_union.filter("""let c = ["{}"].toSet() in (source == "phesant" && c.contains(phenotype) || (phenotype == "51")) || 
                                                           ((source == "icd10" || source == "finngen") && n_cases < 100)""", keep=False)
    kt_remove_phesant = kt_union.filter('let c = ["{}"].toSet() in (source == "phesant") && (c.contains(phenotype) || (phenotype == "51"))'.format('","'.join(to_remove[sex])))
    kt_remove_icd10_finngen = kt_union.filter('(source == "icd10" || source == "finngen") && n_cases < 100')

    #kt_keep.export('gs://ukb31063-mega-gwas/annotations/phenotypes.{}.tsv.bgz'.format(sex))

    phesant_to_copy.append(
        (kt_keep.filter('source == "phesant"')
                .annotate('path = "gs://ukb31063-mega-gwas/export2-results-tsvs/" + phenotype + ".gwas.imputed_v3.{}.tsv.bgz"'.format(sex))
                .select('path'))
    )

    icd10_finngen_to_remove.append(
        (kt_remove_icd10_finngen.annotate('path = "gs://ukbb-gwas-imputed-v3-results/export2/" + phenotype + ".gwas.imputed_v3.{}.tsv.bgz"'.format(sex))
                                .select('path'))
    )

    print 'n_phenotypes, union, {0}: {1:}'.format(sex, kt_union.count())
    print 'n_phenotypes, to remove phesant, {0}: {1:}'.format(sex, kt_remove_phesant.count())
    print 'n_phenotypes, to remove icd10/finngen, {0}: {1:}'.format(sex, kt_remove_icd10_finngen.count())

#KeyTable.union(*phesant_to_copy).export('gs://ukb31063-mega-gwas/phesant_to_copy.txt', header=False)
#KeyTable.union(*icd10_finngen_to_remove).export('gs://ukb31063-mega-gwas/icd10_finngen_to_remove.txt', header=False)
