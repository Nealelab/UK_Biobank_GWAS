
from hail import *
hc = HailContext()

n_samples = {'both_sexes': 361194, 'female': 194174, 'male': 167020}
n_phesant_pipelines = {'both_sexes': 40, 'female': 35, 'male': 31}
n_icd10_pipelines = {'both_sexes': 8, 'female': 6, 'male': 6}
n_finngen_pipelines = {'both_sexes': 7, 'female': 6, 'male': 5}

for sex in ['both_sexes', 'female', 'male']:
    
    kt_phesant = hc.import_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/phesant/ukb31063.{}.phenosummary.*.tsv'.format(sex), key='FieldID', types={'N.cases': TInt(), 'N.controls': TInt()})
    kt_phesant = kt_phesant.rename({'FieldID': 'phenotype', 'Field': 'description', 'N.cases': 'n_cases', 'N.controls': 'n_controls'})
    kt_phesant = kt_phesant.annotate('source = "phesant"')
    kt_phesant = kt_phesant.select(['phenotype', 'description', 'source', 'n_controls', 'n_cases'])

    phesant_count_keytables = []
    for i in xrange(n_phesant_pipelines[sex]):
        kt_phesant_pipeline = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/phesant/ukb31063.{0}.phesant.pipeline.{1:}.kt'.format(sex, i))
        cols = [x for x in kt_phesant_pipeline.columns if x != 's']
        counts = [{'phenotype': x[0], 'n_missing': n_samples[sex] - x[1], 'n_non_missing': x[1]} for x in zip(cols, kt_phesant_pipeline.query(['`{}`.filter(x => isDefined(x)).count()'.format(y) for y in cols]))]
        kt_count = KeyTable.from_py(hc=hc, rows_py=counts, schema=TStruct(['phenotype', 'n_missing', 'n_non_missing'], [TString(), TInt(), TInt()]))
        phesant_count_keytables.append(kt_count)
    kt_phesant_counts = KeyTable.union(*phesant_count_keytables)
    kt_phesant_counts = kt_phesant_counts.key_by('phenotype')
    kt_phesant = kt_phesant.join(kt_phesant_counts, how='inner')

    kt_icd10 = hc.import_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/icd10/coding19.tsv', key='coding')
    kt_icd10 = kt_icd10.rename({'coding': 'phenotype', 'meaning': 'description'})
    kt_icd10 = kt_icd10.select(['phenotype', 'description'])
    kt_icd10 = kt_icd10.annotate('description = "Diagnoses - main ICD10: " + description, source = "icd10"')
    
    icd10_count_keytables = []
    for i in xrange(n_icd10_pipelines[sex]):
        kt_icd10_pipeline = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/icd10/ukb31063.{0}.icd10.pipeline.{1:}.kt'.format(sex, i))
        cols = [x for x in kt_icd10_pipeline.columns if x != 's']
        counts = [{'phenotype': x[0], 'n_controls': n_samples[sex] - x[1], 'n_cases': x[1], 'n_missing': n_samples[sex] - x[2], 'n_non_missing': x[2]} 
                  for x in zip(cols, kt_icd10_pipeline.query(['`{}`.filter(x => x).count()'.format(y) for y in cols]), kt_icd10_pipeline.query(['`{}`.filter(x => isDefined(x)).count()'.format(y) for y in cols]))]
        kt_count = KeyTable.from_py(hc=hc, rows_py=counts, schema=TStruct(['phenotype', 'n_controls', 'n_cases', 'n_missing', 'n_non_missing'], [TString(), TInt(), TInt(), TInt(), TInt()]))
        icd10_count_keytables.append(kt_count)
    kt_icd10_counts = KeyTable.union(*icd10_count_keytables)
    kt_icd10_counts = kt_icd10_counts.key_by('phenotype')
    kt_icd10 = kt_icd10.join(kt_icd10_counts, how='inner')
    
    kt_finngen = hc.import_table('gs://ukb31063-mega-gwas/curated-phenotypes/2018-04-07_ukb-finngen-all-pheno-counts.tsv', key='NAME', missing='')
    kt_finngen = kt_finngen.rename({'NAME': 'phenotype', 'LONGNAME': 'description'})
    kt_finngen = kt_finngen.select(['phenotype', 'description'])
    kt_finngen = kt_finngen.annotate('source = "finngen"')

    finngen_count_keytables = []
    for i in xrange(n_finngen_pipelines[sex]):
        kt_finngen_pipeline = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/finngen/ukb31063.{0}.finngen.pipeline.{1:}.kt'.format(sex, i))
        cols = [x for x in kt_finngen_pipeline.columns if x != 's']
        counts = [{'phenotype': x[0], 'n_controls': n_samples[sex] - x[1], 'n_cases': x[1], 'n_missing': n_samples[sex] - x[2], 'n_non_missing': x[2]} 
                  for x in zip(cols, kt_finngen_pipeline.query(['`{}`.filter(x => x == 1).count()'.format(y) for y in cols]), kt_finngen_pipeline.query(['`{}`.filter(x => isDefined(x)).count()'.format(y) for y in cols]))]
        kt_count = KeyTable.from_py(hc=hc, rows_py=counts, schema=TStruct(['phenotype', 'n_controls', 'n_cases', 'n_missing', 'n_non_missing'], [TString(), TInt(), TInt(), TInt(), TInt()]))
        finngen_count_keytables.append(kt_count)
    kt_finngen_counts = KeyTable.union(*finngen_count_keytables)
    kt_finngen_counts = kt_finngen_counts.key_by('phenotype')
    kt_finngen = kt_finngen.join(kt_finngen_counts, how='inner')

    kt = KeyTable.union(kt_phesant, kt_icd10, kt_finngen)
    print 'nPhenotypes, {0}: {1:,}'.format(sex, kt.count())
    kt.export('gs://ukb31063-mega-gwas/hail-0.1/annotations/phenotypes.{}.tsv.gz'.format(sex))
    