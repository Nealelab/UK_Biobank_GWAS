
from hail import *
hc = HailContext()

for i in xrange(8):
    kt = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/icd10/ukb31063.both_sexes.icd10.pipeline.{:}.kt'.format(i))
    cols = [x for x in kt.columns if x != 's']
    qry = kt.query(['{}.counter()'.format(x) for x in cols])
    counts = zip(cols, qry)
    with hadoop_write('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/icd10/ukb31063.both_sexes.icd10.phenosummary.pipeline.{:}.tsv'.format(i)) as f:
        f.write('\t'.join(['code', 'n_controls', 'n_cases']) + '\n')
        for x in counts:
            f.write('\t'.join([x[0], str(x[1][False]), str(x[1][True])]) + '\n')

for i in xrange(6):
    kt = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/icd10/ukb31063.female.icd10.pipeline.{:}.kt'.format(i))
    cols = [x for x in kt.columns if x != 's']
    qry = kt.query(['{}.counter()'.format(x) for x in cols])
    counts = zip(cols, qry)
    with hadoop_write('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/icd10/ukb31063.female.icd10.phenosummary.pipeline.{:}.tsv'.format(i)) as f:
        f.write('\t'.join(['code', 'n_controls', 'n_cases']) + '\n')
        for x in counts:
            f.write('\t'.join([x[0], str(x[1][False]), str(x[1][True])]) + '\n')

for i in xrange(6):
    kt = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/icd10/ukb31063.male.icd10.pipeline.{:}.kt'.format(i))
    cols = [x for x in kt.columns if x != 's']
    qry = kt.query(['{}.counter()'.format(x) for x in cols])
    counts = zip(cols, qry)
    with hadoop_write('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/icd10/ukb31063.male.icd10.phenosummary.pipeline.{:}.tsv'.format(i)) as f:
        f.write('\t'.join(['code', 'n_controls', 'n_cases']) + '\n')
        for x in counts:
            f.write('\t'.join([x[0], str(x[1][False]), str(x[1][True])]) + '\n')
