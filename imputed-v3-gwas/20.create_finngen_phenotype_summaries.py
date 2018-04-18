
from hail import *
hc = HailContext()

for i in xrange(7):
    kt = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/finngen/ukb31063.both_sexes.finngen.pipeline.{:}.kt'.format(i))
    cols = [x for x in kt.columns if x != 's']
    qry = kt.query(['`{}`.counter()'.format(x) for x in cols])
    counts = zip(cols, qry)
    with hadoop_write('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/finngen/ukb31063.both_sexes.finngen.phenosummary.pipeline.{:}.tsv'.format(i)) as f:
        f.write('\t'.join(['code', 'n_controls', 'n_cases']) + '\n')
        for x in counts:
            f.write('\t'.join([x[0], str(x[1][False]), str(x[1][True])]) + '\n')

for i in xrange(6):
    kt = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/finngen/ukb31063.female.finngen.pipeline.{:}.kt'.format(i))
    cols = [x for x in kt.columns if x != 's']
    qry = kt.query(['`{}`.counter()'.format(x) for x in cols])
    counts = zip(cols, qry)
    with hadoop_write('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/finngen/ukb31063.female.finngen.phenosummary.pipeline.{:}.tsv'.format(i)) as f:
        f.write('\t'.join(['code', 'n_controls', 'n_cases']) + '\n')
        for x in counts:
            f.write('\t'.join([x[0], str(x[1][False]), str(x[1][True])]) + '\n')

for i in xrange(5):
    kt = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/finngen/ukb31063.male.finngen.pipeline.{:}.kt'.format(i))
    cols = [x for x in kt.columns if x != 's']
    qry = kt.query(['`{}`.counter()'.format(x) for x in cols])
    counts = zip(cols, qry)
    with hadoop_write('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/finngen/ukb31063.male.finngen.phenosummary.pipeline.{:}.tsv'.format(i)) as f:
        f.write('\t'.join(['code', 'n_controls', 'n_cases']) + '\n')
        for x in counts:
            f.write('\t'.join([x[0], str(x[1][False]), str(x[1][True])]) + '\n')
