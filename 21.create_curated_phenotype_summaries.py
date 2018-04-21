
from hail import *
hc = HailContext()

sexes = ['both_sexes', 'female', 'male']
for sex in sexes:
    kt = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/curated/ukb31063.{}.curated.kt'.format(sex)).cache()
    cols = [x for x in kt.columns if x != 's']
    queries = kt.query(['{}.collectAsSet()'.format(x) for x in cols])
    export = [['phenotype', 'n_controls', 'n_cases']]
    for i, c in enumerate(cols):
        if queries[i] - set([None, 0, 1]):
            continue
        counts = kt.query('`{}`.counter()'.format(c))
        export.append([c, counts[0], counts[1]])
    with hadoop_write('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/curated/ukb31063.{}.curated.tsv'.format(sex)) as f:
        for x in export:
            f.write('\t'.join([x[0], str(x[1]), str(x[2])]) + '\n')
