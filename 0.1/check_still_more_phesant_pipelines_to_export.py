
import sys
from pprint import pprint
from hail import *
hc = HailContext()

n_still_more = {
    'both_sexes': 8,
    'female': 20,
    'male': 19
}

for sex in ['both_sexes', 'female', 'male']:

    if sex == 'both_sexes':
        plural_sex = 'both_sexes'
    if sex == 'female':
        plural_sex = 'females'
    if sex == 'male':
        plural_sex = 'males'

    with hadoop_read('gs://ukb31063-mega-gwas/phesant_and_still_more.{}.txt'.format(sex)) as f:
        still_more_overwrite_phenotypes = set([x.strip().split('/')[-1].split('.')[0] for x in f.readlines()])
    with hadoop_read('gs://ukb31063-mega-gwas/additional_and_still_more.{}.txt'.format(sex)) as f:
        still_more_overwrite_phenotypes.update(set([x.strip().split('/')[-1].split('.')[0] for x in f.readlines()]))

    pipelines = {}

    for i in xrange(n_still_more[sex]):
        kt = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/still-more-phesant/ukb31063.{0}.still_more_phesant.pipeline.{1:}.kt'.format(plural_sex, i))
        pipelines[i] = [x for x in kt.columns if x in still_more_overwrite_phenotypes]

    with hadoop_write('gs://ukb31063-mega-gwas/still_more_phesant.to_export.{}.tsv'.format(sex)) as f:
        f.write('\t'.join(['pipeline', 'phenotype']) + '\n')
        for i, x in pipelines.iteritems():
            for p in x:
                f.write('\t'.join([str(i), p]) + '\n')
