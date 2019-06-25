
import sys
from pprint import pprint
from hail import *
hc = HailContext()

sex = sys.argv[1]

if sex == 'both_sexes':
    plural_sex = 'both_sexes'
if sex == 'female':
    plural_sex = 'females'
if sex == 'male':
    plural_sex = 'males'
 
n_phesant = {
    'both_sexes': 40,
    'female': 35,
    'male': 31
}
n_additional = {
    'both_sexes': 14,
    'female': 12,
    'male': 11
}
n_still_more = {
    'both_sexes': 8,
    'female': 20,
    'male': 19
}

phesant_phenotypes = set()
for i in xrange(n_phesant[sex]):
    kt = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/phesant/ukb31063.{0}.phesant.pipeline.{1:}.kt'.format(sex, i))
    cols = [x for x in kt.columns if x != 's']
    phesant_phenotypes.update(set(cols))

additional_phenotypes = set()
for i in xrange(n_additional[sex]):
    kt = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/phesant-additional/ukb31063.{0}.phesant_additional.pipeline.{1:}.kt'.format(plural_sex, i))
    cols = [x for x in kt.columns if x != 's']
    additional_phenotypes.update(set(cols))

still_more_phenotypes = set()
for i in xrange(n_still_more[sex]):
    kt = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/still-more-phesant/ukb31063.{0}.still_more_phesant.pipeline.{1:}.kt'.format(plural_sex, i))
    cols = [x for x in kt.columns if x != 's']
    still_more_phenotypes.update(set(cols))  

kt = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/still-more-phesant-cts-irnt/ukb31063.{}.still_more_phesant_cts_irnt.pipeline.0.kt'.format(plural_sex)) 
still_more_irnt_phenotypes = set([x for x in kt.columns if x != 's'])

kt = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/still-more-phesant-cts-raw/ukb31063.{}.still_more_phesant_cts_raw.pipeline.0.kt'.format(plural_sex)) 
still_more_raw_phenotypes = set([x for x in kt.columns if x != 's'])

phesant_and_additional = phesant_phenotypes.intersection(additional_phenotypes)
phesant_and_still_more = phesant_phenotypes.intersection(still_more_phenotypes)
phesant_and_cts_irnt = phesant_phenotypes.intersection(still_more_irnt_phenotypes)
phesant_and_cts_raw = phesant_phenotypes.intersection(still_more_raw_phenotypes)

additional_and_still_more = additional_phenotypes.intersection(still_more_phenotypes)
additional_and_still_more_cts_irnt = additional_phenotypes.intersection(still_more_irnt_phenotypes)
additional_and_still_more_cts_raw = additional_phenotypes.intersection(still_more_raw_phenotypes)

with hadoop_write('gs://ukb31063-mega-gwas/phesant_and_additional.{}.txt'.format(sex)) as f:
    for x in phesant_and_additional:
        f.write('gs://ukb31063-mega-gwas/export2-results-tsvs/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')
with hadoop_write('gs://ukb31063-mega-gwas/phesant_and_still_more.{}.txt'.format(sex)) as f:
    for x in phesant_and_still_more:
        f.write('gs://ukb31063-mega-gwas/export2-results-tsvs/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')
with hadoop_write('gs://ukb31063-mega-gwas/phesant_and_still_more_cts_irnt.{}.txt'.format(sex)) as f:
    for x in phesant_and_cts_irnt:
        f.write('gs://ukb31063-mega-gwas/export2-results-tsvs/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')
with hadoop_write('gs://ukb31063-mega-gwas/phesant_and_still_more_cts_raw.{}.txt'.format(sex)) as f:
    for x in phesant_and_cts_raw:
        f.write('gs://ukb31063-mega-gwas/export2-results-tsvs/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')

with hadoop_write('gs://ukb31063-mega-gwas/additional_and_still_more.{}.txt'.format(sex)) as f:
    for x in additional_and_still_more:
        f.write('gs://ukb31063-mega-gwas/export2-results-tsvs/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')
with hadoop_write('gs://ukb31063-mega-gwas/additional_and_still_more_cts_irnt.{}.txt'.format(sex)) as f:
    for x in additional_and_still_more_cts_irnt:
        f.write('gs://ukb31063-mega-gwas/export2-results-tsvs/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')
with hadoop_write('gs://ukb31063-mega-gwas/additional_and_still_more_cts_raw.{}.txt'.format(sex)) as f:
    for x in additional_and_still_more_cts_raw:
        f.write('gs://ukb31063-mega-gwas/export2-results-tsvs/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')

with hadoop_write('gs://ukb31063-mega-gwas/remove_from_export.{}.txt'.format(sex)) as f:
    for x in phesant_and_cts_irnt:
        f.write('gs://ukbb-gwas-imputed-v3-results/export2/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')
    for x in phesant_and_cts_raw:
        f.write('gs://ukbb-gwas-imputed-v3-results/export2/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')
    for x in additional_and_still_more_cts_irnt:
        f.write('gs://ukbb-gwas-imputed-v3-results/export2/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')
    for x in additional_and_still_more_cts_raw:
        f.write('gs://ukbb-gwas-imputed-v3-results/export2/' + x + '.gwas.imputed_v3.{}.tsv.bgz'.format(sex) + '\n')

