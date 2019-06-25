
import gzip
from hail import *
hc = HailContext()

files = ['gs://ukbb-gwas-imputed-v3-results/export2/variants.tsv.bgz',
         'gs://ukbb-gwas-imputed-v3-results/export2/phenotypes.both_sexes.tsv.bgz',
         'gs://ukbb-gwas-imputed-v3-results/export2/phenotypes.female.tsv.bgz',
         'gs://ukbb-gwas-imputed-v3-results/export2/phenotypes.male.tsv.bgz']

with hadoop_read('gs://ukb31063-mega-gwas/annotations/phenotypes.both_sexes.tsv.bgz') as f:
    f.readline()
    for x in f:
        files.append('gs://ukbb-gwas-imputed-v3-results/export2/' + x.split('\t')[0] + '.gwas.imputed_v3.both_sexes.tsv.bgz')

with hadoop_read('gs://ukb31063-mega-gwas/annotations/phenotypes.female.tsv.bgz') as f:
    f.readline()
    for x in f:
        files.append('gs://ukbb-gwas-imputed-v3-results/export2/' + x.split('\t')[0] + '.gwas.imputed_v3.female.tsv.bgz')

with hadoop_read('gs://ukb31063-mega-gwas/annotations/phenotypes.male.tsv.bgz') as f:
    f.readline()
    for x in f:
        files.append('gs://ukbb-gwas-imputed-v3-results/export2/' + x.split('\t')[0] + '.gwas.imputed_v3.male.tsv.bgz')

with hadoop_write('gs://ukb31063-mega-gwas/master_export_file_list.txt') as f:
    for x in files:
        f.write(x + '\n')
