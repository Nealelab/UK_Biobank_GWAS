
import sys
import hail as hl

sex = sys.argv[1]

"""
ht_hm3 = hl.read_table('gs://ukb31063-mega-gwas/hm3_variants.ht')
mt = hl.read_matrix_table('gs://ukb31063-mega-gwas/results-0.2-tables/gwas.imputed_v3.{}.annotated.mt'.format(sex))
mt = mt.filter_rows(hl.is_defined(ht_hm3[mt.locus, mt.alleles]))
mt.write('gs://ukb31063-mega-gwas/ldsc-export/matrix-tables/{}.mt'.format(sex), overwrite=True)
"""
"""
mt = hl.read_matrix_table('gs://ukb31063-mega-gwas/ldsc-export/matrix-tables/{}.mt'.format(sex))
mt = mt.annotate_rows(SNP=mt.rsid, A1=mt.alleles[1], A2=mt.alleles[0])
mt = mt.annotate_entries(N=mt.n_complete_samples, Z=mt.tstat)
mt = mt.key_rows_by('SNP')
mt = mt.key_cols_by()
mt = mt.select_rows('A1', 'A2')
mt = mt.select_cols('phenotype')
mt = mt.select_entries('N', 'Z')
mt = mt.repartition(36)
mt.write('gs://ukb31063-mega-gwas/ldsc-export/matrix-tables/{}.36_partitions.mt'.format(sex), overwrite=True)
"""
"""
ht_hm3 = hl.import_table('gs://ukb31063-mega-gwas/hm3_variants.tsv.bgz')
ht_hm3 = ht_hm3.annotate(**hl.parse_variant(ht_hm3.v))
ht_hm3 = ht_hm3.select('locus', 'alleles')
ht_hm3 = ht_hm3.key_by('locus', 'alleles')
print(ht_hm3.count())
ht_hm3.write('gs://ukb31063-mega-gwas/hm3_variants.ht', overwrite=True)
"""

from subprocess import check_output
already_exported = set([x.split('/')[-1].split('.')[0] for x in check_output('gsutil ls -l gs://ukb31063-mega-gwas/ldsc-export/sumstats-files/*.ldsc.imputed_v3.both_sexes.tsv.bgz', shell=True).decode('utf-8').split() if x.endswith('.tsv.bgz')])
print(len(already_exported))

mt = hl.read_matrix_table('gs://ukb31063-mega-gwas/ldsc-export/matrix-tables/{}.36_partitions.mt'.format(sex)).cache()
phenotypes = mt.phenotype.collect()

for p in phenotypes:
    if p in already_exported:
        continue
    (mt.filter_cols(mt.phenotype == p)
       .entries()
       .select('A1', 'A2', 'N', 'Z')
       .export('gs://ukb31063-mega-gwas/ldsc-export/sumstats-files/{0}.ldsc.imputed_v3.{1}.tsv.bgz'.format(p, sex)))
