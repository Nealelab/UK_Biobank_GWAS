
import sys
import hail as hl

sex = sys.argv[1]

"""
ht_hm3 = hl.read_table('gs://ukb31063-mega-gwas/hm3_variants.ht')
mt = hl.read_matrix_table('gs://ukb31063-mega-gwas/results-0.2-tables/gwas.imputed_v3.{}.annotated.mt'.format(sex))
mt = mt.filter_rows(hl.is_defined(ht_hm3[mt.locus, mt.alleles]))
mt.write('gs://ukb31063-mega-gwas/ldsc-export/matrix-tables/{}.mt'.format(sex), overwrite=True)
"""


mt = hl.read_matrix_table('gs://ukb31063-mega-gwas/ldsc-export/matrix-tables/{}.mt'.format(sex))
mt = mt.annotate_rows(SNP=mt.rsid, A1=mt.alleles[1], A2=mt.alleles[0])
mt = mt.annotate_entries(N=mt.n_complete_samples, Z=mt.tstat)
mt = mt.select_rows('SNP', 'A1', 'A2')
mt = mt.select_entries('N', 'Z')
mt = mt.select_cols()
mt = mt.repartition(500)
mt.write('gs://ukb31063-mega-gwas/ldsc-export/matrix-tables/{}.2500_partitions.mt'.format(sex), overwrite=True)


"""
ht_hm3 = hl.import_table('gs://ukb31063-mega-gwas/hm3_variants.tsv.bgz')
ht_hm3 = ht_hm3.annotate(**hl.parse_variant(ht_hm3.v))
ht_hm3 = ht_hm3.select('locus', 'alleles')
ht_hm3 = ht_hm3.key_by('locus', 'alleles')
print(ht_hm3.count())
ht_hm3.write('gs://ukb31063-mega-gwas/hm3_variants.ht', overwrite=True)
"""

"""
mt = hl.read_matrix_table('gs://ukb31063-mega-gwas/ldsc-export/matrix-tables/{}.2500_partitions.mt'.format(sex))
phenotypes = mt.phenotype.collect()

mt = mt.cache()
mt.describe()

for p in phenotypes:
    mt = mt.filter_cols(mt.phenotype == p)
    ht = mt.entries().key_by(None).select('SNP', 'A1', 'A2', 'N', 'Z')
    ht.export('gs://ukb31063-mega-gwas/ldsc-export/sumstats-files/{0}.ldsc.imputed_v3.{1}.tsv.bgz'.format(p, sex))
"""