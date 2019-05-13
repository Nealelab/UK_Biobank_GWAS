
import sys
import hail as hl

sex = sys.argv[1]

ht = hl.read_table('gs://ukb31063-mega-gwas/results-0.2-tables/gwas.imputed_v3.{}.ht'.format(sex))

ht_phenotypes = hl.import_table('gs://ukbb-gwas-imputed-v3-results/export2/phenotypes.{}.tsv.bgz'.format(sex),
                                key='phenotype', 
                                types={'n_non_missing': hl.tint,
                                       'n_missing': hl.tint,
                                       'n_controls': hl.tint,
                                       'n_cases': hl.tint})

ht_variants = hl.import_table('gs://ukbb-gwas-imputed-v3-results/export2/variants.tsv.bgz',
                              types={'pos': hl.tint,
                                     'info': hl.tfloat,
                                     'call_rate': hl.tfloat,
                                     'AC': hl.tint,
                                     'AF': hl.tfloat,
                                     'minor_AF': hl.tfloat,
                                     'p_hwe': hl.tfloat,
                                     'n_called': hl.tint,
                                     'n_not_called': hl.tint,
                                     'n_hom_ref': hl.tint,
                                     'n_het': hl.tint,
                                     'n_hom_var': hl.tint,
                                     'n_non_ref': hl.tint,
                                     'r_heterozygosity': hl.tfloat,
                                     'r_het_hom_var': hl.tfloat,
                                     'r_expected_het_frequency': hl.tfloat})
ht_variants = ht_variants.annotate(**hl.parse_variant(ht_variants.variant))
ht_variants = ht_variants.select('locus',
                                 'alleles',
                                 'variant',
                                 'rsid',
                                 'varid',
                                 'consequence',
                                 'consequence_category',
                                 'info',
                                 'call_rate',
                                 'AC',
                                 'AF',
                                 'minor_allele',
                                 'minor_AF',
                                 'p_hwe',
                                 'n_called',
                                 'n_not_called',
                                 'n_hom_ref',
                                 'n_het',
                                 'n_hom_var',
                                 'n_non_ref',
                                 'r_heterozygosity',
                                 'r_het_hom_var',
                                 'r_expected_het_frequency')
ht_variants = ht_variants.rename({'AC': 'alt_AC'})
ht_variants = ht_variants.key_by('locus', 'alleles')

mt = ht.to_matrix_table(row_key=['locus', 'alleles'],
                        col_key=['phenotype'],
                        row_fields=['variant', 'minor_allele', 'minor_AF'],
                        partition_key=['locus'])
mt = mt.annotate_cols(**ht_phenotypes[mt.phenotype])
mt = mt.annotate_rows(**ht_variants[mt.locus, mt.alleles])
mt.write('gs://ukb31063-mega-gwas/results-0.2-tables/gwas.imputed_v3.{}.annotated.mt'.format(sex), overwrite=True)
