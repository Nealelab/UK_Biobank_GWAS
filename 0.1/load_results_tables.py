
import sys
import hail as hl

sex = sys.argv[1]

phenotypes_both_sexes = hl.import_table('gs://ukbb-gwas-imputed-v3-results/export2/phenotypes.both_sexes.tsv.bgz', 
                                types={'n_non_missing': hl.tint,
                                       'n_missing': hl.tint,
                                       'n_controls': hl.tint,
                                       'n_cases': hl.tint})

phenotypes_female = hl.import_table('gs://ukbb-gwas-imputed-v3-results/export2/phenotypes.female.tsv.bgz', 
                                types={'n_non_missing': hl.tint,
                                       'n_missing': hl.tint,
                                       'n_controls': hl.tint,
                                       'n_cases': hl.tint})

phenotypes_male = hl.import_table('gs://ukbb-gwas-imputed-v3-results/export2/phenotypes.male.tsv.bgz', 
                                types={'n_non_missing': hl.tint,
                                       'n_missing': hl.tint,
                                       'n_controls': hl.tint,
                                       'n_cases': hl.tint})

phenotypes = {
    'both_sexes': set(phenotypes_both_sexes.phenotype.collect()),
    'female': set(phenotypes_female.phenotype.collect()),
    'male': set(phenotypes_male.phenotype.collect())  
}

n_phenotypes = len(phenotypes[sex])
counter = 1

hts = []
for phenotype in phenotypes[sex]:
    print('Loading {0}, {1:,}/{2:,}'.format(phenotype, counter, n_phenotypes))
    ht = hl.import_table('gs://ukbb-gwas-imputed-v3-results/export2/{0}.gwas.imputed_v3.{1}.tsv.bgz'.format(phenotype, sex),
                         types={'minor_AF': hl.tfloat,
                                'low_confidence_variant': hl.tbool,
                                'n_complete_samples': hl.tint,
                                'AC': hl.tfloat,
                                'ytx': hl.tfloat,
                                'beta': hl.tfloat,
                                'se': hl.tfloat,
                                'tstat': hl.tfloat,
                                'pval': hl.tfloat})

    cols = [x for x in ht.row]
    ht = ht.annotate(**hl.parse_variant(ht.variant))
    ht = ht.annotate(phenotype=phenotype)
    if 'expected_case_minor_AC' in cols:
        ht = ht.annotate(expected_case_minor_AC=hl.float(ht.expected_case_minor_AC))
        ht = ht.annotate(expected_min_category_minor_AC=hl.or_missing(False, 0.0))
    elif 'expected_min_category_minor_AC' in cols:
        ht = ht.annotate(expected_case_minor_AC=hl.or_missing(False, 0.0))
        ht = ht.annotate(expected_min_category_minor_AC=hl.float(ht.expected_min_category_minor_AC))
    else:
        ht = ht.annotate(expected_case_minor_AC=hl.or_missing(False, 0.0))
        ht = ht.annotate(expected_min_category_minor_AC=hl.or_missing(False, 0.0))

    ht = ht.annotate(phenotype=phenotype)
    ht = ht.select('variant',
                   'minor_allele',
                   'minor_AF',
                   'expected_case_minor_AC',
                   'expected_min_category_minor_AC',
                   'low_confidence_variant',
                   'n_complete_samples',
                   'AC',
                   'ytx',
                   'beta',
                   'se',
                   'tstat',
                   'pval',
                   'locus',
                   'alleles',
                   'phenotype')
    ht = ht.key_by('locus', 'alleles')
    hts.append(ht)
    counter += 1

ht_union = hl.Table.union(*hts)
ht_union.write('gs://ukb31063-mega-gwas/results-0.2-tables/gwas.imputed_v3.{}.ht'.format(sex), overwrite=True)
