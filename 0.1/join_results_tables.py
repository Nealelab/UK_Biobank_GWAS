
import hail as hl

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

union_set = both_sexes_set.union(female_set, male_set)
n_phenotypes = len(union_set)

counter = 1
for phenotype in union_set:
    print('Loading {0}, {1:,}/{2:,}'.format(phenotype, counter, n_phenotypes))
    try:
        ht_both_sexes = hl.import_table('gs://ukbb-gwas-imputed-v3-results/export2/{}.gwas.imputed_v3.both_sexes.tsv.bgz'.format(phenotype),
                                        min_partitions=12,
                                        types={'minor_AF': hl.tfloat,
                                               'expected_case_minor_AC': hl.tfloat,
                                               'low_confidence_variant': hl.tbool,
                                               'n_complete_samples': hl.tint,
                                               'AC': hl.tfloat,
                                               'ytx': hl.tfloat,
                                               'beta': hl.tfloat,
                                               'se': hl.tfloat,
                                               'tstat': hl.tfloat,
                                               'pval': hl.tfloat})

    except:
        ht_both_sexes = False

    else:
        cols = [x for x in ht_both_sexes.row]
        ht_both_sexes = ht_both_sexes.annotate(**hl.parse_variant(ht_both_sexes.variant))
        if 'expected_case_minor_AC' in cols:
            ht_both_sexes = ht_both_sexes.annotate(both_sexes=hl.struct(expected_case_minor_AC=hl.float(ht_both_sexes.expected_case_minor_AC),
                                                                        low_confidence_variant=ht_both_sexes.low_confidence_variant,
                                                                        n_complete_samples=ht_both_sexes.n_complete_samples,
                                                                        AC=ht_both_sexes.AC,
                                                                        ytx=ht_both_sexes.ytx,
                                                                        beta=ht_both_sexes.beta,
                                                                        se=ht_both_sexes.se,
                                                                        tstat=ht_both_sexes.tstat,
                                                                        pval=ht_both_sexes.pval))
        elif 'expected_min_category_minor_AC' in cols:
            ht_both_sexes = ht_both_sexes.annotate(both_sexes=hl.struct(expected_min_category_minor_AC=hl.float(ht_both_sexes.expected_min_category_minor_AC),
                                                                        low_confidence_variant=ht_both_sexes.low_confidence_variant,
                                                                        n_complete_samples=ht_both_sexes.n_complete_samples,
                                                                        AC=ht_both_sexes.AC,
                                                                        ytx=ht_both_sexes.ytx,
                                                                        beta=ht_both_sexes.beta,
                                                                        se=ht_both_sexes.se,
                                                                        tstat=ht_both_sexes.tstat,
                                                                        pval=ht_both_sexes.pval))
        else:
            ht_both_sexes = ht_both_sexes.annotate(both_sexes=hl.struct(low_confidence_variant=ht_both_sexes.low_confidence_variant,
                                                                        n_complete_samples=ht_both_sexes.n_complete_samples,
                                                                        AC=ht_both_sexes.AC,
                                                                        ytx=ht_both_sexes.ytx,
                                                                        beta=ht_both_sexes.beta,
                                                                        se=ht_both_sexes.se,
                                                                        tstat=ht_both_sexes.tstat,
                                                                        pval=ht_both_sexes.pval))
        ht_both_sexes = ht_both_sexes.select('variant', 
                                             'locus', 
                                             'alleles', 
                                             'minor_allele', 
                                             'minor_AF', 
                                             'both_sexes')
        ht_both_sexes = ht_both_sexes.key_by('variant')

    try:
        ht_female = hl.import_table('gs://ukbb-gwas-imputed-v3-results/export2/{}.gwas.imputed_v3.female.tsv.bgz'.format(phenotype),
                                    min_partitions=12,
                                    types={'minor_AF': hl.tfloat,
                                           'expected_case_minor_AC': hl.tfloat,
                                           'low_confidence_variant': hl.tbool,
                                           'n_complete_samples': hl.tint,
                                           'AC': hl.tfloat,
                                           'ytx': hl.tfloat,
                                           'beta': hl.tfloat,
                                           'se': hl.tfloat,
                                           'tstat': hl.tfloat,
                                           'pval': hl.tfloat})

    except:
        ht_female = False

    else:
        cols = [x for x in ht_female.row]
        if 'expected_case_minor_AC' in cols:
            ht_female = ht_female.annotate(female=hl.struct(expected_case_minor_AC=hl.float(ht_female.expected_case_minor_AC),
                                                            low_confidence_variant=ht_female.low_confidence_variant,
                                                            n_complete_samples=ht_female.n_complete_samples,
                                                            AC=ht_female.AC,
                                                            ytx=ht_female.ytx,
                                                            beta=ht_female.beta,
                                                            se=ht_female.se,
                                                            tstat=ht_female.tstat,
                                                            pval=ht_female.pval))
        elif 'expected_min_category_minor_AC' in cols:
            ht_female = ht_female.annotate(female=hl.struct(expected_min_category_minor_AC=hl.float(ht_female.expected_min_category_minor_AC),
                                                            low_confidence_variant=ht_female.low_confidence_variant,
                                                            n_complete_samples=ht_female.n_complete_samples,
                                                            AC=ht_female.AC,
                                                            ytx=ht_female.ytx,
                                                            beta=ht_female.beta,
                                                            se=ht_female.se,
                                                            tstat=ht_female.tstat,
                                                            pval=ht_female.pval))
        else:
            ht_female = ht_female.annotate(female=hl.struct(low_confidence_variant=ht_female.low_confidence_variant,
                                                            n_complete_samples=ht_female.n_complete_samples,
                                                            AC=ht_female.AC,
                                                            ytx=ht_female.ytx,
                                                            beta=ht_female.beta,
                                                            se=ht_female.se,
                                                            tstat=ht_female.tstat,
                                                            pval=ht_female.pval))
        if not ht_both_sexes:
            ht_female = ht_female.annotate(**hl.parse_variant(ht_female.variant))
            ht_female = ht_female.select('variant', 
                                         'locus', 
                                         'alleles', 
                                         'minor_allele', 
                                         'minor_AF', 
                                         'female')
        else:
            ht_female = ht_female.select('variant', 'female')
        ht_female = ht_female.key_by('variant')

    try:
        ht_male = hl.import_table('gs://ukbb-gwas-imputed-v3-results/export2/{}.gwas.imputed_v3.male.tsv.bgz'.format(phenotype),
                                  min_partitions=12,
                                  types={'minor_AF': hl.tfloat,
                                         'expected_case_minor_AC': hl.tfloat,
                                         'low_confidence_variant': hl.tbool,
                                         'n_complete_samples': hl.tint,
                                         'AC': hl.tfloat,
                                         'ytx': hl.tfloat,
                                         'beta': hl.tfloat,
                                         'se': hl.tfloat,
                                         'tstat': hl.tfloat,
                                         'pval': hl.tfloat})

    except:
        ht_male = False

    else:
        cols = [x for x in ht_male.row]
        if 'expected_case_minor_AC' in cols:
            ht_male = ht_male.annotate(male=hl.struct(expected_case_minor_AC=hl.float(ht_male.expected_case_minor_AC),
                                                      low_confidence_variant=ht_male.low_confidence_variant,
                                                      n_complete_samples=ht_male.n_complete_samples,
                                                      AC=ht_male.AC,
                                                      ytx=ht_male.ytx,
                                                      beta=ht_male.beta,
                                                      se=ht_male.se,
                                                      tstat=ht_male.tstat,
                                                      pval=ht_male.pval))
        elif 'expected_min_category_minor_AC' in cols:
            ht_male = ht_male.annotate(male=hl.struct(expected_min_category_minor_AC=hl.float(ht_male.expected_min_category_minor_AC),
                                                      low_confidence_variant=ht_male.low_confidence_variant,
                                                      n_complete_samples=ht_male.n_complete_samples,
                                                      AC=ht_male.AC,
                                                      ytx=ht_male.ytx,
                                                      beta=ht_male.beta,
                                                      se=ht_male.se,
                                                      tstat=ht_male.tstat,
                                                      pval=ht_male.pval))
        else:
            ht_male = ht_male.annotate(male=hl.struct(low_confidence_variant=ht_male.low_confidence_variant,
                                                      n_complete_samples=ht_male.n_complete_samples,
                                                      AC=ht_male.AC,
                                                      ytx=ht_male.ytx,
                                                      beta=ht_male.beta,
                                                      se=ht_male.se,
                                                      tstat=ht_male.tstat,
                                                      pval=ht_male.pval))

        if not ht_both_sexes and not ht_female:
            ht_male = ht_male.annotate(**hl.parse_variant(ht_male.variant))
            ht_male = ht_male.select('variant', 
                                     'locus', 
                                     'alleles', 
                                     'minor_allele', 
                                     'minor_AF', 
                                     'male')
        else:
            ht_male = ht_male.select('variant', 'male')
        ht_male = ht_male.key_by('variant')

    if all([ht_both_sexes, ht_female, ht_male]):
        ht = ht_both_sexes.join(ht_female, how='inner').join(ht_male, how='inner')

    elif all([ht_both_sexes, ht_female, not ht_male]):
        ht = ht_both_sexes.join(ht_female, how='inner')

    elif all([ht_both_sexes, not ht_female, ht_male]):
        ht = ht_both_sexes.join(ht_male, how='inner')

    elif all([not ht_both_sexes, ht_female, ht_male]):
        ht = ht_female.join(ht_male, how='inner')

    elif all([ht_both_sexes, not ht_female, not ht_male]):
        ht = ht_both_sexes

    elif all([not ht_both_sexes, not ht_female, ht_male]):
        ht = ht_male

    elif all([not ht_both_sexes, ht_female, not ht_male]):
        ht = ht_female

    elif all([not ht_both_sexes, not ht_female, not ht_male]):
        pass

    ht = ht.annotate(phenotype=phenotype)
    ht.write('gs://ukb31063-mega-gwas/results-hail-tables/per-sex/{}.{}.ht'.format(phenotype, sex), overwrite=True)

    counter += 1

