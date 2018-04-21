
from hail import *
hc = HailContext()

kt = hc.import_table('gs://phenotype_31063/ukb31063_sample_qc.tsv', key='iid', min_partitions=12)
kt = kt.drop(['fid', 'mid', 'sex', 'phenotype', 'ID_1', 'ID_2'])

kt = kt.rename({'iid': 's',
                'genotyping.array': 'genotyping_array',
                'Plate.Name': 'Plate_Name',
                'Cluster.CR': 'Cluster_CR',
                'Internal.Pico..ng.uL.': 'Internal_Pico_ng_uL',
                'Submitted.Gender': 'Submitted_Gender',
                'Inferred.Gender': 'Inferred_Gender',
                'X.intensity': 'X_intensity',
                'Y.intensity': 'Y_intensity',
                'Submitted.Plate.Name': 'Submitted_Plate_Name',
                'Submitted.Well': 'Submitted_Well',
                'sample.qc.missing.rate': 'sample_qc_missing_rate',
                'heterozygosity.pc.corrected': 'heterozygosity_pc_corrected',
                'het.missing.outliers': 'het_missing_outliers',
                'putative.sex.chromosome.aneuploidy': 'putative_sex_chromosome_aneuploidy',
                'in.kinship.table': 'in_kinship_table',
                'excluded.from.kinship.inference': 'excluded_from_kinship_inference',
                'excess.relatives': 'excess_relatives',
                'in.white.British.ancestry.subset': 'in_white_British_ancestry_subset',
                'used.in.pca.calculation': 'used_in_pca_calculation',
                'in.Phasing.Input.chr1_22': 'in_Phasing_Input_chr1_22',
                'in.Phasing.Input.chrX': 'in_Phasing_Input_chrX',
                'in.Phasing.Input.chrXY': 'in_Phasing_Input_chrXY'})

kt = kt.annotate(['Cluster_CR = Cluster_CR.toFloat()',
                  'dQC = dQC.toFloat()',
                  'Internal_Pico_ng_uL = Internal_Pico_ng_uL.toFloat()',
                  'is_inferred_female = Inferred_Gender == "F"',
                  'is_submitted_female = Submitted_Gender == "F"',
                  'X_intensity = X_intensity.toFloat()',
                  'Y_intensity = Y_intensity.toFloat()',
                  'sample_qc_missing_rate = sample_qc_missing_rate.toFloat()',
                  'heterozygosity = heterozygosity.toFloat()',
                  'heterozygosity_pc_corrected = heterozygosity_pc_corrected.toFloat()',
                  'het_missing_outliers = het_missing_outliers == "1"',
                  'putative_sex_chromosome_aneuploidy = putative_sex_chromosome_aneuploidy == "1"',
                  'in_kinship_table = in_kinship_table == "1"',
                  'excluded_from_kinship_inference = excluded_from_kinship_inference == "1"',
                  'excess_relatives = excess_relatives == "1"',
                  'in_white_British_ancestry_subset = in_white_British_ancestry_subset == "1"',
                  'used_in_pca_calculation = used_in_pca_calculation == "1"',
                  'in_Phasing_Input_chr1_22 = in_Phasing_Input_chr1_22 == "1"',
                  'in_Phasing_Input_chrX = in_Phasing_Input_chrX == "1"',
                  'in_Phasing_Input_chrXY = in_Phasing_Input_chrXY == "1"'])

kt = kt.annotate(['PC{0:d} = PC{0:d}.toFloat()'.format(i) for i in xrange(1, 41)])

kt_ancestry = hc.import_table('gs://ukb31063-mega-gwas/eur_selection/ukb31063_eur_samples.tsv', key='userId', quote='"')
kt_ancestry = kt_ancestry.select(['userId', 'eur_select'])
kt = kt.join(kt_ancestry, how='left')
kt = kt.annotate('is_in_ancestry_subset = eur_select == "TRUE"')

kt_phenotypes = hc.import_table('gs://phenotype_31063/ukb31063.raw.csv', delimiter='","')
kt_phenotypes = kt_phenotypes.annotate(["s = `\"eid`.replace('^\"', '')", "age = `21022-0.0`.toInt()"]) 
kt_phenotypes = kt_phenotypes.select(['s', 'age'])
kt_phenotypes = kt_phenotypes.key_by('s')

kt = kt.join(kt_phenotypes, how='left')
kt = kt.annotate('is_in_phenotype_data = isDefined(age)')
kt = kt.drop(['Inferred_Gender', 'Submitted_Gender', 'eur_select'])

kt.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.sample_qc.kt', overwrite=True)
