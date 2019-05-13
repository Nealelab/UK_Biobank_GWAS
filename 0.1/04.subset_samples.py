
from hail import *
hc = HailContext()

kt = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.sample_qc.kt')

kt = kt.filter("""in_Phasing_Input_chr1_22 &&
                  in_Phasing_Input_chrX &&
                  is_in_phenotype_data &&
                  is_in_ancestry_subset &&
                  used_in_pca_calculation""")

kt_both = kt.select('s')
kt_both.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.kt', overwrite=True)
kt_both.export('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.txt')
print 'nSamples: {:,}'.format(kt_both.count())

kt_female = kt.filter('is_inferred_female')
kt_female = kt_female.select('s')
kt_female.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.females.kt', overwrite=True)
kt_female.export('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.females.txt')
print 'nFemales: {:,}'.format(kt_female.count())

kt_male = kt.filter('!is_inferred_female')
kt_male = kt_male.select('s')
kt_male.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.males.kt', overwrite=True)
kt_male.export('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.males.txt')
print 'nMales: {:,}'.format(kt_male.count())
