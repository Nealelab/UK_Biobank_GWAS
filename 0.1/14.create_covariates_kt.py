
from hail import *
hc = HailContext()

kt_qc = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.sample_qc.kt')
kt_pca = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.pca_scores.kt')

kt_qc = kt_qc.rename({'is_inferred_female': 'isFemale'})
kt_qc = kt_qc.drop(['PC{:d}'.format(i) for i in range(1, 41)])
kt_qc = kt_qc.annotate(['age_squared = pow(age, 2)',
                        'age_isFemale = age * isFemale.toInt()',
                        'age_squared_isFemale = pow(age, 2) * isFemale.toInt()'])
kt_qc = kt_qc.select(['s', 'isFemale', 'age', 'age_squared', 'age_isFemale', 'age_squared_isFemale'])

kt_pca = kt_pca.annotate(['PC{0:d} = sa.pca_scores.PC{0:d}'.format(i) for i in xrange(1, 21)])
kt_pca = kt_pca.select(['s'] + ['PC{:d}'.format(i) for i in xrange(1, 21)])
kt = kt_qc.join(kt_pca, how='inner') 

kt.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_covariates.both_sexes.kt', overwrite=True)

kt_female = kt.filter('isFemale')
kt_female = kt_female.drop(['isFemale', 'age_isFemale', 'age_squared_isFemale'])
kt_female.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_covariates.female.kt', overwrite=True)

kt_male = kt.filter('!isFemale')
kt_male = kt_male.drop(['isFemale', 'age_isFemale', 'age_squared_isFemale'])
kt_male.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_covariates.male.kt', overwrite=True)
