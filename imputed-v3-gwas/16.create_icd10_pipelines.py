
from __future__ import print_function
from hail import *
hc = HailContext()

kt = hc.import_table('gs://phenotype_31063/ukb31063.raw.csv', delimiter='","')

icd10_cols = [x for x in kt.columns if x.startswith('41202')]

kt = kt.select(['\"eid'] + icd10_cols)
kt = kt.annotate("s = `\"eid`.replace('^\"', '')")
kt = kt.drop('\"eid')
kt = kt.key_by('s')

kt = kt.annotate('set_of_codes = [`{}`[:3]].toSet()'.format(icd10_cols[0]))
for x in icd10_cols[1:]:
    kt = kt.annotate('set_of_codes = set_of_codes.add(`{}`[:3])'.format(x))
kt = kt.annotate('set_of_codes = set_of_codes.remove("")')
kt = kt.select(['s', 'set_of_codes'])

kts = {
    'both_sexes': kt.join(hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.kt'), how='inner').cache(),
    'female': kt.join(hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.females.kt'), how='inner').cache(),
    'male': kt.join(hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.males.kt'), how='inner').cache()
}

codes = {
    'both_sexes': [x[0] for x in kts['both_sexes'].query('set_of_codes.flatMap(x => x.toArray()).counter()').iteritems() if x[1] >= 50],
    'female': [x[0] for x in kts['female'].query('set_of_codes.flatMap(x => x.toArray()).counter()').iteritems() if x[1] >= 50],
    'male': [x[0] for x in kts['male'].query('set_of_codes.flatMap(x => x.toArray()).counter()').iteritems() if x[1] >= 50]    
}

for sex, code_list in codes.iteritems():

    code_groups = [code_list[i:(i+110)] for i in range(0, len(code_list), 110)]
    current_pipeline = 0

    for cg in code_groups:

        print('Sex: {}'.format(sex))
        print('Pipeline: {:}'.format(current_pipeline))
        print('Codes: ', cg)

        kt_pipeline = kts[sex].annotate(['{0} = set_of_codes.contains("{0}")'.format(x) for x in cg])
        kt_pipeline = kt_pipeline.drop('set_of_codes')
        kt_pipeline.write('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/icd10/ukb31063.{0}.icd10.pipeline.{1:}.kt'.format(sex, current_pipeline), overwrite=True)

        current_pipeline += 1
