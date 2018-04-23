
from __future__ import print_function
from hail import *
hc = HailContext()

kt = hc.import_table('gs://ukb31063-mega-gwas/curated-phenotypes/2018-04-06_ukb-finngen-pheno-for-analysis.tsv', key='eid', types={'eid': TString()}, impute=True)
kt = kt.rename({'eid': 's'})

kt_counts = hc.import_table('gs://ukb31063-mega-gwas/curated-phenotypes/2018-04-07_ukb-finngen-all-pheno-counts.tsv', key='NAME', 
                            types={'Ncase': TInt(), 'Nctrl': TInt(), 'Prop': TDouble()}, missing='')
exclude = set(kt_counts.query('NAME.filter(x => isDefined(HD_ICD_10)).collect()'))

traits = [x for x in kt.columns if x != 's' and x not in exclude]

kts = {
    'both_sexes': kt.join(hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.kt'), how='inner').cache(),
    'female': kt.join(hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.females.kt'), how='inner').cache(),
    'male': kt.join(hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.males.kt'), how='inner').cache()
}

queries = {
    'both_sexes': kts['both_sexes'].query(['`{}`.filter(s => s != 0).count()'.format(x) for x in traits]),
    'female': kts['female'].query(['`{}`.filter(s => s != 0).count()'.format(x) for x in traits]),
    'male': kts['male'].query(['`{}`.filter(s => s != 0).count()'.format(x) for x in traits])
}

selections = {
    'both_sexes': [x for i, x in enumerate(traits) if queries['both_sexes'][i] >= 50],
    'female': [x for i, x in enumerate(traits) if queries['female'][i] >= 50],
    'male': [x for i, x in enumerate(traits) if queries['male'][i] >= 50]
}

for k, table in kts.iteritems():

    kt_filter = table.select(['s'] + selections[k])
    cols = [x for x in kt_filter.columns if x != 's']
    groups = [cols[i:(i+110)] for i in xrange(0, len(cols), 110)]
    current_pipeline = 0

    for g in groups:

        print('Sex: {}'.format(k))
        print('Pipeline: {:}'.format(current_pipeline))
        print('Traits: ', g)

        (kt_filter.select(['s'] + g)
                  .write('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/finngen/ukb31063.{0}.finngen.pipeline.{1:}.kt'.format(k, current_pipeline), overwrite=True))

        print('Sex: {}'.format(k))
        print('Pipeline: {:,}'.format(current_pipeline))
        print('nSamples: {:,}'.format(hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/finngen/ukb31063.{0}.finngen.pipeline.{1:}.kt'.format(k, current_pipeline)).count()))

        current_pipeline += 1
