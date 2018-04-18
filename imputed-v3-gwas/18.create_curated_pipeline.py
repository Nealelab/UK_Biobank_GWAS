
from hail import *
hc = HailContext()

kt = hc.import_table('gs://ukb31063-mega-gwas/curated-phenotypes/phenos_curated_for_MEGA_GWAS.tsv', types={'ID': TString()}, impute=True, key='ID')
kt = kt.rename({'ID': 's'})

traits = [x for x in kt.columns if x != 's']

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
    (table.select(['s'] + selections[k])
          .write('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/curated/ukb31063.{}.curated.kt'.format(k), overwrite=True))
    print k + ': {:,}'.format(hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/curated/ukb31063.{}.curated.kt'.format(k)).count())
