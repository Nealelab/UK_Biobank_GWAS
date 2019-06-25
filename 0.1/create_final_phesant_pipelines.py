
from hail import *
hc = HailContext()

for sex in ['both_sexes', 'female', 'male']:

    kt = hc.import_table('gs://phenotype_31063/phesant_output_combined_{}.tsv.bgz'.format(sex),
                         key='userID', quote='"', impute=True, missing='',
                         types={'userID': TString()}).cache()

    codes = sorted([x for x in kt.columns if x != 'userID'])

    phenotype_groups = []
    for code in codes:
        prefix = code.split('_')[0]
        try:
            idx = [x[0] for x in phenotype_groups].index(prefix)
        except ValueError:
            phenotype_groups.append([prefix, [code]])
        else:
            phenotype_groups[idx][1].append(code)

    max_phenotypes = 110
    max_linreg = 37
    linreg_calls = 0
    effective_phenotypes_run = 0
    pipelines = []
    next_pipeline = []

    for group, codes in phenotype_groups:
        blocks = [codes[j:(j+max_phenotypes)] for j in range(0, len(codes), max_phenotypes)]
        for block in blocks:
            effective_phenotypes_run += 2 + len(block)
            linreg_calls += 1
            if (effective_phenotypes_run > max_phenotypes) or (linreg_calls > max_linreg):
                pipelines.append(next_pipeline)
                next_pipeline = []
                effective_phenotypes_run = 2 + len(block)
                linreg_calls = 1
            next_pipeline.append(block)
    pipelines.append(next_pipeline)

    current_pipeline = 0
    for pipeline in pipelines:

        pipeline_codes = [code for linreg_call in pipeline for code in linreg_call]

        print('Sex: {}'.format(sex))
        print('Pipeline: {:}'.format(current_pipeline))
        print('Codes: ', pipeline_codes)

        (kt.select(['userID'] + pipeline_codes)
           .rename({'userID': 's'})
           .write('gs://ukb31063-mega-gwas/phenotype-pipelines/final-phesant/ukb31063.{0}.phesant.pipeline.{1:}.kt'.format(sex, current_pipeline), overwrite=True))

        current_pipeline += 1
