from hail import *

hc = HailContext()

phenotypes = 'gs://ukbb_association/ukb1859/ukb1859_output.tsv'
qc = 'gs://ukbb_association/ukb1859/ukb1859_qc.tsv'

kt_phenotypes = (
    hc
    .import_table(
        paths=phenotypes,
        impute=True,
        types={'userId': TString()},
        min_partitions=128
    )
    .rename({'userId': 'ID'})
    .key_by('ID')
    .drop(['age', 'sex'])
)

phenotype_groups = {}
for x in kt_phenotypes.columns[1:]:
    cols = x.split('_')
    try:
        phenotype_groups[cols[0]].append(x)
    except KeyError:
        phenotype_groups[cols[0]] = [x]

max_phenotypes = 110
max_linreg = 37
linreg_calls = 0
effective_phenotypes_run = 0
pipelines = []
next_pipeline = []

for group, fields in phenotype_groups.iteritems():
    blocks = [fields[i:(i+max_phenotypes)] for i in xrange(0, len(fields), max_phenotypes)]
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

kt_qc_path = 'gs://ukbb_association/ukb1859/qc.kt'
(
    hc
    .import_table(
        paths=qc,
        comment='SKIP',
        delimiter='\t',
        impute=True,
        types={'eid': TString()},
        min_partitions=16
    )
    .rename({'eid': 'ID'})
    .key_by('ID')
    .annotate([
        'in_white_British_ancestry_subset = `in.white.British.ancestry.subset` == 1',
        'used_in_pca_calculation = `used.in.pca.calculation` == 1',
        'excess_relatives = `excess.relatives` == 1',
        'putative_sex_chromosome_aneuploidy = `putative.sex.chromosome.aneuploidy` == 1',
        'isFemale = `Inferred.Gender` == "F"'
    ])
    .filter('in_white_British_ancestry_subset && used_in_pca_calculation && !excess_relatives && !putative_sex_chromosome_aneuploidy')
    .select([
        'ID', 
        'isFemale',
        'PC1',
        'PC2',
        'PC3',
        'PC4',
        'PC5',
        'PC6',
        'PC7',
        'PC8',
        'PC9',
        'PC10'
    ])
    .write(kt_qc_path, overwrite=True)
)

def subset_keytable(fields):
    return (
        hc
        .import_table(
            paths=phenotypes,
            impute=True,
            types={'userId': TString()},
            min_partitions=128
        )
        .rename({'userId': 'ID'})
        .key_by('ID')
        .select(['ID'] + fields)
        .join(hc.read_table(kt_qc_path), how='inner')
    )

for i, p in enumerate(pipelines):
    fields = [field for linreg_call in p for field in linreg_call]
    subset_keytable(fields).write('gs://ukbb_association/ukb1859/pipelines/pipeline{}.kt'.format(i), overwrite=True)
