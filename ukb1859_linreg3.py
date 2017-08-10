from hail import *
from pprint import pprint

hc = HailContext(min_block_size=0)

bgens = 'gs://ukbb-bgens/ukb_imp_chr*_v2.bgen'
sample = 'gs://ukbb_association/ukb1859/ukb1859.sample'
variant_list = 'gs://ukbb-cseed-ws/variants.vds'

pipelines = 'gs://ukbb_association/ukb1859/pipelines/'
results = 'gs://ukbb_association/ukb1859/results/'
map_path = 'gs://ukbb_association/ukb1859/results/result_map.tsv'

def reset_vds(pipeline):
    return (
        hc
        .import_bgen(path=bgens, sample_file=sample, tolerance=0.2, min_partitions=None)
        .annotate_samples_table(hc.read_table(pipelines + 'pipeline{}.kt'.format(pipeline)), root='sa')
        .annotate_variants_vds(hc.read(variant_list).annotate_variants_expr('va = true'), expr='va.keep = vds')
        .filter_variants_expr('isDefined(va.keep)')
        .annotate_variants_expr('va = drop(va, varid)') 
    )

def write_keytable(vds, pipeline):
    (
        vds
        .annotate_variants_expr('va = drop(va, linreg)')
        .variants_table()
        .write(results + 'results{}.kt'.format(pipeline), overwrite=True)
    )

def group_phenotypes(pipeline):
    phenotype_groups = []
    for field in hc.read_table(pipelines + 'pipeline{}.kt'.format(pipeline)).columns:
        if field in set(['ID', 'isFemale'] + ['PC{}'.format(i) for i in xrange(1,11)]):
            continue
        prefix = field.split('_', 1)[0]
        try:
            idx = [x[0] for x in phenotype_groups].index(prefix)
        except ValueError:
            phenotype_groups.append([prefix, [field]])
        else:
            phenotype_groups[idx][1].append(field)
            
    return phenotype_groups

def update_result_map(map_path, current_map, phenotypes, group):
    current_map.extend([[field, str(i), str(phenotypes.index(group)), str(group[1].index(field))] for field in group[1]])
    with hadoop_write(map_path) as f:
        f.write('\n'.join(['\t'.join(row) for row in current_map]))

with hadoop_write(map_path) as f:
    current_map = [['phenotype', 'pipeline', 'group', 'element']]
    f.write('\t'.join(current_map[0]) + '\n')

n_pipelines = 21

for i in xrange(n_pipelines):
    phenotypes = group_phenotypes(i)
    vds = reset_vds(i)
    first = True
    for group in phenotypes:
        if first:
            expr = 'va.results = [va.linreg]'
        else:
            expr = 'va.results = va.results.append(va.linreg)'
        vds = (
            vds
            .linreg3(
                ys=['sa.`{}`'.format(field) for field in group[1]],
                covariates=[
                    'sa.isFemale',
                    'sa.PC1',
                    'sa.PC2',
                    'sa.PC3',
                    'sa.PC4',
                    'sa.PC5',
                    'sa.PC6',
                    'sa.PC7',
                    'sa.PC8',
                    'sa.PC9',
                    'sa.PC10'
                ],
                use_dosages=True,
                variant_block_size=8
            )
            .annotate_variants_expr(expr)
        )
        first = False
        update_result_map(map_path, current_map, phenotypes, group)
    write_keytable(vds, i)
