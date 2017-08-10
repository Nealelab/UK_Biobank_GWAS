from hail import *

hc = HailContext()

pipelines = 'gs://ukbb_association/ukb1859/pipelines/'
map_path = 'gs://ukbb-gwas-results/ukb1859/results_map.tsv'

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
	print('phenotypes: ', phenotypes)
	for group in phenotypes:
		update_result_map(map_path, current_map, phenotypes, group)
