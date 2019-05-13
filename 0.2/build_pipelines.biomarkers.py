
import sys
import hail as hl

sex = sys.argv[1]

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError(f'Invalid sex argument "{sex}" - must be one of {{"both_sexes", "female", "male"}}.')

ht_phenotypes = hl.read_table(f'gs://ukb31063/hail/ukb31063.biomarkers_gwas.{sex}.ht')
phenotypes = list(ht_phenotypes.row_value)

chunk_size = 10
groups = [phenotypes[i:(i + chunk_size)] for i in range(0, len(phenotypes), chunk_size)]

i = 0
for group in groups:
    ht = ht_phenotypes.select(*group)
    ht.write(f'gs://ukb31063-mega-gwas/biomarkers/pipelines/ukb31063.biomarkers_gwas.{sex}.pipeline_{i}.ht', overwrite=True)
    i += 1
