from __future__ import print_function
from hail import *
from pprint import pprint

hc = HailContext()

GWAS_VARIANTS = 'gs://ukbb_association/gwas_variants.vds'
VEP_VARIANTS = 'gs://ukbb_association/gwas_variants.vep.vds'

vds = hc.read(GWAS_VARIANTS)
n_variants = vds.count_variants()

print('')
print('nVariants: ', '{:,}'.format(n_variants))
pprint(vds.variant_schema)

(
	vds
	.annotate_variants_db('va.vep')
	.write(VEP_VARIANTS, overwrite=True)
)
