
from hail import *
hc = HailContext()

kt = hc.import_table('gs://phenotype_31063/ukb_snp_qc.txt', delimiter=' ', missing='NA', min_partitions=16, no_header=True, comment='rs_id')

cols = (['rs_id',
         'affymetrix_snp_id',
         'affymetrix_probeset_id',
         'chromosome',
         'position',
         'allele1_ref',
         'allele2_alt',
         'strand',
         'array'] + 
        ['Batch_b{:03d}_qc'.format(i) for i in xrange(1, 96)] +
        ['UKBiLEVEAX_b{:d}_qc'.format(i) for i in xrange(1, 12)] +
        ['in_HetMiss',
         'in_Relatedness',
         'in_PCA'] +
        ['PC{:d}_loading'.format(i) for i in xrange(1, 41)] +
        ['in_Phasing_Input'])

kt = kt.rename({'f{:d}'.format(i): cols[i] for i in xrange(0, len(cols))})

kt = kt.annotate(['chromosome = chromosome.replace("23", "X").replace("24", "Y").replace("25", "X").replace("26", "MT")',
                  'position = position.toInt()',
                  'in_HetMiss = in_HetMiss == "1"',
                  'in_Relatedness = in_Relatedness == "1"',
                  'in_PCA = in_PCA == "1"',
                  'in_Phasing_Input = in_Phasing_Input == "1"'])

kt = kt.annotate(['Batch_b{0:03d}_qc = Batch_b{0:03d}_qc == "1"'.format(i) for i in xrange(1, 96)])
kt = kt.annotate(['UKBiLEVEAX_b{0:d}_qc = UKBiLEVEAX_b{0:d}_qc == "1"'.format(i) for i in xrange(1, 12)])
kt = kt.annotate(['PC{0:d}_loading = PC{0:d}_loading.toFloat()'.format(i) for i in xrange(1, 41)])
kt = kt.annotate('v = Variant(chromosome, position, allele1_ref, allele2_alt)')
kt = kt.key_by('v')

kt.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.genotype_snp_qc.kt', overwrite=True)
from pprint import pprint
pprint(kt.schema)
