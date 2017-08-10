from __future__ import print_function
from pprint import pprint
import sys
from hail import *

hc = HailContext()

BUCKET = '...'

APPLICATION = sys.argv[1]

PREFIX = BUCKET + APPLICATION + '/' + APPLICATION + '_'

WITHDRAWN_SAMPLES = PREFIX + 'withdrawn.csv'
SAMPLE_QC_FILE = PREFIX + 'qc.tsv'
SAMPLE_QC_TABLE = PREFIX + 'qc.kt' 

kt_withdrawn = (
    hc
    .import_table(WITHDRAWN_SAMPLES, delimiter=',', no_header=True)
    .rename({'f0': 'sample'})
    .key_by('sample')
    .annotate('withdrawn = true')
)

(
    hc
    .import_table(
    	SAMPLE_QC_FILE,
    	types={
    		'PC1': TDouble(),
    		'PC2': TDouble(),
    		'PC3': TDouble(),
    		'PC4': TDouble(),
    		'PC5': TDouble(),
    		'PC6': TDouble(),
    		'PC7': TDouble(),
    		'PC8': TDouble(),
    		'PC9': TDouble(),
    		'PC10': TDouble()
    	}
    )
    .annotate([
        'sample = eid',
        'in_white_British_ancestry_subset = `in.white.British.ancestry.subset` == "1"',
        'used_in_pca_calculation = `used.in.pca.calculation` == "1"',
        'excess_relatives = `excess.relatives` == "1"',
        'putative_sex_chromosome_aneuploidy = `putative.sex.chromosome.aneuploidy` == "1"',
        'isFemale = `Inferred.Gender` == "F"'
    ])
    .key_by('sample')
    .join(kt_withdrawn, how='left')
    .annotate('withdrawn = if (isDefined(withdrawn)) true else false')
    .filter("""
        sample != "-1" &&
        sample != "-2" &&
        sample != "-3" &&
    	!withdrawn &&
    	in_white_British_ancestry_subset &&
    	used_in_pca_calculation &&
    	!excess_relatives &&
    	!putative_sex_chromosome_aneuploidy
    """)
    .select([
    	'sample', 
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
    .write(SAMPLE_QC_TABLE, overwrite=True)
)

kt = hc.read_table(SAMPLE_QC_TABLE)
n_samples = kt.count()

print('')
print('Application: ', APPLICATION)
print('nSamples: ', n_samples)
pprint(kt.schema)
