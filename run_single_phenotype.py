from __future__ import print_function
import argparse
import os
from pprint import pprint
from hail import *

parser = argparse.ArgumentParser(description='Run association on single UKBB phenotype.')
parser.add_argument('--bgens', type=str, required=True, help='Google storage path to UKBB .bgen file(s) (can use "*" to regex match).')
parser.add_argument('--sample', type=str, required=True, help='Google storage path to application-specific .sample file.')
parser.add_argument('--fam', type=str, required=True, help='Google storage path to headerless, application-specific .fam file.')
parser.add_argument('--withdrawn', type=str, required=True, help='Google storage path to headerless, application-specific list of sample IDs withdrawn from the study.')
parser.add_argument('--phenotype', type=str, required=True, help='Google storage path to two-column, headerless tsv file where the first column is application-specific sample IDs and the second column is the phenotype (quantitative or 1/0 for case/control; "NA" to denote missing).')
parser.add_argument('--working-directory', type=str, required=True, help='Google storage path to working directory to write results and intermediate files.')
args = parser.parse_args()

GWAS_VARIANTS = 'gs://ukbb_association/gwas_variants.vds'
SAMPLE_QC_TSV = 'gs://ukbb_association/ukb_sqc_v2.tsv'

SAMPLE_KT = os.path.join(args.working_directory, 'sample.kt')
RESULTS_TSV = os.path.join(args.working_directory, 'variants.results.tsv')

hc = HailContext()

kt_fam = (
	hc
	.import_table(args.fam, delimiter=' ', no_header=True)
	.select('f0')
	.rename({'f0': 'id'})
	.indexed()
	.key_by('index')
)

kt_qc = (
	hc
	.import_table(
		SAMPLE_QC_TSV,
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
	.indexed()
	.key_by('index')
)

kt_withdrawn = (
    hc
    .import_table(args.withdrawn, no_header=True)
    .rename({'f0': 'sample'})
    .key_by('sample')
    .annotate('withdrawn = true')
)

kt_phenotype = hc.import_table(args.phenotype, no_header=True, key='f0', types={'f0': TString(), 'f1': TDouble()}).rename({'f1': 'phenotype'})

(
    kt_fam
    .join(kt_qc, how='left')
    .key_by('id')
    .drop('index')
    .annotate([
        'in_white_British_ancestry_subset = `in.white.British.ancestry.subset` == "1"',
        'used_in_pca_calculation = `used.in.pca.calculation` == "1"',
        'excess_relatives = `excess.relatives` == "1"',
        'putative_sex_chromosome_aneuploidy = `putative.sex.chromosome.aneuploidy` == "1"',
        'isFemale = `Inferred.Gender` == "F"'
    ])
    .join(kt_withdrawn, how='left')
    .annotate('withdrawn = if (isDefined(withdrawn)) true else false')
    .filter("""
		id != "-1" &&
		id != "-2" &&
		id != "-3" &&
		!withdrawn &&
		in_white_British_ancestry_subset &&
		used_in_pca_calculation &&
		!excess_relatives &&
		!putative_sex_chromosome_aneuploidy
	""")
	.select([
    	'id', 
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
    .join(kt_phenotype, how='inner')
	.write(SAMPLE_KT, overwrite=True)
)

vds = (
	hc
	.import_bgen(path=args.bgens, sample_file=args.sample, tolerance=0.2, min_partitions=None)
    .annotate_variants_expr('va = {}')
    .annotate_variants_vds(hc.read(GWAS_VARIANTS).annotate_variants_expr('va.keep = true'), root='va')
    .filter_variants_expr('isDefined(va.keep)')
    .annotate_variants_expr('va = select(va, rsid)')
    .annotate_samples_table(hc.read_table(SAMPLE_KT), root='sa')
    .linreg3(
    	ys=['sa.phenotype'], 
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
    .annotate_variants_expr("""
        va.results = {
    		nCompleteSamples: va.linreg.nCompleteSamples,
    		AC: va.linreg.AC,
    		ytx: va.linreg.ytx[0],
    		beta: va.linreg.beta[0],
    		se: va.linreg.se[0],
    		tstat: va.linreg.tstat[0],
    		pval: va.linreg.pval[0]
        }
    """)
)

n_variants = vds.count_variants()
n_samples, n_defined, n_missing = vds.query_samples([
	'samples.count()',
	'samples.filter(s => isDefined(sa.phenotype)).count()',
	'samples.filter(s => isMissing(sa.phenotype)).count()'
])

print('')
print('nVariants: ', '{:,}'.format(n_variants))
print('nSamples: ', '{:,}'.format(n_samples))
print('nSamples, defined phenotype: ', '{:,}'.format(n_defined))
print('nSamples, missing phenotypes: ', '{:,}'.format(n_missing))

vds.export_variants(RESULTS_TSV, 'variant = v, rsid = va.rsid, va.results.*')
