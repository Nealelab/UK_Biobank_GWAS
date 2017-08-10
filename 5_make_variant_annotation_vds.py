from __future__ import print_function
from pprint import pprint
from hail import *

hc = HailContext()

BUCKET = 'gs://ukbb_association/'
APPLICATION = 'test'

BGEN_FILES = 'gs://fc-9a7c5487-04c9-4182-b3ec-13de7f6b409b/imputed/ukb_imp_chr*_v2.bgen'
SAMPLE_FILE = BUCKET + APPLICATION + '/' + APPLICATION + '.sample'

MFI_FILE = BUCKET + 'ukb_mfi_v2.tsv.bgz'
MFI_TABLE = BUCKET + 'mfi.kt'

HRC_FILE = BUCKET + 'HRC.r1-1.GRCh37.wgs.mac5.sites.tab'
HRC_TABLE = BUCKET + 'hrc_autosomes.kt'

QC_TABLE = BUCKET + APPLICATION + '/' + APPLICATION + '_qc.kt'

MFI_HRC_VDS = BUCKET + 'mfi_hrc_annotations.vds'
VARIANT_VDS = BUCKET + 'all_variants.vds'

(
    hc
    .import_table(
        MFI_FILE,
        no_header=True,
        types={
            'f0': TString(),
            'f1': TString(),
            'f2': TInt(),
            'f3': TString(),
            'f4': TString(),
            'f5': TDouble(),
            'f6': TDouble()
        }
    )
    .rename({
        'f0': 'chr',
        'f1': 'rsid',
        'f2': 'pos',
        'f3': 'ref',
        'f4': 'alt',
        'f5': 'maf',
        'f6': 'info'
    })
    .annotate('chr = if (chr.length < 2) ["0",chr].mkString("") else chr')
    .annotate('v = Variant(chr, pos, ref, alt)')
    .key_by('v')
    .select([
    	'v',
    	'rsid',
    	'maf',
    	'info'
    ])
    .write(MFI_TABLE, overwrite=True)
)

(
    hc
    .import_table(
        HRC_FILE,
        types={
            '`#CHROM`': TString(),
            'POS': TInt(),
            'ID': TString(),
            'REF': TString(),
            'ALT': TString(),
            'AC': TInt(),
            'AN': TInt(),
            'AF': TDouble(),
            'AC_EXCLUDING_1000G': TInt(),
            'AN_EXCLUDING_1000G': TInt(),
            'AF_EXCLUDING_1000G': TDouble(),
            'AA': TString()
        }
    )
    .filter('`#CHROM` == "X"', keep=False)
    .annotate('chr = let x = `#CHROM` in if (x.length < 2) ["0",x].mkString("") else x')
    .annotate('v = Variant(chr, POS, REF, ALT)')
    .key_by('v')
    .select('v')
    .write(HRC_TABLE, overwrite=True)
)

(
	VariantDataset
	.from_table(hc.read_table(MFI_TABLE).select(['v', 'rsid', 'info']))
	.annotate_variants_table(hc.read_table(HRC_TABLE), expr='va.isHRC = table')
	.write(MFI_HRC_VDS, overwrite=True)
)

samples = hc.read_table(QC_TABLE).query('sample.collect()')

(
    hc
    .import_bgen(path=BGEN_FILES, sample_file=SAMPLE_FILE, tolerance=0.2, min_partitions=None)
    .annotate_variants_vds(hc.read(MFI_HRC_VDS), root='va')
    .filter_samples_list(samples)
    .variant_qc()
    .drop_samples()
    .write(VARIANT_VDS, overwrite=True)
)

vds = hc.read(VARIANT_VDS)
n_variants = vds.count_variants()

print('')
print('nVariants: ', '{:,}'.format(n_variants))
print('nSamples: ', '{:,}'.format(len(samples)))
pprint(vds.variant_schema)
