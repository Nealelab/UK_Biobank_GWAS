
import sys
import hail as hl

sex = sys.argv[1]
contig = sys.argiv[2]

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError(f'Invalid sex value "{sex}" - must be one of {{"both_sexes", "female", "male"}}.')
if contig not in set(['autosomes', 'chrX', 'chrXY']):
    raise ValueError(f'Invalid contig "{contig}" - must be one of {{"autosomes", "chrX", "chrXY"}}.')