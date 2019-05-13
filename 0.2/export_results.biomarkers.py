ht = hl.read_table(path)
for i, y in enumerate(phenotypes):
    ht_out = ht.select(
        variant=hl.delimit([
            ht['locus'].contig,
            ht['locus'].position,
            ht['alleles'][0],
            ht['alleles'][1]], delimiter=':'),
        )