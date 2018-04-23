
from hail import *
hc = HailContext()

vds_vep = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.vep.vds')

vds_vep = vds_vep.annotate_variants_expr("""va.consequence = va.vep.most_severe_consequence,
                                            va.consequence_category = if (["frameshift_variant",
                                                                        "splice_acceptor_variant",
                                                                        "splice_donor_variant",
                                                                        "stop_gained",
                                                                        "transcript_ablation"].toSet().contains(va.vep.most_severe_consequence)) "ptv"
                                                                   else if (["inframe_deletion",
                                                                             "inframe_insertion",
                                                                             "missense_variant",
                                                                             "protein_altering_variant",
                                                                             "splice_region_variant",
                                                                             "start_lost",
                                                                             "stop_lost",
                                                                             "transcript_amplification"].toSet().contains(va.vep.most_severe_consequence)) "missense"
                                                                   else if (["incomplete_terminal_codon_variant",
                                                                             "stop_retained_variant",
                                                                             "synonymous_variant"].toSet().contains(va.vep.most_severe_consequence)) "synonymous"
                                                                   else if (["3_prime_UTR_variant",
                                                                             "5_prime_UTR_variant",
                                                                             "coding_sequence_variant",
                                                                             "downstream_gene_variant",
                                                                             "feature_elongation",
                                                                             "feature_truncation",
                                                                             "intergenic_variant",
                                                                             "intron_variant",
                                                                             "mature_miRNA_variant",
                                                                             "NMD_transcript_variant",
                                                                             "non_coding_transcript_exon_variant",
                                                                             "non_coding_transcript_variant",
                                                                             "regulatory_region_ablation",
                                                                             "regulatory_region_amplification",
                                                                             "regulatory_region_variant",
                                                                             "TFBS_ablation",
                                                                             "TFBS_amplification",
                                                                             "TF_binding_site_variant",
                                                                             "upstream_gene_variant"].toSet().contains(va.vep.most_severe_consequence)) "non_coding"
                                                                   else "NA" """)


vds_vep.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.imputed_v3.consequences_categorized.vep.vds', overwrite=True)
