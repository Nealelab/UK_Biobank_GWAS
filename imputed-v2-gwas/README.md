# UK Biobank GWAS

# Table of Contents
* [Goals](#goals) 
* [Files](#files)
* [UK Biobank Updates](#uk-biobank-updates)
* [Phenotypes and applications](#phenotypes-and-applications)
  * [Phenotype output](#phenotype-output)
  * [Phenotype to Genotype linking](#phenotype-to-genotype-linking)
* [Sample and Variant QC](#sample-and-variant-qc)
  * [Sample QC](#sample-qc)
  * [Genotype QC](#genotype-qc)
* [Association in Hail](#association-in-hail)
  * [Association Model](#association-model)
  * [Single phenotype association workflow](#single-phenotype-association-workflow)
* [Summary stat output](#summary-stat-output)
* [Quick Manhattan and QQ plot](#quick-manhattan-and-qq-plot)

## Goals:

 * Collect UK Biobank phenotypes from collaborating applications 
 * When necessary, convert phenotypes into clean case/control, true/false, or quantitative values using the PHESANT algorithm
 * Use the same sample and variant QC across all phenotypes
 * Run SNP association on UKBB imputed dosage BGEN files using hail on the google cloud platform
 * Provide per-SNP summary stats for all approved phenotypes
 * GWAS Results available for viewing and download via https://biobankengine.stanford.edu


## Files

  * [1_merge_mfi.sh](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/1_merge_mfi.sh) - merge UKBB mfi files
  * [2_fam_sqc_merge.R](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/fam_sqc_merge.R) - merge ukb_sqc_v2.txt to application-specific .fam file
  * [3_make_sample_qc_table.py](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/3_make_sample_qc_table.py) - make application-specific qc key table
  * [4_build_pipelines.py](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/4_build_pipelines.py) - set up block of phenotypes to be analyzed
  * [5_make_variant_annotation_vds.py](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/5_make_variant_annotation_vds.py) - create sites-only VDS of all SNPs
  * [6_filter_gwas_variants.py](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/6_filter_gwas_variants.py) - filter sites-only VDS to QC'ed SNPs
  * [7_run_linreg3.py](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/7_run_linreg3.py) - run association with linreg3 command
  * [8_export_results.py](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/8_export_results.py) - export results to .tsv file
  * [Manhattan_plot.R](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/Manhattan_plot.R) - create quick Manhattan plot in R from summary stat file
  * [PHESANT_pipeline.pdf](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/PHESANT_pipeline.pdf) - diagram of PHESANT phenotype curation strategy
  * [QQ_plot.R](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/QQ_plot.R) - create quick QQ plot in R from summary stat file
  * [run_single_phenotype.py](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/run_single_phenotype.sh) - run a single phenotype through Hail on the Google cloud platform      

## UK Biobank Updates

**July 27th, 2017 - Errors in imputation identified**
  * Imputation of non-HRC SNPs was mis-mapped. Email details below:

*Dear researcher,  
We have identified a problem with the UK Biobank imputed data and which has come to light following discussion via the UKB-GENETICS mail list.  
This problem relates to the imputed data and does not affect the genotyped data from the Affymetrix array. 
The genetic data was imputed using two different reference panels.  
The Haplotype Reference Consortium (HRC) panel was used as first choice option, but for SNPs not in that reference panel the UK10K + 1000 Genomes panel was used.  
The problem arose in the second set of imputed data from the UK10K + 1000 Genomes panel.  
The genotypes at these SNPs are imputed correctly, but have not been recorded as having the correct genome position in the files.
We have established that the imputed data from the HRC panel is not affected and has the correct positions.
This is about ~40M sites and will include the majority of the common SNPs i.e. sites most likely to show genetic associations.
These sites are readily identified since the HRC site list is public
http://www.haplotype-reference-consortium.org/site  
The problem is not easy to fix post-hoc, so we intend to re-impute the data from the UK10K + 1000 Genomes panel and re-release the imputed data.  
For now we recommend that researchers focus exclusively on SNPs in the HRC panel, or work with the directly genotyped data until the new release is available.  
We will progress the re-imputation as quickly as we can and expect to release a new version of the imputed files ideally in September.
We will send more details about this data release and confirm timelines in due course. 
We can only apologise that this error was not identified during the QA review and do not underestimate the frustration this will cause for the research community.  
If you have any questions about this issue please send them to the UKB-GENETICS mail list to https://www.jiscmail.ac.uk/cgi-bin/webadmin?A0=UKB-GENETICS*

*Kind regards,  
UK Biobank & the Access Team*

**July 26th, 2017 - Samples withdrawn from UK Biobank**
  * 23 samples have withdrawn consent for use of their data
  * 8 samples listed in imputed data sample file
  * 6 samples identified as being part of our QC positive sample set

*Dear Researcher,  
As you are aware, participants are free to withdraw from UK Biobank at any time and request that their data no longer be used.  
Since our last review, some participants involved with Application [XXXX] have requested that their data should no longer be used.  
We have attached a file containing the anonymised IDs of these participants and any others who have withdrawn previously.  
Please remove the corresponding records from your unpublished analyses.  
Note that it is possible that this list contains IDs which you have never received as they may have withdrawn before your final dataset was generated.*
 
*Regards,  
The UK Biobank Access Team*

## Phenotypes and applications

**Phenotype Collection Strategy**
  * Collect as many phenotypes as possible from collaborating UK Biobank applications
  * Ensure all applications permit GWAS analysis and public release of summary statistics
  * Add participating Neale Lab members to all applications

**Phenotype Curation Strategy**
  * Auto-curation of phenotypes using PHESANT: 
    * Source repository: https://github.com/MRCIEU/PHESANT
    * Customized PHESANT repository: https://github.com/astheeggeggs/PHESANT
      * Transformed phenotypes are written to disk
      * Altered .log file for easier creation of downstream plotting and summary files
      * Removed additional analysis and plotting functionality not required for our purposes
      * Altered the variable-list file `variable-info/outcome_info_final.tsv` to restrict to the phenotypes we wish to analyze - see 'EXCLUDED' column
      * Added data codings tables from UK biobank data showcase to allow mapping from original encoding in phenotype file (e.g. 1 = 'bungalow', 2 = 'flat' etc)
      * Added additional scripts for creation of simple test plots and summary tables
      * [PHESANT phenotype curation strategy](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/PHESANT_pipeline.pdf)

**Example PHESANT script**
```
Rscript phenomeScan.r \
--phenofile="../../ukb[XXXX].tsv" \
--variablelistfile="../variable-info/outcome_info_final.tsv" \
--datacodingfile="../variable-info/data-coding-ordinal-info.txt" \
--userId="userId" \
--resDir="../../" \
--log="ukb[XXXX]_output" \
```

generates an output file. Finally, run `get_hists_and_notes` followed by `include_PHESANT_reassignment_names` in `summarise_phenoypes.r` to generate a summary file of the phenotypes for analysis. Note: the sample QC file `ukb[XXXX]_qc.tsv` and participant withdrawal list `w[XXXX]_20170726.csv` for each application is required to ensure that numbers in the summary file reflect the QC pipeline described below.

### Phenotype output
  * Output file: **PHESANT ukb[XXXX]_output.tsv** (full phenotype file with rows=sample ID and columns=phenotype ID)
  * Summary file: **ukb[XXXX]_phenosummary_final.tsv**
    * row number - numeric field matching UK Biobank data showcase field
    * Field - short description of phenotype
    * N.non.missing - number of non-missing QC positive samples
    * N.missing - number of missing QC positive samples
    * N.cases - number of QC positive samples responding affirmatively to phenotype designation (NA if quantitative) 
    * N.controls - number of QC positive samples responding negatively to phenotype designation (NA if quantitative)
    * Notes - extended description of phenotype
    * PHESANT.notes - categorizations of phenotype by PHESANT algorithm
    * PHESANT.reassignments - changes to phenotype values by PHESANT
  * Example phenosummary_final.tsv file
```
                             Field N.non.missing N.missing N.cases N.controls Notes                       PHESANT.notes                                                                         PHESANT.reassignments
46       Hand grip strength (left)        335821      1378      NA         NA Left grip strength                                       46_0|| INTEGER || CONTINUOUS || IRNT ||                  <NA>
47      Hand grip strength (right)        335842      1357      NA         NA Right grip strength                                       47_0|| INTEGER || CONTINUOUS || IRNT ||                  <NA>
20116_1   Smoking status: Previous        336024      1175  118419     217605 current/past smoking status 20116_0 || CAT-SINGLE || CAT-SINGLE-BINARY-VAR: 1  || Inc(>=10): 1(173099) ||                  <NA>
```

### Phenotype to Genotype linking
  * The main link files from application to UK Biobank imputed dataset are the **.fam** and **.sample** file
  * Both files have the same IDs
  	* The .fam file is ordered the same way as the ukb_sqc_v2.txt file
  	* The .sample file is ordered the same way as the .bgen file
  * [R Script](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/2_fam_sqc_merge.R) for linking application-specific .fam file and ukb_sqc_v2.txt file

**Current list of participating UK Biobank Applications**

**Completed associations**
  * 18597 - "Methodological extensions to estimate genetic heritability and shared risk factors for phenotypes of the UK Biobank"
    * Primary investigators: Benjamin Neale / Verneri Anttila
  * 11898 - "Genetic studies of anthropometric traits and methods for analysis of multiple phenotypes"
    * Primary investigator: Joel Hirschhorn

**Ongoing analyses**
  * 24983 - "Generating effective therapeutic hypotheses from genomic and hospital linkage data"
    * Primary investigator: Manuel Rivas
  * 11425 - "The Social Science Genetic Association Consortium"
    * Primary investigator: Daniel Benjamin
  * 32568 - "Phenomewide Heritability Analysis"
    * Primary investigator: Jordan Smoller

## Sample and Variant QC

### Sample QC

**Primary sample QC parameters for GWAS from ukb_sqc_v2.txt file:**
  * in.Phasing.Input.chr1_22==1
  * in.white.British.ancestry.subset==1
  * used.in.pca.calculation==1
  * excess.relatives==0
  * putative.sex.chromosome.aneuploidy==0

**Additional QC parameters**
  * Samples withdrawn un UK Biobank update = 8 
  * Samples redacted = 3 ([-3,-2,-1] in the sample ID) 

**Pre/post QC sample counts**
  * Imputed samples removed from QC file = 151180
  * Imputed samples retained in QC file = 337199
  * NOTE: all samples retained are in the .bgen files
  * NOTE: The ukb_sqc_v2.txt file has more samples than the .bgen files, but the same number of samples as the application specific .sample file

**Description of inclusion parameters:**

```
het.missing.outliers		      (0/1)	(no/yes) Indicates samples identified as outliers in heterozygosity and missing rates, which indicates poor-quality genotypes for these samples.
putative.sex.chromosome.aneuploidy    (0/1)	(no/yes) Indicates samples identified as putatively carrying sex chromosome configurations that are not either XX or XY. These were identified by looking at average log2Ratios for Y and X chromosomes. See genotype QC documentation for details.
in.kinship.table		      (0/1)	(no/yes) Indicates samples which have at least one relative among the set of genotyped individuals. These are exactly the set of samples that appear in the kinship table. See genotype QC documentation for details.
excluded.from.kinship.inference	      (0/1)	(no/yes) Indicates samples which were excluded from the kinship inference procedure. See genotype QC documentation for details.
excess.relatives		      (0/1)	(no/yes) Indicates samples which have more than 10 putative third-degree relatives in the kinship table.
in.white.British.ancestry.subset      (0/1)	(no/yes) Indicates samples who self-reported 'White British' and have very similar genetic ancestry based on a principal components analysis of the genotypes. See genotype QC documentation for details.
used.in.pca.calculation		      (0/1)	(no/yes) Indicates samples which were used to compute principal components. All samples with genotype data were subsequently projected on to each of the 40 computed components. See genotype QC documentation for details.
in.Phasing.Input.chr1_22	      (0/1)	(no/yes) Indicates sample was in the input for phasing of chr1-chr22.
in.Phasing.Input.chrX		      (0/1)	(no/yes) Indicates sample was in the input for phasing of chrX.
in.Phasing.Input.chrXY		      (0/1)	(no/yes) Indicates sample was in the input for phasing of chrXY.

Quick descriptives of categories:

table(in.Phasing.Input.chr1_22)
     0      1 
   969 487409 

table(het.missing.outliers)
     0 
487409 
table(putative.sex.chromosome.aneuploidy)
     0      1
486757    652
table(in.kinship.table)
     0      1 
339678 147731 
table(excluded.from.kinship.inference)
     0      1 
487400      9 
table(excess.relatives)
     0      1 
487221    188 
table(in.white.British.ancestry.subset)
     0      1 
 78437 408972 
table(used.in.pca.calculation)
     0      1 
 80190 407219 

Total samples in imputed BGEN files: 487409
Withdrawn samples (remove): 8
Redacted samples (remove): 3
In white, British ancestry subset (keep): 408972
Used in PCA calculation (keep): 407219
Marked as having excess relatives (remove): 188
Marked as sex chromosome aneuploidies (remove): 652
Samples remaining: 337199
```

### Genotype QC

**Primary genotype QC parameters for inclusion to GWAS from ukb_sqc_v2.txt file**
  * Autosomal SNPs
  * SNPs present in HRC imputation file: `../imputed/resources/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab`
  * pHWE > 1e-10 in 337199 QC positive samples
  * call rate > 0.95 in 337199 QC positive samples
  * UKBB INFO score > 0.8
  * AF > 0.001 and AF < 0.999 (AF = alternate allele frequency in 337199 QC positive samples)
  * NOTE: UKBB INFO score available in per chromosome files: ukb_mfi_chr*_v2.txt

```
Quick descriptive of genotype inclusion counts:

QC positive sample subset (n = 337199):
- Total Variants: 92,693,895
- in HRC imputation: 39,131,578
- w/ pHWE > 1e-10: 92,471,744
- w/ call rate > 95%: 92,693,895
- w/ 0.001 < MAF < 0.999: 16,047,294
- w/ UKBB INFO score > 0.8: 29,447,617
- w/ all of above: 10,894,596
```

## Association in Hail

SNP association is performed on UKBB imputed dosage BGEN files using hail on a google cloud platform.

* [Hail](https://hail.is/) scalable genetic analyses
* [Google Cloud Platform](https://cloud.google.com/) cluster computation and data storage
* [Cloud Tools](https://github.com/Nealelab/cloudtools) script submission to compute clusters

This allows us to take advantage of large cluster computing and parallel processing via Apache Spark

### Association model
 * **To expedite our initial GWAS, we implemented a linear regression model on all phenotypes**
 * Linear regression model covariates were inferred sex and the first 10 PCs (taken directly from ukb_sqc_v2.txt)
 * Hail's [linreg3] command allowed us to run associations on blocks of phenotypes with the same missingness structure at the same time to speed up analysis
 * However, memory constraints limited our blocks to:
    * 110 phenotypes when the missingness structure of all phenotypes is the same
    * 37 phenotypes when each phenotype has a different missingness structure
 * A critical consideration in designing the linreg3 pipeline (step 7 below) was to avoid commands that would invoke shuffles in Spark 
 * The following modifications were made to the Spark configuration while running linreg3 (step 7 below):
    * spark.memory.storageFraction=0.5,spark.memory.fraction=0.4,spark.broadcast.blockSize=128m,spark.dynamicAllocation.enabled=false,spark.executor.memory=38g,spark.executor.cores=8

### Workflow for running associations in Hail for phenotypes across multiple applications:
   1) Merge chromosome-specific `ukb_mfi_chr*_v2.txt` files into one `ukb_mfi_v2.tsv` file with a chromosome field
   2) Merge `ukb_sqc_v2.txt` sample QC file with application-specific sample IDs using application-specific .fam files
   3) Create application-specific Hail keytables with application-specific sample ID, inferred sex (from `ukb_sqc_v2.txt`) and the first 10 PCs (also from `ukb_sqc_v2.txt`), subset to the set of 337,199 samples outlined above
   4) Build application-specific 'pipelines', where each pipeline is defined by a Hail keytable with (for the subset of 337,199 samples):
      * application-specific sample ID
      * inferred sex
      * PCs 1-10
      * A block of phenotypes capable of being run through linreg3 in one pass of reading in BGEN -> writing out results keytable
   5) Create a sites-only VDS for all 92,693,895 sites in the imputed data with:
      * rsid
      * info score (from `ukb_mfi_v2.tsv`)
      * flag indicating if site is in HRC 1.1 release (`HRC.r1-1.GRCh37.wgs.mac5.sites.tab` file from http://www.haplotype-reference-consortium.org/site)
      * variant QC metrics from Hail's `variant_qc()` method, calculated on the subset of 337,199 samples
   6) Create a filtered sites-only VDS with 10,894,596 variants, that satisfy the criteria:
      * HRC site
      * info score > 0.8
      * 0.001 < AF < 0.999
      * pHWE > 1e-10
      * callRate > 0.95
   7) Run associations using linreg3, where for each pipeline (see step 4) in each application, we:
      * Read in the BGENs
      * Use the filtered sites-only VDS created in step (6) to filter variants
      * Annotate samples with the pipeline-specific keytable
      * Call the linreg3 command once for each group of phenotypes in the pipeline with the same missingness structure
      * Write out association results to a Hail keytable
   8) Export a separate tsv results file (described below) and LDSC summary stat file for each phenotype in each pipelines across all applications

 #### Scripts associated with each step in the workflow above:
   1) `1_merge_mfi.sh`
   2) `2_fam_sqc_merge.R`
   3) `3_make_sample_qc_table.py`
   4) `4_build_pipelines.py`
   5) `5_make_variant_annotation_vds.py`
   6) `6_filter_gwas_variants.py`
   7) `7_run_linreg3.py`
   8) `8_export_results.py`

### Single phenotype association workflow
Use the script `run_single_phenotype.py` to run an association on a single phenotype using Hail v0.1, using the same set of 337,199 UK Biobank samples and 10,894,596 variants described above. The script requires 6 parameters:
 * `--bgens` - Google storage path to UKBB .bgen file(s) (can use "\*" to regex match).
 * `--sample` - Google storage path to application-specific .sample file.
 * `--fam` - Google storage path to headerless, application-specific .fam file.
 * `--withdrawn` - Google storage path to headerless, application-specific list of sample IDs withdrawn from the study.
 * `--phenotype` - Google storage path to two-column, headerless tsv file where the first column is application-specific sample IDs and the second column is the phenotype (quantitative or 1/0 for case/control; "NA" to denote missing).
 * `--working-directory` - Google storage path to working directory to write results and intermediate files.

Here's an example of running this script using the [cloudtools](https://github.com/nealelab/cloudtools) interface:
```
cluster submit my-cluster run_single_phenotypes.py --args "
  --bgens gs://<bgens-bucket>/ukb_imp_chr*_v2.bgen
  --sample gs://mybucket/myapplication.sample
  --fam gs://mybucket/myapplication.fam
  --withdrawn gs://mybucket/myapplication_withdrawn.txt
  --phenotype gs://mybucket/myphenotype.tsv
  --working-directory gs://mybucket/
"
```

This script will write a file `variants.results.tsv` to the specified `--working-directory`, with the format outlined below in the "SNP summary stat file" section.

The script expects the `--phenotype` file to be a two-column tsv file with no header, where the first column contains the application-specific sample ID and the second column contains the phenotype of interest (either quantitative or 1/0 for case/control, with missing values coded as "NA").

## Summary stat output

**NOTES:**
 * Variants are listed in the form CHR:POS:REF:ALT, where the ALT allele is the effect allele in the model
 * The ALT allele is NOT always the minor allele, but the non-reference allele as stated in the UK Biobank ukb_mfi_chr*_v2.txt and our variants.tsv file


**QCed SNP information file (variants.tsv)**
  * variant (hg19) [CHROM:POS:REF:ALT]
  * rsid
  * info (UKBB INFO score)
  * AF (QC positive alternate allele frequency)
  * pHWE
  * callRate
 * Example SNP information: 
```	
variant		rsid		info		AF		pHWE		callRate
10:61334:G:A	rs183305313	8.46690e-01	4.82208e-03	4.77873e-03	1.00000e+00
10:69083:C:T	rs35418599	8.01306e-01	7.77212e-01	7.74927e-01	9.99961e-01
10:90127:C:T	rs185642176	9.84897e-01	8.62740e-02	2.77123e-01	1.00000e+00
10:94263:C:A	rs184120752	9.57328e-01	2.51928e-02	8.26055e-02	1.00000e+00
10:94426:C:T	rs10904045	9.90536e-01	3.96792e-01	3.54267e-01	1.00000e+00
10:94538:C:T	rs189409193	8.24138e-01	5.50862e-03	9.37271e-01	1.00000e+00
```

**SNP summary stat file ([FieldID]_assoc.tsv.gz)**
  * variant (hg19) [CHROM:POS:REF:ALT]
  * rsid
  * nCompleteSamples (non-missing samples)
  * AC (non-missing sample allele count)
  * ytx (case/control = dosage weighted alternate allele count in cases; quantitative = dosage weighted mean trait value among alternate allele carriers)
  * beta
  * se
  * pval
 * Example SNP summary stat: 
```
variant		rsid		nCompleteSamples	AC	ytx		beta		se		tstat		pval
5:43888254:C:T	rs13184706	953		4.17176e+01	5.64980e+01	-1.11569e-01	8.01312e-02	-1.39233e+00	1.64152e-01
5:43888493:C:T	rs58824264	953		9.03529e+00	1.30706e+01	-3.42168e-02	1.68596e-01	-2.02951e-01	8.39217e-01
5:43888556:T:C	rs72762387	953		4.86235e+01	7.81804e+01	1.31571e-01	7.44976e-02	1.76611e+00	7.77023e-02
5:43888648:C:T	rs115032754	953		3.77647e+01	5.40039e+01	-5.98780e-02	8.59590e-02	-6.96588e-01	4.86233e-01
5:43888690:C:G	rs147555725	953		5.87843e+00	9.80000e+00	1.98330e-01	2.11226e-01	9.38946e-01	3.48000e-01
5:43888838:G:C	rs13185925	953		7.21765e+01	1.04306e+02	-2.46341e-02	6.00665e-02	-4.10113e-01	6.81816e-01
```


## Quick Manhattan and QQ plot

**R Scripts to generate basic Manhattan and QQ plots on a single summary stat file**

  * [Manhattan_plot.R](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/Manhattan_plot.R)
    * Creates manhattan plot on SNPs with p < .001
    * Requires > 1 Gb memory to read in summary stat file
   * Example usage:
```
## generic example:
Rscript Manhattan_plot.R [summary stat file] [output name]

## actual example:
Rscript Manhattan_plot.R 2443.assoc.tsv.gz /home/user/graphs/2443

## generates: /home/user/graphs/2443_Manhattan.png file
```

  * [QQ_plot.R](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/QQ_plot.R)
    * Creates QQ plot on all SNPs                    
    * Requires > 1 Gb memory to read in summary stat file
   * Example usage:
```
## generic example:
Rscript QQ_plot.R [summary stat file] [output name]                                     

## actual example:
Rscript QQ_plot.R 2443.assoc.tsv.gz /home/user/graphs/2443        

## generates: /home/user/graphs/2443_QQ.png file
```


