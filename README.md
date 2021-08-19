### V3 Summary statistics are [now available for download on Amazon Web Services](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?usp=sharing) 

### [README](https://github.com/Nealelab/UK_Biobank_GWAS/tree/master/imputed-v2-gwas#readme) from the previous round of GWAS (n=337,199) [are available here](https://github.com/Nealelab/UK_Biobank_GWAS/tree/master/imputed-v2-gwas)


# Table of Contents
* [Updates](#updates) 
* [Change Log](#change-log)
* [imputed-v3 Phenotypes](#imputed-v3-phenotypes)
* [imputed-v3 Sample QC](#imputed-v3-sample-qc)
* [imputed-v3 Variant QC](#imputed-v3-variant-qc)
* [imputed-v3 Association model](#imputed-v3-association-model)

### Updates

With the re-release of UK Biobank genotype imputation (which we term imputed-v3), we have generated an updated set of GWAS summary statistics for the genetics community. 
  * Increased the number of phenotypes with application UKB31063 and addtl. custom curated phenotypes (see imputed-v3 Phenotypes)
  * More liberal inclusion of samples (see imputed-v3 Sample QC)
  * Inclusion of more SNPs (see imputed-v3 Variant QC)
  * Updates to our association model (imputed-v3 Association model)
Our largest change is that for all phenotypes, we have run a female-only and male-only GWAS along with the full set.

Information and scripts from the previous round of GWAS are available in the [imputed-v2-gwas](https://github.com/Nealelab/UK_Biobank_GWAS/tree/master/imputed-v2-gwas) subdirectory

Finally, the [0.1](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.1) and [0.2](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.2) script repositories refer to the version of Hail used to run the GWAS 


### Change log

Updates to the Rapid GWAS summary statistics or download Manifest will be recorded here:

 * _August, 2021_
   * AUTHOR: Raymond Walters (with help from Sam Bryant and Daniel Howrigan, and hat tip to Trevor manz for discovering the issue)
   * ISSUE: We recieved a couple reports that the round 2 biomarker GWAS results tsv were sorted in a different order than the rest of the round 2 (aka imputed-v3) gwas (alphabetically on variant instead of in genome order), meaning they don't correctly mesh with the variants.tsv.bgz file for rsids/info/consequence/etc., requiring a join instead of just pasting the columns together as our readme suggests. So the existing files aren't wrong, but they aren't as compatible as intended with the variant annotation file. A total of 186 GWAS results files, all coming from the [addition of biomarker phenotype sumstats](http://www.nealelab.is/blog/2019/9/16/biomarkers-gwas-results) were affected 
   * SOLUTION: All 186 files have been re-ordered to match properly with the variants.tsv.bgz file. The updated files have "varorder" added to their filename. The earlier files are now removed from the manifest and public download, however we will keep a copy in-house for archive purposes if needed. Below is the list of phenotype codes and description that were adjusted, with 6 subsets for each phenotype (IRNT/RAW,male/female/both_sexes). 
```
30600 - Albumin (g/L)
30610 - Alkaline phosphatase (U/L)
30620 - Alanine aminotransferase (U/L)
30630 - Apoliprotein A (g/L)
30640 - Apoliprotein B (g/L)
30650 - Aspartate aminotransferase (U/L)
30660 - Direct bilirubin (umol/L)
30670 - Urea (mmol/L)
30680 - Calcium (mmol/L)
30690 - Cholesterol (mmol/L)
30700 - Creatinine (umol/L)
30710 - C-reactive protein (mg/L)
30720 - Cystatin C (mg/L)
30730 - Gamma glutamyltransferase (U/L)
30740 - Glucose (mmol/L)
30750 - Glycated haemoglobin (mmol/mol)
30760 - HDL cholesterol (mmol/L)
30770 - IGF-1 (nmol/L)
30780 - LDL direct (mmol/L)
30790 - Lipoprotein A (nmol/L)
30800 - Oestradiol (pmol/L)
30810 - Phosphate (mmol/L)
30820 - Rheumatoid factor (IU/ml)
30830 - SHBG (nmol/L)
30840 - Total bilirubin (umol/L)
30850 - Testosterone (nmol/L)
30860 - Total protein (g/L)
30870 - Triglycerides (mmol/L)
30880 - Urate (umol/L)
30890 - Vitamin D (nmol/L)
30897 - Estimated sample dilution factor (factor)
```
   * DETAILS: 

These phenotypes were added to the list for round 2 late in the process, and were GWASed using Hail 0.2 instead of the Hail 0.1 pipeline used for the rest of the round 2 gwas. The code change from this switch Hail versions is responsible for the different sort order of the output files. Specifically, the tsv export scripts key/sort on "variant" in both pipelines, but the [Hail 0.2 version](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.2/export_results.biomarkers.py#L60) sorts on variant as a constructed string (sorting alphabetically) while the [Hail 0.1 version](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.1/23.export_results.py#L57) used the variant type which sorts by genomic location. (I'm pretty sure at the time Hail 0.2 hadn't implemented it's analogous locus type associated with genome build, though I haven't gone back to confirm.) The [script creating variants.tsv.bgz](https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.1/24.create_variant_annotation_file.py) was part of the primary Hail 0.1 pipeline and so matches that location-based sort used by all of the GWAS outside of the biomarkers.

As far as other releases based on the biomarkers, I haven't found any related issues. Specifically...

Results files for the alternative version of the biomarker GWAS that added dilution fraction as a covariate would probably also be affected by this issue, but those weren't part of the public release (they were just discussed in [a blog post](http://www.nealelab.is/blog/2019/9/16/biomarkers-gwas-results)) and as far as I know aren't getting any other use.

The ldsc-formatted results files ([linked from the h2 results site](https://nealelab.github.io/UKBB_ldsc/downloads.html)) are not affected by this issue. They were exported directly from the Hail GWAS results, and so did not involve re-matching with the variants.tsv.bgz in a way that could introduce issues from sort order. The ldsc format isn't sensitive to sort order, and the export scripts end up sorted on rsid the same for both the biomarkers and everything else anyway.

Round 3 / pan-UKB analyses won't be affected by this issue since they're a separate pipeline unrelated to variants.tsv.bgz and without the separate handling of biomarkers.

As far as I can find, the twitter bot hasn't tweeted any of the round 2 biomarker results. The [manifest of manhattan plots](https://docs.google.com/spreadsheets/d/1T_v3ECnmbxM9ejyDOBn4PA1n_2_fyXySUJn-RpNkpbs/edit#gid=0) that had been created from the UKB results that originally fed the bot doesn't have any of the biomarker phenotypes, so I'm guessing they were never in the bot's rotation.

The results being displayed in CTG-View appear unaffected. For example the [results for albumin](https://view.genoma.io/gwas/5e937684f2e746647caa8962) have the expected variant locations as top hits and they are correctly matched to rsids (i.e. matching both the correctly aligned variants.tsv.bgz and canonical dbSNP entry for the locus).

The [IEU Open GWAS Project](http://gwas.mrcieu.ac.uk/) similarly seems to have correct rsids etc in their phewas lookups, dataset summary, and reformatted vcf.

Harder to confirm for other use the GWAS results have gotten/are getting, and of course should make people aware of the issue, but (fingers crossed) it looks like this hasn't caused widespread problems.	   


 * _January, 2021_
   * We are currently experiencing issues with our DropBox account and working on a fix.  We appreciate your patience during this time
   * UPDATE: GWAS summary statistics are no longer available on Dropbox, and [now available for download on Amazon Web Services](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?usp=sharing)   


 * _Oct 17th, 2019_
   * 89 summary stat files affected by mis-applied low confidence filter have been updated and uploaded to the public release (File Manifest Release 20180731)  

 * _Oct 9, 2019_
   * Summary statistics identified where low confidence filter was mis-applied 
   * [Issue details here](https://github.com/Nealelab/UK_Biobank_GWAS/issues/20)
   * List of files affected (111 files): (GWAS_list_low_confidence_filter_update.txt.gz)[https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/GWAS_list_low_confidence_filter_update.txt.gz]
   * Of these 111 files, 89 require updating, as 22 files are unchanged with the application of updated filter 
   * File column description:
     * phenotype = phenotype number
     * description = UK Biobank description of phenotype
     * min_category = smallest category defined by PHESANT
     * max_category = largest category defined by PHESANT
     * category_distribution = sample counts across categories split by '|'
     * additive_tsvs_list_name = GWAS summary statistic filename
     * n_missing = number of samples without phenotype information
     * tsv_requires_update = TRUE/FALSE does the file require updating of low confidence filter? (phenotypes where min_category < 12500 requires updating)
   * R script to update summary statistic files (Rapid_GWAS_low_confidence_filter_update.R)[https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/Rapid_GWAS_low_confidence_filter_update.R]
     * Requires data.table 1.12.2 R package
     * Requires GWAS_list_low_confidence_filter_update.txt.gz file (or subsetted version with only files you want to update) 

 * _Sept 16, 2019_
   * GWAS summary statistics of Biomarkers now available
   * 34 biomarker meaurements tested
   * [Blog details here](http://www.nealelab.is/blog/2019/9/16/biomarkers-gwas-results)



### imputed-v3 Phenotypes

  * Auto-curated phenotypes using PHESANT: 
    * Source repository: https://github.com/MRCIEU/PHESANT
    * Customized PHESANT repository: https://github.com/astheeggeggs/PHESANT
  * ICD10 codes (all non-coded individuals treated as controls)
  * Curated phenotypes in collaboration with the [FinnGen consortium](https://www.finngen.fi/)

  * Phenotypes in both sexes
    * PHESANT: 2891 total (274 continuous / 271 ordinal / 2346 binary)
    * ICD10: 633 binary
    * FinnGen curated: 559

  * Phenotypes in females
    * PHESANT: 2393 total (259 continuous / 257 ordinal / 1877 binary)
    * ICD10: 482 binary
    * FinnGen curated: 412

  * Phenotypes in males   
    * PHESANT: 2305 total (262 continuous / 259 ordinal / 1784 binary)
    * ICD10: 439 binary
    * FinnGen curated: 400

  * Unique PHESANT phenotypes: 3011, of which 274 are continuous
  * 4203 total unique phenotypes: 3011 PHESANT + 559 finngen + 633 ICD10

  * Summary files: 
	* **phenotypes.both_sexes.tsv.gz** 
	* **phenotypes.female.tsv.gz** 
	* **phenotypes.male.tsv.gz**
    * phenotype - phenotype ID
    * description - short description of phenotype
    * source - PHESANT auto-curation, ICD10, or FinnGen
    * n_controls - number of QC positive samples responding negatively to phenotype designation (NA if quantitative)
    * n_cases - number of QC positive samples responding affirmatively to phenotype designation (NA if quantitative)
    * n_missing - number of missing QC positive samples
    * n_non_missing - number of non-missing QC positive samples

### imputed-v3 Sample QC

  * __imputed-v3 parameters__
	* Used.in.pca.calculation filter (unrelated samples)
	* sex chromosome aneuploidy filter
	* Use provided PCs for European sample selection to determine British ancestry
	  * Use 7 standard deviations away from the 1st 6 PCs
	  * Further Filter to self-reported 'white-British' / 'Irish' / 'White'
	* **QCed sample count: 361194 samples** 
  * __imputed-v2 parameters__
    * Used.in.pca.calculation filter (unrelated samples)
    * sex chromosome aneuploidy filter
    * White.british.ancestry filter
    * **QCed sample count: 337199 samples** 

### imputed-v3 Variant QC

  * __imputed-v3 parameters__
    * Autosomes and X chromosome (including pseudo-autosomal region or XY)
    * SNPs from HRC, UK10K, and 1KG imputation (~90 million)
    * INFO score > 0.8
    * MAF > 0.001
	  * Exception: VEP annotated coding (PTV/Missense/Synonymous) MAF > 1e-6
    * HWE p-value > 1e-10
          * Exception: VEP annotated coding w/MAF < 0.001   
	* **QCed SNP count: 13.7 million** 
  * __imputed-v2 parameters__
    * Autosomes only
    * SNPs from HRC imputation (~40 million)
    * INFO score > 0.8
    * MAF > 0.001
    * HWE p-value > 1e-10
	* **QCed SNP count: 10.9 million** 

### imputed-v3 Association model

  * __imputed-v3 model__
    * Linear regression model in Hail (linreg)
    * Three GWAS per phenotype
	    * Both sexes
	    * Female only
	    * Male only
    * Covariates: 1st 20 PCs + sex + age + age^2 + sex\*age + sex\*age2
    * Sex-specific covariates: 1st 20 PCs + age + age^2 
    * Extra column for variant confidence in case/control phenotypes
      * column name: expected_case_minor_AC
      * Used to filter out false-positive SNPs when case count is low
      * [Blog details here](http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas)  
  * __imputed-v2 model__
    * Linear regression model in Hail (linreg)
    * Covariates: 1st 10 PCs + sex

