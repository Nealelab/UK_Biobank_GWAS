# Table of Contents
* [Updates](#updates) 
* [imputed-v3 Phenotypes](#imputed-v3-phenotypes)
* [imputed-v3 Sample QC](#imputed-v3-sample-qc)
* [imputed-v3 Variant QC](#imputed-v3-variant-qc)
* [imputed-v3 Association model](#imputed-v3-association-model)

### Updates

With the re-release of UK Biobank genotype imputation (which we term imputed-v3), we have generated an updated set of GWAS summary statistics for the genetics community. We have increased the number of phenotypes (our current application: UKB31063), inclusion of samples (see Sample QC changes), inclusion of SNPs (see Variant QC changes), and our association model (see Association model changes). Our largest change is the for all phenotypes, we have run a sex-specific GWAS analysis along with the full set.

Information and scripts from the previous  round of GWAS are available in the [imputed-v2-gwas](https://github.com/Nealelab/UK_Biobank_GWAS/tree/master/imputed-v2-gwas) subdirectory

### imputed-v3 Phenotypes

  * 2455 auto-curated phenotypes using PHESANT: 
    * Source repository: https://github.com/MRCIEU/PHESANT
    * Customized PHESANT repository: https://github.com/astheeggeggs/PHESANT
  * 806 ICD10 codes
  * 694 curated phenotypes in collaboration with the [FinnGen consortium](https://www.finngen.fi/)

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

  * imputed-v3 parameters
	* Used.in.pca.calculation filter (unrelated samples)
	* sex chromosome aneuploidy filter
	* Use provided PCs for European sample selection to determine British ancestry
	  * Use 7 standard deviations away from the 1st 6 PCs
	  * Further Filter to self-reported 'white-British' / 'Irish' / 'White'
	* **QCed sample count: 361194 samples** 
  * imputed-v2 parameters
    * Used.in.pca.calculation filter (unrelated samples)
    * sex chromosome aneuploidy filter
    * White.british.ancestry filter
    * **QCed sample count: 337199 samples** 

### imputed-v3 Variant QC

  * imputed-v3 parameters
    * Autosomes and X chromosome (but not pseudo-autosomal region or XY)
    * SNPs from HRC, UK10K, and 1KG imputation (~90 million)
    * INFO score > 0.8
    * MAF > 0.0001
	  * Exception: VEP annotated Missense and PTV MAF > 1e-6
    * HWE p-value > 1e-10  
	* **QCed SNP count: 13.7 million** 
  * imputed-v2 parameters
    * Autosomes only
    * SNPs from HRC imputation (~40 million)
    * INFO score > 0.8
    * MAF > 0.0001
	* **QCed SNP count: 10.9 million** 

### imputed-v3 Association model

  * imputed-v3 model
    * Linear regression model in Hail (linreg)
    * Three GWAS per phenotype
	    * Both sexes
	    * Female only
	    * Male only
    * Covariates: 1st 20 PCs + sex + age + age^2 + sex*age + sex*age2
    * Sex-specific covariates: 1st 20 PCs + age + age^2 
    * Extra column for variant confidence in case/control phenotypes
      * column name: expected_case_minor_AC
      * Used to filter out false-positive SNPs when case count is low
      * [Blog details here](http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas)  
  * imputed-v2 model
    * Linear regression model in Hail (linreg)
    * Covariates: 1st 10 PCs + sex

