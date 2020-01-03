# GWAS-Parsl
GWAS running with Parsl for running anywhere

## Introduction
There are some limiations inherited in classical GWAS:
- Needs large amounts of data for generating the summary data
- Data must be homogeneous (i.e. similar background, European, African, Asian, Hispanic, etc.)
- Data must be of high quality, so one may lose significant variants in the process

Ideally, we should use a Deep Learning (DL) based algorithm to use all data types to learn and improve a model.
For the purposes of completeness, we will first generate a GWAS algorithm that runs the classical GWAS method. Later we will generate a DL model and compare with our GWAS algorithm.
Due to the nature of the data, we will use the python library [Parsl](https://parsl.readthedocs.io/en/stable/), to be able to run these steps on any platform and fully parallelizable. In addition, this will allow us to follow the [FAIR steps](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0213013) by making the process reproducible.
The steps involved in performing GWAS analysis are the following and are performed in multiple scripts below:
- GWAS quality control
- Association analysis
- PRS on distinct data of similar background

Each step will be further explained below.

## Dependencies
The following is a list of dependencies which will need to be installed in your conda environment. We will be generating a Singularity container that users can take to run on any platform.
- [Plink v1.9](https://www.cog-genomics.org/plink/1.9/)
- [R version 3.6.1](https://www.r-project.org)
    - [QQman v0.1.4](https://cran.r-project.org/web/packages/qqman/index.html)
- Python 3.6.7 installed via anaconda3. The following packages need to be installed. These packages have dependencies that also need to be installed.
    - [Rpy2 v2.9](https://rpy2.readthedocs.io/en/latest/)
    - [Parsl v0.9](https://parsl.readthedocs.io/en/stable/)
    - [ldsc v1.0.0](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) - Needed for post-GWAS validation.
    - [PRSice-2](https://github.com/choishingwan/PRSice) - Needed for PRS
    - [LDPred v1.0.8](https://github.com/bvilhjal/ldpred) - Needed for PRS

## Reference Data
The following reference data needs to be included in the analysis to perform stratification quality control.

Instructions to be added later ...

## GWAS Quality Control
The first step in GWAS analysis is to perform a quality control (QC) of the initial data. The QC of the data will limit the amount of low quality variants used in the learning process. Low quality variants increase the chances of having false positives in the downstream polygenetic risk score (PRS) step. Thus, there are several steps in removing low quality variants.
In this script generated the four types of low quality variants that will be filtered out are:

### Missingness
Investigate missingness per individual and per SNP and generate histograms output: plink.imiss and plink.lmiss, these files show respectively the proportion of missing SNPs per individual and the proportion of missing individuals per SNP.

### Sex Discrepancy
Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate. Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.

### Minor allele Frequency
Generate a bfile with autosomal SNPs only and delte SNPs with a low minor allele Frequency (MAF). Select autosomal SNPs only from chr 1 to 22 and then remove the variants with low MAF.

### Hardy-Weinberg equilibrium
Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE). By default the --hwe option in plink only filters for controls. Therefore, two steps, first we use a stringent HWE threshold for controls, followed by a less stringent threshold for the case data. The HWE threshold for the cases filters out only SNPs which deviate extremely from HWE. This second HWE step only focusses on cases because in the controls all SNPs with a HWE p-value < hwe 1e-6 were already removed. Theoretical background for this step is given in our accompanying [article](https://www.ncbi.nlm.nih.gov/pubmed/29484742).

### Heterozygozity control
Remove individuals with a heterozygosity rate deviating more than 3 standard deviations from the mean.
Checks for heterozygosity are performed on a set of SNPs which are not highly correlated. Therefore, to generate a list of non-(highly)correlated SNPs, we exclude high inversion regions (inversion.txt [High LD regions]) and prune the SNPs using the command --indep-pairwise<92>. The parameters <91>50 5 0.2<92> stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.

### Relatedness control
Remove relatedness individuals. It is essential to check datasets you analyse for cryptic relatedness. Assuming a random population sample we are going to exclude all individuals above the pihat threshold of 0.2

The python script [```gwas-qc-parsl.py```]() runs all these quality control checks. To run this, activate the environemnt and run:

```#python gwas-qc-parsl.py --input-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/inputs/HapMap_3_r3_1.bed --output-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/outputs/  --inversion-regions /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/inversion.txt  --step-start relatedness_qc```

## Stratification and Association Analysis


## To read
[A flexible and parallelizable approach to genome‐wide polygenic risk scores ](https://onlinelibrary.wiley.com/doi/epdf/10.1002/gepi.22245)
