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
    
### Building the Singularity image
We built a singularity image containing the necessary tools to run the GWAS application. The recipe can be accessed [here](https://github.com/arodri7/GWAS-Parsl/blob/master/singularity-gwas.recipe).
Steps to build the image are the following:
- Install Vagrant
- Install Singularity
- Initiate a VM using vagrant, bring up the VM and ssh to VM
    ````
    mkdir singularity-VM
    cd singularity-VM
    vagrant init singularityware/singularity-2.4
    vagrant up
    vagrant ssh
    ````
- Download the recipe and build image
    ````
    wget https://github.com/arodri7/GWAS-Parsl/blob/master/singularity-gwas.recipe
    sudo singularity build singularity-gwas.simg singularity-gwas.recipe
    ````

You should now see an image file called ```singularity-gwas.simg```, which can be used to run your commands in the script.

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

The python script [```gwas-qc-parsl.py```](https://github.com/arodri7/GWAS-Parsl/blob/master/gwas-qc-parsl.py) runs all these quality control checks. To run this, activate the environemnt and run:

```
python gwas-qc-parsl.py \
   --input-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/inputs/HapMap_3_r3_1.bed \ 
   --output-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/outputs/  \
   --inversion-regions /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/inversion.txt \
   --step-start relatedness_qc
   ```

## Stratification
Before performing the association analysis we need to make sure the data are all from the same background. This is one of the limitations on GWAS. We perform stratification quality control to separate the data using the reference data material from 1000 Genome material.
The following are in general the steps for stratification:
- Convert from vcf to plink format
- Take plink input and assign unique indentifiers to the SNPs with a missing rs-identifier
- Remove SNPs with a low MAF frequency
- Run the MDS removal part to generate plots for stratification
    - Extract the variants present in HapMap dataset from the 1000 genomes dataset.
    - Extract the variants present in 1000 Genomes dataset from the HapMap dataset.
    - The datasets must have the same build. Change the build 1000 Genomes data build.
    - Prior to merging 1000 Genomes data with the HapMap data we want to make sure that the files are mergeable:
        1) Make sure the reference genome is similar in the HapMap and the 1000 Genomes Project datasets.
        2) Resolve strand issues.
        3) Remove the SNPs which after the previous two steps still differ between datasets.

    - Merge HapMap with 1000 Genomes Data.
    - Perform MDS on HapMap-CEU data anchored by 1000 Genomes data.
    - Using a set of pruned SNPs
    - Convert population codes into superpopulation codes (i.e., AFR,AMR,ASN, and EUR).
    - Create a racefile of your own data.
    - Concatenate racefiles.
    - Exclude ethnic outliers.
    - Select individuals in HapMap data below cut-off thresholds -The cut-off levels are not fixed thresholds but have to be determined based on the  visualization of the first two dimensions. To exclude ethnic outliers, the thresholds need to be set around the cluster of population of interest. For now we will hardcode to exclude non EURopeian. 
    - Create covariates based on MDS.
    - Perform an MDS ONLY on HapMap data without ethnic outliers - The values of the 10 MDS dimensions are subsequently used as covariates in the association analysis in the third tutorial.
    - Plot stratification of data
    
The python script [```reference-qc-gwas-parsl.py```](https://github.com/arodri7/GWAS-Parsl/blob/master/reference-qc-gwas-parsl.py) runs the stratification step. To run this, activate the environemnt and run:

```
python reference-qc-gwas-parsl.py \
    --input-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/outputs/plink_remove_het.bed \
    --reference /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/references/ALL.2of4intersection.20100804.genotypes.vcf.gz \
    --population /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/references/20100804.ALL.panel \
    --output-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/2_Population_stratification/outputs/ \
    --prune /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/outputs/indepSNP.prune.in \
    --step-start plot_qc
```

## Association Analysis
The next step ingests the filtered data and the multiple stratified groups and performs association analysis on each group which will generate multiple summary association files which can then be used in the PRS downstream process. The steps are simple and performed using plink:
- Binary traits association - Note, the --assoc option does not allow to correct covariates such as principal components (PC's)/ MDS components, which makes it less suited for association analyses.
- Logistic - We will be using 10 principal components as covariates in this logistic analysis.
- The results obtained from these GWAS analyses will be visualized in the last step. This will also show if the data set contains any genome-wide significant SNPs.

The python script [```association-parsl.py```](https://github.com/arodri7/GWAS-Parsl/blob/master/association-parsl.py) runs the association step. To run this, activate the environemnt and run:

```
python association-parsl.py \
    --input-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/2_Population_stratification/outputs/HapMap_3_r3_13 \
    --covariates /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/2_Population_stratification/outputs/ALL.2of4intersection.20100804.genotypes.1kG_MDS.covariants.txt \
    --output-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/2_Population_stratification/outputs/
```

## Polygenic Risk Score (PRS)
PRS implementation of GWAS data. We will implement multiple methods so the user can decide which to use
Initial methods will include:
- LDPred
- PRSice-2
We can add more methods later on to compare which performs best.

In order to run PRS for GWAS data, we need the base and target data. The base data is essentially the association summary file from the GWAS study. The target data are the individuals which we wish to generate a polygenic risk score for. Usually, this should include more than 100 individuals. No matter which PRS method we use, the base and target data need to go through QC.
The QC steps performed for both the base and target data are similar to the QC steps followed for the GWAS analysis. Thus, we can call the same functions we have used previously. Most importantly, we need to make sure that samples in the target data are not included in the base data, otherwise the PRS will be skewed and there will be too many false positives.

#### Outline of QC steps for Base and Target data:
- Heritability check - For now we will assume that is higher than 0.05. This is something that needs to be done for each association file.
- Effect allele - This is clearly marked on the association file (col A1)
- Genome build - This check will need to be made for the target qc step
- Standard GWAS QC - Since we will only have an association file, we will perform QC only on the file.
                   - For the target QC, we should have plink files and we can use regular plink qc steps
- Ambiguous SNPs - Ambiguous SNPs can be removed in the base data and then there will be no such SNPs in the subsequent analyses, since analyses are performed only on SNPs that overlap between the base and target data.
- Mismatching genotypes - since we need the target data to know which SNPs have mismatching genotypes across the data sets, then we will perform this 'allele flipping' in the target data.
- Duplicate SNPs - Most PRS software do not allow duplicated SNPs in the base data input and thus they should be removed
- Sex chromosomes - Previously performed QC on these data removed individuals with mismatching (inferred) biological and reported sex. Already performed on base data.
- Sample overlap with target data - users should ensure that the possibility of sample overlap between the base and target data is minimised
- Relatedness - users should ensure that the possibility of closely related individuals between the base and target data is minimised. This is part of the QC performed in generating the association file.
    
The QC steps will be explained throughout this script for both the base and target data.

### Base QC
Performs qc on the base data. There are some assumptions that come from generating the association file and those steps are skipped since these qc steps are performed before.
Biggest assumption is that the GWAS association file was generated locally.
For external GWAS association files, we will need to perform more checks.
This only covers base files created locally!

Checks to be skipped and performed:
- Heritability check - DONE - this is performed after the GWAS association file has been generated. This should already be known and it should be higher than 0.05
- Effect allele - DONE in post-association step.
- Genome build - will be done at the target level
- Standard GWAS qc -
    - Filter the SNPs according to INFO score and MAF since we only have a summary base file. SNPs with low minor allele frequency (MAF) or imputation information score (INFO) are more likely to generate false positive results due to their lower statistical power. Therefore, SNPs with low MAF and INFO are typically removed before performing downstream analyses.
    - Ambiguous SNPs QC - Ambiguous SNPs can be removed in the base data and then there will be no such SNPs in the subsequent analyses, since analyses are performed only on SNPs that overlap between the base and target data
    - Duplicate SNPs QC - Most PRS software do not allow duplicated SNPs in the base data input and thus they should be removed.
- Ambiguous SNPs - YES
- Mismatching genotypes - will be done at the target level
- Duplicate SNPs - YES
- Sex chromosome - DONE
- Sample overlap with target data - will be done at the target level
- Relatedness - DONE

### Target QC
This consists of individual-level genotype-phenotype data, usually generated within your lab/department/collaboration
- Genome build - YES
- Standard GWAS qc - YES
- Mismatching genotypes - YES
- Duplicate SNPs - YES
- Sex chromosome - YES
- Sample overlap with target data - YES
- Relatedness - YES

## To read
[A flexible and parallelizable approach to genome‐wide polygenic risk scores ](https://onlinelibrary.wiley.com/doi/epdf/10.1002/gepi.22245)
