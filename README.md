# GWAS-Parsl
GWAS running with Parsl for running anywhere

## Dependencies
The following is a list of dependencies which will need to be installed in your conda environment. We will be generating a Singularity container that users can take to run on any platform.
- [Plink v1.9](https://www.cog-genomics.org/plink/1.9/)
- [R version 3.6.1](https://www.r-project.org)
    - [QQman v0.1.4](https://cran.r-project.org/web/packages/qqman/index.html)
- Python 3.6.7 installed via anaconda3. The following packages need to be installed. These packages have dependencies that also need to be installed.
    - [Rpy2 v2.9](https://rpy2.readthedocs.io/en/latest/)
    - [Parsl v0.9](https://parsl.readthedocs.io/en/stable/)
    - [LDPred v1.0.8](https://github.com/bvilhjal/ldpred)
    - [ldsc v1.0.0](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)
    - [PRSice-2](https://github.com/choishingwan/PRSice)

## Reference Data
The following reference data needs to be included in the analysis to perform stratification quality control.

Instructions to be added later ...
