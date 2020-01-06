import parsl
import os, argparse, sys
from parsl.configs.local_threads import config
import parsl_funcs as pf

#parsl.set_stream_logger() # <-- log everything to stdout
#print(parsl.__version__)
parsl.load(config)

# PRS implementation of GWAS data
# We will implement multiple methods so the user can decide which to use
# Initial methods will include:
# - LDPred
# - PRSice-2
# We can add more methods later on to compare which performs best.

# In order to run PRS for GWAS data, we need the base and target data.
# The base data is essentially the association summary file from the GWAS study.
# The target data are the individuals which we wish to generate a polygenic risk score for.
# Usually, this should include more than 100 individuals.
# No matter which PRS method we use, the base and target data need to go through QC.
# The QC steps performed for both the base and target data are similar to the QC steps
# followed for the GWAS analysis. Thus, we can call the same functions we have used previously.
# Most importantly, we need to make sure that samples in the target data are not included in the 
# base data, otherwise the PRS will be skewed and there will be too many false positives.
# The QC steps will be explained throughout this script for both the base and target data.

#steps = ["missingness_qc", "sex_discrepancy_qc", "maf_qc", "hwe_qc", "het_qc", "relatedness_qc"]
steps = []
prs_methods = ['LDPred', 'PRSice2']


# base_qc: performs qc on the base data. There are some assumptions that come from generating the 
#          association file and those steps are skipped since these qc steps are performed before.
# Biggest assumption is that the GWAS association file was generated locally.
# For external GWAS association files, we will need to perform more checks.
# This only covers base files created locally!
def base_qc(inputs=[], outputs=[]):
    # Checks to be skipped and performed: 
    #       Heritability check - DONE - this is performed after the GWAS association file has been generated.
    #                            This should already be known and it should be higher than 0.05
    #       Effect allele - DONE beforehand
    #       Genome build - will be done at the target level
    #       Standard GWAS qc - YES
    #       Ambiguous SNPs - YES
    #       Mismatching genotypes - will be done at the target level
    #       Duplicate SNPs - YES
    #       Sex chromosome - DONE 
    #       Sample overlap with target data - will be done at the target level
    #       Relatedness - DONE

    b_outs = []

    # 1. Standrad GWAS QC
    #    Filter the SNPs according to INFO score and MAF since we only have a summary base file.
    #    SNPs with low minor allele frequency (MAF) or imputation information score (INFO) 
    #    are more likely to generate false positive results due to their lower statistical power 
    #    Therefore, SNPs with low MAF and INFO are typically removed before performing downstream analyses.
    bgqo = [args.output_dir + "tmp1.out"]
    base_gwas_qc_out = pf.base_gwas_qc(inputs=[args.input_base], outputs=bgqo)
    base_gwas_qc_out.result()
    b_outs.extend(base_gwas_qc_out.outputs)
    
    # 2. Ambiguous SNPs QC
    #    Ambiguous SNPs can be removed in the base data and then there will be no such SNPs in the 
    #    subsequent analyses, since analyses are performed only on SNPs that overlap between the base and target data
    asqo = [args.output_dir + "tmp2.out"]
    base_amb_snps_qc_out = pf.base_amb_snps_qc(inputs=bgqo, outputs=asqo)
    base_amb_snps_qc_out.result()
    b_outs.extend(base_amb_snps_qc_out.outputs)

    # 3. Duplicate SNPs QC
    #    Most PRS software do not allow duplicated SNPs in the base data input and thus they should be removed
    dsqo = [args.output_dir + "base.out"] 
    base_dup_snps_qc_out = pf.base_dup_snps_qc(inputs=asqo, outputs=dsqo)
    base_dup_snps_qc_out.result()
    b_outs.extend(base_dup_snps_qc_out.outputs)

    return(dsqo)

# Target qc: This consists of individual-level genotype-phenotype data, 
#            usually generated within your lab/department/collaboration
def target_qc(inputs=[], outputs=[]):
    # Checks to be skipped and performed:
    #       Genome build - YES
    #       Standard GWAS qc - YES
    #       Mismatching genotypes - YES
    #       Duplicate SNPs - YES
    #       Sex chromosome - YES
    #       Sample overlap with target data - YES
    #       Relatedness - YES
    
    # 1. Genome build
    gbqo = [args.output_dir + "tmp1.bed"]
    target_genome_qc_out = pf.target_build_qc(inputs=[inputs[0]], outputs=gbqo)
    target_genome_qc_out.result()
    b_outs.extend(target_genome_qc_out.outputs)

    # 2. Standard GWAS qc
    sgqo = [args.output_dir + "tmp2.bed"]
    target_std_qc_out = pf.target_std_qc(inputs=gbqo, outputs=sgqo)
    target_std_qc_out.result()
    b_outs.extend(target_std_qc_out.outputs)

    # 3. Mismatching genotypes qc
    mgqo = [args.output_dir + "tmp3.bed"]
    target_mismatch_qc_out = pf.target_mismatch_qc(inputs=sgqo, outputs=mgqo)
    target_mismatch_qc_out.result()
    b_outs.extend(target_mismatch_qc_out.outputs)

    # 4. Duplicate snps qc
    dsqo = [args.output_dir + "tmp4.bed"]
    target_dup_qc_out = pf.target_dup_qc(inputs=mgqo, outputs=dsqo)
    target_dup_qc_out.result()
    b_outs.extend(target_dup_qc_out.outputs)

    # 5. Sex chromosome qc
    scqo = [args.output_dir + "tmp5.bed"]
    target_sex_qc_out = pf.target_sex_qc(inputs=dsqo, outputs=scqo)
    target_sex_qc_out.result()
    b_outs.extend(target_sex_qc_out.outputs)

    # 6. Sample overlap qc
    soqo = [args.output_dir + "tmp6.bed"]
    target_overlap_qc_out = pf.target_overlap_qc(inputs=scqo, outputs=soqo)
    target_overlap_qc_out.result()
    b_outs.extend(target_overlap_qc_out.outputs)

    # 7. sample Relatedness qc
    srqo = [outputs[0]]
    target_related_qc_out = pf.target_related_qc(inputs=soqo, outputs=srqo)
    target_related_qc_out.result()
    b_outs.extend(target_related_qc_out.outputs)

    return(srqo)

###########################
# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Process GWAS pipeline for different steps.')
    parser.add_argument('--input-base', dest='input_base', 
                        help='Location of the input GWAS association summary file')
    parser.add_argument('--target-path', dest='target_path',
                        help='Location for the target data files. These are typically plink files (bed, fam, bim)')
    parser.add_argument('--output-directory', dest='output_dir',
                        help='Location of the outputs GWAS data')
    parser.add_argument('--step-start', dest='step_start',
                        choices=steps,
                        help='Step to start in. The assumption is that the input must be ready for the step you are starting from.')
    parser.add_argument('--method', dest='prs_method', choices=prs_methods,
                        help='select the PRS method to use')

    return(parser.parse_args())

# Main function controls the inputs and steps to submit
def main():
    # get input arguments
    args = parse_args()
    print(args)
    flag_start = True
    # workflow outputs variable
    wo = []

    #####################################################################
    #####################################################################
    # Start the Base dataset QC
    # outline:
    # - Heritability check - For now we will assume that is higher than 0.05. 
    #                        This is something that needs to be done for each association file.
    # - Effect allele - This is clearly marked on the association file (col A1)
    # - Genome build - This check will need to be made for the target qc step
    # - Standard GWAS QC - Since we will only have an association file, we will perform QC only on the file.
    #                    - For the target QC, we should have plink files and we can use regular plink qc steps
    # - Ambiguous SNPs - Ambiguous SNPs can be removed in the base data and then there will be no such SNPs 
    #                    in the subsequent analyses, since analyses are performed only on SNPs that 
    #                    overlap between the base and target data.
    # - Mismatching genotypes - since we need the target data to know which SNPs have mismatching genotypes 
    #                           across the data sets, then we will perform this 'allele flipping' in the target data.
    # - Duplicate SNPs - Most PRS software do not allow duplicated SNPs in the base data 
    #                    input and thus they should be removed
    # - Sex chromosomes - Previously performed QC on these data removed individuals with mismatching 
    #                     (inferred) biological and reported sex. Already performed on base data.
    # - Sample overlap with target data - users should ensure that the possibility of sample overlap 
    #                                     between the base and target data is minimised
    # - Relatedness - users should ensure that the possibility of closely related individuals between 
    #                 the base and target data is minimised. This is part of the QC performed in generating the
    #                 association file

    ## STEP 1: perform QC on base data
    #if (args.step_start == "missingness_qc" and flag_start == False) or args.step_start is None:
    #    flag_start = True

    #print("FLAG STEP1: %s" % flag_start)    
    output_s1 = [args.output_dir + "assoc_results.assoc",
                 args.output_dir + "assoc_results.log", 
                 args.output_dir + "logistic_results.assoc.logistic",
                 args.output_dir + "logistic_results.log",
                 args.output_dir + "logistic_results.assoc_2.logistic"]

    if flag_start is True:
        # Run the base qc
        s1 = base_qc(inputs=[args.input_base, args.target_path], outputs=output_s1)
        s1.result()
        wo.extend(s1.outputs)

    # STEP 2: perform QC on target data. Steps are similar to base qc
    output_s2 = [args.output_dir + "assoc_results.assoc",
                 args.output_dir + "assoc_results.log",
                 args.output_dir + "logistic_results.assoc.logistic",
                 args.output_dir + "logistic_results.log",
                 args.output_dir + "logistic_results.assoc_2.logistic"]
    if flag_start is True:
        # Run the target qc.
        s2 = target_qc(inputs=[output_s1[0], output_s1[4]], outputs=output_s2)
        s2.result()
        wo.extend(s2.outputs)

    # STEP 3: Perform the PRS using the desired method
    output_s2 = [args.output_dir + "assoc_results.assoc",
                 args.output_dir + "assoc_results.log"
                ]
    if flag_start is True:
        # Run the PRS for method selected
        s3 = run_prs(inputs=[], outputs=[], method=args.prs_method)
        s3.result()
        wo.extend(s3.outputs)


main()
