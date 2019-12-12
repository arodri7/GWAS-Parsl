import parsl
import os, argparse
from parsl.configs.local_threads import config
import sys
import parsl_funcs as pf

#parsl.set_stream_logger() # <-- log everything to stdout

#python gwas-qc-parsl.py --input-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/inputs/HapMap_3_r3_1.bed --output-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/outputs/  --inversion-regions /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/inversion.txt  --step-start relatedness_qc

#print(parsl.__version__)
parsl.load(config)

steps = ["missingness_qc", "sex_discrepancy_qc", "maf_qc", "hwe_qc", "het_qc", "relatedness_qc"]

###########################
# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Process GWAS pipeline for different steps.')
    parser.add_argument('--input-directory', dest='input_dir', 
                        help='Location of the input GWAS data')
    parser.add_argument('--inversion-regions', dest='inversion_regions',
                        help='Location of the input GWAS data')
    parser.add_argument('--output-directory', dest='output_dir',
                        help='Location of the outputs GWAS data')
    parser.add_argument('--step-start', dest='step_start',
                        choices=steps,
                        help='Step to start in. The assumption is that the input must be ready for the step you are starting from.')

    return(parser.parse_args())

# Main function controls the inputs and steps to submit
def main():
    # get input arguments
    args = parse_args()
    print(args)
    flag_start = False
    # workflow outputs variable
    wo = []

    ## STEP 1: Missingness
    if (args.step_start == "missingness_qc" and flag_start == False) or args.step_start is None:
        flag_start = True

    print("FLAG STEP1: %s" % flag_start)    
    output_mq = [args.output_dir + "plink_missingness.imiss", 
                 args.output_dir + "plink_missingness.lmiss"]

    output_mq_plots = [args.output_dir + "histimiss.pdf",
                       args.output_dir + "histlmiss.pdf"]

    rsnps_outputs = [args.output_dir + "plink_missingness_remove_snps.bed",
                     args.output_dir + "plink_missingness_remove_snps.hh",
                     args.output_dir + "plink_missingness_remove_snps.fam",
                     args.output_dir + "plink_missingness_remove_snps.bim",
                     args.output_dir + "plink_missingness_remove_snps.log",
                    ]

    if flag_start is True:
        # Run missingness QC
        mq = pf.missingness_qc(inputs=[args.input_dir], outputs=[args.output_dir + "plink_missingness.imiss"])
        #if flag_start is True: mq.result()
        mq.result()
        wo.extend(mq.outputs)

        # Generate missingness plots
        mq_plots = pf.missingness_qc_plots(inputs=output_mq, outputs=output_mq_plots)
        #if flag_start is True: mq_plots.result()
        mq_plots.result()
        wo.extend(mq_plots.outputs)

        # Delete SNPs and individuals with high levels of missingness
        rsnps = pf.plink_remove_missingness(inputs=[args.input_dir], outputs=rsnps_outputs)
        #if flag_start is True: rsnps.result()
        rsnps.result()
        wo.extend(rsnps.outputs)

    ###########################################################
    ## STEP 2: Sex discrepancy
    if args.step_start == "sex_discrepancy_qc" and flag_start == False:
        flag_start = True
 
    print("FLAG STEP 2: %s" % flag_start)
    output_sq = [args.output_dir + "plink_sex_discrepancy.sexcheck" ]

    output_sq_plots = [args.output_dir + "Gender_check.pdf",
                       args.output_dir + "Men_check.pdf",
                       args.output_dir + "Women_check.pdf"]

    rsex_outputs = [args.output_dir + "plink_remove_sex_discrepancy.bed",
                    args.output_dir + "plink_remove_sex_discrepancy.hh",
                    args.output_dir + "plink_remove_sex_discrepancy.fam",
                    args.output_dir + "plink_remove_sex_discrepancy.bim",
                    args.output_dir + "plink_remove_sex_discrepancy.log",
                    ]

    if flag_start is True:
        # Run sex discrepancy QC
        sq = pf.sexdiscrepancy_qc(inputs=rsnps_outputs, outputs=output_sq)
        #if flag_start is True: sq.result()
        sq.result()
        wo.extend(sq.outputs)

        # Generate sex discrepancy plots
        sq_plots = pf.sexdiscrepancy_qc_plots(inputs=output_sq, outputs=output_sq_plots)
        #if flag_start is True: sq_plots.result()
        sq_plots.result()
        wo.extend(sq_plots.outputs)

        # Delete Sex discrepancy individuals
        rsex = pf.remove_sex_discrepancy(inputs=[output_sq[0], rsnps_outputs[0]], outputs=rsex_outputs)
        #if flag_start is True: rsex.result()
        rsex.result()
        wo.extend(rsex.outputs)

    ###########################################################
    ## STEP 3: Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF)
    if args.step_start == "maf_qc" and flag_start == False:
        flag_start = True

    print("FLAG STEP 3: %s" % flag_start)
    output_autosomal_snps = [args.output_dir + "plink_autosomal_snps.bed",
                             args.output_dir + "plink_autosomal_snps.fam",
                             args.output_dir + "plink_autosomal_snps.bim",
                             args.output_dir + "plink_autosomal_snps.log" ]
    output_maf_check_qc = [args.output_dir + "plink_maf_freq.frq"]
    output_maf_plots = [args.output_dir + "MAF_distribution.pdf"]
    output_delete_maf = [args.output_dir + "plink_remove_maf.bed",
                         args.output_dir + "plink_remove_maf.fam",
                         args.output_dir + "plink_remove_maf.bim",
                         args.output_dir + "plink_remove_maf.log",
                        ]

    # Run autosomal SNP detection
    autosomal = pf.select_autosomal_snps(inputs=[rsex_outputs[3]], outputs=output_autosomal_snps)
    if flag_start is True: autosomal.result()
    #autosomal.result()
    wo.extend(autosomal.outputs)

    if flag_start is True:
        # Generate MAF distribution
        maf_dist = pf.maf_check_qc(inputs=output_autosomal_snps, outputs=output_maf_check_qc)
        #if flag_start is True: maf_dist.result()
        maf_dist.result()
        wo.extend(maf_dist.outputs)

        # Generate maf distribution plots
        maf_plots = pf.maf_distribution_qc_plots(inputs=output_maf_check_qc, outputs=output_maf_plots)
        #if flag_start is True: maf_plots.result()
        maf_plots.result()
        wo.extend(maf_plots.outputs)

        # Delete SNPs with low MAF frequency
        low_maf_remove = pf.remove_snps_low_maf(inputs=output_autosomal_snps, outputs=output_delete_maf)
        #if flag_start is True: low_maf_remove.result()
        low_maf_remove.result()
        wo.extend(low_maf_remove.outputs)

    ################################################
    # Step 4
    # Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE)
    if args.step_start == "hwe_qc" and flag_start == False:
        flag_start = True

    print("FLAG STEP 4: %s" % flag_start)
    output_hwe_qc_plots = [args.output_dir + "histhwe_below_theshold.pdf"]
    output_hwe_qc = [args.output_dir + "plink_hwe_qc.hwe",
                     args.output_dir + "plinkzoomhwe.hwe" ]

    output_hwe = [args.output_dir + "plink_remove_hwe.bed",
                  args.output_dir + "plink_remove_hwe.fam",
                  args.output_dir + "plink_remove_hwe.bim",
                  args.output_dir + "plink_remove_hwe.log" ]

    if flag_start is True:
        # Check distribution of HWE p-values of all SNPs
        hwe_dist = pf.check_hwe_distribution(inputs=output_delete_maf, outputs=output_hwe_qc)
        #if flag_start is True: hwe_dist.result()
        hwe_dist.result()
        wo.extend(hwe_dist.outputs)
   
        # HWE distribution plots
        hwe_dist_plots = pf.check_hwe_distribution_plots(inputs=output_hwe_qc, outputs=output_hwe_qc_plots)
        #if flag_start is True: hwe_dist_plots.result()
        hwe_dist_plots.result()
        wo.extend(hwe_dist_plots.outputs)

        # By default the --hwe option in plink only filters for controls.
        # Therefore, two steps, first we use a stringent HWE threshold for controls, 
        # followed by a less stringent threshold for the case data.
        # The HWE threshold for the cases filters out only SNPs which deviate extremely from HWE.
        # This second HWE step only focusses on cases because in the controls all SNPs 
        # with a HWE p-value < hwe 1e-6 were already removed
        # Theoretical background for this step is given in our accompanying article: 
        # https://www.ncbi.nlm.nih.gov/pubmed/29484742 . 
        remove_hwe = pf.remove_hwe_snps(inputs=output_delete_maf, outputs=output_hwe)
        #if flag_start is True: remove_hwe.result()
        remove_hwe.result()
        wo.extend(remove_hwe.outputs)

    ################################################
    # Step 5
    # remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.
    # Generate a plot of the distribution of the heterozygosity rate of your subjects.
    # And remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.

    # Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
    # Therefore, to generate a list of non-(highly)correlated SNPs, we exclude high inversion regions 
    # (inversion.txt [High LD regions]) and prune the SNPs using the command --indep-pairwise<92>.
    # The parameters <91>50 5 0.2<92> stand respectively for: the window size, the number of SNPs 
    # to shift the window at each step, and the multiple correlation coefficient for a SNP being 
    # regressed on all other SNPs simultaneously.

    if args.step_start == "het_qc" and flag_start == False:
        flag_start = True

    print("FLAG STEP 5: %s" % flag_start)
    output_ld_qc = [args.output_dir + "indepSNP.prune.in",
                    args.output_dir + "indepSNP.prune.out",
                    args.output_dir + "plink_check_het.het"]
    output_het_qc_plots = [args.output_dir + "heterozygosity.pdf",
                           args.output_dir + "fail-het-qc.txt"]
    output_het = [args.output_dir + "plink_remove_het.bed",
                  args.output_dir + "plink_remove_het.fam",
                  args.output_dir + "plink_remove_het.bim",
                  args.output_dir + "plink_remove_het.log" ]

    if flag_start is True:
        # Remove high LD regions
        remove_ld = pf.run_remove_ld(inputs=[output_hwe[0], args.inversion_regions], outputs=output_ld_qc)
        remove_ld.result()
        wo.extend(remove_ld.outputs)

        # Plot heterozygosity rate distribution and generate list of individuals who 
        # deviate more than 3 std from the heterozygosity rate mean 
        het_qc_plots = pf.run_het_qc_plots(inputs=[output_ld_qc[2]], outputs=output_het_qc_plots)
        het_qc_plots.result()
        wo.extend(het_qc_plots.outputs)

        # remove the failed het individuals
        remove_het = pf.run_remove_het(inputs=[output_hwe[0], output_het_qc_plots[1]], outputs=output_het)
        remove_het.result()
        wo.extend(remove_het.outputs)

    ############################################
    # STEP 6: Remove relatedness individuals
    # It is essential to check datasets you analyse for cryptic relatedness.
    # Assuming a random population sample we are going to exclude all individuals 
    # above the pihat threshold of 0.2 in this tutorial.

    # Check for relationships between individuals with a pihat > 0.2.
    if args.step_start == "relatedness_qc" and flag_start == False:
        flag_start = True

    print("FLAG STEP 6: %s" % flag_start)
    output_relatedness_qc = [args.output_dir + "pihat_min0.2.genome",
                             args.output_dir + "zoom_pihat.genome"
                            ]
    output_relatedness_qc_plots = [args.output_dir + "relatedness.pdf",
                                   args.output_dir + "zoom_relatedness.pdf",
                                   args.output_dir + "hist_relatedness.pdf"
                                  ]
    output_relatedness_list = [args.output_dir + "fail_IBD.txt" ]
    output_remove_relatedness = [args.output_dir + "plink_remove_het.bed",
                                 args.output_dir + "plink_remove_het.fam",
                                 args.output_dir + "plink_remove_het.bim",
                                 args.output_dir + "plink_remove_het.log" ]

    if flag_start is True:
        # Find relatedness individuals to get ready to plot and remove (IBD)
        relatedness_qc = pf.run_relatedness_qc(inputs=[output_het[0], output_ld_qc[0]], outputs=output_relatedness_qc)
        relatedness_qc.result()
        wo.extend(relatedness_qc.outputs)

        # get relatedness individuals
        get_remove_relatedness = pf.run_get_remove_relatedness(inputs=[output_relatedness_qc[0],output_mq[0]], outputs=output_relatedness_list, pi_hat_value="0.2")
        get_remove_relatedness.result()
        wo.extend(get_remove_relatedness.outputs)

        # remove relatedness individuals
        remove_relatedness = pf.run_remove_relatedness(inputs=[output_het[0], output_relatedness_list[0]], outputs=output_remove_relatedness)
        remove_relatedness.result()
        wo.extend(remove_relatedness.outputs)

        # plot relatedness of individuals
        relatedness_qc_plot = pf.run_relatedness_qc_plot(inputs=output_relatedness_qc, outputs=output_relatedness_qc_plots)
        relatedness_qc_plot.result()
        wo.extend(relatedness_qc_plot.outputs)

    # Finish the QC section
    #####################################################################
    #####################################################################
    # Start the Stratification section
    # outline:
    # - download reference file (to be done externally)
    # - make sure SNPids are in the reference file
    # - make sure build is same for reference and my gwas samples
    # - perform QC on the reference file
    # - merge reference and my gwas sample
    # - stratisfy!


main()
