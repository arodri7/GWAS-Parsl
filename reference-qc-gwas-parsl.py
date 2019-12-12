import parsl
import os, argparse, sys
from parsl.configs.local_threads import config
import parsl_funcs as pf

#parsl.set_stream_logger() # <-- log everything to stdout

# python reference-qc-gwas-parsl.py --input-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/outputs/plink_remove_het.bed --reference /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/references/ALL.2of4intersection.20100804.genotypes.vcf.gz --population /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/references/20100804.ALL.panel --output-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/2_Population_stratification/outputs/ --prune /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/1_QC_GWAS/outputs/indepSNP.prune.in --step-start plot_qc

#print(parsl.__version__)
parsl.load(config)

steps = ["convert_to_plink", "qc_fill_ids", "missingness_qc", "maf_qc", "merge_qc", "plot_qc" ]

###########################
# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Process GWAS pipeline for different steps.')
    parser.add_argument('--input-directory', dest='input_dir', 
                        help='Location of the input GWAS data')
    parser.add_argument('--reference', dest='reference',
                        help='Location of the reference file')
    parser.add_argument('--population', dest='pop',
                        help='population information of the 1000 genomes dataset')
    parser.add_argument('--prune', dest='prune',
                        help='prune file')
    parser.add_argument('--output-directory', dest='output_dir',
                        help='Location of the outputs GWAS data')
    parser.add_argument('--step-start', dest='step_start',
                        choices=steps,
                        help='Step to start in. The assumption is that the input must be ready for the step you are starting from.')
    return(parser.parse_args())

# Main function controls the inputs and steps to submit
def main():
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

    # get input arguments
    args = parse_args()
    print(args)
    flag_start = False
    # workflow outputs variable
    wo = []

    ## STEP 1: Convert VCF to Plink
    if (args.step_start == "convert_to_plink" and flag_start == False) or args.step_start is None:
        flag_start = True

    #print("FLAG STEP1: %s" % flag_start)    
    name_prefix = "/" + os.path.basename(args.reference.replace(".vcf.gz", ""))
    output_vcf2plink = [args.output_dir + name_prefix + ".bed",
                        args.output_dir + name_prefix + ".fam", 
                        args.output_dir + name_prefix + ".bim",
                        args.output_dir + name_prefix + ".log"]

    if flag_start is True:
        # Run missingness QC
        s1 = pf.convert_vcf2plink(inputs=[args.reference], outputs=output_vcf2plink)
        s1.result()
        wo.extend(s1.outputs)


    ## STEP 2: Insert unique missing ids to plink files
    if (args.step_start == "qc_fill_ids" and flag_start == False) or args.step_start is None:
        flag_start = True
    output_s2 = [args.output_dir + name_prefix + ".genotypes_no_missing_IDs.bed",
                 args.output_dir + name_prefix + ".genotypes_no_missing_IDs.fam",
                 args.output_dir + name_prefix + ".genotypes_no_missing_IDs.bim",
                 args.output_dir + name_prefix + ".genotypes_no_missing_IDs.log"]
    if flag_start is True:
        # Insert unique missing ids to plink files
        s2 = pf.plink_assign_missing_ids(inputs=output_vcf2plink, outputs=output_s2)
        s2.result()
        wo.extend(s2.outputs)


    ## STEP 3: Insert unique missing ids to plink files
    if (args.step_start == "missingness_qc" and flag_start == False) or args.step_start is None:
        flag_start = True
    output_s3 = [args.output_dir + name_prefix + ".1kG_MIS.bed",
                 args.output_dir + name_prefix + ".1kG_MIS.fam",
                 args.output_dir + name_prefix + ".1kG_MIS.bim",
                 args.output_dir + name_prefix + ".1kG_MIS.log"]
    if flag_start is True:
        # Delete SNPs and individuals with high levels of missingness
        s3 = pf.plink_remove_missingness(inputs=output_s2, outputs=output_s3)
        s3.result()
        wo.extend(s3.outputs)

    ## STEP 4: MAF QC
    if (args.step_start == "maf_qc" and flag_start == False) or args.step_start is None:
        flag_start = True
    output_s4 = [args.output_dir + name_prefix + ".1kG_MAF.bed",
                 args.output_dir + name_prefix + ".1kG_MAF.fam",
                 args.output_dir + name_prefix + ".1kG_MAF.bim",
                 args.output_dir + name_prefix + ".1kG_MAF.log",
                ]
    if flag_start is True:
        # Removal based on MAF
        s4 = pf.plink_remove_maf(inputs=output_s3, outputs=output_s4)
        s4.result()
        wo.extend(s4.outputs)

    ## STEP 5: Merge datasets
    if (args.step_start == "merge_qc" and flag_start == False) or args.step_start is None:
        flag_start = True
    output_s5 = [args.output_dir + name_prefix + ".1kG_MDS.mds",
                 args.output_dir + name_prefix + ".1kG_MDS.cluster1",
                 args.output_dir + name_prefix + ".1kG_MDS.cluster2",
                 args.output_dir + name_prefix + ".1kG_MDS.log",
                 args.output_dir + name_prefix + ".1kG_MDS.covariants.txt",
                 args.output_dir + "/MDS_merge2.mds",
                 args.output_dir + "/racefile.txt"
                ]
    if flag_start is True:
        # Merge hapmap and my data
        s5 = pf.plink_mds(inputs=[output_s4[0], args.input_dir, args.pop, args.prune], outputs=output_s5)
        s5.result()
        wo.extend(s5.outputs)

    ## STEP 6: Plot stratification
    if (args.step_start == "plot_qc" and flag_start == False) or args.step_start is None:
        flag_start = True
    output_s6 = [args.output_dir + "/MDS.pdf"]
    if flag_start is True:
        # Generate population stratification plot.
        s6 = pf.plot_mds(inputs=[output_s5[5], output_s5[6]], outputs=output_s6)
        s6.result()
        wo.extend(s6.outputs)


main()
