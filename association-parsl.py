import parsl
import os, argparse, sys
from parsl.configs.local_threads import config
import parsl_funcs as pf

#parsl.set_stream_logger() # <-- log everything to stdout

# python association-parsl.py --input-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/2_Population_stratification/outputs/HapMap_3_r3_13 --covariates /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/2_Population_stratification/outputs/ALL.2of4intersection.20100804.genotypes.1kG_MDS.covariants.txt --output-directory /Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/GWA_tutorial/2_Population_stratification/outputs/

#print(parsl.__version__)
parsl.load(config)

#steps = ["missingness_qc", "sex_discrepancy_qc", "maf_qc", "hwe_qc", "het_qc", "relatedness_qc"]
steps = []

###########################
# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Process GWAS pipeline for different steps.')
    parser.add_argument('--input-directory', dest='input_dir', 
                        help='Location of the input GWAS data')
    parser.add_argument('--covariates', dest='covar',
                        help='population information of the 1000 genomes dataset')
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
    # Start the association section
    # outline:
    # - Binary traits association
    # - Logistic

    # get input arguments
    args = parse_args()
    print(args)
    flag_start = True
    # workflow outputs variable
    wo = []

    ## STEP 1: association and logistic
    #if (args.step_start == "missingness_qc" and flag_start == False) or args.step_start is None:
    #    flag_start = True

    #print("FLAG STEP1: %s" % flag_start)    
    output_s1 = [args.output_dir + "assoc_results.assoc",
                 args.output_dir + "assoc_results.log", 
                 args.output_dir + "logistic_results.assoc.logistic",
                 args.output_dir + "logistic_results.log",
                 args.output_dir + "logistic_results.assoc_2.logistic"]

    if flag_start is True:
        # Run association and logistic
        s1 = pf.plink_assoc(inputs=[args.input_dir, args.covar], outputs=output_s1)
        s1.result()
        wo.extend(s1.outputs)

    output_s2 = [args.output_dir + "/Logistic_manhattan.jpeg",
                 args.output_dir + "/assoc_manhattan.jpeg",
                 args.output_dir + "/QQ-Plot_logistic.jpeg",
                 args.output_dir + "/QQ-Plot_assoc.jpeg"
                ]
    if flag_start is True:
        # Generate association plots.
        s2 = pf.plot_assoc(inputs=[output_s1[0], output_s1[4]], outputs=output_s2)
        s2.result()
        wo.extend(s2.outputs)


main()
