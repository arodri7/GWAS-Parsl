import parsl
import os, argparse, sys
from parsl.app.app import python_app, bash_app
from parsl.configs.local_threads import config

#parsl.set_stream_logger() # <-- log everything to stdout

#print(parsl.__version__)
parsl.load(config)

#plink = "/Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/plink-1.07-mac-intel/plink"
plink = "/Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/plink2-mac-intel/plink2"
#steps = ["missingness_qc", "sex_discrepancy_qc", "maf_qc", "hwe_qc", "het_qc", "relatedness_qc"]
steps = []

# run association part
@bash app
def plink_assoc(inputs=[], outputs=[]):
    out_prefix_assoc = outputs[0].replace(".assoc", "")
    out_prefix_logistic = outputs[9].replace(".assoc_2.logistic", "")
    b_prefix = inputs[0].replace(".bed", "")
    covar_file = inputs[1]
    cmd_line = '%s --bfile %s --assoc --out %s' % (plink, b_prefix, out_prefix_assoc)
    cmd_line += ';%s --bfile %s --covar %s --logistic --hide-covar --out %s' % (plink, b_prefix, covar_file, out_prefix)
    cmd_line += 'awk \'!/\'NA\'/\' %s.assoc.logistic > %s.assoc_2.logistic' % (out_prefix_logistic, out_prefix_logistic)
    return cmd_line

# plot association of data
@python_app
def plot_assoc(inputs=[], outputs=[]):
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri as py2ri

    robjects.r('''library("qqman")''')
    robjects.r('''results_as = read.table("%s", header=T)''' % inputs[0])
    robjects.r('''results_log = read.table("%s", header=T)''' % inputs[1])
    robjects.r('''jpeg("%s")''' % outputs[0])
    robjects.r('''manhattan(results_log,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: logistic")''')
    robjects.r('''dev.off()''')
    robjects.r('''jpeg("%s")''' % outputs[1])
    robjects.r('''manhattan(results_as,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: assoc")''')
    robjects.r('''dev.off()''')

    robjects.r('''jpeg("%s")''' % outputs[2])
    robjects.r('''qq(results_log$P, main = "Q-Q plot of GWAS p-values : log")''')
    robjects.r('''dev.off()''')
    robjects.r('''jpeg("%s")''' % outputs[3])
    robjects.r('''qq(results_as$P, main = "Q-Q plot of GWAS p-values : log")''')
    robjects.r('''dev.off()''')


###########################
# parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Process GWAS pipeline for different steps.')
    parser.add_argument('--input-directory', dest='input_dir', 
                        help='Location of the input GWAS data')
    parser.add_argument('--population', dest='pop',
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
                 args.output_dir + "assoc_results.bed",
                 args.output_dir + "assoc_results.fam", 
                 args.output_dir + "assoc_results.bim",
                 args.output_dir + "assoc_results.log", 
                 args.output_dir + "logistic_results.bed",
                 args.output_dir + "logistic_results.fam",
                 args.output_dir + "logistic_results.bim",
                 args.output_dir + "logistic_results.log",
                 args.output_dir + "logistic_results.assoc_2.logistic"]

    if flag_start is True:
        # Run association and logistic
        s1 = plink_assoc(inputs=[args.input_dir], outputs=output_s1)
        s1.result()
        wo.extend(s1.outputs)

    output_s2 = [args.output_dir + "/Logistic_manhattan.jpeg",
                 args.output_dir + "/assoc_manhattan.jpeg",
                 args.output_dir + "/QQ-Plot_logistic.jpeg",
                 args.output_dir + "/QQ-Plot_assoc.jpeg"
                ]
    if flag_start is True:
        # Generate association plots.
        s2 = plot_assoc(inputs=[output_s1[0], output_s1[9]], outputs=output_s2)
        s2.result()
        wo.extend(s2.outputs)


main()
