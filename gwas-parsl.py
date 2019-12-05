import parsl
import os, argparse
from parsl.app.app import python_app, bash_app
from parsl.configs.local_threads import config
import sys

#parsl.set_stream_logger() # <-- log everything to stdout

#print(parsl.__version__)
parsl.load(config)

plink = "/Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/plink-1.07-mac-intel/plink"
steps = ["missingness_qc", "sex_discrepancy_qc", "maf_qc", "hwe_qc", "het_qc", "relatedness_qc"]


# Investigate missingness per individual and per SNP and generate histograms
# output: plink.imiss and plink.lmiss, these files show respectively the proportion
# of missing SNPs per individual and the proportion of missing individuals per SNP.
@bash_app
def missingness_qc(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".imiss", "")
    b_prefix = inputs[0].replace(".bed", "")
    cmd_line = '%s --bfile %s --missing --out %s' % (plink, b_prefix, out_prefix)
    return cmd_line

# Delete SNPs and individuals with high levels of missingness.
@bash_app
def remove_snps(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".bed", "")
    b_prefix =  inputs[0].replace(".bed", "")
    cmd_line = '%s --bfile %s --geno 0.2 --make-bed --out %s.tmp1' % (plink, b_prefix, out_prefix)
    cmd_line += '; %s --bfile %s.tmp1 --mind 0.2 --make-bed --out %s.tmp2' % (plink, out_prefix, out_prefix)
    cmd_line += '; %s --bfile %s.tmp2 --geno 0.02 --make-bed --out %s.tmp3' % (plink, out_prefix, out_prefix)
    cmd_line += '; %s --bfile %s.tmp3 --mind 0.02 --make-bed --out %s' % (plink, out_prefix, out_prefix)
    return cmd_line

# Generate plots to visualize the missingness results.
@python_app
def missingness_qc_plots(inputs=[], outputs=[]):
    import rpy2.robjects as robjects
    indmiss = robjects.DataFrame.from_csvfile("%s" % inputs[0], header=True, sep="", as_is=True) 
    robjects.r('pdf("%s")' % outputs[0] )
    robjects.r('hist(%s, main="Histogram individual missingness")' % (indmiss.rx(True, 6).r_repr())) #selects column 6, names header of file
    robjects.r('dev.off()')

    lndmiss = robjects.DataFrame.from_csvfile("%s" % inputs[1], header=True, sep="", as_is=True)
    robjects.r('pdf("%s")' % outputs[1] )
    robjects.r('hist(%s, main="Histogram SNP missingness")' % (indmiss.rx(True, 4).r_repr())) #selects column 5, names header of file
    robjects.r('dev.off()')
    return(outputs)

##############

# Check for sex discrepancy
# Subjects who were a priori determined as females must have a F value of <0.2, 
# and subjects who were a priori determined as males must have a F value >0.8. 
# This F value is based on the X chromosome inbreeding (homozygosity) estimate.
# Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.
@bash_app
def sexdiscrepancy_qc(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".sexcheck", "")
    b_prefix = inputs[0].replace(".bed", "")
    cmd_line = '%s --bfile %s --check-sex --out %s' % (plink, b_prefix, out_prefix)
    return cmd_line

# Delete sex discrepancy individuals.
@bash_app
def remove_sex_discrepancy(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".bed", "")
    b_prefix =  inputs[1].replace(".bed", "")
    sex_disc_txt = inputs[0].replace(".sexcheck", ".sex_discrepancy.txt")
    cmd_line = 'grep "PROBLEM" %s| awk \'{print$1,$2}\'> %s; %s --bfile %s --remove %s --make-bed --out %s' % (inputs[1], sex_disc_txt, plink, b_prefix, sex_disc_txt, out_prefix)
    return cmd_line

# Generate plots to visualize the sex discrepancy results.
@python_app
def sexdiscrepancy_qc_plots(inputs=[], outputs=[]):
    import rpy2.robjects as robjects
    infile = robjects.DataFrame.from_csvfile("%s" % inputs[0], header=True, sep="", as_is=True)
    robjects.r('pdf("%s")' % outputs[0] )
    robjects.r('hist(%s, main="Gender", xlab="F")' % (infile.rx(True, 6).r_repr())) #selects column 6, names header of file
    robjects.r('dev.off()')

    robjects.r('pdf("%s")' % outputs[1] )
    male = infile.rx(infile.rx2('PEDSEX').ro == '1', True)
    robjects.r('hist(%s, main="Men", xlab="F")' % (male.rx(True, 6).r_repr())) #selects column 6, names header of file
    robjects.r('dev.off()')

    robjects.r('pdf("%s")' % outputs[2] )
    female = infile.rx(infile.rx2('PEDSEX').ro == '2', True)
    robjects.r('hist(%s, main="Women", xlab="F")' % (female.rx(True, 6).r_repr())) #selects column 6, names header of file
    robjects.r('dev.off()')
    return(outputs)


#################################################
# Generate a bfile with autosomal SNPs only and delte SNPs with a low minor allele Frequency
# Select autosomal SNPs only from chr 1 to 22
@bash_app
def select_autosomal_snps(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".bed", "")
    b_prefix =  inputs[0].replace(".bim", "")
    snps_loc_txt = inputs[0].replace(".bim", ".snp_1_22.txt")
    cmd_line = 'awk \'{ if ($1 >= 1 && $1 <= 22) print $2 }\' %s > %s; %s --bfile %s --extract %s --make-bed --out %s' % (inputs[0], snps_loc_txt, plink, b_prefix, snps_loc_txt, out_prefix)
    return cmd_line

# maf QC
@bash_app
def maf_check_qc(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".frq", "")
    b_prefix = inputs[0].replace(".bed", "")
    cmd_line = '%s --bfile %s --freq --out %s' % (plink, b_prefix, out_prefix)
    return cmd_line

# Remove SNPs with a low MAF frequency
@bash_app
def remove_snps_low_maf(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".bed", "")
    b_prefix = inputs[0].replace(".bed", "")
    cmd_line = '%s --bfile %s --maf 0.05 --make-bed --out %s' % (plink, b_prefix, out_prefix)
    return cmd_line

# Generate plots to visualize the MAF distribution.
@python_app
def maf_distribution_qc_plots(inputs=[], outputs=[]):
    import rpy2.robjects as robjects
    infile = robjects.DataFrame.from_csvfile("%s" % inputs[0], header=True, sep="", as_is=True)
    robjects.r('pdf("%s")' % outputs[0] )
    robjects.r('hist(%s, main="MAF Distribution", xlab="MAF")' % (infile.rx(True, 5).r_repr())) #selects column 5, names header of file
    robjects.r('dev.off()')

##################
# Check HWE distribution and generate hwe SNPs
@bash_app
def check_hwe_distribution(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".hwe", "")
    b_prefix =  inputs[0].replace(".bed", "")
    hwe_loc_txt = outputs[1]
    cmd_line = '%s --bfile %s --hardy --out %s; awk \'{ if ($9 <0.00001) print $0 }\' %s > %s' % (plink, b_prefix, out_prefix, outputs[0], hwe_loc_txt)
    print (cmd_line)
    return cmd_line

# Generate plots to visualize the MAF distribution.
@python_app
def check_hwe_distribution_plots(inputs=[], outputs=[]):
    import rpy2.robjects as robjects
    infile = robjects.DataFrame.from_csvfile("%s" % inputs[0], header=True, sep="", as_is=True)
    robjects.r('pdf("%s")' % outputs[0] )
    robjects.r('hist(%s, main="Histogram HWE: strongly deviating SNPs only")' % (infile.rx(True, 9).r_repr())) #selects column 9, names header of file
    robjects.r('dev.off()')

# Remove HWE SNPs in two steps
@bash_app
def remove_hwe_snps(inputs=[], outputs=[]):
    out_prefix_s1 = outputs[0].replace(".bed", "_s1")
    out_prefix = outputs[0].replace(".bed", "")
    b_prefix =  inputs[0].replace(".bed", "")
    cmd_line = '%s --bfile %s --hwe 1e-6 --make-bed --out %s; %s --bfile %s --hwe 1e-10 --hwe-all --make-bed --out %s' % (plink, b_prefix, out_prefix_s1, plink, out_prefix_s1, out_prefix)
    return cmd_line

##########################

# Generate plots to visualize the HET distribution and generate list of individuals
# who deviate more than 3 std from the heterozygosity rate mean
@python_app
def run_het_qc_plots(inputs=[], outputs=[]):
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    r_base = importr('base')
    stats = importr('stats')
    utils_package = importr("utils")

    infile = robjects.DataFrame.from_csvfile("%s" % inputs[0], header=True, sep="", as_is=True)
    data = ((infile.rx2("N.NM.").ro - infile.rx2("O.HOM.")).ro)/infile.rx2("N.NM.")
    colnames = infile.colnames
    new_colnames = infile.colnames + "HET_RATE"
    infile = robjects.r.cbind(infile, data)
    infile.colnames = new_colnames

    robjects.r('pdf("%s")' % outputs[0] )
    robjects.r('hist(%s, main="Heterozygosity Rate", ylab="Frequency", xlab="Heterozygosity Rate")' % (data.r_repr())) #selects column 9, names header of file
    robjects.r('dev.off()')

    # generate list of individuals that deviate 3std form the mean
    het_fail = infile.rx(infile.rx2('HET_RATE').ro < (r_base.mean(infile.rx2('HET_RATE')).ro-(robjects.r.sd(infile.rx2('HET_RATE')).ro * 3))[0], True) or (infile.rx(infile.rx2('HET_RATE').ro > (r_base.mean(infile.rx2('HET_RATE')).ro+(robjects.r.sd(infile.rx2('HET_RATE')).ro * 3))[0], True))
    het_dst = (het_fail.rx2("HET_RATE").ro - r_base.mean(infile.rx2('HET_RATE'))).ro/(robjects.r.sd(infile.rx2('HET_RATE')))
    new_colnames = het_fail.colnames + "HET_DST"
    het_fail = robjects.r.cbind(het_fail, het_dst)
    het_fail.colnames = new_colnames
    utils_package.write_table(het_fail, outputs[1], sep="\t", row_names=False)

# Remove LD regions and generate pruned data set
@bash_app
def run_remove_ld(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".prune.in", "")
    b_prefix =  inputs[0].replace(".bed", "")
    inversions_regions = inputs[1]
    plink_check_het = outputs[1].replace(".het", "")
    cmd_line = '%s --bfile %s --exclude %s --range --indep-pairwise 50 5 0.2 --out %s; %s --bfile %s --extract %s --het --out %s' % (plink, b_prefix, inversions_regions, out_prefix, plink, b_prefix, outputs[0], plink_check_het)
    print (cmd_line)
    return cmd_line

# Remove regions that failed het qc
@bash_app
def run_remove_het(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".bed", "")
    b_prefix =  inputs[0].replace(".bed", "")
    het_fail = inputs[1]
    het_fail_ind = "%s.mod.txt" % het_fail
    cmd_line = 'sed \'s/"// g\' %s | awk \'{print$1, $2}\'> %s; %s --bfile %s --remove %s --make-bed --out %s' % (het_fail, het_fail_ind, plink, b_prefix, het_fail_ind, out_prefix)
    print (cmd_line)
    return cmd_line

########################
# Identify related individuals to remove
@bash_app
def run_relatedness_qc(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".genome", "")
    b_prefix =  inputs[0].replace(".bed", "")
    pihat_file = outputs[1]
    cmd_line = '%s --bfile %s --extract %s --genome --min 0.2 --out %s; awk \'{ if ($8 >0.9) print $0 }\' %s >%s' % (plink, b_prefix, inputs[1], out_prefix, outputs[0], pihat_file)
    print (cmd_line)
    return cmd_line

# get list of relatedness individuals to remove
@python_app
def run_get_remove_relatedness(inputs=[], outputs=[], pi_hat_value=None):
    import pandas as pd

    EOL=chr(10)

    str_type = {"FID":str, "IID":str, "FID1":str, "IID1":str, "FID2":str, 'IID2':str}

    imissf = pd.read_csv(inputs[1],delim_whitespace=True,dtype=str_type)
    imissf.set_index(["FID","IID"],inplace=True)
    genomef = pd.read_csv(inputs[0],delim_whitespace=True,usecols=["FID1","IID1","FID2","IID2","PI_HAT"],dtype=str_type)

    outf   =open(outputs[0],"w")
    super_pi_hat = float(pi_hat_value)

    def getDegrees(remove):
        elts = set(imissf.index.values)
        degd = {}
        rel  = {}
        for elt in elts: 
            degd[elt]=0
            rel[elt]=[]
        deg = pd.Series(degd)
        for i,row in genomef.iterrows():
            x=tuple(row[["FID1","IID1"]].tolist())
            y=tuple(row[["FID2","IID2"]].tolist())
            if x in remove or y in remove : pass
            try:
                deg[x]=deg[x]+1
                deg[y]=deg[y]+1
                rel[x].append(y)
                rel[y].append(x)
            except:
                print("saw the problem")
        return rel, deg

    remove = set(map (tuple,genomef[genomef['PI_HAT']>super_pi_hat][["FID1","IID1"]].to_records(index=False)))\
             | set(map(tuple,genomef[genomef['PI_HAT']>super_pi_hat][["FID2","IID2"]].to_records(index=False)))

    rel, deg = getDegrees(remove) 

    candidates = deg[deg>=1].sort_values(ascending=False)

    for i,c in candidates.iteritems():
        if deg[i]>0:
            remove.add(i)
            deg[i]=deg[i]-1
            for other in rel[i]:
                deg[other]=deg[other]-1

    remove = sorted(list(remove))
    outf.write(EOL.join(map (lambda x: "%s %s"%(x[0],x[1]),remove)))
    outf.close()


# remove relatedness individuals
@bash_app
def run_remove_relatedness(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".bed", "")
    b_prefix =  inputs[0].replace(".bed", "")
    remove_file = inputs[1]
    cmd_line = '%s --bfile %s --remove %s --make-bed --out %s' % (plink, b_prefix, remove_file, out_prefix)
    print (cmd_line)
    return cmd_line

# plot relatedness individuals
@python_app
def run_relatedness_qc_plot(inputs=[], outputs=[]):
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri as py2ri

    robjects.r('''pdf("%s")''' % outputs[2] )
    robjects.r('''relatedness = read.table("%s", header=T)''' % inputs[0])
    robjects.r('''hist(relatedness[,10],main="Histogram relatedness", xlab= "Pihat")''')
    robjects.r('''dev.off()''')

    #infile = robjects.DataFrame.from_csvfile("%s" % inputs[0], header=True, sep="", as_is=True)
    #robjects.r('hist(%s, main="Histogram Relatedness", xlab="Pihat")' % (infile.rx(True, 10).r_repr()))
    #robjects.r('dev.off()')

    robjects.r('''pdf("%s")''' % outputs[0])
    robjects.r('''relatedness = read.table("%s", header=T)''' % inputs[0])
    robjects.r('''with(relatedness,plot(Z0,Z1, xlim=c(0,1), ylim=c(0,1), type="n"))''')
    robjects.r('''with(subset(relatedness,RT=="PO") , points(Z0,Z1,col=4))''')
    robjects.r('''with(subset(relatedness,RT=="UN") , points(Z0,Z1,col=3))''')
    robjects.r('''legend(1,1, xjust=1, yjust=1, legend=levels(relatedness$RT), pch=16, col=c(4,3))''')
    robjects.r('''dev.off()''')

    robjects.r('''pdf("%s")''' % outputs[1])
    robjects.r('''relatedness_zoom = read.table("%s", header=T)''' % inputs[1])
    robjects.r('''par(pch=16, cex=1)''')
    robjects.r('''with(relatedness_zoom,plot(Z0,Z1, xlim=c(0,0.02), ylim=c(0.98,1), type="n"))''')
    robjects.r('''with(subset(relatedness_zoom,RT=="PO") , points(Z0,Z1,col=4))''')
    robjects.r('''with(subset(relatedness_zoom,RT=="UN") , points(Z0,Z1,col=3))''')
    robjects.r('''legend(0.02,1, xjust=1, yjust=1, legend=levels(relatedness$RT), pch=16, col=c(4,3))''')
    robjects.r('''dev.off()''')


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
        mq = missingness_qc(inputs=[args.input_dir], outputs=[args.output_dir + "plink_missingness.imiss"])
        #if flag_start is True: mq.result()
        mq.result()
        wo.extend(mq.outputs)

        # Generate missingness plots
        mq_plots = missingness_qc_plots(inputs=output_mq, outputs=output_mq_plots)
        #if flag_start is True: mq_plots.result()
        mq_plots.result()
        wo.extend(mq_plots.outputs)

        # Delete SNPs and individuals with high levels of missingness
        rsnps = remove_snps(inputs=[args.input_dir], outputs=rsnps_outputs)
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
        sq = sexdiscrepancy_qc(inputs=rsnps_outputs, outputs=output_sq)
        #if flag_start is True: sq.result()
        sq.result()
        wo.extend(sq.outputs)

        # Generate sex discrepancy plots
        sq_plots = sexdiscrepancy_qc_plots(inputs=output_sq, outputs=output_sq_plots)
        #if flag_start is True: sq_plots.result()
        sq_plots.result()
        wo.extend(sq_plots.outputs)

        # Delete Sex discrepancy individuals
        rsex = remove_sex_discrepancy(inputs=[output_sq[0], rsnps_outputs[0]], outputs=rsex_outputs)
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
    autosomal = select_autosomal_snps(inputs=[rsex_outputs[3]], outputs=output_autosomal_snps)
    if flag_start is True: autosomal.result()
    #autosomal.result()
    wo.extend(autosomal.outputs)

    if flag_start is True:
        # Generate MAF distribution
        maf_dist = maf_check_qc(inputs=output_autosomal_snps, outputs=output_maf_check_qc)
        #if flag_start is True: maf_dist.result()
        maf_dist.result()
        wo.extend(maf_dist.outputs)

        # Generate maf distribution plots
        maf_plots = maf_distribution_qc_plots(inputs=output_maf_check_qc, outputs=output_maf_plots)
        #if flag_start is True: maf_plots.result()
        maf_plots.result()
        wo.extend(maf_plots.outputs)

        # Delete SNPs with low MAF frequency
        low_maf_remove = remove_snps_low_maf(inputs=output_autosomal_snps, outputs=output_delete_maf)
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
        hwe_dist = check_hwe_distribution(inputs=output_delete_maf, outputs=output_hwe_qc)
        #if flag_start is True: hwe_dist.result()
        hwe_dist.result()
        wo.extend(hwe_dist.outputs)
   
        # HWE distribution plots
        hwe_dist_plots = check_hwe_distribution_plots(inputs=output_hwe_qc, outputs=output_hwe_qc_plots)
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
        remove_hwe = remove_hwe_snps(inputs=output_delete_maf, outputs=output_hwe)
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
        remove_ld = run_remove_ld(inputs=[output_hwe[0], args.inversion_regions], outputs=output_ld_qc)
        remove_ld.result()
        wo.extend(remove_ld.outputs)

        # Plot heterozygosity rate distribution and generate list of individuals who 
        # deviate more than 3 std from the heterozygosity rate mean 
        het_qc_plots = run_het_qc_plots(inputs=[output_ld_qc[2]], outputs=output_het_qc_plots)
        het_qc_plots.result()
        wo.extend(het_qc_plots.outputs)

        # remove the failed het individuals
        remove_het = run_remove_het(inputs=[output_hwe[0], output_het_qc_plots[1]], outputs=output_het)
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
        relatedness_qc = run_relatedness_qc(inputs=[output_het[0], output_ld_qc[0]], outputs=output_relatedness_qc)
        relatedness_qc.result()
        wo.extend(relatedness_qc.outputs)

        # get relatedness individuals
        get_remove_relatedness = run_get_remove_relatedness(inputs=[output_relatedness_qc[0],output_mq[0]], outputs=output_relatedness_list, pi_hat_value="0.2")
        get_remove_relatedness.result()
        wo.extend(get_remove_relatedness.outputs)

        # remove relatedness individuals
        remove_relatedness = run_remove_relatedness(inputs=[output_het[0], output_relatedness_list[0]], outputs=output_remove_relatedness)
        remove_relatedness.result()
        wo.extend(remove_relatedness.outputs)

        # plot relatedness of individuals
        relatedness_qc_plot = run_relatedness_qc_plot(inputs=output_relatedness_qc, outputs=output_relatedness_qc_plots)
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
