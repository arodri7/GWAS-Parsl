from parsl.app.app import python_app, bash_app
from parsl.configs.local_threads import config


#plink = "/Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/plink-1.07-mac-intel/plink"
#plink = "/Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/plink2-mac-intel/plink2"
plink = "/Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/plink-1.90-mac-intel/plink"

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
def plink_remove_missingness(inputs=[], outputs=[]):
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

# convert from vcf to plink format
@bash_app
def convert_vcf2plink(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".bed", "")
    cmd_line = "%s --vcf %s --max-alleles 2 --make-bed --out %s" % (plink, inputs[0], out_prefix)
    print(cmd_line)
    return cmd_line

# take plink input and assign unique indentifiers to the SNPs with a missing rs-identifier
@bash_app
def plink_assign_missing_ids(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".bed", "")
    b_prefix = inputs[0].replace(".bed", "")
    cmd_line = "%s --bfile %s --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out %s" % (plink, b_prefix, out_prefix)
    print(cmd_line)
    return cmd_line

# Remove SNPs with a low MAF frequency
@bash_app
def plink_remove_maf(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".bed", "")
    b_prefix = inputs[0].replace(".bed", "")
    cmd_line = '%s --bfile %s --maf 0.05 --allow-no-sex --make-bed --out %s' % (plink, b_prefix, out_prefix)
    return cmd_line

# Run the MDS removal part to generate plots for stratification
@bash_app
def plink_mds(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".mds", "")
    ref_b_prefix = inputs[0].replace(".bed", "")
    outd = os.path.dirname(outputs[0])
    gwas_b_prefix = inputs[1].replace(".bed", "")
    pop_loc_file = inputs[2]
    prune_f = inputs[3]
    out_covs_mds = outputs[4]

    # Extract the variants present in HapMap dataset from the 1000 genomes dataset.
    cmd_line = 'awk \'{print$2}\' %s.bim > %s/HapMap_SNPs.txt' % (gwas_b_prefix, outd)
    cmd_line += ';%s --bfile %s --extract %s/HapMap_SNPs.txt --make-bed --out %s/1kG_MDS6' % (plink, ref_b_prefix, outd, outd)

    # Extract the variants present in 1000 Genomes dataset from the HapMap dataset.
    cmd_line += ';awk \'{print$2}\' %s/1kG_MDS6.bim > %s/1kG_MDS6_SNPs.txt' % (outd,outd)
    cmd_line += ';%s --bfile %s --extract %s/1kG_MDS6_SNPs.txt --recode --make-bed --out %s/HapMap_MDS' % (plink, gwas_b_prefix, outd, outd)

    ## The datasets must have the same build. Change the build 1000 Genomes data build.
    cmd_line += ';awk \'{print$2,$4}\' %s/HapMap_MDS.map > %s/buildhapmap.txt' % (outd, outd)
    cmd_line += ';%s --bfile %s/1kG_MDS6 --update-map %s/buildhapmap.txt --make-bed --out %s/1kG_MDS7' % (plink, outd, outd, outd)
    # 1kG_MDS7 and HapMap_MDS now have the same build.

    # Prior to merging 1000 Genomes data with the HapMap data we want to make sure that the files are mergeable:
    # 1) Make sure the reference genome is similar in the HapMap and the 1000 Genomes Project datasets.
    # 2) Resolve strand issues.
    # 3) Remove the SNPs which after the previous two steps still differ between datasets.

    # Start those steps
    # 1) set reference genome
    cmd_line += ';awk \'{print$2,$5}\' %s/1kG_MDS7.bim > %s/1kg_ref-list.txt' % (outd, outd)
    cmd_line += ';%s --bfile %s/HapMap_MDS --reference-allele %s/1kg_ref-list.txt --make-bed --out %s/HapMap-adj' % (plink, outd, outd, outd)
    # The 1kG_MDS7 and the HapMap-adj have the same reference genome for all SNPs.

    # 2) Resolve strand issues.
    # Check for potential strand issues.
    cmd_line += ';awk \'{print$2,$5,$6}\' %s/1kG_MDS7.bim > %s/1kGMDS7_tmp' % (outd, outd)
    cmd_line += ';awk \'{print$2,$5,$6}\' %s/HapMap-adj.bim > %s/HapMap-adj_tmp' % (outd, outd)
    cmd_line += ';sort %s/1kGMDS7_tmp %s/HapMap-adj_tmp |uniq -u > %s/all_differences.txt' % (outd, outd, outd)

    ## Flip SNPs for resolving strand issues.
    # Print SNP-identifier and remove duplicates.
    cmd_line += ';awk \'{print$1}\' %s/all_differences.txt | sort -u > %s/flip_list.txt' % (outd, outd)
    # Generates a file of SNPs. These are the non-corresponding SNPs between the two files.
    # Flip the non-corresponding SNPs.
    cmd_line += ';%s --bfile %s/HapMap-adj --flip %s/flip_list.txt --reference-allele %s/1kg_ref-list.txt --make-bed --out %s/corrected_hapmap' % (plink, outd, outd, outd, outd)

    # Check for SNPs which are still problematic after they have been flipped.
    cmd_line += ';awk \'{print$2,$5,$6}\' %s/corrected_hapmap.bim > %s/corrected_hapmap_tmp' % (outd, outd)
    cmd_line += ';sort %s/1kGMDS7_tmp %s/corrected_hapmap_tmp |uniq -u  > %s/uncorresponding_SNPs.txt' % (outd, outd, outd)
    # This file demonstrates that there are differences between the files.

    # 3) Remove problematic SNPs from HapMap and 1000 Genomes.
    cmd_line += ';awk \'{print$1}\' %s/uncorresponding_SNPs.txt | sort -u > %s/SNPs_for_exlusion.txt' % (outd, outd)
    # The command above generates a list of the SNPs which caused the differences between
    #the HapMap and the 1000 Genomes data sets after flipping and setting of the reference genome.

    # Remove the problematic SNPs from both datasets.
    cmd_line += ';%s --bfile %s/corrected_hapmap --exclude %s/SNPs_for_exlusion.txt --make-bed --out %s/HapMap_MDS2' % (plink, outd, outd, outd)
    cmd_line += ';%s --bfile %s/1kG_MDS7 --exclude %s/SNPs_for_exlusion.txt --make-bed --out %s/1kG_MDS8' % (plink, outd, outd, outd)

    # Merge HapMap with 1000 Genomes Data.
    cmd_line += ';%s --bfile %s/HapMap_MDS2 --bmerge %s/1kG_MDS8.bed %s/1kG_MDS8.bim %s/1kG_MDS8.fam --allow-no-sex --make-bed --out %s/MDS_merge2' % (plink, outd, outd, outd, outd, outd)

    ## Perform MDS on HapMap-CEU data anchored by 1000 Genomes data.
    # Using a set of pruned SNPs
    cmd_line += ';%s --bfile %s/MDS_merge2 --extract %s --genome --out %s/MDS_merge2' % (plink, outd, prune_f, outd)
    cmd_line += ';%s --bfile %s/MDS_merge2 --read-genome %s/MDS_merge2.genome --cluster --mds-plot 10 --out %s/MDS_merge2' % (plink, outd, outd, outd)

    # Convert population codes into superpopulation codes (i.e., AFR,AMR,ASN, and EUR).
    cmd_line += ';awk \'{print$1,$1,$2}\' %s > %s/race_1kG.txt' % (pop_loc_file, outd)
    cmd_line += ';sed \'s/JPT/ASN/g\' %s/race_1kG.txt>%s/race_1kG2.txt' % (outd, outd)
    cmd_line += ';sed \'s/ASW/AFR/g\' %s/race_1kG2.txt>%s/race_1kG3.txt' % (outd, outd)
    cmd_line += ';sed \'s/CEU/EUR/g\' %s/race_1kG3.txt>%s/race_1kG4.txt' % (outd, outd)
    cmd_line += ';sed \'s/CHB/ASN/g\' %s/race_1kG4.txt>%s/race_1kG5.txt' % (outd, outd)
    cmd_line += ';sed \'s/CHD/ASN/g\' %s/race_1kG5.txt>%s/race_1kG6.txt' % (outd, outd)
    cmd_line += ';sed \'s/YRI/AFR/g\' %s/race_1kG6.txt>%s/race_1kG7.txt' % (outd, outd)
    cmd_line += ';sed \'s/LWK/AFR/g\' %s/race_1kG7.txt>%s/race_1kG8.txt' % (outd, outd)
    cmd_line += ';sed \'s/TSI/EUR/g\' %s/race_1kG8.txt>%s/race_1kG9.txt' % (outd, outd)
    cmd_line += ';sed \'s/MXL/AMR/g\' %s/race_1kG9.txt>%s/race_1kG10.txt' % (outd, outd)
    cmd_line += ';sed \'s/GBR/EUR/g\' %s/race_1kG10.txt>%s/race_1kG11.txt' % (outd, outd)
    cmd_line += ';sed \'s/FIN/EUR/g\' %s/race_1kG11.txt>%s/race_1kG12.txt' % (outd, outd)
    cmd_line += ';sed \'s/CHS/ASN/g\' %s/race_1kG12.txt>%s/race_1kG13.txt' % (outd, outd)
    cmd_line += ';sed \'s/PUR/AMR/g\' %s/race_1kG13.txt>%s/race_1kG14.txt' % (outd, outd)

    # Create a racefile of your own data.
    cmd_line += ';awk \'{print$1,$2,"OWN"}\' %s/HapMap_MDS.fam>%s/racefile_own.txt' % (outd, outd)

    # Concatenate racefiles.
    cmd_line += ';echo "FID IID race" > %s/header_line.txt' % (outd)
    #cmd_line += ';cat %s/race_1kG14.txt %s/racefile_own.txt | sed -e \'1i\FID IID race\' > %s/racefile.txt' % (outd, outd, outd)
    cmd_line += ';cat %s/header_line.txt %s/race_1kG14.txt %s/racefile_own.txt >  %s/racefile.txt' % (outd, outd, outd, outd)

    ## Exclude ethnic outliers.
    # Select individuals in HapMap data below cut-off thresholds.
    # The cut-off levels are not fixed thresholds but have to be determined based on the
    # visualization of the first two dimensions. To exclude ethnic outliers, the thresholds
    # need to be set around the cluster of population of interest.
    # For now we will hardcode to exclude non EURopeian
    cmd_line += ';awk \'{ if ($4 <-0.04 && $5 >0.03) print $1,$2 }\' %s/MDS_merge2.mds > %s/EUR_MDS_merge2' % (outd, outd)

    # Extract these individuals in HapMap data.
    cmd_line += ';%s --bfile %s --keep %s/EUR_MDS_merge2 --make-bed --out %s/HapMap_3_r3_13' % (plink, gwas_b_prefix, outd, outd)

    ## Create covariates based on MDS.
    # Perform an MDS ONLY on HapMap data without ethnic outliers.
    # The values of the 10 MDS dimensions are subsequently used as covariates in the
    # association analysis in the third tutorial.
    cmd_line += ';%s --bfile %s/HapMap_3_r3_13 --extract %s --genome --out %s/HapMap_3_r3_13' % (plink, outd, prune_f, outd)
    cmd_line += ';%s --bfile %s/HapMap_3_r3_13 --read-genome %s/HapMap_3_r3_13.genome --cluster --mds-plot 10 --out %s' % (plink, outd, outd, out_prefix)
    # Change the format of the .mds file into a plink covariate file.
    cmd_line += ';awk \'{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}\' %s.mds > %s' % (out_prefix, out_covs_mds)
    print(cmd_line)
    return(cmd_line)

# plot stratification of data
@python_app
def plot_mds(inputs=[], outputs=[]):
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri as py2ri
    import pandas as pd
    from pandas import read_csv

    robjects.r('''pdf("%s")''' % outputs[0] )
    #robjects.r('''data = read.table("%s", header=T)''' % inputs[0])
    #data = robjects.DataFrame.from_csvfile("%s" % inputs[0], header=True, sep="")
    data = read_csv("%s" % inputs[0], delimiter=r"\s+")
    #race = robjects.DataFrame.from_csvfile("%s" % inputs[1], header=True, sep=" ")
    race = read_csv("%s" % inputs[1], delimiter=r"\s+")
    datafile = pd.concat([data, race], axis=1,ignore_index=False, sort=True)
    for index, row in datafile.iterrows():
        if row['race'] == "EUR":
            robjects.r('''plot(%s,%s,type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="green")''' % (row['C1'], row['C5']))
        robjects.r('''par(new=T)''')
        if row['race'] == "ASN":
            robjects.r('''plot(%s,%s,type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="red")''' % (row['C1'], row['C2']))
        robjects.r('''par(new=T)''')
        if row['race'] == "AMR":
            robjects.r('''plot(%s,%s,type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col=470)''' % (row['C1'], row['C2']))
        robjects.r('''par(new=T)''')
        if row['race'] == "AFR":
            robjects.r('''plot(%s,%s,type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="blue")''' % (row['C1'], row['C2']))
        robjects.r('''par(new=T)''')
        robjects.r('''par(new=T)''')
        if row['race'] == "OWN":
            robjects.r('''plot(%s,%s,type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=3,cex=0.7,col="black")''' % (row['C1'], row['C2']))
        robjects.r('''par(new=T)''')
    robjects.r('''abline(v=-0.035,lty=3)''')
    robjects.r('''abline(h=0.035,lty=3)''')
    robjects.r('''legend("topright", pch=c(1,1,1,1,3),c("EUR","ASN","AMR","AFR","OWN"),col=c("green","red",470,"blue","black"),bty="o",cex=1)''')
    robjects.r('''dev.off()''')

# run association part
@bash_app
def plink_assoc(inputs=[], outputs=[]):
    out_prefix_assoc = outputs[0].replace(".assoc", "")
    out_prefix_logistic = outputs[4].replace(".assoc_2.logistic", "")
    b_prefix = inputs[0].replace(".bed", "")
    covar_file = inputs[1]
    cmd_line = '%s --bfile %s --assoc --out %s' % (plink, b_prefix, out_prefix_assoc)
    cmd_line += ';%s --bfile %s --covar %s --logistic --hide-covar --out %s' % (plink, b_prefix, covar_file, out_prefix_logistic)
    cmd_line += ';awk \'!/\'NA\'/\' %s.assoc.logistic > %s.assoc_2.logistic' % (out_prefix_logistic, out_prefix_logistic)
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
    robjects.r('''manhattan(results_as, ylim = c(0, 10), chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: assoc")''')
    robjects.r('''dev.off()''')
    robjects.r('''jpeg("%s")''' % outputs[2])
    robjects.r('''qq(results_log$P, main = "Q-Q plot of GWAS p-values : log")''')
    robjects.r('''dev.off()''')
    robjects.r('''jpeg("%s")''' % outputs[3])
    robjects.r('''qq(results_as$P, main = "Q-Q plot of GWAS p-values : log")''')
    robjects.r('''dev.off()''')


