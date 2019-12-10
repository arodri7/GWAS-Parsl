import parsl
import os, argparse, sys
from parsl.app.app import python_app, bash_app
from parsl.configs.local_threads import config

#parsl.set_stream_logger() # <-- log everything to stdout

#print(parsl.__version__)
parsl.load(config)

#plink = "/Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/plink-1.07-mac-intel/plink"
#plink = "/Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/plink2-mac-intel/plink2"
plink = "/Users/arodri7/Documents/Work/DOE-MVP/GWAS-VA/plink-1.90-mac-intel/plink"
#steps = ["missingness_qc", "sex_discrepancy_qc", "maf_qc", "hwe_qc", "het_qc", "relatedness_qc"]
steps = ["convert_to_plink", "qc_fill_ids", "missingness_qc", "maf_qc", "merge_qc", "plot_qc" ]

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

# Remove SNPs and individuals with high levels of missingness.
@bash_app
def plink_remove_missingness(inputs=[], outputs=[]):
    out_prefix = outputs[0].replace(".bed", "")
    b_prefix =  inputs[0].replace(".bed", "")
    # Remove variants based on missing genotype data.
    cmd_line = '%s --bfile %s --geno 0.2 --make-bed --out %s.tmp1' % (plink, b_prefix, out_prefix)
    # Remove individuals based on missing genotype data.
    cmd_line += '; %s --bfile %s.tmp1 --mind 0.2 --make-bed --out %s.tmp2' % (plink, out_prefix, out_prefix)
    # Remove variants based on missing genotype data.
    cmd_line += '; %s --bfile %s.tmp2 --geno 0.02 --make-bed --out %s.tmp3' % (plink, out_prefix, out_prefix)
    # Remove individuals based on missing genotype data.
    cmd_line += '; %s --bfile %s.tmp3 --mind 0.02 --make-bed --out %s' % (plink, out_prefix, out_prefix)
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
    #robjects.r('''datafile<- merge(%s,%s,by=c("IID","FID"))''' % (data.r_repr(), race.r_repr()))
    #print(datafile)
    #robjects.r('''for (i in 1:nrow(datafile)){''')
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
        if row['race'] == "OWN":
            robjects.r('''plot(%s,%s,type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=3,cex=0.7,col="black")''' % (row['C1'], row['C2']))
        robjects.r('''par(new=T)''')
    robjects.r('''abline(v=-0.035,lty=3)''')
    robjects.r('''abline(h=0.035,lty=3)''')
    robjects.r('''legend("topright", pch=c(1,1,1,1,3),c("EUR","ASN","AMR","AFR","OWN"),col=c("green","red",470,"blue","black"),bty="o",cex=1)''')
    robjects.r('''dev.off()''')
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
        s1 = convert_vcf2plink(inputs=[args.reference], outputs=output_vcf2plink)
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
        s2 = plink_assign_missing_ids(inputs=output_vcf2plink, outputs=output_s2)
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
        s3 = plink_remove_missingness(inputs=output_s2, outputs=output_s3)
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
        s4 = plink_remove_maf(inputs=output_s3, outputs=output_s4)
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
        s5 = plink_mds(inputs=[output_s4[0], args.input_dir, args.pop, args.prune], outputs=output_s5)
        s5.result()
        wo.extend(s5.outputs)

    ## STEP 6: Plot stratification
    if (args.step_start == "plot_qc" and flag_start == False) or args.step_start is None:
        flag_start = True
    output_s6 = [args.output_dir + "/MDS.pdf"]
    if flag_start is True:
        # Generate population stratification plot.
        s6 = plot_mds(inputs=[output_s5[5], output_s5[6]], outputs=output_s6)
        s6.result()
        wo.extend(s6.outputs)


main()
