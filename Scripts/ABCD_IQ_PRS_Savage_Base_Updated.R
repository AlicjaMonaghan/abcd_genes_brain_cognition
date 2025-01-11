# This script details calculation of polygenic risk scores (PRSs) for 
# intelligence (IQ) in the Adolescent Brain and Cognitive Development (ABCD) 
# Study, using protocol and guidelines developed by Choi, Mak, and O'Reilly 
# (2020), with an  accompanying GitHub tutorial 
# (https://choishingwan.github.io/PRS-Tutorial/), through R version 4.0.3. Our 
# base data has been provided by Savage and colleagues (2018), who conducted a 
# GWAS on ~260,000 participants. 

### PART 1 - Setting up the Work Space ####
# Several packages need to be installed using BioConductor/BiocManager, such as
# rtracklayer...
rm(list = ls())
library(plink)
library(plinkQC)
library(data.table)
library(magrittr)
library(ssh)
library(genio)
library(rtracklayer)
library(dplyr)
library(readxl)
library(tidyverse)
library(stats)
library(remotes)
library(bigsnpr)
library(runonce)
library(snpStats)
library(hexbin)
library(biomaRt)
library(ggpubr)

options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

setwd("//cbsu/data/Imaging/astle/am10/ABCD/Data")

# Load the base data
IQ_Base_Savage = fread('Savage_IQ_2018_sumstats/SavageJansen_2018_intelligence_metaanalysis.txt')

# Load the batch information for the current release! There are 11099 participants
# with genomic data! 
ABCD_Release3_Batch_Info = fread('ABCD_Genomics/ABCD_release3.0_.batch_info.txt')
# Correct formatting for some participant numbers here i.e. ensure there is a _
# after NDAR
ABCD_Release3_Batch_Info$abcd.id_redcap[6184] <- 'NDAR_INVPWLFYWTX'
ABCD_Release3_Batch_Info$abcd.id_redcap[2881] <- 'NDAR_INVF3FYXH1G'

### PART 1B - Pull Demographic Information for ABCD ####
# Some demographic data from the ABCD genomics files are missing, such as sex. 
# We also need to create covariate (sex) and phenotype (cognition) files, and 
# save them within the same directory as the ABCD genomics files.

# Load up the batch information file which accompanied the 4th data release (
# titled ABCD_release3.0_.batch_info)
Batch_Info = fread("ABCD_Genomics/ABCD_release3.0_.batch_info.txt")
# Exclude participants with axiom plate 461, which was identified as problematic
Batch_Info_Retained = Batch_Info[!which(Batch_Info$Axiom_Plate==461)]
colnames(Batch_Info_Retained) = c("FID","SUBJECTKEY","Axiom Plate", "BATCH")

# This is the fam file, containing the family ID (V1), individual ID (V2), 
# paternal ID (V3), maternal ID (V4), sex (V5) and phenotype (V6). We shall 
# update the sex column! This contains 11099 participants. 
ABCD_FAM_File = fread('ABCD_Genomics/genotype_QCed/ABCD_release_3.0_QCed.fam')
colnames(ABCD_FAM_File)[2] = "SUBJECTKEY"
# Now correct those with mis-formatted subject keys
ABCD_FAM_File$SUBJECTKEY[6184] <- 'NDAR_INVPWLFYWTX'
ABCD_FAM_File$SUBJECTKEY[2881] <- 'NDAR_INVF3FYXH1G'

# Load up cognitive and sex information for the entire baseline (bar 2 subjects)
# , and correct mis-formatted participant IDs. The parental demographics data
# frame contains all imputed cognitive and demographic data for all children
# from the second ABCD release, bar 2. 
Participant_Demographics = as.data.frame(fread("ABCD_Participant_Information/Parental_Demographic_Imputed.csv"))
# Now select the demographics for participants with good QC genomic data
# Participant_Demographics_QC = Participant_Demographics[which(ABCD_FAM_File$SUBJECTKEY %in% Participant_Demographics$SUBJECTKEY),]
ABCD_NIH_TB_Total_Cognition_with_Subject_Keys = dplyr::select(Participant_Demographics, SUBJECTKEY, sex, PC1_Scores)
ABCD_NIH_TB_Total_Cognition_with_Subject_Keys$SUBJECTKEY[1] = "NDAR_MC003PZF"

# Merge the ABCD_FAM_File with the above, and correct mis-formatted participant IDs
ABCD_FAM_File_with_Demographics = merge(ABCD_NIH_TB_Total_Cognition_with_Subject_Keys, ABCD_FAM_File,by = c("SUBJECTKEY"))
ABCD_FAM_File_with_Demographics = merge(ABCD_FAM_File_with_Demographics, Batch_Info_Retained, by = c("SUBJECTKEY"))

### PART 1C - Reading in Individual ABCD Genomics Files and Checking ####
# We need to check that the rsid's for ABCD are all properly formatted and, if
# not, pull the correct information.
ABCD_BIM = fread('/imaging/astle/am10/ABCD/Data/ABCD_Genomics/genotype_QCed/ABCD_release_3.0_QCed.bim')
colnames(ABCD_BIM) = c("chr","rsid","coord", "pos","A1","A2")
Misformatted_RSID = ABCD_BIM[which(ABCD_BIM$rsid %like% "chr")]

# Separate queries according to chromosome to allow for easier processing
for (i in 1:length(unique(ABCD_BIM$chr))){
  filename = paste("Misformatted_RSID_Chromosome_", i, ".txt", sep="")
  write.table(Misformatted_RSID$rsid[which(Misformatted_RSID$chr==i)], filename, row.names = F, col.names = F, quote = F )
}
# Search for the rsid's for the relevant misformatted rsid's, which are currently
# in the format chrX:XXXXXXX
snp_mart <- useMart(biomart="ENSEMBL_MART_SNP", host="https://grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")
Misformatted_RSID$'chromosomal_region' = paste(Misformatted_RSID$chr,Misformatted_RSID$pos,sep=":",Misformatted_RSID$pos)
Misformatted_RSID_Query = vector(mode = "list", length = length(Misformatted_RSID$chr))

for (i in 1:length(Misformatted_RSID$chr)){
  Misformatted_RSID_Query[[i]] = getBM(attributes = c("refsnp_id",'allele', 'chrom_start'), filters = "chromosomal_region",values = Misformatted_RSID$chromosomal_region[i], mart = snp_mart)
  print(paste("Processed query ", i, " of 8320", sep = ""))
}
save(Misformatted_RSID_Query, file = "Misformatted_RSID_Query.RData")

# To avoid ambiguous SNPs, we shall only select RSIDs which have been matched to
# two alleles e.g. G/A, but not GACACAA. Note '/' counts as a character!
RSID_Query_Dataframe = do.call(rbind,Misformatted_RSID_Query)
RSID_Query_Dataframe_Unique = RSID_Query_Dataframe[which(nchar(RSID_Query_Dataframe$allele)==3),]
# Add in the previous RSID names!
Misformatted_RSID$'chrom_start' = Misformatted_RSID$pos
RSID_Query_Dataframe_Unique = merge(RSID_Query_Dataframe_Unique, Misformatted_RSID[,c(2,7)], by = "chrom_start")
# Reorder the above file
RSID_Query_Dataframe_Unique = RSID_Query_Dataframe_Unique[,c(4,2,1,3)]
# Save the above
write.table(RSID_Query_Dataframe_Unique, "ABCD_Genomics/RSID_Query_Dataframe_Unique.txt", sep="\t", row.names = F, col.names = F, quote = FALSE)
# remove SNPs which have a missing allele
Missing_Allele = which(RSID_Query_Dataframe_Unique$allele %like% "-")
RSID_Query_Dataframe_Unique = RSID_Query_Dataframe_Unique[-c(Missing_Allele),]
write.table(RSID_Query_Dataframe_Unique[,c(1:2)], "ABCD_Genomics/RSID_Query_Dataframe_Unique_Two_Columns.txt", sep="\t", row.names = F, col.names = F, quote = FALSE)

### PART 1D - Update PLINK files with phenotype, rsid, and sex information ####
# First, note that 120 participants have missing demographic information i.e. all
# respective rows are empty. Therefore, we shall create a new list of 
# participants to keep, and update the phenotypic and sex information for them!
Participants_with_Complete_Data = dplyr::select(ABCD_FAM_File_with_Demographics, V1, SUBJECTKEY)
write.table(Participants_with_Complete_Data, "ABCD_Genomics/Genomic_Participants_with_Complete_Data.txt", sep = "\t", row.names = FALSE, col.names = F, quote = F)

# Create a file with sex information
ABCD_Sex_File = dplyr::select(ABCD_FAM_File_with_Demographics, V1, SUBJECTKEY, sex)
# Re-code females as 2 and males as 1
ABCD_Sex_File$sex = as.factor(ABCD_Sex_File$sex)
levels(ABCD_Sex_File)[1] = "2"
levels(ABCD_Sex_File)[2] = "1"
write.table(ABCD_Sex_File, "ABCD_Genomics/ABCD_Sex_File.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# And now create a file with the phenotype information
ABCD_Pheno_File = dplyr::select(ABCD_FAM_File_with_Demographics, V1, SUBJECTKEY, PC1_Scores)
write.table(ABCD_Pheno_File, "ABCD_Genomics/ABCD_Pheno_File.txt", sep="\t",row.names = F, col.names = F, quote = FALSE)

# Now update the phenotypes and sex information! Note, that for the main PRS 
# analyses and QC, we shall only use the non-autosomal chromosomes i.e. 1:22,24,
# 25,26
session = ssh_connect("am10@login-f01")
Update_Sex_and_Pheno_PLINK = ssh_exec_wait(session, command = "module load plink;
                                            cd /imaging/astle/am10/ABCD/Data/ABCD_Genomics;
                                           plink --bfile genotype_QCed/ABCD_release_3.0_QCed --update-sex ABCD_Sex_File.txt --update-name RSID_Query_Dataframe_Unique_Two_Columns.txt --pheno ABCD_Pheno_File.txt --keep Genomic_Participants_with_Complete_Data.txt --make-bed --out /imaging/projects/external/abcd/analyses/Genomics/ABCD_Updated")

# Thus, our dataset which is UPDATED is now: /imaging/projects/external/abcd/analyses/Genomics/ABCD_Updated

### PART 2 - Quality Control (QC) of Base Data #### 
# Note we shall use base summary GWAS statistics from Savage and colleagues 
# (2018), who conducted a GWAS meta-analysis of 269,867 healthy participants, 
# ranging from children to older adults, and found 205 genomic loci and 1016 
# genes significantly linked to IQ. Note that 9,295,118 SNPs were initially 
# included in the base summary GWAS statistics (i.e. before QC). 

### REQUIREMENT 1: HERITABILITY CHECK ####
# Savage and colleagues (201) found a heritability score of single nucleotide 
# polymorphisms (SNPs) of .19. Choi, Mak and  O'Reilly recommend that base 
# statistics should have a SNP heritability larger than .05, hence we shall 
# continue using summary statistics from Savage and colleagues (2018). 

### REQUIREMENT 2: GENOME BUILD ####
# Both the base and target data appear to have the same genome build (GRCh37)!
### REQUIREMENT 3: EFFECT ALLELE ####
# The data produced from Savage and colleagues (2018) clearly differentiates 
# between the effect allele (A1) and non-effect allele (A2). However, we need to
# ensure that the A1 and A2 alleles are upper-class. 
IQ_Base_Savage$A1 = toupper(IQ_Base_Savage$A1)
IQ_Base_Savage$A2 = toupper(IQ_Base_Savage$A2)

### REQUIREMENT 4: STANDARD GWAS QC ####
# Remove SNPs with low minor allele frequency (MAF/EAF_HRC) or imputation 
# information scores (INFO), with thresholds of .01 and .80 respectively, in 
# order to reduce the likelihood of false positives. 
IQ_Base_Savage = IQ_Base_Savage[IQ_Base_Savage$EAF_HRC>.01]
IQ_Base_Savage = IQ_Base_Savage[IQ_Base_Savage$minINFO>.80]

# Now check for any duplicated SNPs by counting the number of unique SNP IDs, 
# and comparing them with the size of the data frame overall.
if (length(unique(IQ_Base_Savage$SNP))== dim(IQ_Base_Savage)[1]) {
  print("No duplicate SNPs!")
} else {
  print("There are duplicate SNPs!")
  IQ_Base_Savage[!duplicated(IQ_Base_Savage$SNP),]
  print("Duplicated SNPs have been removed!")
}
### REQUIREMENT 5: AMBIGUOUS SNPs ####
# Now check for ambiguous SNPs i.e. those where nucleotides are in the incorrect
# pairs on the incorrect strands - we should always expect A/T and C/G 
# combinations on the effect and non-effect alleles, with only a single 
# nucleotide. Note, once these are removed, they will not be indexed in the 
# target data set, so this step only needs to be carried out once.
IQ_Base_Savage = IQ_Base_Savage[!(IQ_Base_Savage$A2=="A" && IQ_Base_Savage$A1=="T") |
                                      !(IQ_Base_Savage$A2=="T" && IQ_Base_Savage$A1=="A") |
                                      !(IQ_Base_Savage$A2=="C" && IQ_Base_Savage$A1=="G") |
                                      !(IQ_Base_Savage$A2=="G" && IQ_Base_Savage$A1=="C")]


# Remove rows where there is more than 1 nucleotide!
IQ_Base_Savage = IQ_Base_Savage[nchar(A2)==1 & nchar(A1)==1]

# It is unclear whether Savage and colleagues (2018) conducted a sex check on 
# their data. Save the above data frame!  
IQ_Base_Savage_QCed = IQ_Base_Savage

# Now, calculate beta and se beta for the summary data
IQ_Base_Savage_QCed$'beta' = IQ_Base_Savage_QCed$Zscore / sqrt(2 * IQ_Base_Savage_QCed$EAF_HRC * (1 - IQ_Base_Savage_QCed$EAF_HRC) * (IQ_Base_Savage_QCed$N_analyzed + (IQ_Base_Savage_QCed$Zscore)^2))
IQ_Base_Savage_QCed$'beta_se' = 1 / sqrt(2 * IQ_Base_Savage_QCed$EAF_HRC * (1 - IQ_Base_Savage_QCed$EAF_HRC) * (IQ_Base_Savage_QCed$N_analyzed + (IQ_Base_Savage_QCed$Zscore)^2))

# Ensure that the formatting matches that of the tutorial by Choi and colleagues (2020).
IQ_Base_Savage_QCed = IQ_Base_Savage_QCed[, c("CHR", "POS", "SNP", "A1", "A2", "UNIQUE_ID", "N_analyzed", "Zscore", "SE", "P", "stdBeta","minINFO","EAF_HRC","EffectDirection","beta","beta_se")]
colnames(IQ_Base_Savage_QCed)[which(names(IQ_Base_Savage_QCed) == "POS")] <- "BP"
# Now save the QC-ed table!
fwrite(IQ_Base_Savage_QCed,"Savage_IQ_2018_sumstats/IQ_Base_Savage_QCed.txt",sep="\t", row.names = F)

### PART 3 - TARGET (ABCD) DATA QUALITY CONTROL #### 
# We shall be using SNP data from ABCD Release 3, and filtering participants to
# include only those in Release 2. We shall conduct the same QS as above. The 
# first step is to remove SNPs with genotyping rates, low MAF frequency (<.01),
# out of Hardy-Weinberg equilibrium (with a p-threshold of 1e-06), and removing
# individuals with low genotyping rates. 

# Remote connect to a CBU cluster node
session = ssh_connect("am10@login-f01")

### REQUIREMENT 1 - STANDARD GWAS QC #### 
# Now run the PLINK command to control for minor allele frequency (MAF), Hardy-
# Weinberg equilibrium (HWE), gene information etc. The first line of the command
# adds plink to the path, and the second line changes the directory to that 
# containing the PLINK files. The third line runs PLINK to control for MAF, HWE,
# GENO, and MIND. The fourth line extracts SNPs which survived QC, and sorts them
# across 50 SNP intervals, retaining those with linkage disequilibrium (LD) scores
# of less than .25. The final line extracts the pruned QC data and calculates
# heterozygosity for each SNP. 

setwd("//cbsu/data/imaging/projects/external/abcd/analyses/Genomics")

ABCD_Target_QC_PLINK_Commands = ssh_exec_wait(session, command = "module load plink;
                    cd /imaging/projects/external/abcd/analyses/Genomics; 
                    plink --bfile ABCD_Updated --maf 0.01 --hwe 1e-6 --geno 0.01 --mind 0.01 --write-snplist --not-chr XY X Y --make-just-fam --out ABCD_QC_AM;
                    plink --bfile ABCD_Updated --keep ABCD_QC_AM.fam --extract ABCD_QC_AM.snplist --indep-pairwise 200 50 0.25 --not-chr XY X Y --out ABCD_QC_AM;
                    plink --bfile ABCD_Updated --extract ABCD_QC_AM.prune.in --keep ABCD_QC_AM.fam --het --not-chr XY X Y --out ABCD_QC_AM")

# Now remove SNPs with high heterozygosity F coefficients from the .het file. 
HET_data_file = fread("ABCD_QC_AM.het")
# Get samples with F coefficients within 3 SD of the population mean
Valid_HET = HET_data_file[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)]
# Now write the FID and IID for all valid samples as a separate file
fwrite(Valid_HET[,c("FID","IID")],"ABCD_QC_AM_Valid_Sample",sep = "\t")

### REQUIREMENT 2 - Mismatching SNPs #### 
# We need to ensure that alleles are complementary in the target and base data. 
# Therefore, we shall check for any mismatching SNPs, and strand-flip the relevant
# alleles. 
### PART A - Reading In Necessary Data ####
# Again, ensure that these processes are conducted separately for the European
# and Non-European samples!
# Read in the bim file and pipe using magrittr
bim <- fread("ABCD_Updated.bim") %>% 
  # Note: . represents the output from previous step
  # The syntax here means, set names of the data read from the bim file, and 
  # replace the original column names by the new names
  setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
  # And immediately change the alleles to upper cases
  .[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]

# Read in summary statistic data (require data.table v1.12.0+)
IQ_Base_Savage_QCed = fread('/imaging/astle/am10/ABCD/Data/Savage_IQ_2018_sumstats/IQ_Base_Savage_QCed.txt') %>%
  # And immediately change the alleles to upper cases
  .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]

# Read in QCed SNPs
qc <- fread("ABCD_QC_AM.snplist", header=F)

### PART B - Identify Strands Requiring Flipping ####
# Merge QCed summary statistic with target, but make sure column names are 
# compatible!
IQ_Base_Savage_QCed$CHR = as.character(IQ_Base_Savage_QCed$CHR)
bim$CHR = as.character(bim$CHR)
info <- merge(bim, IQ_Base_Savage_QCed, by=c("SNP", "CHR", "BP")) %>%
  # And filter out QCed SNPs
  .[SNP %in% qc[,V1]]

# Function for calculating the complementary allele
complement <- function(x){
  switch (x,
          "A" = "T",
          "C" = "G",
          "T" = "A",
          "G" = "C",
          return(NA)
  )
} 
# Get SNPs that have the same alleles across base and target
info.match <- info[A1 == B.A1 & A2 == B.A2, SNP]
# Identify SNPs that are complementary between base and target
com.snps <- info[sapply(B.A1, complement) == A1 & sapply(B.A2, complement) == A2, SNP]
# Now update the bim file
bim[SNP %in% com.snps, c("B.A1", "B.A2") := list(sapply(B.A1, complement), sapply(B.A2, complement))]

### PART C - Identify SNPs Requiring Re-coding ####
recode.snps <- info[B.A1==A2 & B.A2==A1, SNP]
# Update the bim file
bim[SNP %in% recode.snps, c("B.A1", "B.A2") := list(B.A2, B.A1)]
# Identify SNPs that need recoding & complement
com.recode <- info[sapply(B.A1, complement) == A2 & sapply(B.A2, complement) == A1, SNP]
# Now update the bim file
bim[SNP %in% com.recode, c("B.A1", "B.A2") := list(sapply(B.A2, complement), sapply(B.A1, complement))]
# Write the updated bim file
fwrite(bim[,c("SNP", "B.A1")], "ABCD_QC_AM.A1", col.names=F, sep="\t")  
### PART D - Identify Mismatching SNPs ####
mismatch <- bim[!(SNP %in% info.match|SNP %in% com.snps|SNP %in% recode.snps|
                    SNP %in% com.recode),SNP]
write.table(mismatch, "ABCD_QC_AM.mismatch", quote=F, row.names=F, col.names=F)

### REQUIREMENT 3 - DUPLICATE SNPs #### 
# As we conducted for the base data, we shall remove duplicate SNPs. 
if (length(unique(bim$SNP))== dim(bim)[1]) {
  print("No duplicate SNPs!")
} else {
  print("There are duplicate SNPs!")
}

### REQUIREMENT 4 - SEX CHECK #### 
# A difference between biological and reported sex could be the result of sample
# mislabeling, and hence mislabeled samples will be removed. However, Release 
# Notes for ABCD Release 4 noted that sex had been removed from the PLINK files,
# and that researchers should use phenotypic data instead. Therefore, we shall
# include all participant, and update the required tables/data frames. In this
# case, we shall only include the XY chromosomes. Therefore, we will need to 
# extract the sex chromosomes for the included participants and use this as the
# input for the sex check!
ABCD_Target_Sex_Check_PLINK_Commands = ssh_exec_wait(session, command = "module load plink;
                                                    cd /imaging/projects/external/abcd/analyses/Genomics;
                                                    plink --bfile ABCD_Updated --chr X Y XY --check-sex;
                                                    plink --bfile ABCD_Updated --not-chr X Y XY --extract ABCD_QC_AM.prune.in --keep ABCD_QC_AM_Valid_Sample --out ABCD_QC_AM")

# Now update the sex information for the samples!
valid <- read.table("ABCD_QC_AM_Valid_Sample", header=F)
colnames(valid) = c("FID","IID")
dat <- fread("ABCD_QC_AM.sexcheck")[FID%in%valid$FID]
write.table(dat[STATUS=="OK",c("FID", "IID")], "ABCD_QC_AM_Valid", row.names=F, col.names=F, sep="\t", quote=F)

### REQUIREMENT 5 - SAMPLE OVERLAP ####
# A high degree of overlap in characteristics between the base and target data
# can result in inflated PRS values. We believe this issue has been minimised,
# as Savage and colleagues '18 analysed data from both children and adults, thus 
# in contrast to the US ABCD Study. 
### REQUIREMENT 6 - RELATEDNESS #### 
# Highly related individuals, such as first or second-degree relatives, can 
# also inflate PRS values, and thus such individuals must be removed. Note that 
# we are using PLINK1.9 to conduct this analysis 
ABCD_Target_Relatedness_PLINK_Commands = 
  ssh_exec_wait(session, command = "module load plink;
                cd /imaging/projects/external/abcd/analyses/Alicja/Genomics;
                plink --bfile ABCD_Updated --extract ABCD_QC_AM.prune.in --keep ABCD_QC_AM_Valid --rel-cutoff 0.125 --not-chr XY X Y --out ABCD_QC_AM")

# Create the final QCed ABCD file! #
ABCD_Target_Final_QC_PLINK_Commands = 
  ssh_exec_wait(session, command = "module load plink;
                cd /imaging/projects/external/abcd/analyses/Alicja/Genomics;
                plink --bfile ABCD_Updated --make-bed --keep ABCD_QC_AM.rel.id --out ABCD_QC_AM --extract ABCD_QC_AM.snplist --exclude ABCD_QC_AM.mismatch --not-chr XY X Y --a1-allele ABCD_QC_AM.A1")

# Load and inspect the final QCed ABCD file.
Final_QCed_ABCD_File_BIM = fread('/imaging/projects/external/abcd/analyses/Alicja/Genomics/ABCD_QC_AM.bim')
Final_QCed_ABCD_File_FAM = fread('/imaging/projects/external/abcd/analyses/Alicja/Genomics/ABCD_QC_AM.fam')


# Note, before starting our PRS analysis, it is essential to ensure that the 
# relevant files are available and that they are complete. We should have bed, 
# bim, and fam files, alongside covariate and phenotype files for ABCD. Our base
# data will be the QCed version from Savage and colleagues (2018). 

### PART 1D - Identify European and Non-European Cohorts! ####
# Since the discovery data set is European, we need to ensure that the ancestry
# characteristics of the target data set matches, due to differing populations
# being a major confound in PRSs. Therefore, we shall divide the participants 
# into those of European (white) heritage, and Non-European (remaining) heritage.
# To do so, we ran fastStructure with 4 ancestry groups, and shall visualise 
# the results here...
fs_run_K1 = fread("/imaging/astle/am10/ABCD/Analysis/FastStructure/fs_run_K.1.meanQ")
fs_run_K2 = fread("/imaging/astle/am10/ABCD/Analysis/FastStructure/fs_run_K.2.meanQ")
fs_run_K3 = fread("/imaging/astle/am10/ABCD/Analysis/FastStructure/fs_run_K.3.meanQ")
fs_run_K4 = fread("/imaging/astle/am10/ABCD/Analysis/FastStructure/fs_run_K.4.meanQ")

# To simplify analyses, we shall choose K = 2, suggesting European vs Non-European,
# with a threshold of .9 in cluster 2 denoting European ancestry
ABCD_FAM_fastStructure = fread("/imaging/astle/am10/ABCD/Analysis/FastStructure/ABCD_FAM_Subject_List.txt", header=F)
ABCD_fastStructure_K2 = cbind(ABCD_FAM_fastStructure, fs_run_K2)
colnames(ABCD_fastStructure_K2) = c("SUBJECTKEY","CLUSTER_1","CLUSTER_2")
# Following Loughnan and colleagues (2020), we shall select European participants
# as those with ancestry cluster membership scores higher than 90%
European_Participants = dplyr::filter(ABCD_fastStructure_K2, CLUSTER_2 > .9)
European_Participants = merge(European_Participants, ABCD_FAM_File_with_Demographics, by = c("SUBJECTKEY"))
European_Participants = dplyr::select(European_Participants, V1, SUBJECTKEY)
write.table(European_Participants, file = "European_Participants.txt", quote = F, col.names = F, row.names = F, sep = "\t")

Non_European_Participants = as.data.frame(anti_join(ABCD_fastStructure_K2, European_Participants, by = "SUBJECTKEY"))
Non_European_Participants = merge(Non_European_Participants, ABCD_FAM_File_with_Demographics, by = c("SUBJECTKEY"))
Non_European_Participants = dplyr::select(Non_European_Participants, V1, SUBJECTKEY)
write.table(Non_European_Participants, file = "Non_European_Participants.txt", quote = F, col.names = F, row.names = F, sep = "\t")

# Now create the European and Non-European participant data frame from the QC'ed 
# version above!
Ancestry_Stratification_PLINK_Commands = 
  ssh_exec_wait(session, command = "module load plink;
                cd /imaging/projects/external/abcd/analyses/Genomics;
                plink --bfile ABCD_QC_AM --keep European_Participants.txt --make-bed --out ABCD_QC_AM_European;
                plink --bfile ABCD_QC_AM --keep Non_European_Participants.txt --make-bed --out ABCD_QC_AM_Non_European")

### PART 5 - PGS ANALYSIS USING PLINK! #### 
# The following section will use PLINK, PRSice-2, LDPred-2 and lassosum. Ensure
# that all of these software have been installed and added to the relevant paths
# before starting!

### STEP 1 - Clumping ####
# To control for high linkage disequilibrium, meaning high genetic correlations
# between variants across the genome, we shall exclude variants with high LD with
# variants, but retain those strongly associated with the phenotype (IQ). Then,
# extract the indices of the clumped SNPs in the 4th line of the command. Note
# that we shall repeat all steps for European and Non-European subsets.
session = ssh_connect("am10@login-g01")

ABCD_Clumping_Plink_Commands = 
  ssh_exec_wait(session, command = "module load plink; 
              cd /imaging/projects/external/abcd/analyses/Genomics; 
              plink --bfile ABCD_QC_AM_European --not-chr XY X Y --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump /imaging/astle/am10/ABCD/Data/Savage_IQ_2018_sumstats/IQ_Base_Savage_QCed.txt --clump-snp-field SNP --clump-field P --out ABCD_release_3.0_QCed_European;
              plink --bfile ABCD_QC_AM_Non_European --not-chr XY X Y --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump /imaging/astle/am10/ABCD/Data/Savage_IQ_2018_sumstats/IQ_Base_Savage_QCed.txt --clump-snp-field SNP --clump-field P --out ABCD_release_3.0_QCed_Non_European")

# Now extract the different p-values and SNPs!
ABCD_Clumping_Extraction = 
  ssh_exec_wait(session, command = "cd /imaging/projects/external/abcd/analyses/Genomics;
                awk 'NR!=1{print $3}' ABCD_release_3.0_QCed_European.clumped >  ABCD.European.valid.snp;
                awk 'NR!=1{print $3}' ABCD_release_3.0_QCed_Non_European.clumped >  ABCD.Non_European.valid.snp;
                awk '{print $3,$10}' /imaging/astle/am10/ABCD/Data/Savage_IQ_2018_sumstats/IQ_Base_Savage_QCed > SNP.pvalue;
                echo 0.001 0 0.001 > range_list;
                echo 0.05 0 0.05 >> range_list;
                echo 0.1 0 0.1 >> range_list;
                echo 0.2 0 0.2 >> range_list;
                echo 0.3 0 0.3 >> range_list;
                echo 0.4 0 0.4 >> range_list;
                echo 0.5 0 0.5 >> range_list")

### STEP 2 - Generate PGS! ####
# Now run the actual PGS PLINK command!
ABCD_PRS_PLINK_Commands = 
  ssh_exec_wait(session, command = "module load plink; 
                cd /imaging/projects/external/abcd/analyses/Genomics;
                plink --bfile ABCD_QC_AM_European --not-chr XY X Y --score /imaging/astle/am10/ABCD/Data/Savage_IQ_2018_sumstats/IQ_Base_Savage_QCed 3 4 11 header center --q-score-range range_list SNP.pvalue --extract ABCD.European.valid.snp --out ABCD_release_3.0_QCed_European;
                plink --bfile ABCD_QC_AM_Non_European --not-chr XY X Y --score /imaging/astle/am10/ABCD/Data/Savage_IQ_2018_sumstats/IQ_Base_Savage_QCed 3 4 11 header center --q-score-range range_list SNP.pvalue --extract ABCD.Non_European.valid.snp --out ABCD_release_3.0_QCed_Non_European")

### STEP 3 - Accounting for Population Stratification! ####
# To account for population stratification, which acts as a major confounder, we
# shall prune our data, and then extract 6 principal components to represent 
# different ethnicities.
ABCD_Population_Stratification_PLINK_Commands = 
  ssh_exec_wait(session, command = "module load plink; 
                cd /imaging/projects/external/abcd/analyses/Alicja/Genomics;
                plink --bfile ABCD_QC_AM_European --not-chr XY X Y --indep-pairwise 200 50 0.25 --out ABCD_release_3.0_QCed_European;
                plink --bfile ABCD_QC_AM_European --not-chr XY X Y --extract ABCD_QC_AM.prune.in --pca 6 --out ABCD_release_3.0_QCed_European;
                plink --bfile ABCD_QC_AM_Non_European --not-chr XY X Y --indep-pairwise 200 50 0.25 --out ABCD_release_3.0_QCed_Non_European;
                plink --bfile ABCD_QC_AM_Non_European --not-chr XY X Y --extract ABCD_QC_AM.prune.in --pca 6 --out ABCD_release_3.0_QCed_Non_European")

# We shall now find which p-value threshold provides the best fit for our PRS!
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
phenotype <- fread("/imaging/astle/am10/ABCD/Data/ABCD_Genomics/ABCD_Pheno_File.txt", header = F)
colnames(phenotype) = c("FID","IID","Pheno")
pcs_European <- fread("/imaging/projects/external/abcd/analyses/Alicja/Genomics/ABCD_release_3.0_QCed_European.eigenvec", header=F) %>%
  setnames(., colnames(.), c("FID", "IID", paste0("PC",1:6)) )
pcs_Non_European <- fread("/imaging/projects/external/abcd/analyses/Alicja/Genomics/ABCD_release_3.0_QCed_Non_European.eigenvec", header=F) %>%
  setnames(., colnames(.), c("FID", "IID", paste0("PC",1:6)) )

# Include sex as an additional covariate
Sex_Covariate = fread("/imaging/astle/am10/ABCD/Data/ABCD_Genomics/ABCD_Sex_File.txt", header = F)
colnames(Sex_Covariate) = c("FID","IID","Sex")

pheno_European <- merge(phenotype, Sex_Covariate) %>%
  merge(., pcs_European)
null.r2.European <- summary(lm(Pheno~., data=pheno_European[,-c("FID", "IID")]))$r.squared
prs.result_European <- NULL

pheno_Non_European <- merge(phenotype, Sex_Covariate) %>%
  merge(., pcs_Non_European)
null.r2.Non_European <- summary(lm(Pheno~., data=pheno_Non_European[,-c("FID", "IID")]))$r.squared
prs.result_Non_European <- NULL

# Conduct stratification/population adjustment for European subset and run PRS
for(i in p.threshold){
  pheno_European.prs <- paste0("/imaging/projects/external/abcd/analyses/Alicja/Genomics/ABCD_release_3.0_QCed_European.", i, ".profile") %>%
    fread(.) %>%
    .[,c("FID", "IID", "SCORE")] %>%
    merge(., pheno_European, by=c("FID", "IID"))
  
  model_European <- lm(Pheno~., data=pheno_European.prs[,-c("FID","IID")]) %>%
    summary
  model_European.r2 <- model_European$r.squared
  prs.r2_European <- model_European.r2-null.r2.European
  prs_European.coef <- model_European$coeff["SCORE",]
  prs.result_European %<>% rbind(.,
                        data.frame(Threshold=i, R2=prs.r2_European, 
                                   P=as.numeric(prs_European.coef[4]), 
                                   BETA=as.numeric(prs_European.coef[1]),
                                   SE=as.numeric(prs_European.coef[2])))
}
print(prs.result_European[which.max(prs.result_European$R2),])
# A threshold of p = .10 for inclusion in the clusters provides the best-fit PRS
# for the European subset.

# Now repeat the stratification/population adjustment for the Non-European subset
for(i in p.threshold){
  pheno_Non_European.prs <- paste0("/imaging/projects/external/abcd/analyses/Alicja/Genomics/ABCD_release_3.0_QCed_Non_European.", i, ".profile") %>%
    fread(.) %>%
    .[,c("FID", "IID", "SCORE")] %>%
    merge(., pheno_Non_European, by=c("FID", "IID"))
  
  model_Non_European <- lm(Pheno~., data=pheno_Non_European.prs[,-c("FID","IID")]) %>%
    summary
  model_Non_European.r2 <- model_Non_European$r.squared
  prs.r2_Non_European <- model_Non_European.r2-null.r2.Non_European
  prs_Non_European.coef <- model_Non_European$coeff["SCORE",]
  prs.result_Non_European %<>% rbind(.,
                                 data.frame(Threshold=i, R2=prs.r2_Non_European, 
                                            P=as.numeric(prs_Non_European.coef[4]), 
                                            BETA=as.numeric(prs_Non_European.coef[1]),
                                            SE=as.numeric(prs_Non_European.coef[2])))
}
print(prs.result_Non_European[which.max(prs.result_Non_European$R2),])
# A threshold of p = .2 for inclusion in the clusters provides the best-fit PRS
# for the Non-European subset.

### STEP 4 - Visualize the PGS! - SUPPLEMENTARY FIGURE 1 ####
### Best-Fitting PGS for European and Non-European Participants ###
setwd("//cbsu/data/Imaging/projects/external/abcd/analyses/Alicja/Genomics")

# Create a tibble which holds the results of the polygenic scores for the 
# European and non-European participants, with an ancestry factor.
thresholding_plot_df = tibble(Threshold = factor(rep(p.threshold,2)),
                              R2 = c(prs.result_European$R2, prs.result_Non_European$R2),
                              p = c(prs.result_European$P, prs.result_Non_European$P),
                              Ancestry = factor(rep(c("European", "Non-European"),each=length(p.threshold))))

# Create a bar plot where the thresholds are the x-axis, and the proportion of 
# variance in cognitive ability accounted for is the y-axis. We will facet the 
# plot by ancestry.
ggplot(thresholding_plot_df, aes(x = Threshold, y = R2)) +
  geom_bar(stat = "identity", aes(fill=R2)) + facet_wrap(~Ancestry, scales = "free") +
  xlab(expression(paste("Polygenic Clumping Threshold")~italic(P)[T])) +
  ylab(expression(paste("Polygenic Score Model Fit ",R^2))) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 15, face = "bold", hjust = .5),
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold", hjust = .5),
        legend.position = "bottom") +
  scale_fill_gradient(low = "red", high = "green")



# Plot first for European participants! These plots will be Supplementary
prs.result_European$print.p = round(prs.result_European$P, digits = 3)
prs.result_European$print.p[!is.na(prs.result_European$print.p) &
                     prs.result_European$print.p == 0] <-
  format(prs.result_European$P[!is.na(prs.result_European$print.p) &
                        prs.result_European$print.p == 0], digits = 2)
prs.result_European$print.p <- sub("e", "*x*10^", prs.result_European$print.p)


# Initialize ggplot, requiring the threshold as the x-axis (use factor so that it is uniformly distributed)
IQ_European_PRS_Bar = ggplot(data = prs.result_European, aes(x = factor(Threshold), y = R2)) +
  # Specify that we want to print p-value on top of the bars
  geom_text(
    aes(label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 8,
    parse = T)  +
  # Specify the range of the plot, *1.25 to provide enough space for the p-values
  scale_y_continuous(limits = c(0, max(prs.result_European$R2) * 1.25)) +
  # Specify the axis labels
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PS model fit:  ", R ^ 2))) +
  ggtitle("European Subset") +
  # Draw a bar plot
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  # Specify the colors
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-4,
    name = bquote(atop(-log[10] ~ model, italic(P) - value),)
  ) +
  # Some beautification of the plot
  theme_classic() +  theme(
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(size = 25),
    legend.title = element_text(face = "bold", size = 25),
    legend.text = element_text(size = 25),
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(size = 25, face="bold",hjust=0.5))

print(IQ_European_PRS_Bar)

# save the plot
ggsave("PS_IQ_European_bar.png", height = 12, width = 16)

# Now do the same for the Non-European participants
prs.result_Non_European$print.p = round(prs.result_Non_European$P, digits = 3)
prs.result_Non_European$print.p[!is.na(prs.result_Non_European$print.p) &
                              prs.result_Non_European$print.p == 0] <-
  format(prs.result_Non_European$P[!is.na(prs.result_Non_European$print.p) &
                                 prs.result_Non_European$print.p == 0], digits = 2)
prs.result_Non_European$print.p <- sub("e", "*x*10^", prs.result_Non_European$print.p)

# Initialize ggplot, requiring the threshold as the x-axis (use factor so that it is uniformly distributed)
IQ_Non_European_PRS_Bar = ggplot(data = prs.result_Non_European, aes(x = factor(Threshold), y = R2)) +
  # Specify that we want to print p-value on top of the bars
  geom_text(
    aes(label = paste(print.p)),
    vjust = -1.5,
    hjust = 0,
    angle = 45,
    cex = 8,
    parse = T
  )  +
  # Specify the range of the plot, *1.25 to provide enough space for the p-values
  scale_y_continuous(limits = c(0, max(prs.result_Non_European$R2) * 1.25)) +
  # Specify the axis labels
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PS model fit:  ", R ^ 2))) +
  ggtitle("Non-European Subset") +
  # Draw a bar plot
  geom_bar(aes(fill = -log10(P)), stat = "identity") +
  # Specify the colors
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-4,
    name = bquote(atop(-log[10] ~ model, italic(P) - value),)
  ) +
  # Some beautification of the plot
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 25),
    axis.text = element_text(size = 25),
    legend.title = element_text(face = "bold", size =
                                  25),
    legend.text = element_text(size = 25),
    axis.text.x = element_text(angle = 45, hjust =1, size = 25),
    plot.title = element_text(size=25, face="bold", hjust=0.5)) 
# save the plot
ggsave("PS_IQ_Non_European_bar.png", height = 12, width = 16)

# Now combine both plots onto a single plot and save
pdf(file = "Combined_IQ_PRS_P_Value_Thresholds_Plot.pdf", width = 16, height = 13)
Combined_IQ_PRS_P_Value_Thresholds_Plot = ggarrange(IQ_European_PRS_Bar, IQ_Non_European_PRS_Bar, ncol=2, nrow=1,common.legend = TRUE, legend = "bottom")
dev.off()
ggsave("Combined_IQ_PRS_P_Value_Thresholds_Plot.pdf", width = 15, height = 12, dpi = 700)

### STEP 5 - Select PRS participants with neuro imaging data ####
# Load the stratified participant list
Stratified_PP_List = fread('//cbsu/data/imaging/astle/am10/ABCD/nda-abcd-s3-downloader-master/subject_list_ABCD_stratified_updated_March_2022.txt', header = F)
# Now load up the participant IDs of those included in the PRS
PRS_Participants = c(pcs_European$IID, pcs_Non_European$IID)
# Reformat the participant IDs
Stratified_PP_List$V1 = gsub("sub-","",Stratified_PP_List$V1)
Stratified_PP_List$V1 = sub("^(.{4})", "\\1_", Stratified_PP_List$V1)
colnames(Stratified_PP_List) = "IID"
# Find common participants
Neuroimaging_plus_PRS = intersect(Stratified_PP_List$IID, PRS_Participants)
### STEP 6 - Generate ranked gene lists for each participant ####
# For each participant, we can find the most influential SNP by finding the ABCD
# SNP with the largest effect size in the base data set. Therefore, first load
# the valid SNPs for Europeans and Non-Europeans
setwd("//cbsu/data/imaging/projects/external/abcd/analyses/Alicja/Genomics")
European_SNPs = fread('ABCD.European.valid.snp',header=F)
Non_European_SNPs = fread('ABCD.Non_European.valid.snp',header=F)

# Now load the effect sizes from the target data set
IQ_Base_Savage_QCed = fread('//cbsu/data/imaging/astle/am10/ABCD/Data/Savage_IQ_2018_sumstats/IQ_Base_Savage_QCed.txt')

# Create a new data frame, where the base SNPs are ranked by their effect sizes.
# We shall rank in a descending order, so that SNPs with the largest positive 
# beta values are at the top of the list.
Ranked_SNPs_Base_Data = dplyr::select(IQ_Base_Savage_QCed, SNP, beta)
Ranked_SNPs_Base_Data = Ranked_SNPs_Base_Data[order(-beta),]

# Now find the effect sizes for the SNPs included in the PRS, across both 
# European and Non-European participants
ABCD_European_Ranked_SNPs = Ranked_SNPs_Base_Data[Ranked_SNPs_Base_Data$SNP %in% European_SNPs$V1,]
ABCD_Non_European_Ranked_SNPs = Ranked_SNPs_Base_Data[Ranked_SNPs_Base_Data$SNP %in% Non_European_SNPs$V1,]

# Combine together, select only unique entries, and save!
ABCD_Ranked_SNPs = rbind(ABCD_European_Ranked_SNPs, ABCD_Non_European_Ranked_SNPs)
ABCD_Ranked_SNPs = ABCD_Ranked_SNPs[!duplicated(ABCD_Ranked_SNPs[,1]),]
ABCD_Ranked_SNPs$'Rank'[order(-ABCD_Ranked_SNPs$beta)] = 1:nrow(ABCD_Ranked_SNPs)
write.table(ABCD_Ranked_SNPs, file = "ABCD_Ranked_SNPS.txt", quote = F, col.names = T, row.names = F)

