# This script details selection of participants, using data from DCAN Labs. Our 
# inclusion criteria are that all participants must have good-quality 
# neuroimaging data in both DTI and rs-fMRI modalities (PART 2), from the second
# ABCD data release, as well as genomic (PART 3), and cognitive (PART 4) data. 
# Note that BIDS data from DCAN labs has passed quality control by the DCAN team. 
# Further, all included participants must be singletons (PART 5A) and have at 
# least 90% European ancestry (PART 5B), in order to prevent biases in the PRS.
# Further, we shall impute (PART 5C) missing demographic information, such as 
# race, sex, and 2 measures of socioeconomic status (parental education, and 
# parental income). From this, we shall randomly select 2000 participants (PART 
# 6A), and compare their distributions (PART 6C) of race, sex, and SES to the 
# entire sample of ABCD release 2.0 (PART 6B), including singletons, twins, and 
# siblings. 
### PART 1 - Set Up the Work Space and Load Neuroimaging Data #### 
rm(list=ls())
library(data.table)
library(dplyr)
library(tidyr)
library(rockchalk)
library(ggplot2)
library(ggExtra)
library(ggridges)
library(stringr)
library(cowplot)
library(mice)
library(rsample)
library(splitstackshape)
library(esc)
library(Matching)
library(ggfortify)
library(FactoMineR)
library(factoextra)
library(psych)
library(grid)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(fmsb)
library(viridis)
library(mice)

setwd("//cbsu/data/Imaging/astle/am10/ABCD")

# Now load the DCAN Labs Data Manifest File, which contains all S3 links for all
# participants with good-quality neuroimaging data, as deemed by the ABCD DAIC
# Team (Data Analysis and Informatics Core). 
DCAN_Manifest = fread('nda-abcd-s3-downloader-master/datastructure_manifest.txt')
DCAN_Manifest = DCAN_Manifest[-c(1)]

### PART 2 - Sorting Based on DTI, T1w, and rsfMRI Data ####
# Note that some participants have 2 DWI runs, which together are usually a 
# quarter of the size of a single DWI run file. Therefore, we shall first select
# those participants with only a single DWI run...
DCAN_Manifest$'PARTICIPANT_ID' = substr(DCAN_Manifest$manifest_name, 1, 19)
# Now find participants with 2 runs
DCAN_DWI_2_Runs = unique(DCAN_Manifest$PARTICIPANT_ID[grep("run-02_dwi", DCAN_Manifest$associated_file)])
sprintf('There are %0.f participants with 2 DWI runs...', length(DCAN_DWI_2_Runs))
# Find the DWI participants
DCAN_DWI = unique(DCAN_Manifest$PARTICIPANT_ID[grep("inputs.dwi.dwi", DCAN_Manifest$manifest_name)])
# And remove those with 2 runs
DCAN_DWI = setdiff(DCAN_DWI,DCAN_DWI_2_Runs)
# Now extract participant IDs who have the following attributes: fmap, dwi, T1w, 
# and rs-fMRI...
DCAN_T1w = unique(DCAN_Manifest$PARTICIPANT_ID[grep("inputs.anat.T1w", DCAN_Manifest$manifest_name)])
DCAN_rsfMRI = unique(DCAN_Manifest$PARTICIPANT_ID[grep("inputs.func.task-rest", DCAN_Manifest$manifest_name)])
DCAN_FMAP = unique(DCAN_Manifest$PARTICIPANT_ID[grep("inputs.fmap", DCAN_Manifest$manifest_name)])
# Now intersect these to find participants with common good-quality data across
# all of the modalities listed above
DCAN_Common_Participants = Reduce(intersect, list(DCAN_FMAP, DCAN_DWI, DCAN_T1w, DCAN_rsfMRI))
# Update format of participant IDs i.e. remove "sub-", and add underscore
DCAN_Common_Participants = str_sub(DCAN_Common_Participants, 5)
DCAN_Common_Participants = sub("^(.{4})", "\\1_", DCAN_Common_Participants)
# Update user about number of participants left
sprintf("%.0f participants retained!", length(DCAN_Common_Participants))
# Now format the above data frames!
DCAN_DWI = data.frame(DCAN_DWI)
DCAN_T1w = data.frame(DCAN_T1w)
DCAN_FMAP = data.frame(DCAN_FMAP)
DCAN_rsfMRI = data.frame(DCAN_rsfMRI)
colnames(DCAN_DWI) = c("SUBJECTKEY")
colnames(DCAN_T1w) = c("SUBJECTKEY")
colnames(DCAN_FMAP) = c("SUBJECTKEY")
colnames(DCAN_rsfMRI) = c("SUBJECTKEY")

DCAN_DWI_T1w = merge(DCAN_DWI, DCAN_T1w, by = c("SUBJECTKEY"))
DCAN_DWI_T1w$SUBJECTKEY = str_sub(DCAN_DWI_T1w$SUBJECTKEY, 5)
DCAN_DWI_T1w$SUBJECTKEY = sub("^(.{4})", "\\1_", DCAN_DWI_T1w$SUBJECTKEY)

DCAN_rsfMRI_FMAP = merge(DCAN_rsfMRI, DCAN_FMAP, by = c("SUBJECTKEY"))
DCAN_rsfMRI_FMAP$SUBJECTKEY = str_sub(DCAN_rsfMRI_FMAP$SUBJECTKEY, 5)
DCAN_rsfMRI_FMAP$SUBJECTKEY = sub("^(.{4})", "\\1_", DCAN_rsfMRI_FMAP$SUBJECTKEY)

### PART 3 - Sorting Based on Genomic Data #### 
# Load the genomic data batch info!
ABCD_Genomics_Release3_Batch_Info = fread("Data/ABCD_Genomics/ABCD_release3.0_.batch_info.txt")
colnames(ABCD_Genomics_Release3_Batch_Info)[2] = "SUBJECTKEY"
# As specified in the release notes, plate 461 proved to be especially 
# problematic in quality control, and hence will be removed from our analyses. 
ABCD_Genomics_Release3_Batch_Info = ABCD_Genomics_Release3_Batch_Info[!(ABCD_Genomics_Release3_Batch_Info$Axiom_Plate==461),]
# Now merge with the common_participants data frame
DCAN_Common_Participants = as.data.frame(DCAN_Common_Participants)
colnames(DCAN_Common_Participants) = c("SUBJECTKEY")
Common_Participants_Dataframe = DCAN_Common_Participants$SUBJECTKEY[DCAN_Common_Participants$SUBJECTKEY %in% ABCD_Genomics_Release3_Batch_Info$SUBJECTKEY]
Common_Participants_Dataframe = as.data.frame(Common_Participants_Dataframe)
colnames(Common_Participants_Dataframe)[1] = "SUBJECTKEY"
# Update user about number of participants left
sprintf("%.0f participants retained!", length(Common_Participants_Dataframe$SUBJECTKEY))
# Create additional DCAN data frame
DCAN_Common_Participants = data.frame(DCAN_Common_Participants)
colnames(DCAN_Common_Participants) = "SUBJECTKEY"
DCAN_Common_Participants_with_Genomics = merge(DCAN_Common_Participants, ABCD_Genomics_Release3_Batch_Info, by = "SUBJECTKEY")
DCAN_DWI_T1w_with_Genomics = merge(DCAN_DWI_T1w, ABCD_Genomics_Release3_Batch_Info, by = "SUBJECTKEY")
DCAN_rsfMRI_FMAP_with_Genomics = merge(DCAN_rsfMRI_FMAP, ABCD_Genomics_Release3_Batch_Info, by = "SUBJECTKEY")

### PART 4 - Cognitive Scores #### 
# For our genetic analyses, the phenotypic outcome will be intelligence, as 
# measured by the NIH Toolbox. We chose this measure as it is a broad cognitive
# measure, spanning executive functioning, episodic memory, language, processing
# speed, working memory, and attention. To ensure that our results are not 
# restricted to any one test, we shall conduct Principal component Analysis (PCA)
# on all age-corrected sub-tests in the Toolbox, and extract the first principal
# component (PC) as our IQ/outcome measure. 

# Read in the required data
ABCD_NIH_TB = fread('Data/ABCD_Tabulated_Behavioural_and_NonImaging_Release4/Package_1196999/abcd_tbss01.txt')
ABCD_NIH_TB = ABCD_NIH_TB[-c(1),]
# Select participants corresponding to baseline_year_1_arm_1
ABCD_NIH_TB = ABCD_NIH_TB[ABCD_NIH_TB$eventname=="baseline_year_1_arm_1"]
ABCD_NIH_TB = as.data.frame(ABCD_NIH_TB)
# Now select all of the age-corrected cognitive sub-scale scores
Age_Corrected_Cognitive_Scores = ABCD_NIH_TB[, grepl( "agecorrected" , colnames(ABCD_NIH_TB))]
# Remove the composite scores i.e. total composite, fluid composite, and 
# crystallized composite scores, as these will be highly correlated with the 
# sub-tests. 
Age_Corrected_Cognitive_Scores = Age_Corrected_Cognitive_Scores[, -which(colnames(Age_Corrected_Cognitive_Scores) %in% c("nihtbx_fluidcomp_agecorrected","nihtbx_cryst_agecorrected","nihtbx_totalcomp_agecorrected"))]
# Now add in the subject key
Age_Corrected_Cognitive_Scores$'SUBJECTKEY' = ABCD_NIH_TB$subjectkey
# Convert the first 7 columns to numeric
Age_Corrected_Cognitive_Scores[,1:7] = as.numeric(unlist(Age_Corrected_Cognitive_Scores[,1:7]))
# Report upon degree of missing data
sprintf('There are %s missing data points.', sum(is.na(Age_Corrected_Cognitive_Scores)))
# Impute the missing data using MICE
Age_Corrected_Cognitive_Scores_Imputed = mice(Age_Corrected_Cognitive_Scores, seed = 1234)
Age_Corrected_Cognitive_Scores = complete(Age_Corrected_Cognitive_Scores_Imputed)
# Now conduct PCA! 
Age_Corrected_Cognitive_Scores.PCA = PCA(Age_Corrected_Cognitive_Scores[,1:7], scale.unit = TRUE)
# Extract the eigenvalues
Age_Corrected_Cognitive_Scores.PCA.eigenvalues = get_eigenvalue(Age_Corrected_Cognitive_Scores.PCA)
# The first PC has an eigenvalue of 2.64, explaining 37.76% of variance in IQ
# scores, and therefore shall be retained. Only 1 PC shall be retained, as the
# second PC has an eigenvalue of 1.21, only accounting for less than half of the
# variance of the first PC (17.25%).

# Now evaluate which cognitive tests are most strongly associated with PC1
Var = get_pca_var(Age_Corrected_Cognitive_Scores.PCA)
# We find that PC1 is most strongly associated with the picture_vocabulary task,
# followed by the flanker task, listing, and card-sorting task. We can define
# the contributions of each of these tasks to the first dimension using...
Dim_Description <- dimdesc(Age_Corrected_Cognitive_Scores.PCA, axes = c(1,2), proba = 0.001)
# We find that all cognitive tests load onto the first component, with p < .001

# Now we shall extract the coordinates of the first dimension for all 
# participants - these are their PC1 scores!
Individuals = get_pca_ind(Age_Corrected_Cognitive_Scores.PCA)
Individual_Coordinates_Dim1 = Individuals$coord[,1]

# Update the data frame with these PC1 scores
Age_Corrected_Cognitive_Scores$'PC1_Scores' = Individual_Coordinates_Dim1
Age_Corrected_PC1 = dplyr::select(Age_Corrected_Cognitive_Scores, SUBJECTKEY, PC1_Scores)
Common_Participants_Dataframe = merge(Common_Participants_Dataframe, Age_Corrected_PC1, by = c("SUBJECTKEY"))

# Also update the common participants data frame using the NIHTBX_TOTALCOMP_AGE
# CORRECTED_SCORE!
Total_Comp_Age_Corrected = dplyr::select(ABCD_NIH_TB, subjectkey, nihtbx_totalcomp_agecorrected)
colnames(Total_Comp_Age_Corrected) = toupper(colnames(Total_Comp_Age_Corrected))
Common_Participants_Dataframe = merge(Common_Participants_Dataframe, Total_Comp_Age_Corrected, by = c("SUBJECTKEY"))

### PART 5A - Sorting Based on Singletons ####
# To minimize possible inflation of genetic analyses, we shall only include 
# participants who are singletons, without brothers or sisters in the sample. 

# First, load ABCD ACS Post Stratification Weights to see which participants 
# have twins or siblings also in the ABCD study!
ACS_Weights = fread('Data/ABCD_Tabulated_Behavioural_and_NonImaging_Release4/Package_1197340/acspsw03.txt')
ACS_Weights = ACS_Weights[-c(1),] 
# Select participants in the baseline release
ACS_Weights = ACS_Weights[which(ACS_Weights$eventname %like% "baseline_year_1_arm_1")]
# Now select those with a unique family ID, signifying that they are the only 
# child in the ABCD study belonging to that family. 
ACS_Weights = unique(ACS_Weights, by = c("rel_family_id"))
# Create a new dataframe just for those children without siblings or twins
Common_Participants_Unrelated = as.data.frame(Reduce(intersect, list(ACS_Weights$subjectkey, Common_Participants_Dataframe$SUBJECTKEY)))
colnames(Common_Participants_Unrelated)[1] = "SUBJECTKEY"
# Now update the common participants data frame - there are now 8643 participants
# from the common_participants_dataframe who are singletons. 
Common_Participants_Dataframe = Common_Participants_Dataframe$SUBJECTKEY[which(Common_Participants_Dataframe$SUBJECTKEY %in% Common_Participants_Unrelated$SUBJECTKEY)]
# Save updated version
write.csv(Common_Participants_Dataframe, file = "Data/ABCD_Participant_Information/ABCD_DCAN_Common_Participants_rsfMRI_DWI_FMAP_T1W_Cognition.csv")
# Update user about number of participants left
sprintf("%.0f participants retained!", length(Common_Participants_Dataframe))
# And update the DCAN data frame
colnames(ACS_Weights)[4] <- "SUBJECTKEY"
DCAN_Common_Participants_with_Genomics = merge(DCAN_Common_Participants_with_Genomics, ACS_Weights, by = "SUBJECTKEY")

### PART 5B - Sorting Based on Demographic Information ####
### Including RACE Information ####
# First, update the common participants data frame to include demographic data,
# namely race/ethnicity and highest household education. We can obtain race
# information from the Parental Demographic Dataframe. 
Parental_Demographic = fread('Data/ABCD_Tabulated_Behavioural_and_NonImaging_Release4/Package_1197342/pdem02.txt')
Parental_Demographic = Parental_Demographic[-c(1),]
# Columns 17 to 34 of the Parental_Demographic data frame correspond to different
# races, which the parent either places 0 or 1 into. 
Race_Information = Parental_Demographic[,c(4,7,17:36)]
Race_Information = as.data.frame(Race_Information)
Race_Information[,c(3:18)] = as.numeric(unlist(Race_Information[,c(3:18)]))

# Create a new column which we will place the new race labels into - the groups
# will be: White, Black, Hispanic, Asian, and Other (as per groupings described
# by Garavan and colleagues '18, in their description of ABCD recruitment). For
# the first 4 categories, ensure that each child only has a single value, else
# they will be placed in the 'All Other' category, which includes children of
# mixed racial backgrounds.

# CATEGORY 1 - WHITE
Coords_White = which(Race_Information$demo_race_a_p___10 == 1)
#Coords_White = which(rowSums(Race_Information[Coords_White,3:18])==1)
Race_Information$'Racial_Group'[Coords_White] <- "White"

# CATEGORY 2 - BLACK
Coords_Black = which(Race_Information$demo_race_a_p___11 == 1)
#Coords_Black = which(rowSums(Race_Information[Coords_Black,c(3:18)])==1)
Race_Information$Racial_Group[Coords_Black] <- "Black"

# CATEGORY 3 - HISPANIC
Coords_Hispanic = which(Race_Information$demo_ethn_v2==1)
#Coords_Hispanic = which(rowSums(Race_Information[Coords_Hispanic,c(3:18)])==1)
Race_Information$Racial_Group[Coords_Hispanic] <- "Hispanic"

# CATEGORY 4 - ASIAN
Coords_Asian = which(Race_Information$demo_race_a_p___18 == 1 | Race_Information$demo_race_a_p___19 == 1 | Race_Information$demo_race_a_p___20 == 1 
                     | Race_Information$demo_race_a_p___21 == 1 | Race_Information$demo_race_a_p___22 == 1 | Race_Information$demo_race_a_p___23 == 1
                     | Race_Information$demo_race_a_p___24 == 1 )
#Coords_Asian = which(rowSums(Race_Information[Coords_Asian,c(3:18)])==1)
Race_Information$Racial_Group[Coords_Asian] <- "Asian"

# CATEGORY 5 - ALL OTHER
# Includes children of Native Hawaiian, Pacific Islander, Alaskan Native, American Indian and multiple races.
Coords_Mixed = which(rowSums(Race_Information[,c(3:18)]) > 1,)
Coords_All_Other = which(Race_Information$demo_race_a_p___12==1 | Race_Information$demo_race_a_p___13==1|
                           Race_Information$demo_race_a_p___14==1 | Race_Information$demo_race_a_p___15==1 |
                           Race_Information$demo_race_a_p___16==1 | Race_Information$demo_race_a_p___17==1 )
Coords_All_Other = c(Coords_Mixed, Coords_All_Other)
Race_Information$Racial_Group[Coords_All_Other] <- "All Other"

# The remainder of the values are NA... Check any overlap between categories, 
# and reassign as appropriate
Potential_Overlap = Reduce(intersect, list(Coords_White, Coords_Black, Coords_All_Other, Coords_Hispanic, Coords_Asian))
# All are correctly classified as 'All Other'

# Convert to factor!
Race_Information$Racial_Group = as.factor(Race_Information$Racial_Group)

# Now select those with complete race information i.e. exclude those where the
# parents were unsure or refused to answer.
colnames(Race_Information)[1] = "SUBJECTKEY"
Common_Participants_Dataframe = as.data.frame(Common_Participants_Dataframe)
colnames(Common_Participants_Dataframe)[1] = "SUBJECTKEY"
Common_Participants_Dataframe = merge(Common_Participants_Dataframe, Race_Information, by = c("SUBJECTKEY"))

# Save updated version
write.csv(Common_Participants_Dataframe, file = "Data/ABCD_Participant_Information/ABCD_DCAN_Common_Participants_rsfMRI_DWI_FMAP_T1W_Cognition.csv")
# Update user about number of participants left
sprintf("%.0f participants retained!", length(Common_Participants_Dataframe$SUBJECTKEY))

# Update Parental_Demographic data frame!
colnames(Parental_Demographic)[4] <- "SUBJECTKEY"
Parental_Demographic = merge(Parental_Demographic, Race_Information, by = c("SUBJECTKEY"))

### Including SES Variable 1 - Parental Education ####
# Do the same with parental education and sex!
# Load parental demographic information and find factor corresponding to 
# parental education
Parental_Demographic$demo_prnt_ed_v2 = as.factor(Parental_Demographic$demo_prnt_ed_v2)

# For the parental education factor, levels 1-12 represent 1st-12th grades, 
# level 13 is High School Graduate, level 14 is GED, level 15 is Some College, 
# level 16 is Occupational Associate Degree, level 17 is Academic Associate 
# Degree, level 18 is Bachelor's degree, level 19 is a Master's Degree, level 20
# is a Professional School Degree, and level 21 is a Doctoral Degree. Level 777 
# refers to those who did not respond. Now structure the 21 levels according to 
# <12th grade (level 1), HS/GED (level 2), some college (level 3), associate's 
# degree (level 4), bachelor's degree (level 5), master's/professional degree 
# (level 6), and doctoral degree (level 7). 

# Create new column from one of the obsolete columns to encode education
colnames(Parental_Demographic)[133] <- "Parental_Education"

Coords_Less_than_12th_Grade = which(Parental_Demographic$demo_prnt_ed_v2 == 1 | Parental_Demographic$demo_prnt_ed_v2 == 2 | Parental_Demographic$demo_prnt_ed_v2 == 3
                                    | Parental_Demographic$demo_prnt_ed_v2 == 4 | Parental_Demographic$demo_prnt_ed_v2 == 5 | Parental_Demographic$demo_prnt_ed_v2 == 6 
                                    | Parental_Demographic$demo_prnt_ed_v2 == 7 | Parental_Demographic$demo_prnt_ed_v2 == 8 | Parental_Demographic$demo_prnt_ed_v2 == 9
                                    | Parental_Demographic$demo_prnt_ed_v2 == 10 | Parental_Demographic$demo_prnt_ed_v2 == 11 | Parental_Demographic$demo_prnt_ed_v2 == 12)
Parental_Demographic$Parental_Education[Coords_Less_than_12th_Grade] <- "<12th_Grade"
Coords_HS_GED = which(Parental_Demographic$demo_prnt_ed_v2 == 13 | Parental_Demographic$demo_prnt_ed_v2 == 14)
Parental_Demographic$Parental_Education[Coords_HS_GED] <- "HS/GED"
Coords_Some_College = which(Parental_Demographic$demo_prnt_ed_v2 == 15)
Parental_Demographic$Parental_Education[Coords_Some_College] <- "Some_College"
Coords_Associates_Degree = which(Parental_Demographic$demo_prnt_ed_v2 == 16 | Parental_Demographic$demo_prnt_ed_v2 == 17)
Parental_Demographic$Parental_Education[Coords_Associates_Degree] <- "Associates_Degree"
Coords_Bachelors = which(Parental_Demographic$demo_prnt_ed_v2 == 18)
Parental_Demographic$Parental_Education[Coords_Bachelors] <- "Bachelors_Degree"
Coords_Masters_Professional = which(Parental_Demographic$demo_prnt_ed_v2 == 19 | Parental_Demographic$demo_prnt_ed_v2 == 20)
Parental_Demographic$Parental_Education[Coords_Masters_Professional] <- "Masters_Professional_Degree"
Coords_PhD = which(Parental_Demographic$demo_prnt_ed_v2 == 21)
Parental_Demographic$Parental_Education[Coords_PhD] <- "Doctoral_Degree"
# The remainder are set to NA and imputed!
Parental_Demographic = replace(Parental_Demographic, Parental_Demographic=="", NA)

# Convert to factor!
Parental_Demographic$Parental_Education = as.factor(Parental_Demographic$Parental_Education)
Common_Participants_Dataframe = merge(Common_Participants_Dataframe, Parental_Demographic, by = c("SUBJECTKEY"))

# Save updated version
write.csv(Common_Participants_Dataframe, file = "Data/ABCD_Participant_Information/ABCD_DCAN_Common_Participants_rsfMRI_DWI_FMAP_T1W_Cognition.csv")
# Update user about number of participants left
sprintf("%.0f participants retained!", length(Common_Participants_Dataframe$SUBJECTKEY))

### Including SES Variable 2 - Parental Income ####
# For each child, we shall identify how large their household is, and the 
# parental income. We shall then compute the income-to-needs ratio as: Household
# income / US Federal Poverty Line for Relevant Household Size. The US Federal
# Poverty Line Guidelines were retrieved from: https://aspe.hhs.gov/topics/poverty
# -economic-mobility/poverty-guidelines/prior-hhs-poverty-guidelines-federal-
# register-references/2018-poverty-guidelines

# Create data frame with subject ID, total household income, and household size.
Parental_Demographic_Income_Ratio = data.frame(Parental_Demographic$SUBJECTKEY, Parental_Demographic$demo_comb_income_v2, Parental_Demographic$demo_roster_v2)
colnames(Parental_Demographic_Income_Ratio) = c("SUBJECTKEY", "TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND", "HOUSHOLD_SIZE")
# Format data frame to show household poverty guidelines
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==1] <- 12140
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==2] <- 16461
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==3] <- 20780
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==4] <- 25100
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==5] <- 29420
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==6] <- 33740
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==7] <- 38060
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==8] <- 42380
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==9] <- 46700
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==10] <- 51020
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==11] <- 55340
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==12] <- 59660
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==13] <- 63980
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==14] <- 68300
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==15] <- 72620
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==16] <- 76940
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==17] <- 81260
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==18] <- 85580
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==19] <- 89900
Parental_Demographic_Income_Ratio$'Poverty_Line'[Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE==20] <- 94220
# Adjust each income bin to its median
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND[Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND==1] <- 5000-(5000/2)
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND[Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND==2] <- 5000+((11999-5000)/2)
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND[Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND==3] <- 12000+((15999-12000)/2)
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND[Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND==4] <- 16000+((24999-16000)/2)
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND[Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND==5] <- 25000+((34999-25000)/2)
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND[Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND==6] <- 35000+((49999-35000)/2)
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND[Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND==7] <- 50000+((74999-50000)/2)
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND[Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND==8] <- 75000+((99999-75000)/2)
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND[Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND==9] <- 100000+((199999-100000)/2)
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND[Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND==10] <- 200000+((224999-200000)/2)
# Convert relevant data frame columns to numerical
Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND = as.numeric(Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND)
Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE = as.numeric(Parental_Demographic_Income_Ratio$HOUSHOLD_SIZE)
Parental_Demographic_Income_Ratio$Poverty_Line = as.numeric(Parental_Demographic_Income_Ratio$Poverty_Line)
# Now calculate the income-to-needs ratio!
Parental_Demographic_Income_Ratio$'Ratio' = Parental_Demographic_Income_Ratio$TOTAL_COMBINED_HOUSEHOLD_INCOME_BAND/Parental_Demographic_Income_Ratio$Poverty_Line
# Subset according to the included participant subject keys!
Common_Participants_Dataframe = merge(Common_Participants_Dataframe, Parental_Demographic_Income_Ratio, by = c("SUBJECTKEY"))
sprintf("%.0f participants retained!", length(Common_Participants_Dataframe$SUBJECTKEY))
# Now update the parental demographic spreadsheet with Parental Income
Parental_Demographic = merge(Parental_Demographic, Parental_Demographic_Income_Ratio, by = c("SUBJECTKEY"))

### Including SES Variable 3 - Parental Occupation Status #####
# The demo_prnt_empl_v2 variable from the parental demographic data frame
# contains data about whether the parent is working, whether they are retired, 
# a stay-at-home parent, student etc.
Parental_Demographic$demo_prnt_empl_v2 = as.factor(Parental_Demographic$demo_prnt_empl_v2)
# Recode into 6 new variables: Working, looking for work, retired, stay at home
# parent, student, or other.
Parental_Demographic$demo_prnt_empl_v2 = combineLevels(Parental_Demographic$demo_prnt_empl_v2, levs = c("1"), newLabel = c("Working"))
Parental_Demographic$demo_prnt_empl_v2 = combineLevels(Parental_Demographic$demo_prnt_empl_v2, levs = c("3"), newLabel = c("Looking For Work"))
Parental_Demographic$demo_prnt_empl_v2 = combineLevels(Parental_Demographic$demo_prnt_empl_v2, levs = c("4"), newLabel = c("Retired"))
Parental_Demographic$demo_prnt_empl_v2 = combineLevels(Parental_Demographic$demo_prnt_empl_v2, levs = c("6"), newLabel = c("Stay At Home Parent"))
Parental_Demographic$demo_prnt_empl_v2 = combineLevels(Parental_Demographic$demo_prnt_empl_v2, levs = c("2","5","8","9","10","11","7"), newLabel = c("Other"))
# Set level 777 (refuse to answer) as missing, so that we can impute!
Parental_Demographic$demo_prnt_empl_v2[which(Parental_Demographic$demo_prnt_empl_v2=="777")] = NA
Parental_Demographic$demo_prnt_empl_v2 = droplevels(Parental_Demographic$demo_prnt_empl_v2)
#Parental_Demographic$demo_prnt_empl_v2 = combineLevels(Parental_Demographic$demo_prnt_empl_v2, levs = c("777"), newLabel = c("NA"))
### Including SEX Information ####
Sex_All_PPs = data.frame(dplyr::select(Parental_Demographic, SUBJECTKEY, sex))
Common_Participants_Dataframe$sex = as.factor(Common_Participants_Dataframe$sex)
# Save updated version
write.csv(Common_Participants_Dataframe, file = "Data/ABCD_Participant_Information/ABCD_DCAN_Common_Participants_rsfMRI_DWI_FMAP_T1W_Cognition.csv")
# Update user about number of participants left
sprintf("%.0f participants retained!", length(Common_Participants_Dataframe$SUBJECTKEY))
### Including AGE Information ####
# We shall include age at scan as the variable here...
ABCD_MRI_Info = fread('Data/ABCD_MRI_Info/abcd_mri01.txt')
ABCD_MRI_Info = ABCD_MRI_Info[-c(1)]
ABCD_MRI_Info = ABCD_MRI_Info[ABCD_MRI_Info$eventname=="baseline_year_1_arm_1"]
ABCD_Ages_Subject_Key = data.frame(ABCD_MRI_Info$subjectkey, ABCD_MRI_Info$interview_age)
colnames(ABCD_Ages_Subject_Key) <- c("SUBJECTKEY", "MRI_INTERVIEW_AGE")
# Now merge with the common participants data frame!
Common_Participants_Dataframe = merge(ABCD_Ages_Subject_Key, Common_Participants_Dataframe, by = c("SUBJECTKEY"))
# Now update the parental demographic spreadsheet with MRI Interview Age 
Parental_Demographic = merge(Parental_Demographic, ABCD_Ages_Subject_Key, by = c("SUBJECTKEY"))
### Format Common_Participants_Dataframe! #####
#Common_Participants_Dataframe = Common_Participants_Dataframe[,c(1:4,158,181:185)]
#Common_Participants_Dataframe = Common_Participants_Dataframe[,c(1:2,31,179:)]
# Save updated version
#write.csv(Common_Participants_Dataframe, file = "Data/ABCD_Participant_Information/ABCD_DCAN_Common_Participants_rsfMRI_DWI_FMAP_T1W_Cognition.csv")
# Update user about number of participants left
#sprintf("%.0f participants retained!", length(Common_Participants_Dataframe$SUBJECTKEY))
### Format the Parental_Demographic Data Frame! ####
# First, ensure the Parental_Demographic data frame is updated with all the 
# necessary information.
Parental_Demographic = merge(Parental_Demographic, Age_Corrected_Cognitive_Scores, by = c("SUBJECTKEY"), all = TRUE)
Parental_Demographic = merge(Parental_Demographic, Total_Comp_Age_Corrected, by = c("SUBJECTKEY"), all = TRUE)
# Now select the columns we're interested in - ADD PARENTAL EDUCATION
Parental_Demographic = Parental_Demographic[,c(1,8,88,133,156:161,169:170)]
# Format the relevant columns
Parental_Demographic$MRI_INTERVIEW_AGE = as.numeric(Parental_Demographic$MRI_INTERVIEW_AGE)
Parental_Demographic$sex = as.factor(Parental_Demographic$sex)
Parental_Demographic$Parental_Education = as.factor(Parental_Demographic$Parental_Education)
Parental_Demographic$demo_prnt_empl_v2 = as.factor(Parental_Demographic$demo_prnt_empl_v2)
Parental_Demographic$NIHTBX_TOTALCOMP_AGECORRECTED = as.numeric(Parental_Demographic$NIHTBX_TOTALCOMP_AGECORRECTED)
# Now impute the missing data, and ensure user is updated about missing data
Parental_Demographic = replace(Parental_Demographic, Parental_Demographic=="", NA)
sprintf("Before MICE, %.0f participants had missing IQ data, equal to %.2f percent of the data frame", sum(is.na(Parental_Demographic$NIHTBX_TOTALCOMP_AGECORRECTED)), sum(is.na(Parental_Demographic$NIHTBX_TOTALCOMP_AGECORRECTED))/nrow(Parental_Demographic)*100)
sprintf("Before MICE, %.0f participants had missing age data, equal to %.2f percent of the data frame", sum(is.na(Parental_Demographic$MRI_INTERVIEW_AGE)), sum(is.na(Parental_Demographic$MRI_INTERVIEW_AGE))/nrow(Parental_Demographic)*100)
sprintf("Before MICE, %.0f participants had missing racial data, equal to %.2f percent of the data frame", sum(is.na(Parental_Demographic$Racial_Group)), sum(is.na(Parental_Demographic$Racial_Group))/nrow(Parental_Demographic)*100)
sprintf("Before MICE, %.0f participants had missing income-to-needs ratio data, equal to %.2f percent of the data frame", sum(is.na(Parental_Demographic$Ratio)), sum(is.na(Parental_Demographic$Ratio))/nrow(Parental_Demographic)*100)
sprintf("Before MICE, %.0f participants had missing gender data, equal to %.2f percent of the data frame", sum(is.na(Parental_Demographic$sex)), sum(is.na(Parental_Demographic$sex))/nrow(Parental_Demographic)*100)
sprintf("Before MICE, %.0f participants had missing parental education data, equal to %.2f percent of the data frame", sum(is.na(Parental_Demographic$Parental_Education)), sum(is.na(Parental_Demographic$Parental_Education))/nrow(Parental_Demographic)*100)
sprintf("Before MICE, %.0f participants had missing parental working status data, equal to %.2f percent of the data frame", sum(is.na(Parental_Demographic$demo_prnt_empl_v2)), sum(is.na(Parental_Demographic$demo_prnt_empl_v2))/nrow(Parental_Demographic)*100)

# Now impute the above variables. Poverty_Line was previously identified as 
# collinear, and so will be removed to avoid problems in MICE. We shall also
# remove total_combined_household_income_band from the imputation, as this can 
# only take one of several values (not continuous) and household size. We shall
# also remove PC1 scores from the imputation, as these were generated previously.
Parental_Demographic_Imputed = mice(Parental_Demographic[,-c(8)],seed = 123)
Parental_Demographic = complete(Parental_Demographic_Imputed)
write.csv(Parental_Demographic, file = "Data/ABCD_Participant_Information/Parental_Demographic_Imputed.csv")

### PART 6A - Stratified Selection of ~ 2200 Participants ####
# Load up the participants with known DWI issues
Failed = fread('/imaging/astle/am10/ABCD/Analysis/QSIprep/sorted_subject_failures.txt', header=FALSE)
Failed$V1 = str_remove(Failed$V1, "sub-")
Failed$V1 = sub("^(.{4})", "\\1_", Failed$V1)

# Stratify the parental demographic data frame by the 3 discrete variables which
# we shall report on, namely race, parental education, and sex. First, we shall
# select the participants with good quality neuro-imaging data across all modalities. 
# First, update the DCAN_Common_Participants_with_Genomics data frame to include 
# all necessary demographic information 
DCAN_Common_Participants_with_Genomics = Parental_Demographic[Parental_Demographic$SUBJECTKEY %in% DCAN_Common_Participants_with_Genomics$SUBJECTKEY,]
# Now remove those with known issues
DCAN_Common_Participants_with_Genomics = DCAN_Common_Participants_with_Genomics[-c(which(DCAN_Common_Participants_with_Genomics$SUBJECTKEY %in% Failed$V1)),]

# Create a loop which will stratify participants with random seeds each time, 
# until there are no significant differences between the two in terms of race, 
# sex, parental education, income-to-needs ratio, age, and cognition. 

# Get population and sample sizes
Required_Sample = 2200 
# for possible repeats...
Required_Sample_Label = c("Required_Sample")
Whole_Sample_Size = dim(Parental_Demographic)[1]

# Specify sex proportions
M_Percent = sum(Parental_Demographic$sex=="M", na.rm=T)/Whole_Sample_Size
F_Percent = 1 - M_Percent 
# Ensure that the order of the percentages matches that of the labels
#Parental_Demographic$sex = as.factor(Parental_Demographic$SEX)
Sex_Categories = c(levels(Parental_Demographic$sex))
Sex_Percentages = c(F_Percent,M_Percent)
Sex_Percent_Labels = c("F_Percent","M_Percent")

# Specify race proportions
Asian_Percent = sum(Parental_Demographic$Racial_Group=="Asian", na.rm = T)/Whole_Sample_Size
Black_Percent = sum(Parental_Demographic$Racial_Group=="Black", na.rm = T)/Whole_Sample_Size
Hispanic_Percent = sum(Parental_Demographic$Racial_Group=="Hispanic", na.rm = T)/Whole_Sample_Size
Other_Percent = sum(Parental_Demographic$Racial_Group=="Other", na.rm = T)/Whole_Sample_Size
White_Percent = sum(Parental_Demographic$Racial_Group=="White", na.rm = T)/Whole_Sample_Size
Race_Percentages = c(Asian_Percent,Black_Percent,Hispanic_Percent,Other_Percent,White_Percent)
Race_Categories = c(levels(Parental_Demographic$Racial_Group))
Race_Percent_Labels = c("Asian_Percent","Black_Percent","Hispanic_Percent","Other_Percent","White_Percent")

# Specify age proportions
Parental_Demographic$MRI_INTERVIEW_AGE = as.numeric(Parental_Demographic$MRI_INTERVIEW_AGE)
Parental_Demographic$'Age_Bins' = make_strata(Parental_Demographic$MRI_INTERVIEW_AGE, breaks=2)
DCAN_Common_Participants_with_Genomics$'Age_Bins' = cut(DCAN_Common_Participants_with_Genomics$MRI_INTERVIEW_AGE,breaks=c(107,119,133),labels=c("A","B"), include.lowest = TRUE)
Age_Bin_1_Percent = sum(Parental_Demographic$Age_Bins=="[107,119]",na.rm=T)/Whole_Sample_Size
Age_Bin_2_Percent = sum(Parental_Demographic$Age_Bins=="(119,133]",na.rm=T)/Whole_Sample_Size

# Specify IQ proportions
Parental_Demographic$'IQ_Bins' = make_strata(Parental_Demographic$NIHTBX_TOTALCOMP_AGECORRECTED, breaks=4)
DCAN_Common_Participants_with_Genomics$'IQ_Bins' = cut(DCAN_Common_Participants_with_Genomics$NIHTBX_TOTALCOMP_AGECORRECTED, breaks = c(32,88,100,112,221),labels=c("A","B","C","D"),include.lowest = TRUE)
IQ_Bin_1_Percent = sum(Parental_Demographic$IQ_Bins=="[32,88]",na.rm = T)/Whole_Sample_Size
IQ_Bin_2_Percent = sum(Parental_Demographic$IQ_Bins=="(88,100]",na.rm = T)/Whole_Sample_Size
IQ_Bin_3_Percent = sum(Parental_Demographic$IQ_Bins=="(100,112]",na.rm = T)/Whole_Sample_Size
IQ_Bin_4_Percent = sum(Parental_Demographic$IQ_Bins=="(112,221]",na.rm = T)/Whole_Sample_Size

set.seed(4321) 
Common_Participants_Stratified_Sample = 
  stratified(DCAN_Common_Participants_with_Genomics,
             c("IQ_Bins", "sex", "Age_Bins"),
             # Females from first age bin across range of cognition
             c("A F A" = IQ_Bin_1_Percent  * Required_Sample *F_Percent *Age_Bin_1_Percent,
               "B F A"  = IQ_Bin_2_Percent * Required_Sample *F_Percent *Age_Bin_1_Percent,
               "C F A" = IQ_Bin_3_Percent * Required_Sample *F_Percent *Age_Bin_1_Percent,
               "D F A" = IQ_Bin_4_Percent * Required_Sample *F_Percent *Age_Bin_1_Percent,
              # Males from second age bin across range of cognition
               "A M B" = IQ_Bin_1_Percent  * Required_Sample *M_Percent *Age_Bin_2_Percent,
               "B M B"  = IQ_Bin_2_Percent * Required_Sample *M_Percent *Age_Bin_2_Percent,
               "C M B" = IQ_Bin_3_Percent * Required_Sample *M_Percent *Age_Bin_2_Percent,
               "D M B" = IQ_Bin_4_Percent * Required_Sample *M_Percent *Age_Bin_2_Percent,
              # Females from second age bin across range of cognition
               "A F B" = IQ_Bin_1_Percent  * Required_Sample *F_Percent *Age_Bin_2_Percent,
               "B F B"  = IQ_Bin_2_Percent * Required_Sample *F_Percent *Age_Bin_2_Percent,
               "C F B" = IQ_Bin_3_Percent * Required_Sample *F_Percent *Age_Bin_2_Percent,
               "D F B" = IQ_Bin_4_Percent * Required_Sample *F_Percent *Age_Bin_2_Percent,
              # Males from first age bin across range of cognition
               "A M A" = IQ_Bin_1_Percent  * Required_Sample *M_Percent *Age_Bin_1_Percent,
               "B M A"  = IQ_Bin_2_Percent * Required_Sample *M_Percent *Age_Bin_1_Percent,
               "C M A" = IQ_Bin_3_Percent * Required_Sample *M_Percent *Age_Bin_1_Percent,
               "D M A" = IQ_Bin_4_Percent * Required_Sample *M_Percent *Age_Bin_1_Percent))

## Select the participants who were successfully processed!
Successfully_Processed_Participants = fread("//cbsu/data/imaging/projects/external/abcd/analyses/Alicja/QSIprep/Passed_All_Parcellations.txt")
# Temporarily change the name of the above variable to 
# 'Common_Participants_Stratified_Sample' so that the script below works!
Common_Participants_Stratified_Sample = Successfully_Processed_Participants
colnames(Common_Participants_Stratified_Sample) = "SUBJECTKEY"
Common_Participants_Stratified_Sample$SUBJECTKEY = str_remove(Common_Participants_Stratified_Sample$SUBJECTKEY, "sub-")
Common_Participants_Stratified_Sample$SUBJECTKEY = sub("^(.{4})","\\1_",Common_Participants_Stratified_Sample$SUBJECTKEY)
Common_Participants_Stratified_Sample = Parental_Demographic[Parental_Demographic$SUBJECTKEY %in% Common_Participants_Stratified_Sample$SUBJECTKEY,]

### PART 6B - Create Out-Of-Sample Demographic Data Frame! ####
# First, find the difference between the sub sample and the larger ABCD baseline
# release! Note, the out-of-sample composition will be made up of ALL participants,
# irrespective of whether they provided good quality MRI data, and whether or not
# they had missing demographic data. 
#Common_Participants_Stratified_Sample = merge(Parental_Demographic_DCAN, Common_Participants_Dataframe, by = "SUBJECTKEY")
Larger_ABCD_Release_Without_Subsample = anti_join(Parental_Demographic, Common_Participants_Stratified_Sample, by = c("SUBJECTKEY"))
# Save the above data frame!
write.csv(Larger_ABCD_Release_Without_Subsample, file = "Data/ABCD_Participant_Information/DCAN_Larger_ABCD_Release_Without_Subsample.csv")

### PART 6B(ii) - Compare Distribution of Participants with Individual GNMs to the Larger Sample####
# After conducting our neuroimaging analyses, we need to compare the distribution
# of the final sample (N = 2193) with the larger sample. First, load in the participants
# with individual GNMs.
Individual_GNMs_Participant_IDs = 
  read.table("//cbsu/data/Imaging/projects/external/abcd/analyses/Alicja/Schaefer_100_Node_Analyses/Individual_Structural_GNMs/Nodal_Parameterised_and_Simulated_Statistics/Individual_Nodal_Parameterised_and_Nodal_Statistics_IDs.txt")
# Remove the multiple headers
Individual_GNMs_Participant_IDs = data.frame(Individual_GNMs_Participant_IDs[-c(1:3),])
colnames(Individual_GNMs_Participant_IDs) = "ID"
# Format participant IDs as above
Individual_GNMs_Participant_IDs = Individual_GNMs_Participant_IDs %>%
  transform(ID = str_replace(ID, "sub-", "")) %>%
  transform(ID = str_replace(ID, "^(.{4})", "\\1_"))
# Extract the associated demographic information
Individual_GNMs_Participant_Demographics = 
  Parental_Demographic[Parental_Demographic$SUBJECTKEY %in% Individual_GNMs_Participant_IDs$ID,]

# !!!! Now we are changing the name of Individual_GNMs_Participant_Demographics to 
# !!!! "Common_Participants_Stratified_Sample", so that the code below will work.
Common_Participants_Stratified_Sample = Individual_GNMs_Participant_Demographics

# And create the out-of-sample dataset.
Larger_ABCD_Release_Without_Subsample = anti_join(Parental_Demographic, Common_Participants_Stratified_Sample, by = c("SUBJECTKEY"))

### PART 6C - Compare In- and Out-Of-Sample Distributions of Demographic Variables ####
### CONTINUOUS DISTRIBUTIONS ####
# Compare age distribution using a Kolmogorov-Smirnov Test:
Common_Participants_Stratified_Sample$MRI_INTERVIEW_AGE = as.numeric(Common_Participants_Stratified_Sample$MRI_INTERVIEW_AGE)
Larger_ABCD_Release_Without_Subsample$MRI_INTERVIEW_AGE = as.numeric(Larger_ABCD_Release_Without_Subsample$MRI_INTERVIEW_AGE)
Age_KS = ks.boot(Common_Participants_Stratified_Sample$MRI_INTERVIEW_AGE, Larger_ABCD_Release_Without_Subsample$MRI_INTERVIEW_AGE, alternative = c("two.sided"))
if (Age_KS$ks.boot.pvalue > .05) {
  print("No significant difference in AGE distribution in- and out- of sample...")
} else {
  print("Significant difference in AGE distribution in- and out- of sample...")
}

# Compare IQ distribution using a Kolmogorov-Smirnov Test...
Common_Participants_Stratified_Sample$NIHTBX_TOTALCOMP_AGECORRECTED = as.numeric(Common_Participants_Stratified_Sample$NIHTBX_TOTALCOMP_AGECORRECTED)
Larger_ABCD_Release_Without_Subsample$NIHTBX_TOTALCOMP_AGECORRECTED = as.numeric(Larger_ABCD_Release_Without_Subsample$NIHTBX_TOTALCOMP_AGECORRECTED)
NIH_TB_KS = ks.boot(Common_Participants_Stratified_Sample$NIHTBX_TOTALCOMP_AGECORRECTED, Larger_ABCD_Release_Without_Subsample$NIHTBX_TOTALCOMP_AGECORRECTED, alternative = c("two.sided"))
if (NIH_TB_KS$ks.boot.pvalue > .05) {
  print("No significant difference in IQ distribution in- and out- of sample...")
} else {
  print("Significant difference in IQ distribution in- and out- of sample...")
}

# Compare income-to-needs ratio using a Kolmogorov-Smirnov Test...
Common_Participants_Stratified_Sample$Ratio = as.numeric(Common_Participants_Stratified_Sample$Ratio)
Larger_ABCD_Release_Without_Subsample$Ratio = as.numeric(Larger_ABCD_Release_Without_Subsample$Ratio)
Income_to_Needs_Ratio_KS = ks.boot(Common_Participants_Stratified_Sample$Ratio, Larger_ABCD_Release_Without_Subsample$Ratio, alternative = c("two.sided"))
if (Income_to_Needs_Ratio_KS$ks.boot.pvalue > .05) {
  print("No significant difference in Income-to-Need Ratio distribution in- and out- of sample...")
} else {
  print("Significant difference in Income-to-Need Ratio distribution in- and out- of sample...")
}

### DISCRETE DISTRIBUTIONS #####
# Compare distribution of sex with a Chi-Square Test
Sex_Table_Stratified_Subsample = table(Common_Participants_Stratified_Sample$sex)
Sex_Table_Out_of_Sample = table(Larger_ABCD_Release_Without_Subsample$sex)
Sex_Comparison = chisq.test(rbind(Sex_Table_Stratified_Subsample, Sex_Table_Out_of_Sample))
if (Sex_Comparison$p.value > .05) {
  print("No significant difference in GENDER distribution in- and out- of sample...")
} else {
  print("Significant difference in GENDER distribution in- and out- of sample...")
}
Sex_Effect_Size = esc_chisq(chisq = Sex_Comparison$statistic, totaln = 2153, es.type = "cox.or")

# Compare distribution of race with a Chi-Square Test. Note that, to avoid errors
# in calculating the Chi-Squared statistic, we shall combine the parents_unsure
# and refuse_to_answer levels to Parents_Unsure/Refuse_to_Answer
Race_Table_Stratified_Subsample = table(Common_Participants_Stratified_Sample$Racial_Group)
Race_Table_Out_of_Sample = table(Larger_ABCD_Release_Without_Subsample$Racial_Group)
Race_Comparison = chisq.test(rbind(Race_Table_Stratified_Subsample, Race_Table_Out_of_Sample))
if (Race_Comparison$p.value > .05) {
  print("No significant difference in RACE distribution in- and out- of sample...")
} else {
  print("Significant difference in RACE distribution in- and out- of sample...")
}
Race_Effect_Size = esc_chisq(chisq = Race_Comparison$statistic, totaln= 2154, es.type = "cox.or")

# Compare distribution of parental education with a Chi-Square Test
Parental_Education_Table_Stratified_Subsample = table(Common_Participants_Stratified_Sample$Parental_Education)
Parental_Education_Table_Out_of_Sample = table(Larger_ABCD_Release_Without_Subsample$Parental_Education)
Parental_Education_Comparison = chisq.test(rbind(Parental_Education_Table_Stratified_Subsample, Parental_Education_Table_Out_of_Sample))
if (Parental_Education_Comparison$p.value > .05) {
  print("No significant difference in PARENTAL EDUCATION distribution in- and out- of sample...")
} else {
  print("Significant difference in PARENTAL EDUCATION distribution in- and out- of sample...")
}

Parental_Education_Effect_Size = esc_chisq(Parental_Education_Comparison$statistic, totaln=2154, es.type = "cox.or")

# Compare distribution of parental working status with a Chi-Square Test
Parental_Working_Table_Stratified_Subsample = table(Common_Participants_Stratified_Sample$demo_prnt_empl_v2)
Parental_Working_Table_Out_of_Sample = table(Larger_ABCD_Release_Without_Subsample$demo_prnt_empl_v2)
Parental_Working_Comparison = chisq.test(rbind(Parental_Working_Table_Stratified_Subsample, Parental_Working_Table_Out_of_Sample))
if (Parental_Working_Comparison$p.value > 0.5){
  print("No significant difference in PARENTAL WORK STATUS distribution in- and out- of sample...")
} else {
  print("Significant difference in PARENTAL WORK STATUS distribution in- and out- of sample...")
}

