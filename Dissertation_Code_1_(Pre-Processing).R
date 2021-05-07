############################################################# PRE-PROCESSING ##############################################################
#
# Refer to subsection 2.2.1 and Figure 1 of the Dissertation
# 
# Code to import and pre-process the ASV dataset
# Revision 05/2021 akjuni@outlook.com
#
# In this code we will do the following things:
# 1) Import the phyloseq object (ASV Dataset) https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217
# 2) Remove potential contaminants such as, chloroplast and mitochondria derived sequences
# 3) Further filter samples that were poorly sequenced and/or low count observations (interfering with the analysis)
#
############################################################# Section 1: Start new Session ################################################

rm(list=ls())
dev.off()

############################################################# Section 2: Load Libraries ###################################################

#required packages
library("phyloseq")
library("dplyr")
library("ggplot2")
library("svglite")

#If not installed, use this code
#install.packages("phyloseq")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("svglite")

#retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

#Set working directory
setwd("")

############################################################# Section 3: Import files #####################################################

#Import the .RDS file
dat_info                                <- readRDS("JH18_silva138_dada2.rds")

#inspect the file
class(dat_info)

#You should get a "phyloseq" result

############################################################# Section 4: Step 1 - Chloroplast & Mitochondria removal ######################

#Getting rid of sequencing reads assigned to Chloroplast and Mitochondria
dat_info_step1                          <- dat_info                             %>%
#Remove reads assigned to the order Chloroplast or any reads assigned to NA at order level
  subset_taxa((Order != "Chloroplast") | is.na(Order))                          %>%
#Remove reads assigned to the family Mitochondria or any reads assigned to NA at family level
  subset_taxa((Family !="Mitochondria") | is.na(Family))

dat_info_step1

#Note that we use the | (or) operator instead of & (and) operator as a read cannot be
#assigned a value AND be NA at the same time

############################################################# Section 5: Step 2 - Prune putative contaminant ASVs #########################

#Import the list of contaminant ASVs from JH06 library
Contaminant_ASVs                        <- read.delim(
  "JH06_contaminant_ASVs_ids.txt", header = FALSE)

#identify the proportion of putative contaminants in the merged object
contaminant_proportion                  <- intersect(taxa_names(dat_info_step1) ,
                                                     Contaminant_ASVs)
contaminant_proportion

#We have no contaminants, so we are good to proceed to the next step
#Increase step counter of dat_info_step1 to reflect this
dat_info_step2                          <- dat_info_step1

############################################################# Section 6: Step 3 - Remove ASVs assigned to NA at phylum level ##############

dat_info_step3                          <- dat_info_step2                       %>%
  subset_taxa(Phylum!="NA")
dat_info_step3

############################################################# Section 7: Step 4 - Abundance threshold #####################################

#Check to see if any samples need removing
sort(sample_sums(dat_info_step3))

#No low counts, no need to remove any samples

#Moving on
#Inspect the experimental design, to work out the minimum amount of
#samples for which a bacterium could be present

sample_data(dat_info_step3)
levels(sample_data(dat_info_step3)$Description)

#Note: We have 10 reps of 6 different conditions

#Abundance threshold: 20 reads, set as 16% the minimum number of samples
#This proportion is given as we have 10 reps for 6 treatments:
#if we have one bacterium in just one treatment this ~16% of the samples

dat_info_step4                          = dat_info_step3                        %>%
  filter_taxa(function(x) sum(x > 20) > (0.16 * length(x)), TRUE)

dat_info_step4
sort(sample_sums(dat_info_step4))

#filtering is very severe for water samples, because we only have 2 water
#samples instead of ten each

#ratio filtered reads/total reads
ratio                                   <- sum(sample_sums(dat_info_step4))     /
    sum(sample_sums(dat_info_step3)) * 100
ratio

#We have retained 87% of our reads

############################################################# Section 8: Step 5 - Aggregate samples at genus level ########################

#(Note: NArm set to false as Phylum NA ASVs pruned above)
dat_info_step5                          <- dat_info_step4                       %>%
  tax_glom(taxrank = "Genus", NArm = FALSE, bad_empty = c(NA, "", " ", "\t"))
dat_info_step5

#lets look at the distribution of reads
shapiro.test(sample_sums(dat_info_step5))

#It's normally distributed
#We can visualize it as well

#Determine binwidth, 60 data points
sqrt(60)
#Find lowest and highest points
summary(sample_sums(dat_info_step5))
#18777 and 203734
#Range divided by our sqrt value
(203734 - 18777)/sqrt(60)
#This is our binwidth

#Code to create a histogram
hist_samplesums                         <-  dat_info_step5@sam_data             %>%
  ggplot(aes(x = sample_sums(dat_info_step5)))                                  +
  geom_histogram(binwidth = 23877.85, fill = "#97878c", col = "#2E282A")        +
  xlab("Number of Reads")                                                       +
  ylab("Number of Samples")                                                     +
  theme_classic()                                                               +
  theme(axis.text.x = element_text(color = c("#2E282A"), size = 10))            +
  theme(axis.text.y = element_text(color = c("#2E282A")))                       +
  theme(axis.title.x = element_text(color = c("#2E282A"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#2E282A"), face = "bold"))       +
  ggtitle(label = "Read Distribution Histogram"                                 ,
          subtitle = "Shapiro-Wilk Test, p-value > 0.05")
hist_samplesums

#Great, we can save this
#hist_samplesums                                                                 %>%
#  ggsave(filename = "Histogram Sample Sums.svg", path =  "Plots/Script 1"       ,
#         device = "svg")

#We will create two phyloseq objects, one with water to see its effect on the
#overall bacterial composition, and another one with plant only samples for
#a more detailed investigation

############################################################# Section 9: Saving the datasets ##############################################

sort(sample_sums(dat_info_step5))
#We'll save the dataset with water at 18k and the dataset w/o water at 42K

#For water
#Rarefy at (18,000) and "freeze" these objects for downstream analyses
#ASVs : ignore the warnings, the object will be saved right after

#!!!!!!!!!!!!!!!!!!!!REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE!!!!!!!!!!!!!!!!!!!!
#dat_info_rare_18K                       <- dat_info_step5                       %>%
#  rarefy_even_depth(18000)
#dat_info_rare_18K                                                               %>%
#  saveRDS(file ="dat_rare_18K.rds")
#!!!!!!!!!!!!!!!!!!!!REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE!!!!!!!!!!!!!!!!!!!!
#Sample file with water created

#Remove water samples prior to creating a w/o water sample file
dat_info_nowater                        <- prune_samples(
  sample_sums(dat_info_step5) > 42000, dat_info_step5)
dat_info_nowater
sort(sample_sums(dat_info_nowater))

#For w/o water
#Rarefy at 42,000 and "freeze" these objects for downstream analyses
#ASVs : ignore the warnings, the object will be saved right after

#!!!!!!!!!!!!!!!!!!!!REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE!!!!!!!!!!!!!!!!!!!!
#dat_info_rare_42K                       <- dat_info_nowater                     %>%
#  rarefy_even_depth(42000)
#dat_info_rare_42K                                                               %>%
#  saveRDS(file ="dat_rare_42K.rds")
#!!!!!!!!!!!!!!!!!!!!REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE!!!!!!!!!!!!!!!!!!!!
#Sample file without water created

#END OF SCRIPT

###Note:
#Further analysis deemed the trade of quality (rarefication depth)
#for data (inclusion of water samples) unnecessary, and as such water
#is no longer included in the analysis