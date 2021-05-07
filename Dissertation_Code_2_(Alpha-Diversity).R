############################################################# ALPHA DIVERSITY ##############################################################
#
# Refer to subsection 2.2.2, Figure 2, and Supplementary Figure 2 of the Dissertation
#
# Code to perform alpha diversity calculations
# Revision 05/2021 akjuni@outlook.com
#
# In this code we will do the following things:
# 1) Visualize alpha diversity indexes http://www.metagenomics.wiki/pdf/definition/alpha-beta-diversity
# 2) Compute statistical analysis to test the hypothesis that field influences alpha diversity
#
############################################################# Section 1: Start new Session #################################################

rm(list=ls())
dev.off()

############################################################# Section 2: Load Libraries ####################################################

#Required packages
library("phyloseq")
library("dplyr")
library("ggplot2")
library("PMCMR")
library("PMCMRplus")
library("hrbrthemes")
library("viridis")

#If not installed, use this code
#install.packages("phyloseq")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("PMCMR")
#install.packages("PMCMRplus")
#install.packages("hrbrthemes")
#install.packages("viridis")
#install.packages("PMCMR")

#Retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

#Set working directory
setwd("")

############################################################# Section 3: Import files ######################################################

#Import the .RDS file
dat_info_42k                            <- readRDS("JH18_rare_nowater_42K_1020.rds")

############################################################# Section 4: Creating the richness plot ########################################

#Create a richness plot for the dataset (https://joey711.github.io/phyloseq/plot_richness-examples.html)

plot_42k                                <- dat_info_42k                         %>%
  plot_richness(x="Description", color="Microhabitat", shape="Field"            ,
                          measures=c("Observed", "Shannon"))                    +
  geom_point(size=4.6, alpha=0.8)                                               +
  xlab("Sample Description")                                                    +
  theme_classic()                                                               +
  scale_fill_manual(values=c("#E4572E", "#76B041", "#17BEBB"))                  +
  scale_color_manual(values=c("#E4572E", "#76B041", "#17BEBB"))                 +
  theme(strip.background =element_rect(fill="#2E282A"))                         +
  theme(axis.text.x = element_text(color = c("#2E282A"), angle = 45             ,
                                   vjust = 0.5, size = 10))                     +
  theme(axis.text.y = element_text(color = c("#2E282A")))                       +
  theme(axis.title.x = element_text(color = c("#2E282A"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#2E282A"), face = "bold"))       +
  theme(strip.text = element_text(colour = 'white', size = 13))

#Generate the plot
plot_42k

#Note that Figure 2 underwent further processing outside of R to make it easier to understand

############################################################# Section 5: Statistical Analysis Preparation ##################################

#First we need to create a new dataset which will be compatible with the statistical analysis
dat_alpha_rare_42k                      <-  estimate_richness(dat_info_42k      ,
                                            measures = c("Observed", "Shannon"))

#Generate a new dataframe for the statistical analysis
design_alpha_42k                        <-  as.data.frame(
                                            as.matrix(
                                            sample_data(
                                                        dat_info_42k)[
                                            rownames(dat_alpha_rare_42k), ]))

dat_alpha_rare_info_42k                 <- cbind(design_alpha_42k               ,
                                                 dat_alpha_rare_42k)

#Check the new dataset: it contains both description of the samples and alpha
dat_alpha_rare_info_42k


#Test the normality of the observed data
shapiro.test(dat_alpha_rare_info_42k$Observed)
shapiro.test(dat_alpha_rare_info_42k$Shannon)

#Okay lets plot it now

###FOR OBSERVED:
#Calculate binwidth, 59 data points
sqrt(59)
#Round up to 8
#Find lowest and highest value points
summary(dat_alpha_rare_info_42k$Observed)
#154 and 263
#Range divided by our sqrt value
(263-154)/8
#This is our binwidth for observed

hist_observed                           <-  dat_alpha_rare_info_42k             %>%
  ggplot(aes(x = Observed))                                                     +
  geom_histogram(binwidth = 13.625, fill = "#97878c", col = "#2E282A")          +
  xlab("Number of ASVs")                                                        +
  ylab("Number of Samples")                                                     +
  theme_classic()                                                               +
  theme(axis.text.x = element_text(color = c("#2E282A"), size = 10))            +
  theme(axis.text.y = element_text(color = c("#2E282A")))                       +
  theme(axis.title.x = element_text(color = c("#2E282A"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#2E282A"), face = "bold"))       +
  ggtitle(label = "Observed Histogram"                                          ,
          subtitle = "Shapiro-Wilk Test, p-value < 0.05")

#Generate the plot
hist_observed

###FOR SHANNON
#Calculate binwidth, 59 data points
sqrt(59)
#Round up to 8
#Find lowest and highest value points
summary(dat_alpha_rare_info_42k$Shannon)
#2.141 and 4.275
#Range divided by our sqrt value
(4.275-2.141)/8
#This is our binwidth for Shannon

hist_shannon                            <-  dat_alpha_rare_info_42k             %>%
  ggplot(aes(x = Shannon))                                                      +
  geom_histogram(binwidth = 0.26675, fill = "#97878c", col = "#2E282A")         +
  xlab("Shannon Index")                                                         +
  ylab("Number of Samples")                                                     +
  theme_classic()                                                               +
  theme(axis.text.x = element_text(color = c("#2E282A"), size = 10))            +
  theme(axis.text.y = element_text(color = c("#2E282A")))                       +
  theme(axis.title.x = element_text(color = c("#2E282A"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#2E282A"), face = "bold"))       +
  ggtitle(label = "Shannon Histogram"                                           ,
          subtitle = "Shapiro-Wilk Test, p-value < 0.05")

#Generate the plot
hist_shannon

#Neither data is normally distributed
#Need to run non-parametric tests on both

############################################################# Section 6: Statistical Analysis ##############################################

#We run kruskal tests and not wilcox tests, because kruskal test is
#an extension of the wilcox test. While wilcox tests only work if
#the variable has 2 factors, kruskal tests can work in the presence
#of 2 or more variables

#First check whether the microhabitat has an effect on observed
kruskal.test(Observed ~ Microhabitat, data = dat_alpha_rare_info_42k)

#Signifcant effect of the microhabitat: we can now consider field
kruskal.test(Observed ~ Field, data = dat_alpha_rare_info_42k)

#No significant effect between fields

#Post-hoc Dunn's test, to look at this variance pairwise
#https://www.rdocumentation.org/packages/PMCMR/versions/4.3/topics/posthoc.kruskal.dunn.test


#Turn the microhabitat into a factor
microhabitat_character_42k              <- as.factor(
                                           dat_alpha_rare_info_42k$Microhabitat)

posthoc.kruskal.dunn.test (x=dat_alpha_rare_info_42k$Observed                   ,
                           g=microhabitat_character_42k, p.adjust.method="BH")

#Repeat for shannon
kruskal.test(Shannon ~ Microhabitat, data= dat_alpha_rare_info_42k)
kruskal.test(Shannon ~ Field, data = dat_alpha_rare_info_42k)

posthoc.kruskal.dunn.test (x=dat_alpha_rare_info_42k$Shannon                    ,
                           g=microhabitat_character_42k, p.adjust.method="BH")

#Same results for shannon as well

#END OF SCRIPT

###Note:
#Comments on the analysis:
#We can now conclude that, in the tested conditions, the fields do
#not have a significant difference in  alpha diversity. Microhabitat
#has a significant effect on the alpha diversity as expected