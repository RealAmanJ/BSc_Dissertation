############################################################# BETA DIVERSITY ###############################################################
#
# Refer to subsection 2.2.3 and Figure 3 of Dissertation
#
# Code to perform beta diversity calculations
# Revision 05/2021 akjuni@outlook.com
#
# In this code we will do the following things:
# 1) Visualize beta diversity indexes http://www.metagenomics.wiki/pdf/definition/alpha-beta-diversity
# 2) Compute statistical analysis to test the hypothesis that fields have a difference in beta diversity
#
############################################################# Section 1: Start new Session #################################################

rm(list=ls())
dev.off()

############################################################# Section 2: Load Libraries ####################################################

#Required packages
library("phyloseq")
library("dplyr")
library("ggplot2")
library("vegan")
library("hrbrthemes")
library("viridis")

#If not installed, use this code
#install.packages("phyloseq")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("vegan")
#install.packages("hrbrthemes")
#install.packages("viridis")

#retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

#Set your working directory
setwd("")

############################################################# Section 3: Import files ######################################################

#Import the .RDS file
dat_info_42k                            <- readRDS("dat_rare_42K.rds")

############################################################# Section 4: Data Visualization ################################################

###Data Visualization (https://joey711.github.io/phyloseq/distance.html)
#Create the bray-distance matrix for the 42k dataset
BC_42k                                  <- dat_info_42k                         %>%
  phyloseq::distance("bray")

#let's inspect the file
BC_42k


#Create the file used for statistical analysis later
dat_info_CAP_BC_42k                     <- dat_info_42k                         %>%
  ordinate("CAP","bray", ~ Field)

#Create the canonical analysis of principal coordinates plot
p_42k                                   <- dat_info_42k                         %>%
  plot_ordination(dat_info_CAP_BC_42k, shape = "Field", color = "Microhabitat") +
  geom_point(size=4.6, alpha=0.8)                                               +
  theme_classic()                                                               +
  scale_fill_manual(values=c("#E4572E", "#76B041", "#17BEBB"))                  +
  scale_color_manual(values=c("#E4572E", "#76B041", "#17BEBB"))                 +
  theme(axis.title.x = element_text(color = c("#2E282A"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#2E282A"), face = "bold"))       +
 ggtitle(label = "Beta Diversity"                                               ,
         subtitle = "Adonis Test, p-value < 0.05, Microhabitat R2 = 0.779, Field R2 = 0.03, Microhabitat * Field R2 = 0.02") +
  xlab("Field [3%]")                                                            +
  ylab("Microhabitat [78%]")

#Generate the plot
p_42k

############################################################# Section 5: Statistical Analysis ##############################################

#ANOVA on the axis (this is to assess the robustness of the ordination)
dat_info_CAP_BC_42k                                                             %>%
  anova(permutations = 5000)

#If this anova result is insignificant, we are not allowed to present this ordination
#because the CAP forces the graph to look nicer, the anova checks if any random presentation
#will be better than the one we have generated.

#Permutational analysis of variance (also called adonis)
adonis(BC_42k ~ Microhabitat * Field                                            ,
       data= as.data.frame(as.matrix(sample_data(dat_info_42k)))                ,
       permutations = 5000)

#R2 = the amount of variation explained by that variable (convert to percentages)
#Residuals is everything not listed as an independent variable, it may be the variation between PCR machines,
#and/or people carrying it out, if residual goes above 50% it becomes a problem


#END OF SCRIPT
###Comments on the analysis
#We can now conclude that our CAP is solid, the microhabitat is the major effect
#accounting for 80% of the variance.
