############################################################# GRAPHICAL OUTPUTS ##############################################
#
# Refer to subsection 2.3, Figures 4-7, Tables 1-3, and Supplementary Figure 3 of the Dissertation
#
# Code to generate the UpsetR plot, the phyloseq bar plots, the lolipop plots and the genera tables
# Revision 05/2021 akjuni@outlook.com
#
# In this code we will do the following things:
# 1) Generate an UpsetR plot to visualize intersections between field microhabitats and their enriched bacteria
# 2) Use phyloseq to visualize the different enriched bacterial communities
# 3) Create lollipop plots to help visualize the bacterial families
# 4) Make tables to help visualize the bacterial genera
#
############################################################# Section 1: Start new Session ###################################

rm(list=ls())
dev.off()

############################################################# Section 2: Load Libraries ######################################

#required packages
library("phyloseq")
library ("DESeq2")
library("dplyr")
library("ggplot2")
library("UpSetR")
library("gt")
library("webshot")
library("hrbrthemes")
library("viridis")

#If not installed, use this code
#install.packages("phyloseq")
#install.packages("DESeq2")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("gt")
#install.packages("webshot")
#install.packages("UpSetR")
#install.packages("hrbrthemes")
#install.packages("viridis")


#Retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

#Set your working directory
setwd("")

############################################################# Section 3: Import files and subset by field ####################

#Import the .RDS file
dat_info_42k                            <- readRDS("JH18_rare_nowater_42K_1020.rds")
dat_info_42k                            <- readRDS("JH18_rare_nowater_42K_1020.rds")
#Inspect the files
class(dat_info_42k)
dat_info_42k

dat_info_42k_dragoni                    <- dat_info_42k                         %>%
  subset_samples(Field == "Dragoni")

dat_info_42k_dragoni

dat_info_42k_minzoni                    <- dat_info_42k                         %>%
  subset_samples(Field == "Minzoni")

#UpSetR plot files
dat_spinach_ASVs                        <- readRDS("JH18_dat_spinach_ASVs.rds")
Dragoni_enriched_FDR005_root            <- readRDS("JH18_dragoni_enriched_root_all.rds")
Dragoni_enriched_FDR005_rhizo           <- readRDS("JH18_dragoni_enriched_rhizo_all.rds")
Minzoni_enriched_FDR005_root            <- readRDS("JH18_minzoni_enriched_root_all.rds")
Minzoni_enriched_FDR005_rhizo           <- readRDS("JH18_minzoni_enriched_rhizo_all.rds")

root_intersect_enriched                 <- readRDS("JH18_enriched_root_intersects.rds")
root_enriched_dragoni                   <- readRDS("JH18_only_dragoni_enriched_root.rds")
root_enriched_minzoni                   <- readRDS("JH18_only_minzoni_enriched_root.rds")

rhizo_intersect_enriched                <- readRDS("JH18_enriched_rhizo_intersect.rds")
rhizo_enriched_dragoni                  <- readRDS("JH18_only_dragoni_enriched_rhizo.rds")
rhizo_enriched_minzoni                  <- readRDS("JH18_only_minzoni_enriched_rhizo.rds")

core_enriched                           <- readRDS("Spinach_core_microbiota.rds")

dragoni_only                            <- readRDS("JH18_dragoni_only.rds")
minzoni_only                            <- readRDS("JH18_minzoni_only.rds")

#Phyloseq analysis files
JH18_root_enriched_intersect            <- readRDS("JH18_root_enriched_intersect_phyloseq.rds")
JH18_root_enriched_dragoni              <- readRDS("JH18_root_enriched_dragoni_phyloseq.rds")
JH18_root_enriched_minzoni              <- readRDS("JH18_root_enriched_minzoni_phyloseq.rds")

JH18_root_enriched_intersect            <- readRDS("JH18_rhizo_enriched_intersect_phyloseq.rds")
JH18_rhizo_enriched_dragoni             <- readRDS("JH18_rhizo_enriched_dragoni_phyloseq.rds")
JH18_rhizo_enriched_minzoni             <- readRDS("JH18_rhizo_enriched_minzoni_phyloseq.rds")

JH18_core_enriched                      <- readRDS("JH18_core_enriched_phyloseq.rds")

JH18_dragoni_only                       <- readRDS("JH18_dragoni_only_phyloseq.rds")
JH18_minzoni_only                       <- readRDS("JH18_minzoni_only_phyloseq.rds")

############################################################# Section 4: UpSetR Plot #########################################

#Generating the UpsetR plot (https://academic.oup.com/bioinformatics/article/33/18/2938/3884387)

######Creating the lists

#All fields and microhabitats
list_all_intersection                   <- list(
                                           "Dragoni_root"                       ,
                                           "Minzoni_root"                       ,
                                           "Dragoni_rhizo"                      ,
                                           "Minzoni_rhizo")

#Dragoni microhabitats
list_dragoni_shared                     <- list(
                                           "Dragoni_root"                       ,
                                           "Dragoni_rhizo")

#Minzoni microhabitats
list_minzoni_shared                     <- list(
                                           "Minzoni_root"                       ,
                                           "Minzoni_rhizo")

#Individual fields and microhabitats
list_dragoni_rhizo                      <- list("Dragoni_rhizo")
list_minzoni_rhizo                      <- list("Minzoni_rhizo")
list_dragoni_root                       <- list("Dragoni_root")
list_minzoni_root                       <- list("Minzoni_root")

#Lists that denotes what lists to intersect
lists_merged                            <- list(
                                           list_all_intersection                ,
                                           list_dragoni_shared                  ,
                                           list_dragoni_rhizo                   ,
                                           list_dragoni_root                    ,
                                           list_minzoni_shared                  ,
                                           list_minzoni_rhizo                   ,
                                           list_minzoni_root)

#We've done all the pre-processing in the previous code, so this is really just one
#line of code to get a graphical output

#Create the plot
upsetr_plot                             <- dat_spinach_ASVs                     %>%
  upset(          sets = c("Dragoni_root"                                       ,
                           "Minzoni_root"                                       ,
                           "Dragoni_rhizo"                                      ,
                           "Minzoni_rhizo")                                     ,
        sets.bar.color = "#56B4E9"                                              ,
        order.by = "freq"                                                       ,
        intersections = lists_merged                                            ,
        matrix.color = "#000000"                                                ,
        main.bar.color = "#000000"                                              ,
        mainbar.y.label = "Number of Shared Bacterial Families"                 ,
        text.scale = 1.3)

#Generate the plot
upsetr_plot

############################################################# Section 5: Bar plots (Phyloseq) ################################

#Create the 29 colors, color-blind friendly palette, courtesy of https://coolors.co/ for generating
#and https://davidmathlogic.com/colorblind/ for double-checking
cbPalette                               <- c("#F633F8", "#A9BED7", "#E4CCA2"    ,
                                             "#8A705F", "#634649", "#675B59"    ,
                                             "#2F4F5C", "#12227A", "#5B4B5F"    ,
                                             "#2EC4B6", "#E71D36", "#FF9F1C"    ,
                                             "#840032", "#E59500", "#E08DAC"    ,
                                             "#6A7FDB", "#1ED7F0", "#AF3B6E"    ,
                                             "#21FA90", "#478978", "#B6EEA6"    ,
                                             "#F2FF49", "#355834", "#6E633D"    ,
                                             "#EEF5DB", "#D5C5C8", "#AA3A84"    ,
                                             "#953AE0", "#47E066")

#Please note that this palette is made with Protoanomaly (the type of colorblindness I have)
#in mind. It may not be suitable for the other types of colorblindness.

###Creating the plots
#Dragoni only enriched
plot_dragoni_only                       <- JH18_dragoni_only                    %>%
  plot_bar("Microhabitat", fill = "Family"                                      ,
           facet_grid= "Field", title= "Dragoni Enriched Both Microhabitats")   +
  scale_fill_manual(values = cbPalette)                                         +
  theme_classic()                                                               +
  theme(strip.background = element_rect(fill = "black"))                        +
  theme(strip.text = element_text(colour = 'white', size = 13))                 +
  theme(axis.title.x = element_text(color = c("#2E282A"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#2E282A"), face = "bold"))


#Minzoni only enriched
plot_minzoni_only                       <- JH18_minzoni_only                    %>%
  plot_bar("Microhabitat", fill = "Family"                                      ,
           facet_grid= "Field", title = "Minzoni Enriched Both Microhabitats")  +
  scale_fill_manual(values = cbPalette)                                         +
  theme_classic()                                                               +
  theme(strip.background =element_rect(fill="black"))                           +
  theme(strip.text = element_text(colour = 'white', size = 13))                 +
  theme(axis.title.x = element_text(color = c("#2E282A"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#2E282A"), face = "bold"))

###Generating the plots
plot_dragoni_only
plot_minzoni_only

############################################################# Section 6: Count Analysis Preparation ##########################

######For the lollipop plot
#We want to order them by the genera count (descending)

###Core enriched
#Create dataframe of families and genera
dat_core_enriched                       <- as.data.frame(
                                           JH18_core_enriched@tax_table@.Data)

#Remove families with NAs at genus level and save as a list
CE_genus_list                           <- na.omit(dat_core_enriched$Genus)

#Create a dataframe of the table
CE_table                                <- as.data.frame(table(
                                           dat_core_enriched$Family))

#Rename the columns and order by decreasing
colnames(CE_table)                      <- c("Family", "Genera Count")

CE_table                                <- CE_table                             %>%
  arrange(-`Genera Count`)                                                      %>%
  mutate(Family = factor(Family, levels = Family))


###Dragoni Enriched
#Dataframe fam+gen
dat_dragoni_enriched                    <- as.data.frame(
                                           JH18_dragoni_only@tax_table@.Data)
#Clear NAs
dragoni_genus_list                      <- na.omit(dat_dragoni_enriched$Genus)

#Table dataframe
dragoni_table                           <- as.data.frame(table(
                                           dat_dragoni_enriched$Family))

#Renaming and ordering
colnames(dragoni_table)                 <- c("Family", "Genera Count")

dragoni_table                           <- dragoni_table                        %>%
  arrange(-`Genera Count`)                                                      %>%
  mutate(Family = factor(Family, levels = Family))

###Minzoni Enriched
#Dataframe fam+gen
dat_minzoni_enriched                    <- as.data.frame(
                                           JH18_minzoni_only@tax_table@.Data)

#Clear NAs
minzoni_genus_list                      <- na.omit(dat_minzoni_enriched$Genus)

#Table Dataframe
minzoni_table                           <- as.data.frame(table(
                                           dat_minzoni_enriched$Family))

#Renaming and ordering
colnames(minzoni_table)                 <- c("Family", "Genera Count")

minzoni_table                           <- minzoni_table                        %>%
  arrange(-`Genera Count`)                                                      %>%
  mutate(Family = factor(Family, levels = Family))

######For the tables
#We want to order them by family (alphabetically)

###Core enriched
#Remove NAs
ce_family_genera                          <- na.omit(dat_core_enriched[,5:6])

#Sort families
ce_family_genera                          <- ce_family_genera %>%
  arrange(Family)

###Dragoni
#Remove NAs
dragoni_family_genera                     <- na.omit(dat_dragoni_enriched[,5:6])

#Sort families
dragoni_family_genera                     <- dragoni_family_genera              %>%
  arrange(Family)

###Minzoni
#Remove NAs
minzoni_family_genera                     <- na.omit(dat_minzoni_enriched[,5:6])

#Sort families
minzoni_family_genera                     <- minzoni_family_genera              %>%
  arrange(Family)

############################################################# Section 7: Lolipop Plots #######################################

###Core Enriched
length(CE_genus_list)

#Create the plot
lolipop_CE                              <- CE_table                             %>%
  ggplot(aes(x = Family, y = `Genera Count`))                                   +
  geom_point(size=3, alpha=1, color = "#0b090a")                                +
  geom_segment(aes(x = Family, xend = Family, y = 0, yend = `Genera Count` ))   +
  coord_flip()                                                                  +
  theme_classic()                                                               +
  theme(axis.title.x = element_text(color = c("#0b090a"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#0b090a"), face = "bold"))       +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11))                     +
  annotate(geom = "text", label = "Total Genera Identified = 38", x = 23, y = 9.6)

#Generate the plot
lolipop_CE

###Dragoni
length(dragoni_genus_list)

#Create the plot
lolipop_dragoni                         <- dragoni_table                        %>%
  ggplot(aes(x = Family, y = `Genera Count`))                                   +
  geom_point(size=3, alpha=1, color = "#0b090a")                                +
  geom_segment(aes(x = Family, xend = Family, y = 0, yend = `Genera Count` ))   +
  coord_flip()                                                                  +
  theme_classic()                                                               +
  theme(axis.title.x = element_text(color = c("#0b090a"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#0b090a"), face = "bold"))       +
  annotate(geom = "text", label = "Total Genera Identified = 17", x = 13, y = 3.5)

#Generate the plot
lolipop_dragoni

###Minzoni
length(minzoni_genus_list)

#Create the plot
lolipop_minzoni                         <- minzoni_table                        %>%
  ggplot(aes(x = Family, y = `Genera Count`))                                   +
  geom_point(size=3, alpha=1, color = "#0b090a")                                +
  geom_segment(aes(x = Family, xend = Family, y = 0, yend = `Genera Count` ))   +
  coord_flip()                                                                  +
  theme_classic()                                                               +
  theme(axis.title.x = element_text(color = c("#0b090a"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#0b090a"), face = "bold"))       +
  scale_y_continuous(breaks = c(0,1,2))                                         +
  annotate(geom = "text", label = "Total Genera Identified = 9", x = 9, y = 1.8)

#Generate the plot
lolipop_minzoni

############################################################# Section 8: Genera Tables #######################################

###Core Enriched
#Create the table
genera_gt_CE                            <- gt(ce_family_genera)                 %>%
  cols_label(Family = md("**Family**"), `Genus` = md("**Genus**"))              %>%
  cols_align(align = c("center"))

#Generate the table
genera_gt_CE

#Save the table
#gtsave(data = genera_gt_CE, filename = "Genera Core Enriched.png")


###Dragoni
#Create the table
genera_gt_dragoni                       <- gt(dragoni_family_genera)            %>%
  cols_label(Family = md("**Family**"), `Genus` = md("**Genus**"))              %>%
  cols_align(align = c("center"))

#Generate the table
genera_gt_dragoni

#Save the table
#gtsave(data = genera_gt_dragoni, filename = "Genera Dragoni Enriched.png")


###Minzoni
#Create the table
genera_gt_minzoni                       <- gt(minzoni_family_genera)            %>%
  cols_label(Family = md("**Family**"), `Genus` = md("**Genus**"))              %>%
  cols_align(align = c("center"))

#Generate the table
genera_gt_minzoni

#Save the table
#gtsave(data = genera_gt_minzoni, filename = "Genera Minzoni Enriched.png")

#END OF SCRIPT