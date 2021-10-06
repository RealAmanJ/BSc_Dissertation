############################################################# GRAPHICAL OUTPUT PRE-PROCESSING ##############################################
#
# Refer to subsection 2.2.4 and Supplementary Figure 1 of the Dissertation
#
# Code to split dataset by microhabitat corrected for field and also to prepare files for graphical outputting
# Revision 05/2021 akjuni@outlook.com
#
# In this code we will do the following things:
# 1) Compute statistical analysis to identify individual bacteria differentiating between microhabitats (corrected for field effect)
# 2) We will also create venn diagrams for the two separate fields
# 3) We will create and save lists for the next code (creating the UpsetR plot)
# 4) We will create and save phyloseq objects for the next code (visualizing enriched bacterial communities in phyloseq)
#
############################################################# Section 1: Start new Session #################################################

rm(list=ls())
dev.off()

############################################################# Section 2: Load Libraries ####################################################

#Required packages
library("phyloseq")
library("DESeq2")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("hrbrthemes")
library("viridis")

#If not installed, use this code
#install.packages("phyloseq")
#install.packages("DESeq2")

# If above line does not work for installing DESeq2, try the lines below
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("hrbrthemes")
#install.packages("viridis")

#Retrieve R and package versions and compare to the uploaded file in GitHub for the reproducibility of the code
sessionInfo()

#Set your working directory
setwd("")

############################################################# Section 3: Import files and subset by field ##################################

#Import the .RDS file
dat_info_42k                            <- readRDS("dat_rare_42k.rds")

#Inspect the file
class(dat_info_42k)
dat_info_42k

#Subset by field
dat_info_42k_dragoni                    <- dat_info_42k                         %>%
  subset_samples(Field == "Dragoni")

dat_info_42k_dragoni

dat_info_42k_minzoni                    <- dat_info_42k                         %>%
  subset_samples(Field == "Minzoni")

dat_info_42k_minzoni

############################################################# Section 4: Construct the DESeq object (Dragoni) ##############################

#Extract count data (usually no issues in exporting this piece of information)
JH18_42k_counts_integer_dragoni         <- dat_info_42k_dragoni                 %>%
  otu_table()

countData_dragoni                       <- JH18_42k_counts_integer_dragoni      %>%
  as.data.frame()

class(countData)
colnames(JH18_42k_counts_integer_dragoni)

#If extracted directly from phyloseq objects (e.g., using sample_data) it will
#return a phyloseq-class object not a data frame and DESeq won't work workaround:
#subset the sample data, first convert it to matrix and then to data frame
colData_dragoni                         <- as.data.frame(
                                           as.matrix(
                                           sample_data(dat_info_42k_dragoni)[
                                           colnames(JH18_42k_counts_integer_dragoni),
                                           ]))

rownames(colData_dragoni)
class(colData_dragoni)

#Construct a DESeq dataset combining count data and sample information
JH18_42k_cds_dragoni                    <- DESeqDataSetFromMatrix(
  countData = countData_dragoni                                                 ,
  colData = colData_dragoni                                                     ,
  design= ~Microhabitat)

#Execute the differential count analysis with the function DESeq
JH18_42k_cds_test_dragoni               <- JH18_42k_cds_dragoni                 %>%
  DESeq(fitType="local", betaPrior=FALSE)

############################################################# Section 5: DESeq calculation (Dragoni: Root vs Bulk) #########################

#Set the contrast
JH18_42k_contrast_dragoni_root          <- JH18_42k_cds_test_dragoni            %>%
  results(contrast = c("Microhabitat",  "Root", "Bulk"))

#Create the MA-plot
rootvbulk_dragoni_maplot                <- JH18_42k_contrast_dragoni_root       %>%
  ggmaplot(size = 3, alpha = 0.8, font.label = c(0.1, "plain", "white"))        +
  ggtitle(label = "Root vs Bulk (Dragoni)", subtitle = "Up = Root, Down = Bulk")+
  theme_classic()                                                               +
  theme(axis.title.x = element_text(color = c("#0b090a"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#0b090a"), face = "bold"))       +
  scale_fill_manual(values=c("#EC0927", "#2D3142", "#BFC0C0"))                  +
  scale_color_manual(values=c("#EC0927", "#2D3142", "#BFC0C0"))                 +
  xlab("Log2 Mean Expression")                                                  +
  ylab("Log2 Fold Change")

#Generate the plot
rootvbulk_dragoni_maplot

#Let's inspect the result file
colnames(JH18_42k_contrast_dragoni_root)

#Extract genera whose adjusted p.value in a given comparison is below 0.05
Dragoni_FDR005_root                    <- JH18_42k_contrast_dragoni_root[(
  rownames(JH18_42k_contrast_dragoni_root)[
  which(JH18_42k_contrast_dragoni_root$padj <0.05)])
  , ]

rownames(Dragoni_FDR005_root)

#Identify Genera enriched in dragoni root (first term of the comparison, positive fold change)
Dragoni_enriched_root                   <-  JH18_42k_contrast_dragoni_root[(
  rownames(JH18_42k_contrast_dragoni_root)[
  which(JH18_42k_contrast_dragoni_root$log2FoldChange > 0)])
  , ]

#Intersection: this is the list of genera significantly enriched in the root field
Dragoni_enriched_FDR005_root            <- intersect(
                                           rownames(Dragoni_FDR005_root)        ,
                                           rownames(Dragoni_enriched_root))
length(Dragoni_enriched_FDR005_root)

############################################################# Section 6: DESeq calculation (Dragoni: Rhizo vs Bulk) ########################

#Repeat what was done for Dragoni: Root vs Bulk
#Contrast
JH18_42k_contrast_dragoni_rhizo         <- JH18_42k_cds_test_dragoni            %>%
  results(contrast = c("Microhabitat",  "Rhizosphere", "Bulk"))

#Create MA-Plot
rhizovbulk_dragoni_maplot               <- JH18_42k_contrast_dragoni_rhizo      %>%
  ggmaplot(size = 3, alpha = 0.8, font.label = c(0.1, "plain", "white"))        +
  ggtitle(label = "Rhizo vs Bulk (Dragoni)", subtitle = "Up = Rhizo, Down = Bulk")+
  theme_classic()                                                               +
  theme(axis.title.x = element_text(color = c("#0b090a"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#0b090a"), face = "bold"))       +
  scale_fill_manual(values=c("#EC0927", "#2D3142", "#BFC0C0"))                  +
  scale_color_manual(values=c("#EC0927", "#2D3142", "#BFC0C0"))                 +
  xlab("Log2 Mean Expression")                                                  +
  ylab("Log2 Fold Change")

#Generate Plot
rhizovbulk_dragoni_maplot

#Inspect the result file
colnames(JH18_42k_contrast_dragoni_rhizo)

#Extract Genera w/ adjusted p.value < 0.05
Dragoni_FDR005_rhizo                    <- JH18_42k_contrast_dragoni_rhizo[(
  rownames(JH18_42k_contrast_dragoni_rhizo)[
  which(JH18_42k_contrast_dragoni_rhizo$padj <0.05)])
  , ]

rownames(Dragoni_FDR005_rhizo)

#Genera enriched in dragoni rhizo
Dragoni_enriched_rhizo                  <- JH18_42k_contrast_dragoni_rhizo[(
  rownames(JH18_42k_contrast_dragoni_rhizo)[
  which(JH18_42k_contrast_dragoni_rhizo$log2FoldChange > 0)])
  , ]

#Intersection
Dragoni_enriched_FDR005_rhizo           <- intersect(
                                           rownames(Dragoni_FDR005_rhizo)       ,
                                           rownames(Dragoni_enriched_rhizo))
length(Dragoni_enriched_FDR005_rhizo)

############################################################# Section 7: Construct the DESeq object (Minzoni) ##############################

#Repeat what was done for Dragoni for Minzoni
JH18_42k_counts_integer_minzoni         <- dat_info_42k_minzoni                 %>%
  otu_table()

#Count Data
countData_minzoni                       <- JH18_42k_counts_integer_minzoni      %>%
  as.data.frame()

class(countData_minzoni)
colnames(JH18_42k_counts_integer_minzoni)

#Column Data
colData_minzoni                         <- as.data.frame(
                                           as.matrix(
                                           sample_data(dat_info_42k_minzoni)[
                                           colnames(JH18_42k_counts_integer_minzoni),
                                           ]))
rownames(colData_minzoni)
class(colData_minzoni)


#DESeq dataset construction
JH18_42k_cds_minzoni                    <- DESeqDataSetFromMatrix(
  countData = countData_minzoni          ,
  colData = colData_minzoni              ,
  design= ~Microhabitat)

#Differential count analysis 
JH18_42k_cds_test_minzoni               <- JH18_42k_cds_minzoni                 %>%
  DESeq(fitType="local", betaPrior=FALSE)

############################################################# Section 8: DESeq calculation (Minzoni: Root vs Bulk) #########################

#Contrast
JH18_42k_contrast_minzoni_root        <- JH18_42k_cds_test_minzoni              %>%
  results(contrast = c("Microhabitat",  "Root", "Bulk"))

#MA-plot creation
rootvbulk_minzoni_maplot              <- JH18_42k_contrast_minzoni_root         %>%
  ggmaplot(size = 3, alpha = 0.8, font.label = c(0.1, "plain", "white"))        +
  ggtitle(label = "Root vs Bulk (Minzoni)", subtitle = "Up = Root, Down = Bulk")+
  theme_classic()                                                               +
  theme(axis.title.x = element_text(color = c("#0b090a"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#0b090a"), face = "bold"))       +
  scale_fill_manual(values=c("#EC0927", "#2D3142", "#BFC0C0"))                  +
  scale_color_manual(values=c("#EC0927", "#2D3142", "#BFC0C0"))                 +
  xlab("Log2 Mean Expression")                                                  +
  ylab("Log2 Fold Change")

#Plot generation
rootvbulk_minzoni_maplot

#Inspect the result file
colnames(JH18_42k_contrast_minzoni_root)

#Extract Genera w/ adjusted p.value < 0.05
Minzoni_FDR005_root                       <- JH18_42k_contrast_minzoni_root[(
  rownames(JH18_42k_contrast_minzoni_root)[
    which(JH18_42k_contrast_minzoni_root$padj <0.05)])
  , ]

rownames(Minzoni_FDR005_root)

#Genera enriched in minzoni root
Minzoni_enriched_root                     <-  JH18_42k_contrast_minzoni_root[(
  rownames(JH18_42k_contrast_minzoni_root)[
    which(JH18_42k_contrast_minzoni_root$log2FoldChange > 0)])
  , ]


#Intersection
Minzoni_enriched_FDR005_root              <- intersect(
                                             rownames(Minzoni_FDR005_root)      ,
                                             rownames(Minzoni_enriched_root))
length(Minzoni_enriched_FDR005_root)

############################################################# Section 9: DESeq calculation (Minzoni: Rhizo vs Bulk) ########################

#Contrast
JH18_42k_contrast_minzoni_rhizo           <- JH18_42k_cds_test_minzoni          %>%
  results(contrast = c("Microhabitat",  "Rhizosphere", "Bulk"))

#MA-plot creation
rhizovbulk_minzoni_maplot                 <- JH18_42k_contrast_minzoni_rhizo    %>%
  ggmaplot(size = 3, alpha = 0.8, font.label = c(0.1, "plain", "white"))        +
  ggtitle(label = "Rhizo vs Bulk (Minzoni)", subtitle = "Up = Rhizo, Down = Bulk")+
  theme_classic()                                                               +
  theme(axis.title.x = element_text(color = c("#0b090a"), face = "bold"))       +
  theme(axis.title.y = element_text(color = c("#0b090a"), face = "bold"))       +
  scale_fill_manual(values=c("#EC0927", "#2D3142", "#BFC0C0"))                  +
  scale_color_manual(values=c("#EC0927", "#2D3142", "#BFC0C0"))                 +
  xlab("Log2 Mean Expression")                                                  +
  ylab("Log2 Fold Change")

#Plot generation
rhizovbulk_minzoni_maplot

#Inspect the result file
colnames(JH18_42k_contrast_minzoni_rhizo)

#Extract Genera w/ adjusted p.value < 0.05
Minzoni_FDR005_rhizo                      <- JH18_42k_contrast_minzoni_rhizo[(
  rownames(JH18_42k_contrast_minzoni_rhizo)[
    which(JH18_42k_contrast_minzoni_rhizo$padj <0.05)])
  , ]

rownames(Minzoni_FDR005_rhizo)

#Genera enriched in minzoni bulk
Minzoni_enriched_rhizo                    <-  JH18_42k_contrast_minzoni_rhizo[(
  rownames(JH18_42k_contrast_minzoni_rhizo)[
    which(JH18_42k_contrast_minzoni_rhizo$log2FoldChange > 0)])
  , ]

#Intersection
Minzoni_enriched_FDR005_rhizo             <- intersect(
                                             rownames(Minzoni_FDR005_rhizo)     ,
                                             rownames(Minzoni_enriched_rhizo))
length(Minzoni_enriched_FDR005_rhizo)



############################################################# Section 10: UpSetR Output Files Preparation ##################################

#Total amount (use this for sanity checks)
length(Dragoni_enriched_FDR005_root)
length(Dragoni_enriched_FDR005_rhizo)
length(Minzoni_enriched_FDR005_root)
length(Minzoni_enriched_FDR005_rhizo)

########Creating the data frames
#Code to turn our lists into dataframes

#Dragoni first
Dragoni_enriched_FDR005_root_df         <- Dragoni_FDR005_root[
                                           Dragoni_enriched_FDR005_root
                                           , ]

Dragoni_enriched_FDR005_rhizo_df        <- Dragoni_FDR005_rhizo[
                                           Dragoni_enriched_FDR005_rhizo
                                           , ]


#Now Minzoni
Minzoni_enriched_FDR005_root_df         <- Minzoni_FDR005_root[
                                           Minzoni_enriched_FDR005_root
                                           , ]

Minzoni_enriched_FDR005_rhizo_df        <- Minzoni_FDR005_rhizo[
                                           Minzoni_enriched_FDR005_rhizo
                                           , ]

#######Modifying the data frames
#We are adding the basemean columns using this code below,
#note the step counter

#Dragoni
Dragoni_enriched_FDR005_root_df_1       <- as.data.frame(
                                           Dragoni_enriched_FDR005_root_df[
                                           ,1])

Dragoni_enriched_FDR005_rhizo_df_1      <- as.data.frame(
                                           Dragoni_enriched_FDR005_rhizo_df[
                                           ,1])


#Minzoni
Minzoni_enriched_FDR005_root_df_1       <- as.data.frame(
                                           Minzoni_enriched_FDR005_root_df[
                                           ,1])

Minzoni_enriched_FDR005_rhizo_df_1      <- as.data.frame(
                                           Minzoni_enriched_FDR005_rhizo_df[
                                           ,1])

######Renaming the rows and columns
#We are renaming the rownames to the ASVs they correspond to

#Dragoni
rownames(Dragoni_enriched_FDR005_root_df_1)   <- Dragoni_enriched_FDR005_root
rownames(Dragoni_enriched_FDR005_rhizo_df_1)  <- Dragoni_enriched_FDR005_rhizo

#Minzoni
rownames(Minzoni_enriched_FDR005_root_df_1)   <- Minzoni_enriched_FDR005_root
rownames(Minzoni_enriched_FDR005_rhizo_df_1)  <- Minzoni_enriched_FDR005_rhizo

#Rename the columns

#Dragoni
colnames(Dragoni_enriched_FDR005_root_df_1)   <- c("Dragoni_root")
colnames(Dragoni_enriched_FDR005_rhizo_df_1)  <- c("Dragoni_rhizo")

#Now Minzoni
colnames(Minzoni_enriched_FDR005_root_df_1)   <- c("Minzoni_root")
colnames(Minzoni_enriched_FDR005_rhizo_df_1)  <- c("Minzoni_rhizo")

######Turning the column from abundances (basemean) to boolean

#Dragoni
Dragoni_enriched_FDR005_root_df_1[Dragoni_enriched_FDR005_root_df_1 > 1]    <- 1
Dragoni_enriched_FDR005_rhizo_df_1[Dragoni_enriched_FDR005_rhizo_df_1 > 1]  <- 1

#Now Minzoni
Minzoni_enriched_FDR005_root_df_1[Minzoni_enriched_FDR005_root_df_1 > 1]    <- 1
Minzoni_enriched_FDR005_rhizo_df_1[Minzoni_enriched_FDR005_rhizo_df_1 > 1]  <- 1

######Creating the list of total unique enriched ASVs

#Union of unique roots first
Enriched_root_ASV                       <- unique(
                                           union(
                                           Dragoni_enriched_FDR005_root         ,
                                           Minzoni_enriched_FDR005_root))
#Union of unique rhizo next
Enriched_rhizo_ASV                      <- unique(
                                           union(
                                           Dragoni_enriched_FDR005_rhizo        ,
                                           Minzoni_enriched_FDR005_rhizo))
#Final union of unique all
Enriched_spinach_ASV                    <- unique(
                                           union(
                                           Enriched_root_ASV                    ,
                                           Enriched_rhizo_ASV))
length(Enriched_spinach_ASV)

######Unifying the dataset with the list of Spinach ASVs

#Dragoni
Dragoni_root_merging                    <- as.data.frame(
                                           Dragoni_enriched_FDR005_root_df_1[
                                           Enriched_spinach_ASV
                                           , ])

Dragoni_rhizo_merging                   <- as.data.frame(
                                           Dragoni_enriched_FDR005_rhizo_df_1[
                                           Enriched_spinach_ASV
                                           , ])


#Minzoni
Minzoni_root_merging                    <- as.data.frame(
                                           Minzoni_enriched_FDR005_root_df_1[
                                           Enriched_spinach_ASV
                                           , ])

Minzoni_rhizo_merging                   <- as.data.frame(
                                           Minzoni_enriched_FDR005_rhizo_df_1[
                                           Enriched_spinach_ASV
                                           , ])

######Renaming the row and column names

#Now we need to rename the row names to the
#Enriched_Spinach_ASV list (order is conserved)

#Dragoni
rownames(Dragoni_root_merging)          <- as.vector(Enriched_spinach_ASV)
rownames(Dragoni_rhizo_merging)         <- as.vector(Enriched_spinach_ASV)

#Minzoni
rownames(Minzoni_root_merging)          <- as.vector(Enriched_spinach_ASV)
rownames(Minzoni_rhizo_merging)         <- as.vector(Enriched_spinach_ASV)

#Rename column names

#Dragoni
colnames(Dragoni_root_merging)          <- c("Dragoni_root")
colnames(Dragoni_rhizo_merging)         <- c("Dragoni_rhizo")

#Minzoni
colnames(Minzoni_root_merging)          <- c("Minzoni_root")
colnames(Minzoni_rhizo_merging)         <- c("Minzoni_rhizo")

######Merging the data frames

dat_spinach_ASVs                        <- cbind(
                                           Dragoni_root_merging                 ,
                                           Dragoni_rhizo_merging)

dat_spinach_ASVs                        <- cbind(
                                           dat_spinach_ASVs                     ,
                                           Minzoni_root_merging)

dat_spinach_ASVs                        <- cbind(
                                           dat_spinach_ASVs                     ,
                                           Minzoni_rhizo_merging)

#Rename the N/A to 0
dat_spinach_ASVs[is.na(dat_spinach_ASVs)]   <- 0


######Intersections

#Bacteria enriched in both roots
root_intersect_enriched                 <- intersect(
                                           Dragoni_enriched_FDR005_root         ,
                                           Minzoni_enriched_FDR005_root)
length(root_intersect_enriched)

#Uniquely enriched bacteria
#You can just subtract from the intersect from the total to find uniquely enriched bacteria,
#but for the sake of simplicity and reproducibility I have included the commands in this code

#Bacteria enriched only in Dragoni root
root_enriched_dragoni                   <- setdiff(
                                           Dragoni_enriched_FDR005_root         ,
                                           Minzoni_enriched_FDR005_root)
length(root_enriched_dragoni)

#Bacteria enriched only in Minzoni root
root_enriched_minzoni                   <- setdiff(
                                           Minzoni_enriched_FDR005_root         ,
                                           Dragoni_enriched_FDR005_root)
length(root_enriched_minzoni)


#Bacteria enriched in both rhizos
rhizo_intersect_enriched                <- intersect(
                                           Dragoni_enriched_FDR005_rhizo        ,
                                           Minzoni_enriched_FDR005_rhizo)
length(rhizo_intersect_enriched)

#Bacteria enriched only in Dragoni rhizo
rhizo_enriched_dragoni                  <- setdiff(
                                           Dragoni_enriched_FDR005_rhizo        ,
                                           Minzoni_enriched_FDR005_rhizo)
length(rhizo_enriched_dragoni)

#Bacteria enriched only in Minzoni rhizo
rhizo_enriched_minzoni                  <- setdiff(
                                           Minzoni_enriched_FDR005_rhizo        ,
                                           Dragoni_enriched_FDR005_rhizo)
length(rhizo_enriched_minzoni)

#Bacteria enriched in the core Spinach microbiome
core_enriched                           <- intersect(
                                           root_intersect_enriched              ,
                                           rhizo_intersect_enriched)
length(core_enriched)

#Bacteria enriched in individual fields
#Dragoni
dragoni_only                            <- setdiff(
                                           unique(
                                           union(
                                           Dragoni_enriched_FDR005_root         ,
                                           Dragoni_enriched_FDR005_rhizo))      ,(
                                           unique(
                                           union(
                                           Minzoni_enriched_FDR005_root         ,
                                           Minzoni_enriched_FDR005_rhizo))))

minzoni_only                            <- setdiff(
                                           unique(
                                           union(
                                           Minzoni_enriched_FDR005_root         ,
                                           Minzoni_enriched_FDR005_rhizo))      ,(
                                           unique(
                                           union(
                                           Dragoni_enriched_FDR005_root         ,
                                           Dragoni_enriched_FDR005_rhizo))))

############################################################# Section 11: UpSetR Output Files Saving #######################################

#Save files for UpsetR analysis in subsequent code

#!!!!!!!!!!!!!!!!!!!!REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE!!!!!!!!!!!!!!!!!!!!
#dat_spinach_ASVs                                                                %>%
#  saveRDS(file ="JH18_dat_spinach_ASVs.rds")

#Dragoni_enriched_FDR005_root                                                    %>%
#  saveRDS(file ="JH18_dragoni_enriched_root_all.rds")
#Dragoni_enriched_FDR005_rhizo                                                   %>%
#  saveRDS(file ="JH18_dragoni_enriched_rhizo_all.rds")
#Minzoni_enriched_FDR005_root                                                    %>%
#  saveRDS(file ="JH18_minzoni_enriched_root_all.rds")
#Minzoni_enriched_FDR005_rhizo                                                   %>%
#  saveRDS(file ="JH18_minzoni_enriched_rhizo_all.rds")

#root_intersect_enriched                                                         %>%
#  saveRDS(file ="JH18_enriched_root_intersects.rds")
#root_enriched_dragoni                                                           %>%
#  saveRDS(file ="JH18_only_dragoni_enriched_root.rds")
#root_enriched_minzoni                                                           %>%
#  saveRDS(file ="JH18_only_minzoni_enriched_root.rds")

#rhizo_intersect_enriched                                                        %>%
#  saveRDS(file ="JH18_enriched_rhizo_intersect.rds")
#rhizo_enriched_dragoni                                                          %>%
#  saveRDS(file ="JH18_only_dragoni_enriched_rhizo.rds")
#rhizo_enriched_minzoni                                                          %>%
#  saveRDS(file ="JH18_only_minzoni_enriched_rhizo.rds")

#core_enriched                                                                   %>%
#  saveRDS(file ="Spinach_core_microbiota.rds")

#dragoni_only                                                                    %>%
#  saveRDS(file ="JH18_dragoni_only.rds")
#minzoni_only                                                                    %>%
#  saveRDS(file ="JH18_minzoni_only.rds")
#!!!!!!!!!!!!!!!!!!!!REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE!!!!!!!!!!!!!!!!!!!!

############################################################# Section 12: Phyloseq Output File Preparation #################################

######Roots
#Root intersect
JH18_root_enriched_intersect_phyloseq   <- prune_taxa(
                                           root_intersect_enriched              ,
                                           dat_info_42k)
#Dragoni root
JH18_root_enriched_dragoni_phyloseq     <- prune_taxa(
                                           root_enriched_dragoni                ,
                                           dat_info_42k)
#Minzoni root
JH18_root_enriched_minzoni_phyloseq     <- prune_taxa(
                                           root_enriched_minzoni                ,
                                           dat_info_42k)

#######Rhizo
#Rhizo intersect
JH18_rhizo_enriched_intersect_phyloseq  <- prune_taxa(
                                           rhizo_intersect_enriched             ,
                                           dat_info_42k)
#Dragoni Rhizo
JH18_rhizo_enriched_dragoni_phyloseq    <- prune_taxa(
                                           rhizo_enriched_dragoni               ,
                                           dat_info_42k)
#Minzoni Rhizo
JH18_rhizo_enriched_minzoni_phyloseq    <- prune_taxa(
                                           rhizo_enriched_minzoni               ,
                                           dat_info_42k)

######Core microbiomes
#Spinach
JH18_core_enriched_phyloseq             <- prune_taxa(
                                           core_enriched                        ,
                                           dat_info_42k)
#Dragoni
JH18_dragoni_only                       <- prune_taxa(
                                           dragoni_only                         ,
                                           dat_info_42k)
#Minzoni
JH18_minzoni_only                       <- prune_taxa(
                                           minzoni_only                         ,
                                           dat_info_42k)

############################################################# Section 13: Phyloseq Output Files Saving #####################################

#!!!!!!!!!!!!!!!!!!!!REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE!!!!!!!!!!!!!!!!!!!!
#JH18_root_enriched_intersect_phyloseq                                           %>%
#  saveRDS(file ="JH18_root_enriched_intersect_phyloseq.rds")
#JH18_root_enriched_dragoni_phyloseq                                             %>%
#  saveRDS(file ="JH18_root_enriched_dragoni_phyloseq.rds")
#JH18_root_enriched_minzoni_phyloseq                                             %>%
#  saveRDS(file ="JH18_root_enriched_minzoni_phyloseq.rds")

#JH18_root_enriched_intersect_phyloseq                                           %>%
#  saveRDS(file ="JH18_rhizo_enriched_intersect_phyloseq.rds")
#JH18_rhizo_enriched_dragoni_phyloseq                                            %>%
#  saveRDS(file ="JH18_rhizo_enriched_dragoni_phyloseq.rds")
#JH18_rhizo_enriched_minzoni_phyloseq                                            %>%
#  saveRDS(file ="JH18_rhizo_enriched_minzoni_phyloseq.rds")

#JH18_core_enriched_phyloseq                                                     %>%
#  saveRDS(file ="JH18_core_enriched_phyloseq.rds")

#JH18_dragoni_only                                                               %>%
#  saveRDS(file ="JH18_dragoni_only_phyloseq.rds")
#JH18_minzoni_only                                                               %>%
#  saveRDS(file ="JH18_minzoni_only_phyloseq.rds")
#!!!!!!!!!!!!!!!!!!!!REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE!!!!!!!!!!!!!!!!!!!!

#In next code we will use these pre-generated objects to visualize the bacteria being differentially enriched

#END OF SCRIPT
