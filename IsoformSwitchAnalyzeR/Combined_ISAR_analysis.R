#Combined IsoformSwitchAnalyzeR Analysis
#Written by: Caleb Embree
#Last Updated: 2/24/2024
#Changelog: 9/1/2023: Created
# 10/1/2023: Modified the graphs for the 2023 RRM
# 2/27/2024: Look at the length of ORF for various classes of transcript
# 3/20/2024: Changed the length plot to match other boxplots made

####Load Packages####
library(readr)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(ggbeeswarm)
library(eulerr)
library(biomaRt)
library(tidyverse)
setwd("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis")

####Analyze the TPMs of PTC+, MANE, and Novel isoforms as determined by ISAR####
##load in data##
#EFTUD2
EFTUD2_TPM <- read_csv("EFTUD2_TPM.csv")
EFTUD2_TPM = EFTUD2_TPM %>% mutate(Sample = "EFTUD2")
#AQR
AQR_TPM <- read_csv("AQR_TPM.csv")
AQR_TPM = AQR_TPM %>% mutate(Sample = "AQR")
#SF3B1
SF3B1_TPM <- read_csv("SF3B1_TPM.csv")
SF3B1_TPM = SF3B1_TPM %>% mutate(Sample = "SF3B1")
#PRPF4
PRPF4_TPM <- read_csv("PRPF4_TPM.csv")
PRPF4_TPM = PRPF4_TPM %>% mutate(Sample = "PRPF4")
#EIF4A3
EIF4A3_TPM <- read_csv("EIF4A3_TPM.csv")
EIF4A3_TPM = EIF4A3_TPM %>% mutate(Sample = "EIF4A3")
#UPF1
UPF1_TPM <- read_csv("UPF1_TPM.csv")
UPF1_TPM = UPF1_TPM %>% mutate(Sample = "UPF1")
#PUM1
# PUM1_TPM <- read_csv("PUM1_TPM.csv")
# PUM1_TPM = PUM1_TPM %>% mutate(Sample = "PUM1")
#ACTN4
# ACTN4_TPM <- read_csv("ACTN4_TPM.csv")
# ACTN4_TPM = ACTN4_TPM %>% mutate(Sample = "ACTN4")
#SNRNP200
SNRNP200_TPM <- read_csv("SNRNP200_TPM.csv")
SNRNP200_TPM = SNRNP200_TPM %>% mutate(Sample = "SNRNP200")
#U2AF1
U2AF1_TPM <- read_csv("U2AF1_TPM.csv")
U2AF1_TPM = U2AF1_TPM %>% mutate(Sample = "U2AF1")
#SF3A3
SF3A3_TPM <- read_csv("SF3A3_TPM.csv")
SF3A3_TPM = SF3A3_TPM %>% mutate(Sample = "SF3A3")
#RBM22
RBM22_TPM <- read_csv("RBM22_TPM.csv")
RBM22_TPM = RBM22_TPM %>% mutate(Sample = "RBM22")
#CDC40
CDC40_TPM <- read_csv("CDC40_TPM.csv")
CDC40_TPM = CDC40_TPM %>% mutate(Sample = "CDC40")
#PRPF3
PRPF3_TPM <- read_csv("PRPF3_TPM.csv")
PRPF3_TPM = PRPF3_TPM %>% mutate(Sample = "PRPF3")
#AGO1
AGO1_TPM <- read_csv("AGO1_TPM.csv")
AGO1_TPM = AGO1_TPM %>% mutate(Sample = "AGO1")
#SF3A1
SF3A1_TPM <- read_csv("SF3A1_TPM.csv")
SF3A1_TPM = SF3A1_TPM %>% mutate(Sample = "SF3A1")
#SF3B3
SF3B3_TPM <- read_csv("SF3B3_TPM.csv")
SF3B3_TPM = SF3B3_TPM %>% mutate(Sample = "SF3B3")
#SNRNP70
SNRNP70_TPM <- read_csv("SNRNP70_TPM.csv")
SNRNP70_TPM = SNRNP70_TPM %>% mutate(Sample = "SNRNP70")
#SNRPC
SNRPC_TPM <- read_csv("SNRPC_TPM.csv")
SNRPC_TPM = SNRPC_TPM %>% mutate(Sample = "SNRPC")
#GNB2L1
GNB2L1_TPM <- read_csv("GNB2L1_TPM.csv")
GNB2L1_TPM = GNB2L1_TPM %>% mutate(Sample = "GNB2L1")
#CDC5L
CDC5L_TPM <- read_csv("CDC5L_TPM.csv")
CDC5L_TPM = CDC5L_TPM %>% mutate(Sample = "CDC5L")
#HNRNPK
HNRNPK_TPM <- read_csv("HNRNPK_TPM.csv")
HNRNPK_TPM = HNRNPK_TPM %>% mutate(Sample = "HNRNPK")
#MAGOH
MAGOH_TPM <- read_csv("MAGOH_TPM.csv")
MAGOH_TPM = MAGOH_TPM %>% mutate(Sample = "MAGOH")

####Determine which PTC+ transcripts are overlapped####
EFTUD2_PTC = EFTUD2_TPM %>% filter(type == "PTC")
AQR_PTC = AQR_TPM %>% filter(type == "PTC")
SF3B1_PTC = SF3B1_TPM %>% filter(type == "PTC")
PRPF4_PTC = PRPF4_TPM %>% filter(type == "PTC")
UPF1_PTC = UPF1_TPM %>% filter(type == "PTC")
EIF4A3_PTC = EIF4A3_TPM %>% filter(type == "PTC")
PTC_iso = list("EFTUD2 PTC isoforms" = EFTUD2_PTC$isoform_id,
               "AQR PTC isoforms" = AQR_PTC$isoform_id,
               "SF3B1 PTC isoforms" = SF3B1_PTC$isoform_id,
               "PRPF4 PTC isoforms" = PRPF4_PTC$isoform_id,
               "UPF1 PTC isoforms" = UPF1_PTC$isoform_id,
               "EIF4A3 PTC isoforms" = EIF4A3_PTC$isoform_id)

PTC_venn = plot(euler(PTC_iso), quantities = TRUE, main = "Shared PTC+ isoforms with TPM>1")
PTC_venn
ggsave("PTC_Isoform_Overlap.pdf",
       plot = PTC_venn,
       device = "pdf",
       width = 8,
       height = 8,
       units = "in",
       dpi = 300)
Shared_PTC = EFTUD2_PTC %>% inner_join(AQR_PTC,by = "isoform_id") %>% inner_join(SF3B1_PTC, by = "isoform_id") %>% 
  inner_join(PRPF4_PTC, by = "isoform_id") %>% inner_join(UPF1_PTC, by = "isoform_id") %>% inner_join(EIF4A3_PTC, by = "isoform_id")
Shared_PTC = Shared_PTC %>% dplyr::select(1,9:12) %>% rename(geneID = ensembl_gene_id.x,gene = external_gene_name.x,
                                                      biotype = transcript_biotype.x, isoformID = isoform_id,
                                                      MANE = transcript_mane_select.x)
Shared_PTC_iso = Shared_PTC %>% dplyr::select(1) #108 shared PTC isoforms
Shared_PTC_genes = Shared_PTC %>% dplyr::select(5) %>% distinct() #101 shared genes
Stringent_PTC_geneID <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/spliceosomeKD_analysis/Stringent_PTC_geneID.csv")
Shared_sPTC = Stringent_PTC_geneID %>% filter(PTC == TRUE) %>% filter(transID %in% Shared_PTC$isoformID)
Shared_sPTC_neg = Stringent_PTC_geneID %>% filter(PTC == FALSE) %>% filter(geneID %in% Shared_PTC$geneID)
Shared_sPTC = Shared_sPTC %>% bind_rows(Shared_sPTC_neg)
write_csv(Shared_sPTC, "shared_stringent_PTC.csv", col_names = TRUE)


#### Pull the TPMs of the shared transcripts#### COMMENTED OUT BECAUSE IT HASN't BEEN UPDATED
# EFTUD2_sPTC = EFTUD2_TPM %>% inner_join(Shared_PTC_iso) #Filter to just the shared PTC isoforms
# EFTUD2_sMANE = EFTUD2_TPM %>% filter(type == "MANE") %>% inner_join(Shared_PTC_genes, by = c("ensembl_gene_id" = "geneID")) #Pull just the MANE isoforms from the shared PTC genes
# EFTUD2_shared_TPM = EFTUD2_TPM %>% filter(type == "Novel") %>% full_join(EFTUD2_sPTC) %>% full_join(EFTUD2_sMANE) #Combine them together
# EFTUD2_shared_sum = EFTUD2_shared_TPM %>% group_by(type) %>% summarise(total_TPM = sum(kdTPM), count = n()) %>% mutate(sample = "EFTUD2")
# #AQR
# AQR_sPTC = AQR_TPM %>% inner_join(Shared_PTC_iso) #Filter to just the shared PTC isoforms
# AQR_sMANE = AQR_TPM %>% filter(type == "MANE") %>% inner_join(Shared_PTC_genes, by = c("ensembl_gene_id" = "geneID")) #Pull just the MANE isoforms from the shared PTC genes
# AQR_shared_TPM = AQR_TPM %>% filter(type == "Novel") %>% full_join(AQR_sPTC) %>% full_join(AQR_sMANE) #Combine them together
# AQR_shared_sum = AQR_shared_TPM %>% group_by(type) %>% summarise(total_TPM = sum(kdTPM), count = n()) %>% mutate(sample = "AQR")
# #SF3B1
# SF3B1_sPTC = SF3B1_TPM %>% inner_join(Shared_PTC_iso) #Filter to just the shared PTC isoforms
# SF3B1_sMANE = SF3B1_TPM %>% filter(type == "MANE") %>% inner_join(Shared_PTC_genes, by = c("ensembl_gene_id" = "geneID")) #Pull just the MANE isoforms from the shared PTC genes
# SF3B1_shared_TPM = SF3B1_TPM %>% filter(type == "Novel") %>% full_join(SF3B1_sPTC) %>% full_join(SF3B1_sMANE) #Combine them together
# SF3B1_shared_sum = SF3B1_shared_TPM %>% group_by(type) %>% summarise(total_TPM = sum(kdTPM), count = n()) %>% mutate(sample = "SF3B1")
# #PRPF4
# PRPF4_sPTC = PRPF4_TPM %>% inner_join(Shared_PTC_iso) #Filter to just the shared PTC isoforms
# PRPF4_sMANE = PRPF4_TPM %>% filter(type == "MANE") %>% inner_join(Shared_PTC_genes, by = c("ensembl_gene_id" = "geneID")) #Pull just the MANE isoforms from the shared PTC genes
# PRPF4_shared_TPM = PRPF4_TPM %>% filter(type == "Novel") %>% full_join(PRPF4_sPTC) %>% full_join(PRPF4_sMANE) #Combine them together
# PRPF4_shared_sum = PRPF4_shared_TPM %>% group_by(type) %>% summarise(total_TPM = sum(kdTPM), count = n()) %>% mutate(sample = "PRPF4")
# #UPF1
# UPF1_sPTC = UPF1_TPM %>% inner_join(Shared_PTC_iso) #Filter to just the shared PTC isoforms
# UPF1_sMANE = UPF1_TPM %>% filter(type == "MANE") %>% inner_join(Shared_PTC_genes, by = c("ensembl_gene_id" = "geneID")) #Pull just the MANE isoforms from the shared PTC genes
# UPF1_shared_TPM = UPF1_TPM %>% filter(type == "Novel") %>% full_join(UPF1_sPTC) %>% full_join(UPF1_sMANE) #Combine them together
# UPF1_shared_sum = UPF1_shared_TPM %>% group_by(type) %>% summarise(total_TPM = sum(kdTPM), count = n()) %>% mutate(sample = "UPF1")
# #EIF4A3
# EIF4A3_sPTC = EIF4A3_TPM %>% inner_join(Shared_PTC_iso) #Filter to just the shared PTC isoforms
# EIF4A3_sMANE = EIF4A3_TPM %>% filter(type == "MANE") %>% inner_join(Shared_PTC_genes, by = c("ensembl_gene_id" = "geneID")) #Pull just the MANE isoforms from the shared PTC genes
# EIF4A3_shared_TPM = EIF4A3_TPM %>% filter(type == "Novel") %>% full_join(EIF4A3_sPTC) %>% full_join(EIF4A3_sMANE) #Combine them together
# EIF4A3_shared_sum = EIF4A3_shared_TPM %>% group_by(type) %>% summarise(total_TPM = sum(kdTPM), count = n()) %>% mutate(sample = "EIF4A3")
# 
# combined_shared = EFTUD2_shared_sum %>% full_join(AQR_shared_sum) %>% full_join(SF3B1_shared_sum) %>% full_join(PRPF4_shared_sum) %>% 
#   full_join(UPF1_shared_sum) %>% full_join(EIF4A3_shared_sum)
# 
# ####Make a graph of the TPMs of the shared PTC isoforms####
# shared_TPMplot = ggplot(data = combined_shared)
# shared_TPMplot = shared_TPMplot + geom_col(aes(y = factor(sample, level = c("EIF4A3","UPF1","AQR","SF3B1","EFTUD2","PRPF4")), 
#                                  x = total_TPM, fill = type), 
#                              position = position_dodge(),
#                              orientation = "y")
# shared_TPMplot
# ggsave("shared_PTC_isoform_TPM.pdf",
#        plot = shared_TPMplot,
#        device = "pdf",
#        width = 8,
#        height = 8,
#        units = "in",
#        dpi = 300)
# 
# ####Determine if the shared PTC list all have MANE isoforms####
# listMarts()
# ensembl <- useMart("ensembl")
# ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
# attributes <- listAttributes(ensembl)
# listgenes <- getBM(attributes = c("external_gene_name","ensembl_gene_id","ensembl_transcript_id_version","transcript_mane_select"),
#                    filters = "ensembl_gene_id",
#                    values = Shared_PTC$ensembl_gene_id.x,
#                    mart = ensembl)
# shared_PTC_sum = Shared_PTC %>% summarise(count = n_distinct(ensembl_gene_id.x)) %>% mutate(type = "shared_PTC") #101 genes
# shared_PTC_sum
# listgenes_sum = listgenes %>% summarise(count = n_distinct(ensembl_gene_id)) %>% mutate(type = "enembl_genes") #101 genes
# listgenes_sum
# mane_sum = listgenes %>% summarise(count = n_distinct(transcript_mane_select)) %>% mutate(type = "MANE")#101 genes
# mane_sum
# unique_gene_summary = shared_PTC_sum %>% full_join(listgenes_sum) %>% full_join(mane_sum)

####Analyize the TPM of all novel transcripts ####
##Import Datasets##
#AQR
AQR_novel_TPM <- read_csv("AQR_novel_TPM.csv")
#EFTUD2
EFTUD2_novel_TPM <- read_csv("EFTUD2_novel_TPM.csv")
#SF3B1
SF3B1_novel_TPM <- read_csv("SF3B1_novel_TPM.csv")
#PRPF4
PRPF4_novel_TPM <- read_csv("PRPF4_novel_TPM.csv")
#UPF1
UPF1_novel_TPM <- read_csv("UPF1_novel_TPM.csv")
#EIF4A3
EIF4A3_novel_TPM <- read_csv("EIF4A3_novel_TPM.csv")
#PUM1
# PUM1_novel_TPM <- read_csv("PUM1_novel_TPM.csv")
#ACTN4
#ACTN4_novel_TPM <- read_csv("ACTN4_novel_TPM.csv")
#SNRNP200
SNRNP200_novel_TPM <- read_csv("SNRNP200_novel_TPM.csv")
#PRPF3
PRPF3_novel_TPM <- read_csv("PRPF3_novel_TPM.csv")
#U2AF1
U2AF1_novel_TPM <- read_csv("U2AF1_novel_TPM.csv")
#CDC40
CDC40_novel_TPM <- read_csv("CDC40_novel_TPM.csv")
#SF3A3
SF3A3_novel_TPM <- read_csv("SF3A3_novel_TPM.csv")
#RBM22
RBM22_novel_TPM <- read_csv("RBM22_novel_TPM.csv")
#AGO1
AGO1_novel_TPM <- read_csv("AGO1_novel_TPM.csv")
#SF3B3
SF3B3_novel_TPM <- read_csv("SF3B3_novel_TPM.csv")
#SF3A1
SF3A1_novel_TPM <- read_csv("SF3A1_novel_TPM.csv")
#SNRPC
SNRPC_novel_TPM <- read_csv("SNRPC_novel_TPM.csv")
#SNRNP70
SNRNP70_novel_TPM <- read_csv("SNRNP70_novel_TPM.csv")
#GNB2L1
GNB2L1_novel_TPM <- read_csv("GNB2L1_novel_TPM.csv")
#CDC5L
CDC5L_novel_TPM <- read_csv("CDC5L_novel_TPM.csv")
#MAGOH
MAGOH_novel_TPM <- read_csv("MAGOH_novel_TPM.csv")
#HNRNPK
HNRNPK_novel_TPM <- read_csv("HNRNPK_novel_TPM.csv")


#Combine the datasets
full_novel_TPM = AQR_novel_TPM %>% full_join(EFTUD2_novel_TPM) %>% full_join(PRPF4_novel_TPM) %>% 
  full_join(SF3B1_novel_TPM) %>% full_join(EIF4A3_novel_TPM) %>% full_join(UPF1_novel_TPM) %>%
  full_join(SNRNP200_novel_TPM)  %>% full_join(PRPF3_novel_TPM) %>% 
  full_join(U2AF1_novel_TPM) %>% full_join(CDC40_novel_TPM) %>% full_join(SF3A3_novel_TPM) %>% full_join(RBM22_novel_TPM) %>% 
  full_join(AGO1_novel_TPM) %>% full_join(SF3B3_novel_TPM) %>% full_join(SF3A1_novel_TPM) %>% full_join(SNRNP70_novel_TPM) %>% 
  full_join(SNRPC_novel_TPM) %>% full_join(GNB2L1_novel_TPM) %>% full_join(CDC5L_novel_TPM) %>% full_join(HNRNPK_novel_TPM) %>% 
  full_join(MAGOH_novel_TPM) %>% 
  rename(Switch_Cons = switchConsequencesGene) #ACTN4 and PUM1 are not included



####Create a range plot for the TPM of different conditions####
EFTUD2_TPM_only = EFTUD2_TPM %>% dplyr::select(1,6:7,9:15) #8330 isoforms
AQR_TPM_only = AQR_TPM %>% dplyr::select(1,6:7,9:15) #8111 isoforms
PRPF4_TPM_only = PRPF4_TPM %>% dplyr::select(1,6:7,9:15) #7810 isoforms
SF3B1_TPM_only = SF3B1_TPM %>% dplyr::select(1,6:7,9:15) #8339 isoforms
EIF4A3_TPM_only = EIF4A3_TPM %>% dplyr::select(1,6:7,9:15) #8558
UPF1_TPM_only = UPF1_TPM %>% dplyr::select(1,6:7,9:15) #8245 isoforms
SNRNP200_TPM_only = SNRNP200_TPM %>% dplyr::select(1,6:7,9:15) #7591 isoforms
# ACTN4_TPM_only = ACTN4_TPM %>% dplyr::select(1,6:7,9:15)#8064 isoforms
# PUM1_TPM_only = PUM1_TPM %>% dplyr::select(1,6:7,9:15) #7759 isoforms
SF3A3_TPM_only = SF3A3_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
U2AF1_TPM_only = U2AF1_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
PRPF3_TPM_only = PRPF3_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
CDC40_TPM_only = CDC40_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
RBM22_TPM_only = RBM22_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
AGO1_TPM_only = AGO1_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
SF3B3_TPM_only = SF3B3_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
SF3A1_TPM_only = SF3A1_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
SNRPC_TPM_only = SNRPC_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
SNRNP70_TPM_only = SNRNP70_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
GNB2L1_TPM_only = GNB2L1_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
CDC5L_TPM_only = CDC5L_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
MAGOH_TPM_only = MAGOH_TPM %>% dplyr::select(1,6:7,9:15) # isoforms
HNRNPK_TPM_only = HNRNPK_TPM %>% dplyr::select(1,6:7,9:15) # isoforms

all_TPM = EFTUD2_TPM_only %>% full_join(AQR_TPM_only) %>% full_join(PRPF4_TPM_only) %>% full_join(SF3B1_TPM_only) %>% 
  full_join(EIF4A3_TPM_only) %>% full_join(UPF1_TPM_only) %>% full_join(SNRNP200_TPM_only) %>%  full_join(SF3A3_TPM_only) %>%
  full_join(U2AF1_TPM_only) %>% full_join(PRPF3_TPM_only) %>% 
  full_join(CDC40_TPM_only) %>% full_join(RBM22_TPM_only) %>% full_join(AGO1_TPM_only) %>% full_join(SF3A1_TPM_only) %>%
  full_join(SF3B3_TPM_only) %>% full_join(SNRPC_TPM_only) %>% full_join(SNRNP70_TPM_only) %>% full_join(GNB2L1_TPM_only) %>% 
  full_join(HNRNPK_TPM_only) %>% full_join(CDC5L_TPM_only) %>% full_join(MAGOH_TPM_only)
all_TPM = all_TPM %>% pivot_longer(2:3,names_to = "TPM_type", values_to = "TPM")
all_TPM = all_TPM %>% mutate(log_10 = log10(TPM))
all_TPM = all_TPM %>% mutate(Sample = fct_relevel(Sample,"SNRPC","GNB2L1","SNRNP70","PRPF3","SNRNP200","HNRNPK","SF3A3","PRPF4","CDC40",
                                                  "EFTUD2","CDC5L","U2AF1","MAGOH","RBM22","SF3B1","SF3A1","SF3B3","AQR","EIF4A3",
                                                  "UPF1","AGO1"))

##Making all of the plots##
range_colors = c("kdTPM" = "#E75A7C", "wtTPM" = "#FCDFA6",
                 "TRUE" = "#17BEBB", "FALSE" = "#E5B6F6")
#PTC plot
all_PTC_range = ggplot(data = dplyr::filter(all_TPM,type == "PTC"))
all_PTC_range = all_PTC_range + geom_violin(aes(x = Sample,
                                            y = log_10, fill = TPM_type),
                                            adjust = 0.5,
                                            draw_quantiles = 0.5) +
  scale_fill_manual(values = range_colors, labels = c("kdTPM" = "KD", "wtTPM" = "WT"))
all_PTC_range = all_PTC_range + labs(title = "TPM of PTC+ isoforms",
                                     fill = "Sample Type",
                                     x = "KD",
                                     y = "Log(10) TPM")
all_PTC_range = all_PTC_range + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
all_PTC_range
ggsave("PTC_range.pdf",
       plot = all_PTC_range,
       dpi = 300,
       device = pdf,
       width = 10,
       height = 8,
       units = "in")

all_PTC_range_sum = all_TPM %>% group_by(Sample,TPM_type) %>% filter(type == "PTC") %>%
  summarise(count = n(), totalTPM = sum(TPM))
all_PTC_range_sum = all_PTC_range_sum %>% mutate(Sample = fct_relevel(Sample,"SNRPC","GNB2L1","SNRNP70","PRPF3","SNRNP200",
                                                                      "HNRNPK","SF3A3","PRPF4","CDC40",
                                                                      "EFTUD2","CDC5L","U2AF1","MAGOH","RBM22","SF3B1",
                                                                      "SF3A1","SF3B3","AQR","EIF4A3",
                                                                      "UPF1","AGO1"))
all_ptc_bar = ggplot(data = all_PTC_range_sum)
all_ptc_bar = all_ptc_bar + geom_col(aes(x = totalTPM, y = Sample, fill = TPM_type),
                                     position = position_dodge(),
                                     color = "black") + 
  scale_fill_manual(values = range_colors, labels = c("kdTPM" = "KD", "wtTPM" = "WT"))
all_ptc_bar = all_ptc_bar + labs(title = "TPM of PTC+ isoforms",
                                fill = "Sample Type",
                                x = "Total TPM",
                                y = "Sample")
all_ptc_bar = all_ptc_bar + theme_bw() + 
  theme(panel.grid.major.y = element_blank())
all_ptc_bar
ggsave("PTC_TPM.pdf",
       plot = all_ptc_bar,
       dpi = 300,
       width = 10,
       height = 4,
       units = "in")

#All MANE transcripts plot
all_MANE_range = ggplot(data = dplyr::filter(all_TPM,type == "MANE"))
all_MANE_range = all_MANE_range + geom_violin(aes(x = Sample, y = log_10, fill = TPM_type),
                                            adjust = 0.5,
                                            draw_quantiles = 0.5) +
  scale_fill_manual(values = range_colors, labels = c("kdTPM" = "KD", "wtTPM" = "WT"))

all_MANE_range = all_MANE_range + labs(title = "TPM of MANE isoforms",
                                     fill = "Sample Type",
                                     x = "KD",
                                     y = "Log(10) TPM")
all_MANE_range = all_MANE_range + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
all_MANE_range
ggsave("MANE_range.pdf",
       plot = all_MANE_range,
       dpi = 300,
       device = pdf,
       width = 10,
       height = 8,
       units = "in")

all_MANE_range_sum = all_TPM %>% group_by(Sample,TPM_type) %>%
  filter(type == "MANE") %>%
  summarise(count = n(), totalTPM = sum(TPM)) %>% 
  mutate(Sample = fct_relevel(Sample,"SNRPC","GNB2L1","SNRNP70","PRPF3","SNRNP200","HNRNPK","SF3A3","PRPF4","CDC40",
                              "EFTUD2","CDC5L","U2AF1","MAGOH","RBM22","SF3B1","SF3A1","SF3B3","AQR","EIF4A3",
                              "UPF1","AGO1"))
all_MANE_bar = ggplot(data = all_MANE_range_sum)
all_MANE_bar = all_MANE_bar + geom_col(aes(x = totalTPM, y = Sample, fill = TPM_type),
                                     position = position_dodge(),
                                     color = "black") + 
  scale_fill_manual(values = range_colors, labels = c("kdTPM" = "KD", "wtTPM" = "WT"))
all_MANE_bar = all_MANE_bar + labs(title = "TPM of MANE+ isoforms",
                                 fill = "Sample Type",
                                 x = "Total TPM",
                                 y = "Sample")
all_MANE_bar = all_MANE_bar + theme_bw() + 
  theme(panel.grid.major.y = element_blank())
all_MANE_bar
ggsave("MANE_TPM.pdf",
       plot = all_MANE_bar,
       dpi = 300,
       width = 10,
       height = 4,
       units = "in")

#MANE PTC plots
ptc_MANE_range = ggplot(data = dplyr::filter(all_TPM,type == "PTC_MANE")) 
ptc_MANE_range = ptc_MANE_range + geom_violin(aes(x = Sample, y = log_10, fill = TPM_type),
                                              adjust = 0.5,
                                              draw_quantiles = 0.5) +
  scale_fill_manual(values = range_colors, labels = c("kdTPM" = "KD", "wtTPM" = "WT"))

ptc_MANE_range = ptc_MANE_range + labs(title = "TPM of MANE isoforms of PTC+ genes",
                                       fill = "Sample Type",
                                       x = "KD",
                                       y = "Log(10) TPM")
ptc_MANE_range = ptc_MANE_range + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
ptc_MANE_range
ggsave("ptc_MANE_range.pdf",
       plot = ptc_MANE_range,
       dpi = 300,
       device = pdf,
       width = 12,
       height = 10,
       units = "in")

ptc_MANE_range_sum = all_TPM %>% group_by(Sample,TPM_type) %>%
  filter(type == "MANE") %>%
  summarise(count = n(), totalTPM = sum(TPM)) %>% 
  mutate(Sample = fct_relevel(Sample,"SNRPC","GNB2L1","SNRNP70","PRPF3","SNRNP200","HNRNPK","SF3A3","PRPF4","CDC40",
                              "EFTUD2","CDC5L","U2AF1","MAGOH","RBM22","SF3B1","SF3A1","SF3B3","AQR","EIF4A3",
                              "UPF1","AGO1"))
ptc_MANE_bar = ggplot(data = ptc_MANE_range_sum)
ptc_MANE_bar = ptc_MANE_bar + geom_col(aes(x = totalTPM, y = Sample, fill = TPM_type),
                                       position = position_dodge(),
                                       color = "black") + 
  scale_fill_manual(values = range_colors, labels = c("kdTPM" = "KD", "wtTPM" = "WT"))
ptc_MANE_bar = ptc_MANE_bar + labs(title = "TPM of MANE in PTC+ genes isoforms",
                                   fill = "Sample Type",
                                   x = "Total TPM",
                                   y = "Sample")
ptc_MANE_bar = ptc_MANE_bar + theme_bw() + 
  theme(panel.grid.major.y = element_blank())
ptc_MANE_bar
ggsave("ptc_MANE_TPM.pdf",
       plot = ptc_MANE_bar,
       dpi = 300,
       width = 10,
       height = 6,
       units = "in")

#Novel TPM plot
full_novel_TPM = full_novel_TPM %>% mutate(sample = fct_relevel(sample,"SNRPC","GNB2L1","SNRNP70","PRPF3","SNRNP200","HNRNPK","SF3A3","PRPF4","CDC40",
                                                                "EFTUD2","CDC5L","U2AF1","MAGOH","RBM22","SF3B1","SF3A1","SF3B3","AQR","EIF4A3",
                                                                "UPF1","AGO1"))
full_novel_TPM = full_novel_TPM %>% mutate(log10 = log10(kdTPM))
novel_switch_range = ggplot(data = full_novel_TPM) 
novel_switch_range = novel_switch_range + geom_violin(aes(x = sample, y = log10, fill = Switch_Cons),
                                                      adjust = 0.5,
                                                      draw_quantiles = 0.5) +
  scale_fill_manual(values = range_colors, labels = c("TRUE" = "NMD", "FALSE" = "No effect on NMD"))
novel_switch_range = novel_switch_range + labs(title = "TPM of Novel isoforms in KD samples",
                          fill = "Isoform Switch Consequence",
                          x = "Sample",
                          y = "Log(10) TPM")
novel_switch_range = novel_switch_range + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
novel_switch_range
ggsave("novel_Switch_range.pdf",
       plot = novel_switch_range,
       dpi = 300,
       width = 10,
       height = 8,
       units = "in")

novel_tpm_summary = full_novel_TPM %>% group_by(sample,Switch_Cons) %>%
  summarise(count = n(), totalTPM_KD = sum(kdTPM)) %>% 
  mutate(Sample = fct_relevel(sample,"SNRPC","GNB2L1","SNRNP70","PRPF3","SNRNP200","HNRNPK","SF3A3","PRPF4","CDC40",
                              "EFTUD2","CDC5L","U2AF1","MAGOH","RBM22","SF3B1","SF3A1","SF3B3","AQR","EIF4A3",
                              "UPF1","AGO1"))
nmd_tpm_sum = novel_tpm_summary %>% filter(Switch_Cons == TRUE)
write_csv(nmd_tpm_sum,"novel_NMD_TPM.csv")
write_csv(nmd_tpm_sum, "C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/novel_NMD_TPM.csv")
novel_tpm_sum = ggplot(data = novel_tpm_summary)
novel_tpm_sum = novel_tpm_sum + geom_col(aes(x = totalTPM_KD, y = Sample, fill = Switch_Cons),
                                       position = position_dodge(),
                                       color = "black") + 
  scale_fill_manual(values = range_colors, labels = c("TRUE" = "NMD", "FALSE" = "No effect on NMD"))
novel_tpm_sum = novel_tpm_sum + labs(title = "TPM of Novel isoforms in KD samples",
                                   fill = "Isoform Switch Consequence",
                                   x = "Total TPM",
                                   y = "Sample")
novel_tpm_sum = novel_tpm_sum + theme_bw() + 
  theme(panel.grid.major.y = element_blank())
novel_tpm_sum
ggsave("Novel_TPM_bar.pdf",
       plot = novel_tpm_sum,
       dpi = 300,
       width = 10,
       height = 4,
       units = "in")

novel_tpm_facet = ggplot(data = novel_tpm_summary)
novel_tpm_facet = novel_tpm_facet + geom_col(aes(y = totalTPM, x = Sample, fill = Switch_Cons),
                                         position = position_dodge(),
                                         color = "black") + 
  scale_fill_manual(values = range_colors,
                    labels = c("TRUE" = "NMD", "FALSE" = "No effect on NMD")) +
  facet_wrap(vars(Switch_Cons))
novel_tpm_facet = novel_tpm_facet + labs(title = "TPM of Novel isoforms in KD samples",
                                     fill = "Isoform Switch Consequence",
                                     x = "Sample",
                                     y = "Total TPM") 
novel_tpm_facet = novel_tpm_facet + theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  guides(x = guide_axis(angle = 45))
novel_tpm_facet
ggsave("Novel_TPM_facet.pdf",
       plot = novel_tpm_facet,
       dpi = 300,
       width = 12,
       height = 6,
       units = "in")

novel_prop = novel_tpm_summary %>% ungroup() %>% dplyr::select(2,4:5) %>% 
  pivot_wider(names_from = Switch_Cons, values_from = totalTPM) %>% 
  rename(No_NMD = "FALSE", NMD = "TRUE") %>% 
  mutate(TotalTPM = No_NMD + NMD) %>% mutate(NMDprop = NMD/TotalTPM)

novel_proportion = ggplot(data = novel_prop)
novel_proportion = novel_proportion + geom_col(aes(x = Sample, y = NMDprop),
                                               fill = "#17BEBB",
                                               color = "black") +
  theme_bw() +
  guides(x = guide_axis(angle = 45)) +
  labs(title = "Proportion of Novel isoforms targeted to NMD",
       x = "Sample",
       y = "Proportion of novel isoform TPM") +
  theme(panel.grid.major.x = element_blank())
novel_proportion
ggsave("Novel_TPM_proportion.pdf",
       plot = novel_proportion,
       dpi = 300,
       width = 12,
       height = 6,
       units = "in")
# Poster for Rustbelt RNA meeting 2023

RRM_novel_TPM = full_novel_TPM %>% filter(sample == "SNRPC" | sample == "AQR" | sample == "SF3B1" | sample == "EFTUD2" |
                                            sample == "CDC40" | sample == "AGO1" | sample == "UPF1" | sample == "EIF4A3")
RRM_novel_TPM = RRM_novel_TPM %>% mutate(sample = fct_relevel(sample,"AGO1","UPF1","EIF4A3","AQR",
                                                              "CDC40","SF3B1",
                                                              "EFTUD2"))
RRM_switch_range = ggplot(data = RRM_novel_TPM) 
RRM_switch_range = RRM_switch_range + geom_violin(aes(x = sample, y = log10, fill = Switch_Cons),
                                                      adjust = 0.5,
                                                      draw_quantiles = 0.5) +
  scale_fill_manual(values = range_colors, labels = c("TRUE" = "NMD", "FALSE" = "No effect on NMD"))
RRM_switch_range = RRM_switch_range + labs(title = "TPM of Novel isoforms in KD samples",
                                               fill = "Isoform Switch Consequence",
                                               x = "Sample",
                                               y = "Log(10) TPM")
RRM_switch_range = RRM_switch_range + theme_bw() + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
RRM_switch_range
ggsave("RRM_Switch_range.pdf",
       plot = RRM_switch_range,
       dpi = 300,
       width = 10,
       height = 8,
       units = "in")

RRM_tpm_summary = RRM_novel_TPM %>% group_by(sample,Switch_Cons) %>%
  summarise(count = n(), totalTPM = sum(kdTPM)) %>% 
  mutate(Sample = fct_relevel(sample,"AGO1","UPF1","EIF4A3","AQR",
                              "CDC40","SF3B1",
                              "EFTUD2"))
RRM_tpm_sum = ggplot(data = RRM_tpm_summary)
RRM_tpm_sum = RRM_tpm_sum + geom_col(aes(y = totalTPM, x = Sample, fill = Switch_Cons),
                                         position = position_dodge(),
                                         color = "black") + 
  scale_fill_manual(values = range_colors, labels = c("TRUE" = "NMD", "FALSE" = "No effect on NMD"))
RRM_tpm_sum = RRM_tpm_sum + labs(title = "TPM of Novel isoforms in KD samples",
                                     fill = "Isoform Switch Consequence",
                                     x = "KD",
                                     y = "Total TPM")
RRM_tpm_sum = RRM_tpm_sum + theme_bw() + 
  theme(panel.grid.major.y = element_blank())
RRM_tpm_sum
ggsave("RRM_TPM_bar.pdf",
       plot = RRM_tpm_sum,
       dpi = 300,
       width = 12,
       height = 10,
       units = "in")
RRM_facet = RRM_tpm_sum + facet_wrap(vars(Switch_Cons))
RRM_facet
ggsave("RRM_TPM_bar_facet.pdf",
       plot = RRM_facet,
       dpi = 300,
       width = 12,
       height = 10,
       units = "in")

#### Look at length of the novel vs MANE isoforms ####
main_fig_list = c("UPF1","EIF4A3",
                  "AQR","EFTUD2","SF3B1",
                  "SNRPC")
all_GOI = c("UPF1","EIF4A3","MAGOH",
            "AQR","RBM22","CDC5L",
            "EFTUD2","SNRNP200",
            "SF3A1","SF3B1","SF3B3","SF3A3","U2AF1",
            "CDC40","HNRNPK",
            "PRPF3","PRPF4","PRPF31_IPSC",
            "GNB2L1",
            "SNRNP70","SNRPC")
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart = ensembl)
attributes = listAttributes(mart = ensembl)
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_ORF_analysis"), read_csv(paste0(i,"_switch_ORF_length.csv")))
  assign(paste0(i,"_annotated_ORF"), getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","external_gene_name",
                                                          "cds_length"),
                                           filters = "ensembl_gene_id",
                                           values = eval(parse(text = paste0(i,"_ORF_analysis$gene_id"))),
                                           mart = ensembl))
  assign(paste0(i,"_both_ORF"), eval(parse(text = paste0(i,"_ORF_analysis"))) %>% 
           mutate(Sample = i) %>% left_join(eval(parse(text=paste0(i,"_annotated_ORF"))), 
                                            by = c("isoform_id" = "ensembl_transcript_id",
                                                   "gene_id" = "ensembl_gene_id",
                                                   "gene_name" = "external_gene_name")))
  assign(paste0(i,"_MANE"), getBM(attributes = c("ensembl_transcript_id","transcript_mane_select"),
                                           filters = "ensembl_transcript_id",
                                           values = eval(parse(text = paste0(i,"_both_ORF$isoform_id"))),
                                           mart = ensembl))
  assign(paste0(i,"_both_ORF"), eval(parse(text = paste0(i,"_both_ORF"))) %>% 
           left_join(eval(parse(text=paste0(i,"_MANE"))), 
                     by = c("isoform_id" = "ensembl_transcript_id")))
  assign(paste0(i,"_both_ORF"), eval(parse(text = paste0(i,"_both_ORF"))) %>% 
           mutate(Novel = if_else(IF2 == 0, "TRUE", "FALSE")))
}
all_both_ORF = UPF1_both_ORF %>% full_join(EIF4A3_both_ORF) %>% full_join(MAGOH_both_ORF) %>%
full_join(AQR_both_ORF) %>% full_join(RBM22_both_ORF) %>% full_join(CDC5L_both_ORF) %>%
full_join(EFTUD2_both_ORF) %>% full_join(SNRNP200_both_ORF) %>%
full_join(SF3A1_both_ORF) %>% full_join(SF3B1_both_ORF) %>% full_join(SF3B3_both_ORF) %>%full_join(SF3A3_both_ORF) %>%full_join(U2AF1_both_ORF) %>%
full_join(CDC40_both_ORF) %>% full_join(HNRNPK_both_ORF) %>%
  full_join(PRPF3_both_ORF) %>% full_join(PRPF4_both_ORF) %>% full_join(PRPF31_IPSC_both_ORF) %>%
  full_join(GNB2L1_both_ORF) %>%
  full_join(SNRNP70_both_ORF) %>% full_join(SNRPC_both_ORF)
#Maybe import the isoform TPM list and then get their length?

ORF_length_plot = ggplot(data = all_both_ORF %>% filter(Sample %in% main_fig_list))
ORF_length_plot = ORF_length_plot +
  geom_boxplot(aes(x = Novel, y = log10(orfTransciptLength),
                   fill = Novel),
               alpha = 0.5,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,) +
  facet_wrap(facets = vars(factor(Sample, levels = main_fig_list)),
             ncol = 3)
ORF_length_plot
ggsave("switcing_length_plot.pdf",
       plot = ORF_length_plot,
       device = pdf,
       height = 10,
       width = 10,
       units = "in",
       dpi = 300)

## Look only at unannotated novel isoforms
novel_ORF = all_both_ORF %>% filter(Novel == "TRUE")
UAnovel_ORF = novel_ORF %>% mutate(unannotated = str_detect(isoform_id,"MSTRG"))%>% 
  filter(unannotated == "TRUE") %>% 
  select(1:26)
unnanotated_genes = UAnovel_ORF %>% distinct(gene_id)
write_csv(unnanotated_genes,"genes_with_novel_isoforms.csv")
unnanotated_genes_ORF = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","cds_length"),
                              filters = "ensembl_gene_id",
                              values = unnanotated_genes,
                              mart = ensembl)
unnanotated_genes_MANE = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","transcript_mane_select"),
                              filters = "ensembl_gene_id",
                              values = unnanotated_genes,
                              mart = ensembl)
unnanotated_genes_ORF = unnanotated_genes_ORF %>% full_join(unnanotated_genes_MANE) %>% 
  rename(MANE = transcript_mane_select,
         MANE_transcript_id = ensembl_transcript_id,
         MANE_length = cds_length) %>% 
  filter(str_detect(MANE,"NM") == T)
UAnovel_ORF = UAnovel_ORF %>% left_join(unnanotated_genes_ORF, by = c("gene_id" = "ensembl_gene_id"))
UAnovel_ORF = UAnovel_ORF %>% pivot_longer(cols = c(orfTransciptLength,MANE_length),names_to = "transcript",
                                           values_to = "Length") %>% 
  filter(!is.na(Length)) %>% 
  mutate(log_length = log10(Length))

UA_sum = UAnovel_ORF %>% group_by(Sample,transcript) %>%
  summarise(n = n(),
            genes = n_distinct(gene_id),
            Med = median(Length, na.rm = T),
            quart = quantile(Length, 0.5,na.rm = T),
            log_Med = median(log_length, na.rm = T),
            log_quart = quantile(log_length, 0.5,na.rm = T)) #calculate the summary stats for labeling the plots
UA_sum = UA_sum %>% group_by(Sample) %>% mutate(Total = sum(n),
                                                Total_genes = sum(genes))

## Calculates the P-value for each plot
UA_res = tibble(Sample = as.character(), P = as.numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_UA"), UAnovel_ORF %>% filter(Sample == i))
  assign(paste0(i,"_res"), wilcox.test(Length ~ transcript, data = eval(parse(text = paste0(i,"_UA"))),
                                       exact = FALSE, alternative = "greater")) #use greater because we expect MANE to be longer than PE
  UA_res = UA_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_res$p.value"))))
}
UA_sum = UA_sum %>% left_join(UA_res, by = "Sample")
#Plot the lengths
UA_boxplot = ggplot(data = UAnovel_ORF)
UA_boxplot = UA_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log_length,
                   fill = transcript),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = c("MANE_length" = "#c270b7",
                               "orfTransciptLength" = "#eda058"),
                    labels = c("MANE_length" = "MANE",
                               "orfTransciptLength" = "Novel")) +
  scale_color_manual(values = c("MANE_length" = "#c270b7",
                                "orfTransciptLength" = "#eda058"),
                     labels = c("MANE_length" = "MANE",
                                "orfTransciptLength" = "Novel"))+
  geom_label(data = UA_sum,
             aes(x = Sample,
                 y = 4.5,
                 label = paste0("p-value = ",signif(P,digits = 3))))+
  geom_label(data = UA_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = log_Med,
                 color = transcript,
                 label = round(log_Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8)) +
  geom_label(data = UA_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = 1.75,
                 color = transcript,
                 label = paste0("n =",n)),
             show.legend = F,
             position = position_dodge2(width = 0.8)) +
  theme_bw() +
  labs(x = "Sample",
       y = "Log10(Transcript Length)",
       fill = "Transcript Type") +
  coord_cartesian(y = c(1.5,4.75))
UA_boxplot
ggsave("Novel_transcript_length_boxplot.pdf",
       device = pdf,
       plot = UA_boxplot,
       width = 35,
       height = 10,
       units = "in",
       dpi = 300)

#Plot the lengths
UA_main_fig = ggplot(data = UAnovel_ORF %>% filter(Sample %in% main_fig_list))
UA_main_fig = UA_main_fig + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log_length,
                   fill = transcript),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = c("MANE_length" = "#c270b7",
                               "orfTransciptLength" = "#eda058"),
                    labels = c("MANE_length" = "MANE",
                               "orfTransciptLength" = "Novel")) +
  scale_color_manual(values = c("MANE_length" = "#c270b7",
                                "orfTransciptLength" = "#eda058"),
                     labels = c("MANE_length" = "MANE",
                                "orfTransciptLength" = "Novel"))+
  geom_label(data = UA_sum %>% filter(Sample %in% main_fig_list),
             aes(x = Sample,
                 y = 4.5,
                 label = paste0("p = ",signif(P,digits = 3))),
             size = 5,
             alpha = 0.75)+
  geom_label(data = UA_sum %>% filter(Sample %in% main_fig_list),
             aes(x = factor(Sample, levels = all_GOI),
                 y = log_Med,
                 color = transcript,
                 label = round(log_Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label(data = UA_sum %>% filter(Sample %in% main_fig_list),
             aes(x = factor(Sample, levels = all_GOI),
                 y = 1.75,
                 color = transcript,
                 label = paste0("n =",n)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5,
             alpha = 0.75) +
  theme_bw() +
  labs(x = "Sample",
       y = "Log10(Transcript Length)",
       fill = "Isoform Type") +
  coord_cartesian(y = c(1.5,4.75))
UA_main_fig
ggsave("Novel_transcript_length__main_fig.pdf",
       device = pdf,
       plot = UA_main_fig,
       width = 16,
       height = 10,
       units = "in",
       dpi = 300)
