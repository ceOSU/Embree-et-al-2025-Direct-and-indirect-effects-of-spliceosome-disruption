#Compare the abundance of transcripts in different KDs
#And make the CDF plots into density plots
#Written by: Caleb Embree
#Written on: 1/10/2024
#Modified on: 3/62024
#Changelog: 2/1/2024- Added density plots for PTC+- 
#2/12/2024- Create TPM scatterplot
#2/14/2024- Made additional scatterplots
#2/16/2024: Added in the no NMD data for the density plots
#3/4/2024: Added in the side by side plots
#3/5/2024: Continued working on the side by side boxplots
#3/6/2024: Modified the median TPM plot to add in more datasets
#3/14/2024: Created the three different types of graphs for comparing log2FC
#3/18/2024: Created the new plots for PTC+/-, the PE, and the disease samples
#4/23/2024: Added in the PRPF8 mutant dataset
#4/29/2024: Updated the PE data with the new PE list combining both Saltzman et al and Thomas et al list of genes
#5/20/204: Made full versions of all the plots, and added in PRPF8 and PRPF6 ENCODE datasets
#5/22/2024: Made PE datasests using just the Saltzman paper
#5/23/2024: Added in the no NMD boxplot and the PTC heatmap
#5/28/2024: Changed and added more options for the no-NMD datasets
#6/3/2024: Added in the AS/noAS isoform boxplots and the noAS+noNMD boxplots
#6/5/2024-6/8/2024: Added in other AS subplots
#6/10/2024: Added in the disease AS dataplot
#6/12-15/2024: Added the no PTC NMD plots
#6/24/2024: Made the plot for the PTC+,MANE,and Other isoforms from the same genes. 
#6/27/2024: Made sure the disease datasets are using the right data. 
#6/28/2024: Filtered the PTC and other isoform plot to just TSL 1/2
#7/3/2024: Removed some excess code and fixed some of the TPM code. Also added in "main figure" plots with reduced number of samples

library(readr)
library(readxl)
library(tidyverse)
library(MetBrewer)
library(ggh4x)
library(ggrepel)
library(ggpmisc)
library(biomaRt)
library(patchwork)
library(pheatmap)
library(ggVennDiagram)
library(MoMAColors)

all_GOI = c("UPF1","EIF4A3","MAGOH",
            "AQR","RBM22","CDC5L",
            "EFTUD2","SNRNP200","PRPF8","PRPF6",
            "SF3A1","SF3B1","SF3B3","SF3A3","U2AF1",
            "CDC40",
            "PRPF3","PRPF4",
            "GNB2L1",
            "SNRNP70","SNRPC")
Pres_KD = c("UPF1","EIF4A3","MAGOH","AQR","RBM22","EFTUD2","SNRNP200","SF3B1","SF3B3","CDC40","SNRPC","SNRNP70")

#Set up biomart
listMarts()
ensembl <- useMart("ensembl",host = "https://feb2023.archive.ensembl.org") #use version 109 of ensembl
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)

#### Histogram of Log2FC from WT filtered cells####
hist_all_GOI = c("UPF1","AGO1","AQR","SF3B3","EIF4A3")
for (i in hist_all_GOI) {
  print(i)
  assign(paste0(i,"_alltrans"), read_csv(paste0(i,"_WTfilt_alltrans.csv")))
  assign(paste0(i,"_alltrans"), eval(parse(text = paste0(i,"_alltrans"))) %>% mutate(Sample = paste0(i)))
}
Hist_combined = UPF1_alltrans %>% full_join(AGO1_alltrans) %>% full_join(AQR_alltrans) %>% full_join(SF3B3_alltrans) %>% 
  full_join(EIF4A3_alltrans)
Hist_combined = Hist_combined %>% mutate(Sig = if_else(padj <=0.1, "S","NS"))

Hist = ggplot(data = Hist_combined)
Hist = Hist + geom_freqpoly(aes(x = log2FoldChange,y = after_stat(density), color = Sig),
                             binwidth = 0.1,
                            linewidth= 1) +
  facet_wrap(vars(Sample)) +
  scale_color_manual(values = c("S" = "blue", "NS" = "grey"),
                     labels = c("S" = "Significant", "NS"="Not Significant"))
Hist = Hist + theme_bw() +
  theme(strip.background = element_blank()) +
  labs(caption = "Filtered from WT only",
       color = "P adjust value",
       title = "Frequency of Log2FC",
       x = "Isoform Log2FC")
Hist
ggsave("WTfilt_Log2FC_freq.pdf",
       plot = Hist,
       units = "in",
       width = 12,
       height = 8,
       device = pdf,
       dpi = 300)
#Import the double filter dataset
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_full_alltrans"), read_csv(paste0(i,"_DBfilt_alltrans.csv")))
  assign(paste0(i,"_full_alltrans"), eval(parse(text = paste0(i,"_DBfilt"))) %>% mutate(Sample = paste0(i)))
}
DBfilt_combined = UPF1_DBfilt %>% full_join(AGO1_DBfilt) %>% full_join(AQR_DBfilt) %>% full_join(SF3B3_DBfilt) %>% 
  full_join(EIF4A3_DBfilt)
DBfilt_combined = DBfilt_combined %>% mutate(Sig = if_else(padj <=0.1, "S","NS"))

DB_hist = ggplot(data = DBfilt_combined)
DB_hist = DB_hist + geom_freqpoly(aes(x = log2FoldChange,y = after_stat(density), color = Sig),
                            binwidth = 0.1,
                            linewidth= 1) +
  facet_wrap(vars(Sample)) +
  scale_color_manual(values = c("S" = "blue", "NS" = "grey"),
                     labels = c("S" = "Significant", "NS"="Not Significant"))
DB_hist = DB_hist + theme_bw() +
  theme(strip.background = element_blank()) +
  labs(caption = "Filtered in both WT and KD",
       color = "P adjust value",
       title = "Frequency of Log2FC",
       x = "Isoform Log2FC")
DB_hist
ggsave("DBfilt_Log2FC_freq.pdf",
       plot = DB_hist,
       units = "in",
       width = 12,
       height = 8,
       device = pdf,
       dpi = 300)




#### Look at the Log2FC of PTC+- gene list####
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_PTC_alltrans"), read_csv(paste0(i,"_PTC_alltrans.csv")))
  assign(paste0(i,"_PTC_alltrans"), eval(parse(text = paste0(i,"_PTC_alltrans"))) %>% mutate(Sample = paste0(i)) %>% 
           dplyr::select(2:8,10:12)) 
  assign(paste0(i,"_full_alltrans_sig"), eval(parse(text = paste0(i,"_full_alltrans"))) %>% mutate(Sample = paste0(i)) %>% 
           filter(padj <= 0.05))
}

sig_combined = MAGOH_full_alltrans_sig %>% full_join(EIF4A3_full_alltrans_sig) %>% full_join(UPF1_full_alltrans_sig) %>% 
  full_join(RBM22_full_alltrans_sig) %>% full_join(AQR_full_alltrans_sig) %>% full_join(SNRNP200_full_alltrans_sig) %>% 
  full_join(EFTUD2_full_alltrans_sig) %>% full_join(SF3B1_full_alltrans_sig) %>% full_join(SF3B3_full_alltrans_sig) %>% 
  full_join(SNRPC_full_alltrans_sig) %>% full_join(SNRNP70_full_alltrans_sig) %>% full_join(PRPF8_full_alltrans_sig) %>%
  full_join(PRPF6_full_alltrans_sig) %>% full_join(CDC5L_full_alltrans_sig) %>% full_join(SF3A1_full_alltrans_sig) %>% 
  full_join(SF3A3_full_alltrans_sig) %>% full_join(U2AF1_full_alltrans_sig) %>% full_join(CDC40_full_alltrans_sig) %>% 
  full_join(PRPF3_full_alltrans_sig) %>% full_join(PRPF4_full_alltrans_sig) %>% full_join(GNB2L1_full_alltrans_sig) %>% 
  mutate(Type = "Sig")
PTC_combined = MAGOH_PTC_alltrans %>% full_join(EIF4A3_PTC_alltrans) %>% full_join(UPF1_PTC_alltrans) %>% 
  full_join(RBM22_PTC_alltrans) %>% full_join(AQR_PTC_alltrans) %>% full_join(SNRNP200_PTC_alltrans) %>% 
  full_join(EFTUD2_PTC_alltrans) %>% full_join(SF3B1_PTC_alltrans) %>% full_join(SF3B3_PTC_alltrans) %>% 
  full_join(SNRPC_PTC_alltrans) %>% full_join(SNRNP70_PTC_alltrans) %>% full_join(PRPF8_PTC_alltrans) %>% 
  full_join(PRPF6_PTC_alltrans) %>% full_join(CDC5L_PTC_alltrans) %>% full_join(SF3A1_PTC_alltrans) %>% 
  full_join(SF3A3_PTC_alltrans) %>% full_join(U2AF1_PTC_alltrans) %>% full_join(CDC40_PTC_alltrans) %>% 
  full_join(PRPF3_PTC_alltrans) %>% full_join(PRPF4_PTC_alltrans) %>% full_join(GNB2L1_PTC_alltrans) %>% 
  mutate(PTC = as.character(PTC),
         Type = "PTC")
sig_and_PTC = sig_combined %>% full_join(PTC_combined) %>%
  mutate(Subcomplex = case_when(Sample == "UPF1" | Sample == "EIF4A3" | Sample == "MAGOH" ~ "NMD",
                                Sample == "SNRNP200" | Sample == "EFTUD2" | Sample == "PRPF8" | Sample == "PRPF6" ~ "U5",
                                Sample == "AQR" | Sample == "RBM22" | Sample == "CDC5L" ~ "NTC related",
                                Sample == "SF3B1" | Sample == "SF3B3" | Sample == "SF3A1" | Sample == "SF3A3" | Sample == "U2AF1" ~ "U2",
                                Sample == "SNRPC" | Sample == "SNRNP70" ~ "U1",
                                Sample == "CDC40" ~ "Bact",
                                Sample == "PRPF3" | Sample == "PRPF4" ~ "U4/U6",
                                Sample == "GNB2L1" ~ "C2"))
PTC_combined = PTC_combined %>%   mutate(Subcomplex = case_when(Sample == "UPF1" | Sample == "EIF4A3" | Sample == "MAGOH" ~ "NMD",
                                                                Sample == "SNRNP200" | Sample == "EFTUD2" | Sample == "PRPF8" | Sample == "PRPF6" ~ "U5",
                                                                Sample == "AQR" | Sample == "RBM22" | Sample == "CDC5L" ~ "NTC related",
                                                                Sample == "SF3B1" | Sample == "SF3B3" | Sample == "SF3A1" | Sample == "SF3A3" | Sample == "U2AF1"  ~ "U2",
                                                                Sample == "SNRPC" | Sample == "SNRNP70" ~ "U1",
                                                                Sample == "CDC40" ~ "Bact",
                                                                Sample == "PRPF3" | Sample == "PRPF4" ~ "U4/U6",
                                                                Sample == "GNB2L1" ~ "C2"))
PTC_sum = PTC_combined %>% group_by(PTC,Sample) %>%
  summarise(med = median(log2FoldChange),
            mean = mean(log2FoldChange),
            n = n()) %>% 
  mutate(PTC = str_replace(PTC,"FALSE","MANE"),
         PTC = str_replace(PTC, "TRUE", "PTC"))
all_PTCres = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_PTC_res"),
         wilcox.test(log2FoldChange ~ PTC, data = PTC_combined %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  all_PTCres = all_PTCres %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_PTC_res$p.value"))))
}
PTC_sum = PTC_sum %>% left_join(all_PTCres)


colors = c("MANE" = "#663171","PTC" = "#ea7428", "Sig" = "grey",
           "TRUE"= "#ea7428", "FALSE" = "#663171",
           "None" = "#B4436C")


PTC_boxplot = ggplot(data = PTC_combined %>% filter(Sample %in% Pres_KD))
PTC_boxplot = PTC_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = Pres_KD),
                   y = log2FoldChange,
                   fill = PTC),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1)+
  scale_fill_manual(values = colors, labels = c("FALSE" = "MANE",
                                                   "TRUE" = "PTC")) +
  scale_color_manual(values = colors, labels = c("FALSE" = "MANE",
                                                    "TRUE" = "PTC"))+
  geom_text(data = PTC_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = Pres_KD),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             color = "black",
            show.legend = F)+
  geom_label(data = PTC_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = Pres_KD),
                 y = med,
                 color = PTC,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = PTC_sum %>% filter(Sample %in% Pres_KD),
                  aes(x = factor(Sample, levels = Pres_KD),
                      y = -3,
                      color = PTC,
                      label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Transcript")+
  coord_cartesian(y = c(-4,4))
PTC_boxplot
#Last saved and modified on 10/7/2024
ggsave("PTC_boxplot_FC.pdf",
       device = pdf,
       plot = PTC_boxplot,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)


full_PTC_boxplot = ggplot(data = PTC_combined)
full_PTC_boxplot = full_PTC_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = PTC),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1)+
  scale_fill_manual(values = colors, labels = c("FALSE" = "MANE",
                                                "TRUE" = "PTC")) +
  scale_color_manual(values = colors, labels = c("FALSE" = "MANE",
                                                 "TRUE" = "PTC"))+
  geom_text(data = PTC_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             color = "black")+
  geom_label(data = PTC_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = PTC,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = PTC_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3,
                 color = PTC,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Transcript")+
  coord_cartesian(y = c(-4,4))
full_PTC_boxplot
#Last saved and modified on 10/7/2024
ggsave("PTC_boxplot_full.pdf",
       device = pdf,
       plot = full_PTC_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)


#import and manipulate data (Manu's no NMD gene list")####
noNMD_transcripts <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/noNMD_isoforms.csv")
all_noNMD = tibble(baseMean = as.numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noNMD_alltrans"), eval(parse(text = paste0(i,"_DBfilt"))) %>% 
           inner_join(noNMD_transcripts, by = c("ENST.ID" = "ensembl_transcript_id"))) %>%
    select(2:9)
  all_noNMD = all_noNMD %>% full_join(eval(parse(text = paste0(i,"_noNMD_alltrans"))))
}


noNMD_sum = all_noNMD %>% group_by(isoform,Sample) %>%
  summarise(med = median(log2FoldChange),
            mean = mean(log2FoldChange),
            n = n())
noNMD_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noNMD_res"),
         wilcox.test(log2FoldChange ~ isoform, data = all_noNMD %>% filter(Sample == i),
                     exact = FALSE, alternative = "less")) 
  noNMD_res = noNMD_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_noNMD_res$p.value"))))
}
noNMD_sum = noNMD_sum %>% left_join(noNMD_res)


noNMD_colors = c("MANE" = "#93B497","Other" = "#577854")


noNMD_boxplot = ggplot(data = all_noNMD)
noNMD_boxplot = noNMD_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1)+
  scale_fill_manual(values = noNMD_colors) +
  scale_color_manual(values = noNMD_colors)+
  geom_text(data = noNMD_sum %>% filter(isoform == "MANE"),
             aes(x = Sample,
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             color = "black")+
  geom_label(data = noNMD_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = noNMD_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3,
                 color = isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Transcript",
       title = "No NMD biotype")+
  coord_cartesian(y = c(-4,4))
noNMD_boxplot
# Last saved and modified on 10/7/2024
ggsave("noNMD_boxplot_FC.pdf",
       device = pdf,
       plot = noNMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

mainFig_noNMD_boxplot = ggplot(data = all_noNMD %>% filter(Sample %in% Pres_KD))
mainFig_noNMD_boxplot = mainFig_noNMD_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1)+
  scale_fill_manual(values = noNMD_colors) +
  scale_color_manual(values = noNMD_colors)+
  geom_text(data = noNMD_sum %>% filter(Sample %in% Pres_KD),
             aes(x = Sample,
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             color = "black")+
  geom_label(data = noNMD_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = noNMD_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3,
                 color = isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Transcript",
       title = "No NMD biotype")+
  coord_cartesian(y = c(-4,4))
mainFig_noNMD_boxplot
# Last saved and modified on 10/7/2024
ggsave("mainFig_noNMD_boxplot_FC.pdf",
       device = pdf,
       plot = mainFig_noNMD_boxplot,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)

#Repeat with the full list of transcripts from Manu's noNMD genes
MS_noNMD_transcripts <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/MS_noNMD_transcripts.csv")
MS_noNMD_transcripts = MS_noNMD_transcripts %>% mutate(isoform = if_else(!is.na(transcript_mane_select),
                                                                         "MANE",
                                                                         "Other"))
all_MSnoNMD = tibble(baseMean = as.numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_MSnoNMD_alltrans"), eval(parse(text = paste0(i,"_DBfilt"))) %>% 
           inner_join(MS_noNMD_transcripts, by = c("ENST.ID" = "ensembl_transcript_id"))) %>%
    select(2:9)
  all_MSnoNMD = all_MSnoNMD %>% full_join(eval(parse(text = paste0(i,"_MSnoNMD_alltrans"))))
}


MSnoNMD_sum = all_MSnoNMD %>% group_by(isoform,Sample) %>%
  summarise(med = median(log2FoldChange),
            mean = mean(log2FoldChange),
            n = n())
MSnoNMD_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_MSnoNMD_res"),
         wilcox.test(log2FoldChange ~ isoform, data = all_MSnoNMD %>% filter(Sample == i),
                     exact = FALSE, alternative = "less")) 
  MSnoNMD_res = MSnoNMD_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_MSnoNMD_res$p.value"))))
}
MSnoNMD_sum = MSnoNMD_sum %>% left_join(MSnoNMD_res)



MSnoNMD_boxplot = ggplot(data = all_MSnoNMD)
MSnoNMD_boxplot = MSnoNMD_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1)+
  scale_fill_manual(values = noNMD_colors) +
  scale_color_manual(values = noNMD_colors)+
  geom_text(data = MSnoNMD_sum %>% filter(isoform == "MANE"),
             aes(x = Sample,
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7)+
  geom_label(data = MSnoNMD_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = MSnoNMD_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3.5,
                 color = isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Transcript",
       title = "MS full list")+
  coord_cartesian(y = c(-4,4))
MSnoNMD_boxplot
# Last saved and modified on 10/7/2024
ggsave("MSnoNMD_boxplot_FC.pdf",
       device = pdf,
       plot = MSnoNMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)


##Look at the non-NMD transcripts from the Lykke-Anderson genes and development paper####
LA_nonNMD = read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/Lykke_non_NMD.csv") #244 isoforms, 132 genes
LA_nonNMD_MANE = LA_nonNMD_MANE %>% 
  mutate(isoform = if_else(str_detect(transcript_mane_select,"NM"),
                           "MANE",
                           "Other"))
LA_summary = LA_nonNMD_MANE %>% group_by(isoform) %>% summarise(n = n(),
                                                                genes = n_distinct(ensembl_gene_id),
                                                                isoforms = n_distinct(ensembl_transcript_id))

all_LAnoNMD = tibble(baseMean = as.numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_LAnoNMD_alltrans"), eval(parse(text = paste0(i,"_DBfilt"))) %>% 
           inner_join(LA_nonNMD_MANE, by = c("ENST.ID" = "ensembl_transcript_id"))) %>%
    select(2:9)
  all_LAnoNMD = all_LAnoNMD %>% full_join(eval(parse(text = paste0(i,"_LAnoNMD_alltrans"))))
}


LAnoNMD_sum = all_LAnoNMD %>% group_by(isoform,Sample) %>%
  summarise(med = median(log2FoldChange),
            mean = mean(log2FoldChange),
            n = n())
LAnoNMD_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_LAnoNMD_res"),
         wilcox.test(log2FoldChange ~ isoform, data = all_LAnoNMD %>% filter(Sample == i),
                     exact = FALSE, alternative = "less")) 
  LAnoNMD_res = LAnoNMD_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_LAnoNMD_res$p.value"))))
}
LAnoNMD_sum = LAnoNMD_sum %>% left_join(LAnoNMD_res)

LA_colors = c("MANE" = "#B27092","Other" = "#512D38")

LAnoNMD_boxplot = ggplot(data = all_LAnoNMD)
LAnoNMD_boxplot = LAnoNMD_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1)+
  scale_fill_manual(values = LA_colors) +
  scale_color_manual(values = LA_colors)+
  geom_text(data = LAnoNMD_sum,
             aes(x = Sample,
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 8)+
  geom_label(data = LAnoNMD_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = LAnoNMD_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3.5,
                 color = isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 8,
             alpha = 0.75) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Transcript",
       title = "Lykke-Anderson 2014 non-NMD Reference",
       caption = "Other = Appris Principal 1-5")+
  coord_cartesian(y = c(-4,4))
LAnoNMD_boxplot
# Last saved and modified on 5/10/2024
ggsave("LAnoNMD_boxplot_FC.pdf",
       device = pdf,
       plot = LAnoNMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)



#### Make side by side TPM box plot ####
all_TPM = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_TPM"),
         read_csv(paste0("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/",
                         i,
                         "_TPM.csv"))) #Import the TPM dataset
  assign(paste0(i,"_TPM"),
         eval(parse(text = paste0(i,"_TPM"))) %>% mutate(Sample = i) %>% 
           select(1,6:14))
  all_TPM = all_TPM %>% full_join(eval(parse(text = paste0(i,"_TPM"))))
}

all_TPM = all_TPM %>% pivot_longer(c(kdTPM,wtTPM),names_to = "Treatment",values_to = "TPM") %>% 
  mutate(Treatment = str_remove(Treatment,"TPM"))


SBS_color = c("wt" = "#EABA0B",
              "kd" = "#ea7428",
              "TRUE" = "#990D35",
              "FALSE" = "#6ACDD8")

PTC_box = ggplot(data = all_TPM %>% filter(type == "PTC"))
PTC_box = PTC_box + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = TPM,
                   fill = factor(Treatment, levels = c("wt","kd"))),
               outlier.shape = 21,
               outlier.alpha = 0.75,
               outlier.colour = NA,
               position =  position_dodge2(width = 0.9),
               width = 0.8) +
  scale_fill_manual(values = SBS_color,
                    labels = c("kd" = "KD",
                                "wt" = "WT"),
                    breaks = c("wt","kd"))+
  scale_color_manual(values = SBS_color,
                    labels = c("kd" = "KD",
                               "wt" = "WT"),
                    breaks = c("wt","KDmean")) +
  theme_bw() +
  labs(x = "Sample",
       y = "TPM",
       fill = "Treatment") +
  coord_cartesian(y = c(0,15))
PTC_box
ggsave("PTC_TPM_boxplot.pdf",
       plot = PTC_box,
       device = pdf,
       units = "in",
       height = 10,
       width = 12,
       dpi = 300)

novel_box = ggplot(data = all_TPM %>% filter(type == "Novel" & Treatment == "kd"))
novel_box = novel_box +
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
               y = TPM,
               fill = NMD_Causing),
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
           position = position_dodge2(width = 0.9),
           width = 0.8) +
  scale_fill_manual(values = SBS_color, labels = c("TRUE" = "NMD",
                                                   "FALSE" = "non-NMD")) +
  scale_color_manual(values = SBS_color, labels = c("TRUE" = "NMD",
                                                   "FALSE" = "non-NMD")) +
  theme_bw() +
  labs(x = "Sample",
       y = "TPM",
       fill = "Novel Isoform") +
   coord_cartesian(y = c(0,15))
novel_box
ggsave("Novel_TPM_boxplot.pdf",
       plot = novel_box,
       device = pdf,
       units = "in",
       width = 12,
       height = 10,
       dpi = 300)

SBS = (PTC_box | novel_box) + plot_layout(guides = "collect",
                                          axes = "collect")
SBS
ggsave("Side-by-side_TPM.pdf",
       plot = SBS,
       device = pdf,
       height = 10,
       width = 40,
       units = "in",
       dpi = 300)



#Look at the total TPM of each category
all_TPM_restructure = all_TPM %>% mutate(NMD_Causing = str_replace(NMD_Causing,"TRUE","NMD")) %>% 
  unite(type,NMD_Causing,col = "Type",sep = "_",na.rm = TRUE) %>% 
  mutate(Type = str_replace(Type,"FALSE","Stable")) #Combines the two novel annotation columns into one

all_TPM_summary = all_TPM_restructure %>% group_by(Sample, Type, Treatment) %>% 
  summarise(med = median(TPM),
            total = sum(TPM))

TPM_summary_colors = c("MANE" = "#663171",
                       "Novel_Stable" = "#6ACDD8",
                       "Novel_NMD" = "#990D35",
                       "PTC" = "#ea7428")
total_TPM_PTC = ggplot(data = all_TPM_summary %>% filter(Treatment == "kd") %>%
                               filter(Type == "MANE" | Type == "PTC"))
total_TPM_PTC = total_TPM_PTC + geom_col(aes(x = factor(Sample, levels = all_GOI),
                                             y = total,
                                             fill = Type),
                                         position = "dodge") +
  scale_fill_manual(values = TPM_summary_colors) +
  labs(x = "Depletion",
       y = "Total TPM",
       title = "Annotated isoforms",
       fill = "Isoform") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust=1))+
  coord_cartesian(ylim = c(0,4000))
total_TPM_PTC
ggsave("total_TPM_annotated.pdf",
       plot = total_TPM_PTC,
       width = 16,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

total_TPM_novel = ggplot(data = all_TPM_summary %>% filter(Treatment == "kd") %>%
                         filter(Type == "Novel_NMD" | Type == "Novel_Stable"))
total_TPM_novel = total_TPM_novel + geom_col(aes(x = factor(Sample, levels = all_GOI),
                                               y = total,
                                               fill = factor(Type, levels = c("Novel_Stable","Novel_NMD"))),
                                             position = "dodge") +
  scale_fill_manual(values = TPM_summary_colors, labels = c("Novel_NMD" = "Novel NMD",
                                                             "Novel_Stable" = "Novel Stable"),
                    breaks = c("Novel_Stable","Novel_NMD")) +
  labs(x = "Depletion",
       y = "Total TPM",
       title = "Novel Isoforms",
       fill = "Isoform") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  coord_cartesian(ylim = c(0,4000))
total_TPM_novel
ggsave("total_TPM_novel.pdf",
       plot = total_TPM_novel,
       width = 16,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

total_TPM_SBS = (total_TPM_PTC | total_TPM_novel) +
  plot_layout(guides = "collect",
              axes = "collect") +
  labs(caption = "TPM from KD samples") 
total_TPM_SBS
ggsave("Total_TPM_SBS.pdf",
       plot = total_TPM_SBS,
       device = pdf,
       height = 10,
       width = 22,
       units = "in",
       dpi = 300)

#make the version just with the main figures
main_fig_PTC_box = ggplot(data = all_TPM %>% filter(type == "PTC" & Sample %in% Pres_KD))
main_fig_PTC_box = main_fig_PTC_box + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = TPM,
                   fill = factor(Treatment, levels = c("wt","kd"))),
               outlier.shape = 21,
               outlier.alpha = 0.75,
               outlier.colour = NA,
               position =  position_dodge2(width = 0.9),
               width = 0.8) +
  scale_fill_manual(values = SBS_color,
                    labels = c("kd" = "KD",
                               "wt" = "WT"),
                    breaks = c("wt","kd"))+
  scale_color_manual(values = SBS_color,
                     labels = c("kd" = "KD",
                                "wt" = "WT"),
                     breaks = c("wt","KDmean")) +
  theme_bw() +
  labs(x = "Sample",
       y = "TPM",
       fill = "Treatment") +
  coord_cartesian(y = c(0,15))
main_fig_PTC_box
ggsave("main_fig_PTC_TPM_boxplot.pdf",
       plot = main_fig_PTC_box,
       device = pdf,
       units = "in",
       height = 10,
       width = 8,
       dpi = 300)

main_fig_novel_box = ggplot(data = all_TPM %>% filter(type == "Novel" & Treatment == "kd" & Sample %in% Pres_KD))
main_fig_novel_box = main_fig_novel_box +
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = TPM,
                   fill = NMD_Causing),
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               position = position_dodge2(width = 0.9),
               width = 0.8) +
  scale_fill_manual(values = SBS_color, labels = c("TRUE" = "NMD",
                                                   "FALSE" = "non-NMD")) +
  scale_color_manual(values = SBS_color, labels = c("TRUE" = "NMD",
                                                    "FALSE" = "non-NMD")) +
  theme_bw() +
  labs(x = "Sample",
       y = "TPM",
       fill = "Novel Isoform") +
  coord_cartesian(y = c(0,15))
main_fig_novel_box
ggsave("main_fig_Novel_TPM_boxplot.pdf",
       plot = main_fig_novel_box,
       device = pdf,
       units = "in",
       width = 8,
       height = 10,
       dpi = 300)

main_fig_SBS = (main_fig_PTC_box | main_fig_novel_box) + plot_layout(guides = "collect",
                                                                     axes = "collect")
main_fig_SBS
ggsave("main_fig_SBS_TPM.pdf",
       plot = main_fig_SBS,
       device = pdf,
       height = 10,
       width = 18,
       units = "in",
       dpi = 300)


main_fig_total_TPM_PTC = ggplot(data = all_TPM_summary %>% filter(Treatment == "kd" & Sample %in% Pres_KD) %>%
                         filter(Type == "MANE" | Type == "PTC"))
main_fig_total_TPM_PTC = main_fig_total_TPM_PTC + geom_col(aes(x = factor(Sample, levels = all_GOI),
                                                               y = total,
                                                               fill = factor(Type,levels = c("MANE","PTC"))),
                                                           position = "dodge") +
  scale_fill_manual(values = TPM_summary_colors) +
  labs(x = "Depletion",
       y = "Total TPM",
       title = "Annotated isoforms",
       fill = "Isoform") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  coord_cartesian(ylim = c(0,4000))
main_fig_total_TPM_PTC
ggsave("main_fig_total_TPM_annotated.pdf",
       plot = main_fig_total_TPM_PTC,
       width = 10,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

main_fig_total_TPM_novel = ggplot(data = all_TPM_summary %>% filter(Treatment == "kd" & Sample %in% Pres_KD) %>%
                           filter(Type == "Novel_NMD" | Type == "Novel_Stable"))
main_fig_total_TPM_novel = main_fig_total_TPM_novel + geom_col(aes(x = factor(Sample, levels = all_GOI),
                                                                   y = total,
                                                                   fill = factor(Type,levels = c("Novel_Stable","Novel_NMD"))),
                                                               position = "dodge") +
  scale_fill_manual(values = TPM_summary_colors, labels = c("Novel_NMD" = "Novel NMD",
                                                            "Novel_Stable" = "Novel Stable")) +
  labs(x = "Depletion",
       y = "Total TPM",
       title = "Novel Isoforms",
       fill = "Isoform") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  coord_cartesian(ylim = c(0,4000))
main_fig_total_TPM_novel
ggsave("main_fig_total_TPM_novel.pdf",
       plot = main_fig_total_TPM_novel,
       width = 8,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

main_fig_total_TPM_SBS = (main_fig_total_TPM_PTC | main_fig_total_TPM_novel) +
  plot_layout(guides = "collect",
              axes = "collect") +
  labs(caption = "TPM from KD samples") 
main_fig_total_TPM_SBS
ggsave("main_fig_total_TPM_SBS.pdf",
       plot = main_fig_total_TPM_SBS,
       device = pdf,
       height = 10,
       width = 18,
       units = "in",
       dpi = 300)


#### Look at FC of PE NMD isoforms ####
#Combines all tables
MF_alltrans = tibble(Sample = character())
for (i in all_GOI) {
  MF_alltrans = MF_alltrans %>% full_join(eval(parse(text = paste0(i,"_full_alltrans"))))
}
write_csv(MF_alltrans, file = "comined_non_disease_alltrans.csv")

## Load PE data
PE_transcripts <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Poison_Exons/PE_NMD_transcripts_both_lists.csv")
PE_full_summary = PE_transcripts %>% group_by(isoform) %>% summarise(n = n(),
                                                                     genes = n_distinct(ensembl_gene_id),
                                                                     NMD = sum(str_detect(transcript_biotype,
                                                                                          "nonsense"))) #Both categories have 379 genes

## Filter the log2FC data to only transcripts in the PE/MANE list
PE_FC = MF_alltrans %>% filter(ENST.ID %in% PE_transcripts$ensembl_transcript_id)
PE_FC = PE_FC %>% left_join(PE_transcripts, by = c("ENST.ID" = "ensembl_transcript_id"))
PE_FC_sum = PE_FC %>% group_by(Sample,isoform) %>%
  summarise(n = n(),
            Med = median(log2FoldChange, na.rm = T),
            quart = quantile(log2FoldChange, 0.5)) #calculate the summary stats for labeling the plots
PE_FC_sum = PE_FC_sum %>% group_by(Sample) %>% mutate(Total = sum(n)) %>% ungroup()

## Calculates the P-value for each plot
PE_res = tibble(Sample = as.character(), P = as.numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_PE"), PE_FC %>% filter(Sample == i))
  assign(paste0(i,"_PE_res"), wilcox.test(log2FoldChange ~ isoform, data = eval(parse(text = paste0(i,"_PE"))),
                                          exact = FALSE, alternative = "less")) #use less because we expect MANE to be lower than PE
  PE_res = PE_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_PE_res$p.value"))))
}
PE_FC_sum = PE_FC_sum %>% left_join(PE_res, by = "Sample")
PE_pres_sum = PE_FC_sum %>% filter(Sample %in% Pres_KD)

PE_Colors = c("NMD" = "#e94220",
              "MANE" = "#5a439d")



PE_boxplot = ggplot(data = PE_FC %>% filter(Sample %in% Pres_KD))
PE_boxplot = PE_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = PE_Colors) +
  scale_color_manual(values = PE_Colors)+
  geom_text(data = PE_pres_sum %>% filter(Sample %in% Pres_KD & isoform == "NMD"),
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             color = "black",
            show.legend = F)+
  geom_label(data = PE_pres_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = Med,
                 color = isoform,
                 label = round(Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = PE_pres_sum %>% filter(Sample %in% Pres_KD),
                  aes(x = factor(Sample, levels = all_GOI),
                      y = -3,
                      color = isoform,
                      label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "isoform Type")+
  coord_cartesian(y = c(-4,4))
PE_boxplot
ggsave("PE_boxplot_FC.pdf",
       device = pdf,
       plot = PE_boxplot,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)

full_PE_boxplot = ggplot(data = PE_FC)
full_PE_boxplot = full_PE_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = PE_Colors) +
  scale_color_manual(values = PE_Colors)+
  geom_text(data = PE_FC_sum %>% filter(isoform == "NMD"),
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             color = "black")+
  geom_label(data = PE_FC_sum %>% filter(Sample %in% all_GOI),
             aes(x = factor(Sample, levels = all_GOI),
                 y = Med,
                 color = isoform,
                 label = round(Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = PE_FC_sum %>% filter(Sample %in% all_GOI),
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3,
                 color = isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "isoform Type")+
  coord_cartesian(y = c(-4,4))
full_PE_boxplot
ggsave("full_PE_boxplot_Full.pdf",
       device = pdf,
       plot = full_PE_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

##Look at the PE only from the Saltzman inclusion list##
#Load the data
Saltzman_PE = read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Poison_Exons/Saltzman_inclusion_NMD_transcripts.csv")
Saltzman_PE = Saltzman_PE %>% mutate(Isoform = "NMD")
Saltzman_MANE = getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id","external_gene_name",
                                     "transcript_biotype","transcript_mane_select"),
                      filters = "ensembl_gene_id",
                      values = Saltzman_PE$ensembl_gene_id,
                      mart = ensembl)
Saltzman_MANE = Saltzman_MANE %>% filter(str_detect(transcript_mane_select,"NM")) %>%  #filter to only MANE transcripts
  mutate(Isoform = "MANE")
Saltzman_PE = Saltzman_PE %>% full_join(Saltzman_MANE) #combine the tables

#Make the transcript table
Saltzman_alltrans = MF_alltrans %>% filter(ENST.ID %in% Saltzman_PE$ensembl_transcript_id) %>% 
  select(2:9)
Saltzman_alltrans = Saltzman_alltrans %>% left_join(Saltzman_PE, by = c("ENST.ID" = "ensembl_transcript_id"),
                                                    relationship = "many-to-many")
Saltzman_alltrans_sum = Saltzman_alltrans %>% group_by(Sample,Isoform) %>%
  summarise(n = n(),
            Med = median(log2FoldChange, na.rm = T),
            quart = quantile(log2FoldChange, 0.5)) #calculate the summary stats for labeling the plots
Saltzman_alltrans_sum = Saltzman_alltrans_sum %>% group_by(Sample) %>% mutate(Total = sum(n)) %>% ungroup()

## Calculates the P-value for each plot
Saltzman_PE_res = tibble(Sample = as.character(), P = as.numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_PE"), Saltzman_alltrans %>% filter(Sample == i))
  assign(paste0(i,"_PE_res"), wilcox.test(log2FoldChange ~ Isoform, data = eval(parse(text = paste0(i,"_PE"))),
                                          exact = FALSE, alternative = "less")) #use less because we expect MANE to be lower than PE
  Saltzman_PE_res = Saltzman_PE_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_PE_res$p.value"))))
}
Saltzman_alltrans_sum = Saltzman_alltrans_sum %>% left_join(Saltzman_PE_res, by = "Sample")



#Plot the poison exons
Saltzman_PE_boxplot = ggplot(data = Saltzman_alltrans)
Saltzman_PE_boxplot = Saltzman_PE_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = Isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = PE_Colors) +
  scale_color_manual(values = PE_Colors)+
  geom_text(data = Saltzman_alltrans_sum %>% filter(Isoform == "NMD"),
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
            color = "black")+
  geom_label(data = Saltzman_alltrans_sum %>% filter(Sample %in% all_GOI),
             aes(x = factor(Sample, levels = all_GOI),
                 y = Med,
                 color = Isoform,
                 label = round(Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = Saltzman_alltrans_sum %>% filter(Sample %in% all_GOI),
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type")+
  coord_cartesian(y = c(-4,4))
Saltzman_PE_boxplot
ggsave("Saltzman_PE_boxplot_Full.pdf",
       device = pdf,
       plot = Saltzman_PE_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

Saltzman_PE_boxplot_short  = ggplot(data = Saltzman_alltrans %>% filter(Sample %in% Pres_KD))
Saltzman_PE_boxplot_short = Saltzman_PE_boxplot_short + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = Isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = PE_Colors) +
  scale_color_manual(values = PE_Colors)+
  geom_text(data = Saltzman_alltrans_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
            color = "black")+
  geom_label(data = Saltzman_alltrans_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = Med,
                 color = Isoform,
                 label = round(Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = Saltzman_alltrans_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type")+
  coord_cartesian(y = c(-4,4))
Saltzman_PE_boxplot_short
ggsave("Saltzman_PE_boxplot_main_fig.pdf",
       device = pdf,
       plot = Saltzman_PE_boxplot_short,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)



#### Make FC plot for AS genes vs not####
Novel_splicing = read_csv("all_novel_splice_events.csv")
AF_AS_combined = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_annotated_AS"),
         read_csv(paste0(i,"_annotated_AS.csv")))
  
  assign(paste0(i,"_sig_AS_genes"),
         Novel_splicing %>%
           filter(Sample == i) %>% 
           full_join(eval(parse(text = paste0(i,"_annotated_AS")))) %>% 
           distinct(GeneID))
  
  assign(paste0(i,"_alltrans_genes"),
         getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","transcript_mane_select",
                              "transcript_biotype","external_gene_name"),
               filters = "ensembl_transcript_id",
               values = eval(parse(text = paste0(i,"_full_alltrans$ENST.ID"))),
               mart = ensembl))
  
  assign(paste0(i,"_AS_alltrans"),
         eval(parse(text = paste0(i,"_full_alltrans"))) %>% 
           full_join(eval(parse(text = paste0(i,"_alltrans_genes"))),
                     by = c("ENST.ID" = "ensembl_transcript_id")) %>% 
           mutate(AS_gene = if_else(ensembl_gene_id %in% eval(parse(text = paste0(i,"_sig_AS_genes$GeneID"))),
                                    TRUE,
                                    FALSE)))
  AF_AS_combined = AF_AS_combined %>% full_join(eval(parse(text = paste0(i,"_AS_alltrans"))))
}


AF_AS_MANE = AF_AS_combined %>% filter(str_detect(transcript_mane_select,"NM"))
AS_MANE_sum = AF_AS_MANE %>% group_by(Sample,AS_gene) %>% summarise(n = n(),
                                                                    med = median(log2FoldChange),
                                                                    genes = n_distinct(ensembl_gene_id))
all_ASres = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_AS_res"),
         wilcox.test(log2FoldChange ~ AS_gene, data = AF_AS_MANE %>% filter(Sample == i),
                     exact = FALSE, alternative = "greater"))
  all_ASres = all_ASres %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_AS_res$p.value"))))
}
AS_MANE_sum = AS_MANE_sum %>% left_join(all_ASres, by = "Sample")

AS_colors = c("TRUE" = "#6DC5B9",
              "FALSE" = "#DE4D86")



AS_boxplot = ggplot(data = AF_AS_MANE %>% filter(Sample %in% Pres_KD))
AS_boxplot = AS_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = AS_gene),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = AS_colors, labels = c("FALSE" = "no AS",
                                                   "TRUE" = "AS genes")) +
  geom_text(data = AS_MANE_sum %>% filter(Sample %in% Pres_KD),
            aes(x = factor(Sample,
                           levels = all_GOI),
                y = 3,
                label = signif(P,digits = 3)),
             size = 7,
            show.legend = F,
            color = "black")+
  geom_label(data = AS_MANE_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = AS_gene,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = AS_MANE_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -3,
                 color = AS_gene,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Gene Type")+
  coord_cartesian(y = c(-4,4))
AS_boxplot
ggsave("AS_boxplot_FC_plot.pdf",
       device = pdf,
       plot = AS_boxplot,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)

full_AS_boxplot = ggplot(data = AF_AS_MANE)
full_AS_boxplot = full_AS_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = AS_gene),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = AS_colors, labels = c("FALSE" = "no AS",
                                                   "TRUE" = "AS genes")) +
  geom_text(data = AS_MANE_sum,
            aes(x = factor(Sample,
                           levels = all_GOI),
                y = 3,
                label = signif(P,digits = 3)),
            size = 7,
            color = "black",
            show.legend = F)+
  geom_label(data = AS_MANE_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = AS_gene,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = AS_MANE_sum,
                  aes(x = factor(Sample,
                                 levels = all_GOI),
                      y = -3,
                      color = AS_gene,
                      label = n),
                  show.legend = F,
                  position = position_dodge2(width = 0.9),
                  size = 7,
                  direction = "y",
                  segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Gene Type")+
  coord_cartesian(y = c(-4,4))
full_AS_boxplot
ggsave("full_AS_boxplot_FC_plot.pdf",
       device = pdf,
       plot = full_AS_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)


#Look at MANE vs other transcripts for AS-NMD genes
AS_only = AF_AS_combined %>% filter(AS_gene == "TRUE") %>%
  mutate(Isoform = if_else(str_detect(transcript_mane_select,"NM"),
                           "MANE",
                           "Other")) %>% 
  filter(!is.na(ensembl_gene_id))
No_AS_only = AF_AS_combined %>% filter(AS_gene == "FALSE") %>%
  mutate(Isoform = if_else(str_detect(transcript_mane_select,"NM"),
                           "MANE",
                           "Other")) %>% 
  filter(!is.na(ensembl_gene_id))
AS_only_sum = AS_only %>% group_by(Sample,Isoform) %>% summarise(n = n(),
                                                                    med = median(log2FoldChange),
                                                                    genes = n_distinct(ensembl_gene_id))
no_AS_sum = No_AS_only %>% group_by(Sample,Isoform) %>% summarise(n = n(),
                                                                 med = median(log2FoldChange),
                                                                 genes = n_distinct(ensembl_gene_id))
AS_only_res = tibble(Sample = character(),P = numeric())
no_AS_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_ASonly_res"),
         wilcox.test(log2FoldChange ~ Isoform, data = AS_only %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  assign(paste0(i,"_noAS_res"),
         wilcox.test(log2FoldChange ~ Isoform, data = No_AS_only %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  AS_only_res = AS_only_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_ASonly_res$p.value"))))
  no_AS_res = no_AS_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_noAS_res$p.value"))))
}
AS_only_sum = AS_only_sum %>% left_join(AS_only_res, by = "Sample")
no_AS_sum = no_AS_sum %>% left_join(no_AS_res, by = "Sample")


noAS_isoform_colors = c("MANE" = "#DE4D86","Other" = "#9B1C4F")
onlyAS_isoform_colors = c("MANE" = "#6DC5B9","Other" = "#2F756B")

no_AS_boxplot = ggplot(data = No_AS_only)
no_AS_boxplot = no_AS_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = Isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = noAS_isoform_colors) +
  scale_color_manual(values = noAS_isoform_colors) +
  geom_text(data = no_AS_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
            color = "black")+
  geom_label(data = no_AS_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = no_AS_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -3,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type",
       title = "Genes with no AS events")+
  coord_cartesian(y = c(-4,4))
no_AS_boxplot
ggsave("no_AS_boxplot_FC_plot.pdf",
       device = pdf,
       plot = no_AS_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

no_AS_main_boxplot = ggplot(data = No_AS_only %>% filter(Sample %in% Pres_KD))
no_AS_main_boxplot = no_AS_main_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = Isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = noAS_isoform_colors) +
  scale_color_manual(values = noAS_isoform_colors) +
  geom_text(data = no_AS_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
            color = "black")+
  geom_label(data = no_AS_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = no_AS_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -3,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type",
       title = "Genes with no AS events")+
  coord_cartesian(y = c(-4,4))
no_AS_main_boxplot
ggsave("no_AS_main_boxplot_FC_plot.pdf",
       device = pdf,
       plot = no_AS_main_boxplot,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)

only_AS_boxplot = ggplot(data = AS_only)
only_AS_boxplot = only_AS_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = Isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = onlyAS_isoform_colors) +
  scale_color_manual(values = onlyAS_isoform_colors) +
  geom_text(data = AS_only_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
            color = "black")+
  geom_label(data = AS_only_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = AS_only_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -3,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type",
       title = "Genes with AS events")+
  coord_cartesian(y = c(-4,4))
only_AS_boxplot
ggsave("only_AS_boxplot_FC_plot.pdf",
       device = pdf,
       plot = only_AS_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

only_AS_main_boxplot = ggplot(data = AS_only %>% filter(Sample %in% Pres_KD))
only_AS_main_boxplot = only_AS_main_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = Isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = onlyAS_isoform_colors) +
  scale_color_manual(values = onlyAS_isoform_colors) +
  geom_text(data = AS_only_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
            color = "black")+
  geom_label(data = AS_only_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = AS_only_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type",
       title = "Genes with AS events")+
  coord_cartesian(y = c(-4,4))
only_AS_main_boxplot
ggsave("only_AS_main_boxplot_FC_plot.pdf",
       device = pdf,
       plot = only_AS_main_boxplot,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)


####Compare no-AS genes to no-NMD genes####
MS_noNMD_genes = MS_noNMD_transcripts %>% distinct(ensembl_gene_id)
all_noAS_noNMD = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noNMD_noAS"),
         No_AS_only %>% filter(Sample == i) %>% 
           filter(ensembl_gene_id %in% MS_noNMD_genes$ensembl_gene_id)) 
  all_noAS_noNMD = all_noAS_noNMD %>% full_join(eval(parse(text = paste0(i,"_noNMD_noAS"))))
}
noAS_noNMD_sum = all_noAS_noNMD %>% group_by(Sample, Isoform) %>% 
  summarise(n = n(),
            genes = n_distinct(ensembl_gene_id),
            med = median(log2FoldChange))
no_ASorNMD_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noASorNMD_res"),
         wilcox.test(log2FoldChange ~ Isoform, data = all_noAS_noNMD %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  no_ASorNMD_res = no_ASorNMD_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_noASorNMD_res$p.value"))))
}
noAS_noNMD_sum = noAS_noNMD_sum %>% left_join(no_ASorNMD_res, by = "Sample")

noAS_noNMD_colors = c("MANE" = "#F69F8E","Other" = "#ED6145")
noAS_noNMD_boxplot = ggplot(data = all_noAS_noNMD)
noAS_noNMD_boxplot = noAS_noNMD_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = Isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = noAS_noNMD_colors) +
  scale_color_manual(values = noAS_noNMD_colors) +
  geom_text(data = noAS_noNMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
            color = "black")+
  geom_label(data = noAS_noNMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = noAS_noNMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type",
       title = "Genes with no AS events or NMD isoforms")+
  coord_cartesian(y = c(-4,4))
noAS_noNMD_boxplot
ggsave("noAS_noNMD_boxplot_FC_plot.pdf",
       device = pdf,
       plot = noAS_noNMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

noAS_noNMD_main_boxplot = ggplot(data = all_noAS_noNMD %>% filter(Sample %in% Pres_KD))
noAS_noNMD_main_boxplot = noAS_noNMD_main_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = Isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = noAS_noNMD_colors) +
  scale_color_manual(values = noAS_noNMD_colors) +
  geom_text(data = noAS_noNMD_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
            color = "black")+
  geom_label(data = noAS_noNMD_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = noAS_noNMD_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type",
       title = "Genes with no AS events or NMD isoforms")+
  coord_cartesian(y = c(-4,4))
noAS_noNMD_main_boxplot
ggsave("noAS_noNMD_main_boxplot_FC_plot.pdf",
       device = pdf,
       plot = noAS_noNMD_main_boxplot,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)

#Compare AS noNMD datasets
all_AS_noNMD = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noNMD_AS"),
         AS_only %>% filter(Sample == i) %>% 
           filter(ensembl_gene_id %in% MS_noNMD_genes$ensembl_gene_id)) 
  all_AS_noNMD = all_AS_noNMD %>% full_join(eval(parse(text = paste0(i,"_noNMD_AS"))))
}
AS_noNMD_sum = all_AS_noNMD %>% group_by(Sample, Isoform) %>% 
  summarise(n = n(),
            genes = n_distinct(ensembl_gene_id),
            med = median(log2FoldChange))
AS_noNMD_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_AS_noNMD_res"),
         wilcox.test(log2FoldChange ~ Isoform, data = all_AS_noNMD %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  AS_noNMD_res = AS_noNMD_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_AS_noNMD_res$p.value"))))
}
AS_noNMD_sum = AS_noNMD_sum %>% left_join(AS_noNMD_res, by = "Sample")

AS_noNMD_colors = c("MANE" = "#86C1B9","Other" = "#3A9286")
AS_noNMD_boxplot = ggplot(data = all_AS_noNMD)
AS_noNMD_boxplot = AS_noNMD_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = Isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = AS_noNMD_colors) +
  scale_color_manual(values = AS_noNMD_colors) +
  geom_label(data = AS_noNMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
             color = "black")+
  geom_label(data = AS_noNMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label_repel(data = AS_noNMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type",
       title = "Genes with AS events but no NMD isoforms")+
  coord_cartesian(y = c(-4,4))
AS_noNMD_boxplot
ggsave("AS_noNMD_boxplot_FC_plot.pdf",
       device = pdf,
       plot = AS_noNMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)


####Look at PTC+ transcripts that are not AS####
sPTC_list = read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/Stringent_PTC_MANE_CE.csv")
sPTC_list = sPTC_list %>% select(2:4)
noAS_NMD = No_AS_only %>% inner_join(sPTC_list,by = c("ENST.ID" = "transID"))
noAS_NMD_sum = noAS_NMD %>% group_by(Sample,PTC) %>% 
  summarise(n = n(),
            med = median(log2FoldChange))
noAS_NMD_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noAS_NMD_res"),
         wilcox.test(log2FoldChange ~ Isoform, data = noAS_NMD %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  noAS_NMD_res = noAS_NMD_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_noAS_NMD_res$p.value"))))
}
noAS_NMD_sum = noAS_NMD_sum %>% left_join(noAS_NMD_res, by = "Sample")

noAS_NMD_colors = c("FALSE" = "#663171","TRUE" = "#ea7428")
noAS_NMD_boxplot = ggplot(data = noAS_NMD)
noAS_NMD_boxplot = noAS_NMD_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = PTC),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = noAS_NMD_colors, labels = c("FALSE"="MANE",
                                                         "TRUE"="NMD")) +
  scale_color_manual(values = noAS_NMD_colors, labels = c("FALSE"="MANE",
                                                          "TRUE"="NMD")) +
  geom_text(data = noAS_NMD_sum %>% filter(PTC == "TRUE"),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
            color = "black") +
  geom_label(data = noAS_NMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = PTC,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = noAS_NMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = PTC,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type",
       title = "Genes with no AS events on the PTC list")+
  coord_cartesian(y = c(-4,4))
noAS_NMD_boxplot
ggsave("noAS_NMD_boxplot_FC_plot.pdf",
       device = pdf,
       plot = noAS_NMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)


####Determine if NMD genes are AS####
NMD_genes = c("EIF4A3","RBM8A","MAGOH","RNPS1","CASC3","ICE1","PYM1",
              "UPF1","UPF2","UPF3B","UPF3A","ETF1","GSPT1","NCBP1","NCBP2","EIF4E","SMG1","SMG8","SMG9","DHX34","RUVBL1","RUVBL2",
              "SMG5","SMG7","SMG6","CNOT8","DCP1A","PNRC2","DCP2","MOV10","PP2CA","PPP2R1A","XRN1","DIS3L","DIS3L2","EXOSC10","PARN",
              "NBAS",
              "PABPC1","EIF4G1","FMR1","EIF3E","SRSF1","SEC13","GNL2")
NMD_AS = AF_AS_MANE %>% filter(external_gene_name %in% NMD_genes)
NMD_AS_sum = NMD_AS %>% group_by(Sample, AS_gene) %>% summarise(n = n(),
                                                                med = median(log2FoldChange))
  
NMD_AS_boxplot = ggplot(data = NMD_AS)
NMD_AS_boxplot = NMD_AS_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = AS_gene),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = AS_colors, labels = c("FALSE" = "no AS",
                                                   "TRUE" = "AS genes")) +
  geom_label(data = NMD_AS_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = AS_gene,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = NMD_AS_sum,
            aes(x = factor(Sample,
                           levels = all_GOI),
                y = -3,
                color = AS_gene,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Gene Type",
       title = "NMD Genes AS",
       caption = "full NMD factor list")+
  coord_cartesian(y = c(-4,4))
NMD_AS_boxplot
ggsave("NMD_AS_boxplot.pdf",
       device = pdf,
       plot = NMD_AS_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

NMD_AS_boxplot_main = ggplot(data = NMD_AS %>% filter(Sample %in% Pres_KD))
NMD_AS_boxplot_main = NMD_AS_boxplot_main + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = AS_gene),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = AS_colors, labels = c("FALSE" = "no AS",
                                                   "TRUE" = "AS genes")) +
  geom_label(data = NMD_AS_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = AS_gene,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = NMD_AS_sum %>% filter(Sample %in% Pres_KD),
            aes(x = factor(Sample,
                           levels = all_GOI),
                y = -2.5,
                color = AS_gene,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Gene Type",
       title = "NMD Genes AS",
       caption = "full NMD factor list")+
  coord_cartesian(y = c(-4,4))
NMD_AS_boxplot_main
ggsave("NMD_AS_boxplot_mainFig.pdf",
       device = pdf,
       plot = NMD_AS_boxplot_main,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)

filt_NMD_genes = c("EIF4A3","RBM8A","MAGOH","RNPS1","CASC3","ICE1","PYM1",
                   "UPF1","UPF2","UPF3B","UPF3A","ETF1","GSPT1","NCBP1","NCBP2","EIF4E","SMG1","SMG8","SMG9","DHX34","RUVBL1","RUVBL2",
                   "SMG5","SMG7","SMG6","MOV10")
NMD_AS_filt = AF_AS_MANE %>% filter(external_gene_name %in% filt_NMD_genes)
NMD_AS_filt_sum = NMD_AS_filt %>% group_by(Sample, AS_gene) %>% summarise(n = n(),
                                                                med = median(log2FoldChange))

NMD_AS_filt_boxplot = ggplot(data = NMD_AS_filt)
NMD_AS_filt_boxplot = NMD_AS_filt_boxplot + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = AS_gene),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = AS_colors, labels = c("FALSE" = "no AS",
                                                   "TRUE" = "AS genes")) +
  geom_label(data = NMD_AS_filt_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = AS_gene,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = NMD_AS_filt_sum,
            aes(x = factor(Sample,
                           levels = all_GOI),
                y = -3,
                color = AS_gene,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Gene Type",
       title = "NMD Genes AS",
       caption = "core NMD factor list")+
  coord_cartesian(y = c(-4,4))
NMD_AS_filt_boxplot
ggsave("NMD_AS_filt_boxplot.pdf",
       device = pdf,
       plot = NMD_AS_filt_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

NMD_AS_filt_boxplot_main = ggplot(data = NMD_AS_filt %>% filter(Sample %in% Pres_KD))
NMD_AS_filt_boxplot_main = NMD_AS_filt_boxplot_main + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = AS_gene),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = AS_colors, labels = c("FALSE" = "no AS",
                                                   "TRUE" = "AS genes")) +
  geom_label(data = NMD_AS_filt_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = AS_gene,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = NMD_AS_filt_sum %>% filter(Sample %in% Pres_KD),
            aes(x = factor(Sample,
                           levels = all_GOI),
                y = -3,
                color = AS_gene,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Gene Type",
       title = "NMD Genes AS",
       caption = "core NMD factor list")+
  coord_cartesian(y = c(-4,4))
NMD_AS_filt_boxplot_main
ggsave("NMD_AS_filt_boxplot_mainFig.pdf",
       device = pdf,
       plot = NMD_AS_filt_boxplot_main,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)


####Analyze PRPF31 RP IPSC and RPE####
Disease_sample = c("PRPF31_RP","PRPF8_RP")
for (i in Disease_sample) {
  print(i)
  assign(paste0(i,"_PTC_alltrans"), read_csv(paste0(i,"_PTC_alltrans.csv")))
  assign(paste0(i,"_PTC_alltrans"), eval(parse(text = paste0(i,"_PTC_alltrans"))) %>% mutate(Sample = paste0(i)) %>% 
           dplyr::select(2:8,10:12)) 
}
Disease_PTC = PRPF31_RP_PTC_alltrans %>% full_join(PRPF8_RP_PTC_alltrans)
Disease_PTC_sum = Disease_PTC %>% group_by(PTC,Sample) %>%
  summarise(med = median(log2FoldChange),
            mean = mean(log2FoldChange),
            n = n()) %>% 
  mutate(PTC = str_replace(PTC,"FALSE","MANE"),
         PTC = str_replace(PTC, "TRUE", "PTC"))
all_DIS_res = tibble(Sample = character(),P = numeric())
for (i in Disease_sample) {
  print(i)
  assign(paste0(i,"_PTC_res"),
         wilcox.test(log2FoldChange ~ PTC, data = Disease_PTC %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  all_DIS_res = all_DIS_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_PTC_res$p.value"))))
}
Disease_PTC_sum = Disease_PTC_sum %>% left_join(all_DIS_res)
write_csv(Disease_PTC,"Disease_PTC_transcripts.csv")

DIS_boxplot = ggplot(data = Disease_PTC)
DIS_boxplot = DIS_boxplot + 
  geom_boxplot(aes(x = Sample,
                   y = log2FoldChange,
                   fill = PTC),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = colors, labels = c("FALSE" = "MANE",
                                                "TRUE" = "PTC")) +
  scale_color_manual(values = colors) +
  geom_text(data = Disease_PTC_sum,
             aes(x = Sample,
                 y = 3.5,
                 label = signif(P,digits = 3)),
             size = 7,
             show.legend = F,
            color = "black")+
  geom_label(data = Disease_PTC_sum,
             aes(x = Sample,
                 y = med,
                 color = PTC,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = Disease_PTC_sum,
             aes(x = Sample,
                 y = -3,
                 color = PTC,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Gene Type")+
  coord_cartesian(y = c(-4,4))
DIS_boxplot
ggsave("Disease_boxplot_FC_plot.pdf",
       device = pdf,
       plot = DIS_boxplot,
       width = 10,
       height = 10,
       units = "in",
       dpi = 300)

##Volcano plots of Disease PTC
Disease_PTC_vol_data = Disease_PTC %>% mutate(Sig = case_when(log2FoldChange > 0.58 & padj < 0.05 ~ "UP",
                                                                    log2FoldChange < -0.58 & padj < 0.05 ~ "DOWN",
                                                                    TRUE ~ "NS"))
Disease_PTC_vol_labs = Disease_PTC_vol_data %>% group_by(Sig,Sample) %>% slice_min(padj, n = 10) %>% filter(Sig != "NS")
write_csv(Disease_PTC_vol_labs,"top10_sig_Disease_PTC.csv")

vol_colors = c("UP" = "#F5BC00",
               "DOWN" = "#840B2D",
               "NS" = "#9D9C9D")
Disease_PTC_vol = ggplot(data = Disease_PTC_vol_data)
Disease_PTC_vol = Disease_PTC_vol + geom_vline(xintercept = c(-0.58,0.58),
                                               color = "grey",
                                               linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             color = "grey",
             linetype = "dashed") +
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj),
                 color = Sig),
             alpha = 0.5) +
  geom_label_repel(data = Disease_PTC_vol_labs,
                   aes(x = log2FoldChange,
                       y = -log10(padj),
                       color = Sig,
                       label = Gene),
                   show.legend = FALSE,
                   max.overlaps = NA) +
  scale_color_manual(values = vol_colors)+
  labs(title = "Disease vs WT",
       y = "-Log10(padj)",
       x = "Log2(Fold Change)",
       color = "Significance",
       caption = "PTC gene transcripts") +
  facet_wrap(facets = vars(Sample)) +
  theme_bw()
Disease_PTC_vol
ggsave("Disease_PTC_volcano.pdf",
       Disease_PTC_vol,
       width = 20,
       height = 10,
       units = "in",
       dpi = 300,
       device = pdf)

##Look at the disease effect on PE##
for (i in Disease_sample) {
  print(i)
  assign(paste0(i,"_full_alltrans"),
         read_csv(paste0(i,"_DBfilt_alltrans.csv")))
  assign(paste0(i,"_full_alltrans"),
         eval(parse(text = paste0(i,"_full_alltrans"))) %>% 
           mutate(Sample = i))
}
DIS_alltrans = PRPF31_RP_full_alltrans %>% full_join(PRPF8_RP_full_alltrans)

## Filter the log2FC data to only transcripts in the PE/MANE list
PE_DIS = DIS_alltrans %>% filter(ENST.ID %in% PE_transcripts$ensembl_transcript_id)
PE_DIS = PE_DIS %>% left_join(PE_transcripts, by = c("ENST.ID" = "ensembl_transcript_id"))
PE_DIS_sum = PE_DIS %>% group_by(Sample,isoform) %>%
  summarise(n = n(),
            genes = n_distinct(ensembl_gene_id),
            Med = median(log2FoldChange, na.rm = T),
            quart = quantile(log2FoldChange, 0.5)) #calculate the summary stats for labeling the plots
PE_DIS_sum = PE_DIS_sum %>% group_by(Sample) %>% mutate(Total = sum(n),
                                                        Total_genes = sum(genes))
write_csv(PE_DIS,"Disease_PE_transcripts.csv")

## Calculates the P-value for each plot
DIS_PE_res = tibble(Sample = as.character(), P = as.numeric())
for (i in Disease_sample) {
  print(i)
  assign(paste0(i,"_PE"), PE_DIS %>% filter(Sample == i))
  assign(paste0(i,"_res"), wilcox.test(log2FoldChange ~ isoform, data = eval(parse(text = paste0(i,"_PE"))),
                                       exact = FALSE, alternative = "less")) #use less because we expect MANE to be lower than PE
  DIS_PE_res = DIS_PE_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_res$p.value"))))
}
PE_DIS_sum = PE_DIS_sum %>% left_join(DIS_PE_res, by = "Sample")


DIS_PE_boxplot = ggplot(data = PE_DIS)
DIS_PE_boxplot = DIS_PE_boxplot + 
  geom_boxplot(aes(x =Sample,
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = PE_Colors) +
  scale_color_manual(values = PE_Colors)+
  geom_text(data = PE_DIS_sum %>% filter(isoform == "MANE"),
             aes(x = Sample,
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 7,
            show.legend = F,
            color = "black")+
  geom_label(data = PE_DIS_sum,
             aes(x = Sample,
                 y = Med,
                 color = isoform,
                 label = round(Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = PE_DIS_sum,
             aes(x = Sample,
                 y = -2.5,
                 color = isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 7,
             direction = "y",
             segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type")+
  coord_cartesian(y = c(-4,4))
DIS_PE_boxplot
ggsave("Disease_PE_boxplot.pdf",
       device = pdf,
       plot = DIS_PE_boxplot,
       width = 10,
       height = 10,
       units = "in",
       dpi = 300)

##Volcano plots of Disease PE
Disease_PE_vol_data = PE_DIS %>% mutate(Sig = case_when(log2FoldChange > 0.58 & padj < 0.05 ~ "UP",
                                                              log2FoldChange < -0.58 & padj < 0.05 ~ "DOWN",
                                                              TRUE ~ "NS"))
Disease_PE_vol_labs = Disease_PE_vol_data %>% group_by(Sig,Sample) %>% slice_min(padj, n = 10) %>% filter(Sig != "NS")
write_csv(Disease_PE_vol_labs,"top10_sig_Disease_PE.csv")

Disease_PE_vol = ggplot(data = Disease_PE_vol_data)
Disease_PE_vol = Disease_PE_vol + geom_vline(xintercept = c(-0.58,0.58),
                                               color = "grey",
                                               linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             color = "grey",
             linetype = "dashed") +
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj),
                 color = Sig),
             alpha = 0.5) +
  geom_label_repel(data = Disease_PE_vol_labs,
                   aes(x = log2FoldChange,
                       y = -log10(padj),
                       color = Sig,
                       label = external_gene_name),
                   show.legend = FALSE,
                   max.overlaps = NA) +
  scale_color_manual(values = vol_colors)+
  labs(title = "Disease vs WT",
       y = "-Log10(padj)",
       x = "Log2(Fold Change)",
       color = "Significance",
       caption = "PE gene transcripts") +
  facet_wrap(facets = vars(Sample)) +
  theme_bw()
Disease_PE_vol
ggsave("Disease_PE_volcano.pdf",
       Disease_PE_vol,
       width = 20,
       height = 10,
       units = "in",
       dpi = 300,
       device = pdf)

#Volcano plot of all transcripts#
Disease_all_vol_data = DIS_alltrans %>% mutate(Sig = case_when(log2FoldChange > 0.58 & padj < 0.05 ~ "UP",
                                                        log2FoldChange < -0.58 & padj < 0.05 ~ "DOWN",
                                                        TRUE ~ "NS"))
Dis_all_genes = getBM(attributes = c("ensembl_gene_id","external_gene_name","ensembl_transcript_id","transcript_biotype",
                                     "transcript_mane_select"),
                      filters = "ensembl_transcript_id",
                      values = Disease_all_vol_data$ENST.ID,
                      mart = ensembl)
Dis_all_genes = Dis_all_genes %>% distinct(ensembl_transcript_id, .keep_all = TRUE)
Disease_all_vol_data = Disease_all_vol_data %>% left_join(Dis_all_genes, by = c("ENST.ID" = "ensembl_transcript_id"))
Disease_all_vol_labs = Disease_all_vol_data %>% group_by(Sig,Sample) %>% slice_min(padj, n = 10) %>% filter(Sig != "NS")
write_csv(Disease_all_vol_labs,"top10_sig_Disease_all_transcripts.csv")

Disease_all_vol = ggplot(data = Disease_all_vol_data)
Disease_all_vol = Disease_all_vol + geom_vline(xintercept = c(-0.58,0.58),
                                             color = "grey",
                                             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             color = "grey",
             linetype = "dashed") +
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj),
                 color = Sig),
             alpha = 0.5) +
  geom_label_repel(data = Disease_PTC_vol_labs,
                   aes(x = log2FoldChange,
                       y = -log10(padj),
                       label = Gene),
                   show.legend = FALSE,
                   color = "black",
                   max.overlaps = NA) +
  geom_label_repel(data = Disease_all_vol_labs,
                   aes(x = log2FoldChange,
                       y = -log10(padj),
                       label = external_gene_name,
                       color = Sig),
                   show.legend = FALSE,
                   max.overlaps = NA) +
  scale_color_manual(values = vol_colors)+
  labs(title = "Disease vs WT",
       y = "-Log10(padj)",
       x = "Log2(Fold Change)",
       color = "Significance",
       caption = "all transcripts") +
  facet_wrap(facets = vars(Sample)) +
  theme_bw()
Disease_all_vol
ggsave("Disease_all_volcano.pdf",
       Disease_all_vol,
       width = 30,
       height = 10,
       units = "in",
       dpi = 300,
       device = pdf)

## Look at the no NMD transcripts in the disease samples##
Dis_no_NMD = DIS_alltrans %>% inner_join(MS_noNMD_transcripts, by = c("ENST.ID" = "ensembl_transcript_id"))
Dis_no_NMD_sum = Dis_no_NMD %>% group_by(Sample,isoform) %>% 
  summarise(n = n(),
            med = median(log2FoldChange))
Dis_no_NMD_res = tibble(Sample = character(), P = numeric())
for (i in Disease_sample) {
  print(i)
  assign(paste0(i,"_noNMD_res"), wilcox.test(log2FoldChange ~ isoform, data = Dis_no_NMD %>% filter(Sample == i),
                                       exact = FALSE, alternative = "less")) #use less because we expect MANE to be lower than PE
  Dis_no_NMD_res = Dis_no_NMD_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_noNMD_res$p.value"))))
}
Dis_no_NMD_sum = Dis_no_NMD_sum %>% left_join(Dis_no_NMD_res, by = "Sample")

DIS_no_NMD_boxplot = ggplot(data = Dis_no_NMD)
DIS_no_NMD_boxplot = DIS_no_NMD_boxplot + 
  geom_boxplot(aes(x =Sample,
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = noNMD_colors) +
  scale_color_manual(values = noNMD_colors)+
  geom_text(data = Dis_no_NMD_sum %>% filter(isoform == "MANE"),
            aes(x = Sample,
                y = 3,
                label = signif(P,digits = 3)),
            size = 7,
            show.legend = F,
            color = "black")+
  geom_label(data = Dis_no_NMD_sum,
             aes(x = Sample,
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = Dis_no_NMD_sum,
            aes(x = Sample,
                y = -2.5,
                color = isoform,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type")+
  coord_cartesian(y = c(-4,4))
DIS_no_NMD_boxplot
ggsave("Disease_no_NMD_boxplot.pdf",
       device = pdf,
       plot = DIS_no_NMD_boxplot,
       width = 10,
       height = 10,
       units = "in",
       dpi = 300)

#Look at the AS genes in the disease sample
all_disease_AS = tibble(Sample = character())
for (i in Disease_sample) {
  print(i)
  assign(paste0(i,"_AS_events"),
         read_csv(paste0(i,"_annotated_AS.csv")))
  assign(paste0(i,"_alltrans_genes"),
         getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","transcript_mane_select"),
               filters = "ensembl_transcript_id",
               values = eval(parse(text = paste0(i,"_full_alltrans$ENST.ID"))),
               mart = ensembl))
  assign(paste0(i,"_AS_alltrans"),
         eval(parse(text = paste0(i,"_full_alltrans"))) %>% 
           left_join(eval(parse(text = paste0(i,"_alltrans_genes"))),
                     by = c("ENST.ID" = "ensembl_transcript_id")))
  assign(paste0(i,"_AS_alltrans"),
         eval(parse(text = paste0(i,"_AS_alltrans"))) %>% 
           mutate(Gene_type = if_else(ensembl_gene_id %in% eval(parse(text = paste0(i,"_AS_events$GeneID"))),
                                      "AS",
                                      "non-AS"),
                  Transcript_type = if_else(str_detect(transcript_mane_select,"NM"),
                                            "MANE",
                                            "Other")))
  all_disease_AS = all_disease_AS %>% full_join(eval(parse(text = paste0(i,"_AS_alltrans"))))
}
disease_AS_sum = all_disease_AS %>% group_by(Sample,Gene_type) %>% 
  filter(Transcript_type == "MANE") %>% 
  summarise(n = n(),
            med = median(log2FoldChange))
disease_AS_res = tibble(Sample = character(), P = numeric())
for (i in Disease_sample) {
  print(i)
  assign(paste0(i,"_AS_res"), wilcox.test(log2FoldChange ~ Gene_type, data = all_disease_AS %>% filter(Sample == i & Transcript_type == "MANE"),
                                             exact = FALSE, alternative = "less")) #use less because we expect MANE to be lower than PE
  disease_AS_res = disease_AS_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_AS_res$p.value"))))
}
disease_AS_sum = disease_AS_sum %>% left_join(disease_AS_res, by = "Sample")

disease_AS_colors = c("AS" = "#6DC5B9",
                      "non-AS" = "#DE4D86")
disease_AS_boxplot = ggplot(data = all_disease_AS %>% filter(Transcript_type == "MANE"))
disease_AS_boxplot = disease_AS_boxplot + 
  geom_boxplot(aes(x =Sample,
                   y = log2FoldChange,
                   fill = factor(Gene_type,levels = c("non-AS","AS"))),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = disease_AS_colors,labels = c("non-AS" = "no AS",
                                                          "AS" = "AS genes")) +
  scale_color_manual(values = disease_AS_colors)+
  geom_text(data = disease_AS_sum %>% filter(Gene_type == "AS"),
            aes(x = Sample,
                y = 3,
                label = signif(P,digits = 3)),
            size = 7,
            show.legend = F,
            color = "black")+
  geom_label(data = disease_AS_sum,
             aes(x = Sample,
                 y = med,
                 color = factor(Gene_type,levels = c("non-AS","AS")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = disease_AS_sum,
            aes(x = Sample,
                y = -2.5,
                color = factor(Gene_type,levels = c("non-AS","AS")),
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Gene Type")+
  coord_cartesian(y = c(-4,4))
disease_AS_boxplot
ggsave("Disease_AS_boxplot.pdf",
       device = pdf,
       plot = disease_AS_boxplot,
       width = 10,
       height = 10,
       units = "in",
       dpi = 300)


####Look at fold change of non-PTC NMD targets####
noPTC_NMD_MANE <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/NMD_MANE_noPTC.csv")

#Make the datatables
all_noPTC_NMD = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noPTC_NMD"),
         eval(parse(text = paste0(i,"_full_alltrans"))) %>% 
           inner_join(noPTC_NMD_MANE, by = c("ENST.ID" = "transID")))
  all_noPTC_NMD = all_noPTC_NMD %>% full_join(eval(parse(text = paste0(i,"_noPTC_NMD"))))
} #make full datatable


noPTC_NMD_summary = all_noPTC_NMD %>% group_by(Sample,isoform) %>% 
  summarise(n = n(),
            med = median(log2FoldChange)) #Make summary datatable
all_noPTC_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noPTC_res"),
         wilcox.test(log2FoldChange ~ isoform, data = all_noPTC_NMD %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  all_noPTC_res = all_noPTC_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_noPTC_res$p.value"))))
} #Do p-value calculation
noPTC_NMD_summary = noPTC_NMD_summary %>% left_join(all_noPTC_res)

noPTC_NMD_colors = c("NMD" = "#EBB028",
                      "MANE" = "#1F2947")
noPTC_NMD_boxplot = ggplot(data = all_noPTC_NMD)
noPTC_NMD_boxplot = noPTC_NMD_boxplot + 
  geom_boxplot(aes(x =factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = noPTC_NMD_colors) +
  scale_color_manual(values = noPTC_NMD_colors)+
  geom_text(data = noPTC_NMD_summary %>% filter(isoform == "MANE"),
            aes(x = factor(Sample, levels = all_GOI),
                y = 3,
                label = signif(P,digits = 3)),
            size = 7,
            color = "black")+
  geom_label(data = noPTC_NMD_summary,
             aes(x =factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = noPTC_NMD_summary,
            aes(x = factor(Sample, levels = all_GOI),
                y = -3,
                color = isoform,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "isoform Type")+
  coord_cartesian(y = c(-4,4))
noPTC_NMD_boxplot
ggsave("noPTC_NMD_boxplot.pdf",
       device = pdf,
       plot = noPTC_NMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

#Look only at those with no AS genes
#Make the datatables
all_noPTC_noAS_NMD = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noPTC_noAS_NMD"),
         No_AS_only %>% filter(Sample == i) %>% 
           filter(ENST.ID %in% noPTC_NMD_MANE$transID) %>% 
           mutate(isoform = if_else(str_detect(transcript_mane_select,"NM"),
                                    "MANE",
                                    "NMD")))
  all_noPTC_noAS_NMD = all_noPTC_noAS_NMD %>% full_join(eval(parse(text = paste0(i,"_noPTC_noAS_NMD"))))
} #make full datatable


noPTC_noAS_NMD_summary = all_noPTC_noAS_NMD %>% group_by(Sample,isoform) %>% 
  summarise(n = n(),
            med = median(log2FoldChange)) #Make summary datatable
all_noPTC_noAS_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noPTC_noAS_res"),
         wilcox.test(log2FoldChange ~ isoform, data = all_noPTC_noAS_NMD %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  all_noPTC_noAS_res = all_noPTC_noAS_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_noPTC_noAS_res$p.value"))))
} #Do p-value calculation
noPTC_noAS_NMD_summary = noPTC_noAS_NMD_summary %>% left_join(all_noPTC_noAS_res)

noPTC_noAS_NMD_colors = c("NMD" = "#F1C86A",
                     "MANE" = "#445A9C")
noPTC_noAS_NMD_boxplot = ggplot(data = all_noPTC_noAS_NMD)
noPTC_noAS_NMD_boxplot = noPTC_noAS_NMD_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = noPTC_noAS_NMD_colors) +
  scale_color_manual(values = noPTC_noAS_NMD_colors)+
  geom_text(data = noPTC_noAS_NMD_summary %>% filter(isoform == "MANE"),
            aes(x = factor(Sample, levels = all_GOI),
                y = 3,
                label = signif(P,digits = 3)),
            size = 7,
            color = "black")+
  geom_label(data = noPTC_noAS_NMD_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text_repel(data = noPTC_noAS_NMD_summary,
            aes(x = factor(Sample, levels = all_GOI),
                y = -2.5,
                color = isoform,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "isoform Type")+
  coord_cartesian(y = c(-4,4))
noPTC_noAS_NMD_boxplot
ggsave("noPTC_noAS_NMD_boxplot.pdf",
       device = pdf,
       plot = noPTC_noAS_NMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

#Look only at those with AS genes
#Make the datatables
all_noPTC_AS_NMD = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noPTC_AS_NMD"),
         AS_only %>% filter(Sample == i) %>% 
           filter(ENST.ID %in% noPTC_NMD_MANE$transID) %>% 
           mutate(isoform = if_else(str_detect(transcript_mane_select,"NM"),
                                    "MANE",
                                    "NMD")))
  all_noPTC_AS_NMD = all_noPTC_AS_NMD %>% full_join(eval(parse(text = paste0(i,"_noPTC_AS_NMD"))))
} #make full datatable


noPTC_AS_NMD_summary = all_noPTC_AS_NMD %>% group_by(Sample,isoform) %>% 
  summarise(n = n(),
            med = median(log2FoldChange)) #Make summary datatable
all_noPTC_AS_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noPTC_AS_res"),
         wilcox.test(log2FoldChange ~ isoform, data = all_noPTC_AS_NMD %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  all_noPTC_AS_res = all_noPTC_AS_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_noPTC_AS_res$p.value"))))
} #Do p-value calculation
noPTC_AS_NMD_summary = noPTC_AS_NMD_summary %>% left_join(all_noPTC_AS_res)

noPTC_AS_NMD_colors = c("NMD" = "#CC9C00",
                          "MANE" = "#8EABE1")
noPTC_AS_NMD_boxplot = ggplot(data = all_noPTC_AS_NMD)
noPTC_AS_NMD_boxplot = noPTC_AS_NMD_boxplot + 
  geom_boxplot(aes(x =factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = noPTC_AS_NMD_colors) +
  scale_color_manual(values = noPTC_AS_NMD_colors)+
  geom_text(data = noPTC_AS_NMD_summary %>% filter(isoform == "MANE"),
            aes(x = factor(Sample, levels = all_GOI),
                y = 3,
                label = signif(P,digits = 3)),
            size = 7,
            color = "black")+
  geom_label(data = noPTC_AS_NMD_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = noPTC_AS_NMD_summary,
            aes(x = factor(Sample, levels = all_GOI),
                y = -2.5,
                color = isoform,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "isoform Type")+
  coord_cartesian(y = c(-4,4))
noPTC_AS_NMD_boxplot
ggsave("noPTC_AS_NMD_boxplot.pdf",
       device = pdf,
       plot = noPTC_AS_NMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

#Look at MANE vs other transcripts for genes on the no PTC NMD list that are undergoing AS
all_noPTC_AS_MANE = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noPTC_AS_MANE"),
         AS_only %>% filter(Sample == i) %>% 
           filter(ensembl_gene_id %in% noPTC_NMD_MANE$ensembl_gene_id))
  all_noPTC_AS_MANE = all_noPTC_AS_MANE %>% full_join(eval(parse(text = paste0(i,"_noPTC_AS_MANE"))))
}
noPTC_AS_TSL = getBM(attributes = c("ensembl_transcript_id","transcript_tsl"),
                     filters = "ensembl_transcript_id",
                     values = all_noPTC_AS_MANE$ENST.ID,
                     mart = ensembl)
noPTC_AS_TSL = noPTC_AS_TSL %>% mutate(TSL = str_extract(transcript_tsl,"tsl[:digit:]")) #Fixes the tsl annotation to just be the tsl number
noPTC_AS_TSL = noPTC_AS_TSL %>% filter(TSL == "tsl1") %>% select("ensembl_transcript_id","TSL") %>% 
  distinct(ensembl_transcript_id,.keep_all = TRUE)
all_noPTC_AS_tsl = all_noPTC_AS_MANE %>% inner_join(noPTC_AS_TSL, by = c("ENST.ID" = "ensembl_transcript_id"))

noPTC_AS_tsl_summary = all_noPTC_AS_tsl %>% group_by(Sample,Isoform) %>% 
  summarise(n = n(),
            med = median(log2FoldChange)) #Make summary datatable
all_noPTC_AS_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_noPTC_AS_res"),
         wilcox.test(log2FoldChange ~ Isoform, data = all_noPTC_AS_tsl %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  all_noPTC_AS_res = all_noPTC_AS_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_noPTC_AS_res$p.value"))))
} #Do p-value calculation
noPTC_AS_tsl_summary = noPTC_AS_tsl_summary %>% left_join(all_noPTC_AS_res)

noPTC_AS_tsl_colors = c("Other" = "#274C91",
                        "MANE" = "#8EABE1")
noPTC_AS_tsl_boxplot = ggplot(data = all_noPTC_AS_tsl)
noPTC_AS_tsl_boxplot = noPTC_AS_tsl_boxplot + 
  geom_boxplot(aes(x =factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = Isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = noPTC_AS_tsl_colors) +
  scale_color_manual(values = noPTC_AS_tsl_colors)+
  geom_text(data = noPTC_AS_tsl_summary %>% filter(Isoform == "MANE"),
            aes(x = factor(Sample, levels = all_GOI),
                y = 3,
                label = signif(P,digits = 3)),
            size = 7,
            color = "black")+
  geom_label(data = noPTC_AS_tsl_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = noPTC_AS_tsl_summary,
            aes(x = factor(Sample, levels = all_GOI),
                y = -2.5,
                color = Isoform,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type")+
  coord_cartesian(y = c(-4,4))
noPTC_AS_tsl_boxplot
ggsave("noPTC_AS_tsl_boxplot.pdf",
       device = pdf,
       plot = noPTC_AS_tsl_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)


####Look at NMD transcripts that are MANE themselves####
MANE_NMD = read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/NMD_MANE_and_TSL1.csv")
all_MANE_NMD = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_MANE_NMD"),
         eval(parse(text = paste0(i,"_full_alltrans"))) %>% 
           inner_join(MANE_NMD, by = c("ENST.ID" = "transID")))
  all_MANE_NMD = all_MANE_NMD %>% full_join(eval(parse(text = paste0(i,"_MANE_NMD"))))
}

all_MANE_NMD_summary = all_MANE_NMD %>% group_by(Sample,isoform) %>% 
  summarise(n = n(),
            med = median(log2FoldChange))

all_MANE_NMD_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_MANE_NMD_res"),
         wilcox.test(log2FoldChange ~ isoform, data = all_MANE_NMD %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  all_MANE_NMD_res = all_MANE_NMD_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_MANE_NMD_res$p.value"))))
} #Do p-value calculation
all_MANE_NMD_summary = all_MANE_NMD_summary %>% left_join(all_MANE_NMD_res)

MANE_NMD_colors = c("Other" = "#542344",
                    "MANE" = "#D33F49")
MANE_NMD_boxplot = ggplot(data = all_MANE_NMD)
MANE_NMD_boxplot = MANE_NMD_boxplot + 
  geom_boxplot(aes(x =factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = isoform),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = MANE_NMD_colors) +
  scale_color_manual(values = MANE_NMD_colors)+
  geom_text(data = all_MANE_NMD_summary %>% filter(isoform == "MANE"),
            aes(x = factor(Sample, levels = all_GOI),
                y = 3,
                label = signif(P,digits = 3)),
            size = 8,
            color = "black")+
  geom_label(data = all_MANE_NMD_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = all_MANE_NMD_summary,
            aes(x = factor(Sample, levels = all_GOI),
                y = -2.5,
                color = isoform,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 8) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type")+
  coord_cartesian(y = c(-4,4))
MANE_NMD_boxplot
ggsave("MANE_NMD_boxplot.pdf",
       device = pdf,
       plot = MANE_NMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

#Look at MANE NMD targets vs other MANE transcripts
all_MANE = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_alltrans_MANE"),
         eval(parse(text = paste0(i,"_full_alltrans"))) %>% left_join(eval(parse(text = paste0(i,"_alltrans_genes"))),
                                                                           by = c("ENST.ID" = "ensembl_transcript_id")) %>% 
           filter(str_detect(transcript_mane_select,"NM")) %>% 
    mutate(NMD = if_else(ENST.ID %in% MANE_NMD$transID,"NMD","Stable")))
  all_MANE = all_MANE %>% full_join(eval(parse(text = paste0(i,"_alltrans_MANE"))))
}

all_MANE_summary = all_MANE %>% group_by(Sample, NMD) %>% summarise(n = n(),
                                                                    med = median(log2FoldChange))
all_MANE_res = tibble(Sample = character(), P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_MANE_res"),
         wilcox.test(log2FoldChange ~ NMD, data = all_MANE %>% filter(Sample == i),
                     exact = FALSE, alternative = "greater"))
  all_MANE_res = all_MANE_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_MANE_res$p.value"))))
} #Do p-value calculation
all_MANE_summary = all_MANE_summary %>% left_join(all_MANE_res)

all_MANE_colors = c("NMD" = "#F27D2E",
                    "Stable" = "#542344")
MANE_boxplot = ggplot(data = all_MANE)
MANE_boxplot = MANE_boxplot + 
  geom_boxplot(aes(x =factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = NMD),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = all_MANE_colors) +
  scale_color_manual(values = all_MANE_colors)+
  geom_text(data = all_MANE_summary %>% filter(NMD == "NMD"),
            aes(x = factor(Sample, levels = all_GOI),
                y = 3,
                label = signif(P,digits = 3)),
            size = 8,
            color = "black")+
  geom_label(data = all_MANE_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = NMD,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = all_MANE_summary,
            aes(x = factor(Sample, levels = all_GOI),
                y = -2.5,
                color = NMD,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 8) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Isoform Type")+
  coord_cartesian(y = c(-4,4))
MANE_NMD_boxplot
ggsave("MANE_boxplot.pdf",
       device = pdf,
       plot = MANE_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)


#### Look at MANE vs PTC vs Other genes from the same set of genes ####
PTC_genes = getBM(attributes = "ensembl_gene_id",
                  filters = "ensembl_transcript_id",
                  values = sPTC_list$transID,
                  mart = ensembl)
PTC_isoforms = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","transcript_biotype",
                                    "transcript_mane_select"),
                     filters = "ensembl_gene_id",
                     values = PTC_genes$ensembl_gene_id,
                     mart = ensembl)
all_PTCgenes = sPTC_list %>% 
  full_join(PTC_isoforms,by = c("transID" = "ensembl_transcript_id")) %>% 
  mutate(Type = case_when(str_detect(transcript_mane_select,"NM") ~ "MANE",
                          PTC == "TRUE" ~ "PTC",
                          TRUE ~ "Other")) %>% 
  select(transID,ensembl_gene_id,transcript_biotype,transcript_mane_select,Type)

all_PTC_foldchange = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_PTC_all_isoforms"),
         eval(parse(text = paste0(i,"_full_alltrans"))) %>% inner_join(all_PTCgenes, by = c("ENST.ID" = "transID")))
  all_PTC_foldchange = all_PTC_foldchange %>% full_join(eval(parse(text = paste0(i,"_PTC_all_isoforms"))))
}
all_PTC_FC_summary = all_PTC_foldchange %>% group_by(Sample,Type) %>% 
  summarise(n = n(),
            med = median(log2FoldChange))

all_PTC_colors = c("MANE" = "#663171",
                   "PTC" = "#EA7428",
                   "NMD" = "#EA7428",
                   "Other" = "#6FC37D")
all_PTC_boxplot = ggplot(all_PTC_foldchange)
all_PTC_boxplot = all_PTC_boxplot + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                                     y = log2FoldChange,
                                                     fill = Type),
                                                 position = position_dodge2(width = 0.9),
                                                 width = 0.8,
                                                 outlier.shape = 21,
                                                 outlier.alpha = 0.5,
                                                 outlier.colour = NA,
                                                 linewidth = 1) +
  scale_fill_manual(values = all_PTC_colors) +
  scale_color_manual(values = all_PTC_colors) +
  geom_text_repel(data = all_PTC_FC_summary,
            aes(x = factor(Sample, levels = all_GOI),
                y = -3,
                color = Type,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7,
            direction = "y",
            segment.color = NA) +
  geom_label(data = all_PTC_FC_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = Type,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 3) +
  labs(x = "Depletion",
       y = "Log2(Fold Change)",
       fill = "Isoform Type") +
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()
all_PTC_boxplot
ggsave("PTC_genes_all_isoforms.pdf",
       plot = all_PTC_boxplot,
       width = 40,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

#Filter to only the TSL 1 or 2 isoforms
PTC_isoforms_tsl = getBM(attributes = c("ensembl_transcript_id","transcript_tsl"),
                         filters = "ensembl_transcript_id",
                         values = all_PTC_foldchange$ENST.ID,
                         mart = ensembl)
PTC_isoforms_tsl = PTC_isoforms_tsl %>% distinct(ensembl_transcript_id, .keep_all = TRUE)
PTC_isoforms_tsl = all_PTC_foldchange %>% left_join(PTC_isoforms_tsl, by = c("ENST.ID" = "ensembl_transcript_id"), relationship = "many-to-many")
PTC_isoforms_tsl = PTC_isoforms_tsl %>% mutate(TSL = str_extract(transcript_tsl,"tsl."))
PTC_isoforms_tsl = PTC_isoforms_tsl %>% filter(TSL == "tsl1" | TSL == "tsl2")
PTC_tsl_summary = PTC_isoforms_tsl %>% group_by(Sample, Type) %>% summarise(med = median(log2FoldChange),
                                                                            n = n())
PTC_TSL_boxplot = ggplot(PTC_isoforms_tsl)
PTC_TSL_boxplot = PTC_TSL_boxplot + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                                     y = log2FoldChange,
                                                     fill = Type),
                                                 position = position_dodge2(width = 0.9),
                                                 width = 0.8,
                                                 outlier.shape = 21,
                                                 outlier.alpha = 0.5,
                                                 outlier.colour = NA,
                                                 linewidth = 1) +
  scale_fill_manual(values = all_PTC_colors) +
  scale_color_manual(values = all_PTC_colors) +
  geom_text_repel(data = PTC_tsl_summary,
                  aes(x = factor(Sample, levels = all_GOI),
                      y = -3,
                      color = Type,
                      label = n),
                  show.legend = F,
                  position = position_dodge2(width = 0.9),
                  size = 7,
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = PTC_tsl_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = Type,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 3) +
  labs(x = "Depletion",
       y = "Log2(Fold Change)",
       fill = "Isoform Type",
       caption = "TSL 1 or 2") +
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()
PTC_TSL_boxplot
ggsave("PTC_genes_TSL12_all_isoforms.pdf",
       plot = PTC_TSL_boxplot,
       width = 40,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

mainfig_PTC_TSL_boxplot = ggplot(PTC_isoforms_tsl %>% filter(Sample %in% Pres_KD))
mainfig_PTC_TSL_boxplot = mainfig_PTC_TSL_boxplot + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                                     y = log2FoldChange,
                                                     fill = Type),
                                                 position = position_dodge2(width = 0.9),
                                                 width = 0.8,
                                                 outlier.shape = 21,
                                                 outlier.alpha = 0.5,
                                                 outlier.colour = NA,
                                                 linewidth = 1) +
  scale_fill_manual(values = all_PTC_colors) +
  scale_color_manual(values = all_PTC_colors) +
  geom_text_repel(data = PTC_tsl_summary %>% filter(Sample %in% Pres_KD),
                  aes(x = factor(Sample, levels = all_GOI),
                      y = -3,
                      color = Type,
                      label = n),
                  show.legend = F,
                  position = position_dodge2(width = 0.9),
                  size = 7,
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = PTC_tsl_summary %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = Type,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 3) +
  labs(x = "Depletion",
       y = "Log2(Fold Change)",
       fill = "Isoform Type",
       caption = "TSL 1 or 2") +
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()
mainfig_PTC_TSL_boxplot
ggsave("main_fig_PTC_genes_TSL12_all_isoforms.pdf",
       plot = mainfig_PTC_TSL_boxplot,
       width = 20,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)



####ID genes with AS events in multiple samples####
shared_AS = MAGOH_sig_AS_genes %>% inner_join(EIF4A3_sig_AS_genes) %>% inner_join(UPF1_sig_AS_genes) %>% 
  inner_join(RBM22_sig_AS_genes) %>% inner_join(AQR_sig_AS_genes) %>% inner_join(SNRNP200_sig_AS_genes) %>% 
  inner_join(EFTUD2_sig_AS_genes) %>% inner_join(SF3B1_sig_AS_genes) %>% inner_join(SF3B3_sig_AS_genes) %>% 
  inner_join(SNRPC_sig_AS_genes) %>% inner_join(SNRNP70_sig_AS_genes) %>% inner_join(PRPF8_sig_AS_genes) %>% 
  inner_join(PRPF6_sig_AS_genes) %>% inner_join(CDC5L_sig_AS_genes) %>% inner_join(SF3A1_sig_AS_genes) %>% 
  inner_join(SF3A3_sig_AS_genes) %>% inner_join(U2AF1_sig_AS_genes) %>% inner_join(CDC40_sig_AS_genes) %>% 
  inner_join(PRPF3_sig_AS_genes) %>% inner_join(PRPF4_sig_AS_genes) %>% inner_join(GNB2L1_sig_AS_genes) %>%
  distinct(GeneID)
shared_AS_annotation = getBM(attributes = c("ensembl_gene_id","external_gene_name","ensembl_transcript_id",
                                            "transcript_biotype","cds_length"),
                             filters = "ensembl_gene_id",
                             values = shared_AS$GeneID,
                             mart = ensembl)
Shared_AS_NMD_iso = shared_AS_annotation %>% filter(transcript_biotype == "nonsense_mediated_decay")
genes_with_novel_isoforms <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/genes_with_novel_isoforms.csv")
shared_AS_novel = Shared_AS_NMD_iso %>% filter(ensembl_gene_id %in% genes_with_novel_isoforms$gene_id)
shared_AS_novel = shared_AS_novel %>% filter(ensembl_transcript_id %in% Saltzman_inclusion_NMD_transcripts$ensembl_transcript_id &
                                               cds_length < 1000)




####Make a master list of transcripts in the KDs####
master_list = tibble(Sample = character())
for (i in all_GOI) {
  print(i)
  master_list = master_list %>% full_join(eval(parse(text = paste0(i,"_full_alltrans"))))
}
master_list = master_list %>% select(Sample,log2FoldChange,ENST.ID,baseMean,padj)
master_list = master_list %>% group_by(Sample, ENST.ID) %>% distinct(ENST.ID, .keep_all = TRUE) %>% ungroup() #removes some rows that got duplicated when making the ENST.ID column
master_list = master_list %>% pivot_wider(names_from = Sample, values_from = c(log2FoldChange,baseMean,padj),
                                          names_vary = "slowest",values_fill = NA)
master_list = master_list %>% left_join(sPTC_list %>% select(transID,PTC), by = c("ENST.ID" = "transID")) %>% 
  left_join(Saltzman_PE %>% select(ensembl_transcript_id,Isoform),
                                   by = c("ENST.ID" = "ensembl_transcript_id"))
master_list = master_list %>% 
  mutate(PTC = str_replace(PTC,"FALSE","MANE"),
         PTC = str_replace(PTC,"TRUE","PTC")) %>% 
  rename(PE = Isoform)
master_annotations = getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id","external_gene_name",
                                          "transcript_mane_select"),
                           filters = "ensembl_transcript_id",
                           values = master_list$ENST.ID,
                           mart = ensembl)
master_list = master_list %>% left_join(master_annotations, by = c("ENST.ID" = "ensembl_transcript_id"))
write_csv(master_list,"KD_DEseq_masterlist.csv")



###Look at highly upregulated transcripts####
upreg_test = c("EIF4A3","EFTUD2","CDC40","SF3B1","AQR")
for (i in upreg_test) {
  print(i)
  assign(paste0(i,"_Upregulated"),
         eval(parse(text = paste0(i,"_AS_alltrans"))) %>% 
           filter(baseMean > 100 & log2FoldChange >= 0.58 & padj < 0.05 & transcript_biotype == "nonsense_mediated_decay"))
  assign(paste0(i,"_genes_up"),
         getBM(attributes = c("external_gene_name","ensembl_gene_id"),
               filters = "ensembl_gene_id",
               values = eval(parse(text = paste0(i,"_Upregulated$ensembl_gene_id"))),
               mart = ensembl))
  assign(paste0(i,"_Upregulated"),
         eval(parse(text = paste0(i,"_Upregulated"))) %>% left_join(eval(parse(text = paste0(i,"_genes_up")))))
  write_csv(eval(parse(text = paste0(i,"_Upregulated"))),
            paste0(i,"_Upregulated_transcripts.csv"))
  assign(paste0(i,"_up_geneset"),
         eval(parse(text = paste0(i,"_Upregulated"))) %>% distinct(ENST.ID))
}
set.seed(1)
s <- list("EIF4A3" = EIF4A3_up_geneset$ENST.ID,
          "EFTUD2" = EFTUD2_up_geneset$ENST.ID,
          "CDC40" = CDC40_up_geneset$ENST.ID,
          "SF3B1" = SF3B1_up_geneset$ENST.ID,
          "AQR" = AQR_up_geneset$ENST.ID)
ggvenn = ggVennDiagram(s,
                       label_alpha = 0.5,
                       label = "count") +
  scale_fill_moma_c("Ernst") +
  labs(title = "Upregulated NMD Biotype Transcripts",
       caption = "Base Mean > 100, 1.5 fold up, padj < 0.05")
ggvenn
ggsave("Upregulated_overlap.pdf",
       plot = ggvenn,
       device = pdf,
       width = 10,
       height = 10,
       units = "in",
       dpi =300)


####Make a scatter plot of KDs compared to EIF4A3 KD####
for (i in upreg_test) {
  print(i)
  assign(paste0(i,"_high_mean"),
         eval(parse(text = paste0(i,"_AS_alltrans"))) %>% 
           filter(baseMean > 100 & transcript_biotype == "nonsense_mediated_decay") %>% 
           select(baseMean,log2FoldChange,padj,ENST.ID,Sample,ensembl_gene_id,transcript_mane_select,
                  transcript_biotype))
}

all_high_mean = AQR_high_mean %>% full_join(EFTUD2_high_mean) %>% 
  full_join(CDC40_high_mean) %>% full_join(SF3B1_high_mean)
EIF4A3_high_mean = EIF4A3_high_mean %>% rename(EIF4A3_FC = log2FoldChange) %>% select(ENST.ID,EIF4A3_FC)
all_high_mean = all_high_mean %>% inner_join(EIF4A3_high_mean, relationship = "many-to-one")

base_mean_colors = c("EFTUD2" = "#472965",
                     "CDC40" = "#E0A500",
                     "AQR" = "#3891A6",
                     "SF3B1" = "#D84654")
base_mean_scatter = ggplot(data = all_high_mean)
base_mean_scatter = base_mean_scatter + geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0,color = "black") +
  geom_point(aes(x = log2FoldChange,
                 y = EIF4A3_FC,
                 color = Sample),
             alpha = 0.5) +
  geom_smooth(aes(x = log2FoldChange,
                  y = EIF4A3_FC),
              color = "darkgrey",
              method = "lm") +
  facet_wrap(vars(Sample),ncol = 2) +
  theme_bw() +
  labs(title = "NMD Biotype Comparison",
       caption = "Base Mean > 100",
       y = "EIF4A3 log2(FC)")
base_mean_scatter
ggsave("NMD_biotype_comparison.pdf",
       plot = base_mean_scatter,
       device = pdf,
       width = 12,
       height = 10,
       units = "in",
       dpi = 300)

up_mean = all_high_mean %>% filter(log2FoldChange > 0.58)
up_mean_scatter = ggplot(data = up_mean)
up_mean_scatter = up_mean_scatter + geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0,color = "black") +
  geom_point(aes(x = log2FoldChange,
                 y = EIF4A3_FC,
                 color = Sample),
               alpha = 0.5) +
  geom_smooth(aes(x = log2FoldChange,
                  y = EIF4A3_FC),
              color = "darkgrey",
              method = "lm") +
  facet_wrap(vars(Sample),ncol = 2) +
  theme_bw() +
  labs(title = "NMD Biotype Comparison",
       caption = "Base Mean > 100, non-EIF4A3 FC>1.5",
       y = "EIF4A3 log2(FC)")
up_mean_scatter
ggsave("NMD_biotype_comparison_upreg.pdf",
       plot = up_mean_scatter,
       device = pdf,
       width = 12,
       height = 10,
       units = "in",
       dpi = 300)


####Pull the transcripts that are up and down in each sample ####
full_upregulated = tibble(Sample = character())
for (i in all_GOI) {
  assign(paste0(i,"_upregulated"),
         eval(parse(text = paste0(i,"_AS_alltrans"))) %>% 
           filter(log2FoldChange > 0.58 & baseMean > 100 & padj < 0.05) %>% 
           select(baseMean,log2FoldChange,padj,ENST.ID,Sample,ensembl_gene_id))
  full_upregulated = full_upregulated %>% full_join(eval(parse(text = paste0(i,"_upregulated"))))
}
full_upregulated_count = full_upregulated %>% filter(Sample != "UPF1" | Sample != "EIF4A3" | Sample != "MAGOH") %>% 
  group_by(ENST.ID) %>% 
  summarise(count = n(),
            med_BM = median(baseMean),
            med_FC = median(log2FoldChange),
            avg_BM = mean(baseMean),
            avg_FC = mean(log2FoldChange))
shared_upregulated = full_upregulated_count %>% filter(count >= 10)
shared_genes = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),
                     filters = "ensembl_transcript_id",
                     values = shared_upregulated$ENST.ID,
                     mart = ensembl)
shared_upregulated = shared_upregulated %>% full_join(shared_genes,by = c("ENST.ID" = "ensembl_transcript_id"))
write_csv(shared_upregulated,"upregulated_transcripts_10ormore.csv")

full_downregulated = tibble(Sample = character())
for (i in all_GOI) {
  assign(paste0(i,"_downregulated"),
         eval(parse(text = paste0(i,"_AS_alltrans"))) %>% 
           filter(log2FoldChange < 0.58 & baseMean > 100 & padj < 0.05) %>% 
           select(baseMean,log2FoldChange,padj,ENST.ID,Sample,ensembl_gene_id))
  full_downregulated = full_downregulated %>% full_join(eval(parse(text = paste0(i,"_downregulated"))))
}
full_downregulated_count = full_downregulated %>% filter(Sample != "UPF1" | Sample != "EIF4A3" | Sample != "MAGOH") %>% 
  group_by(ENST.ID) %>% 
  summarise(count = n(),
            med_BM = median(baseMean),
            med_FC = median(log2FoldChange),
            avg_BM = mean(baseMean),
            avg_FC = mean(log2FoldChange))
shared_downregulated = full_downregulated_count %>% filter(count >= 10)
shared_genes_down = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),
                     filters = "ensembl_transcript_id",
                     values = shared_downregulated$ENST.ID,
                     mart = ensembl)
shared_downregulated = shared_downregulated %>% full_join(shared_genes_down,by = c("ENST.ID" = "ensembl_transcript_id"))
write_csv(shared_downregulated,"downregulated_transcripts_10ormore.csv")


#Look at the PTC+ isoforms that are upregulated
all_PTC_list = read_delim("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/ENST_PTC-EPI-TFG.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

full_sPTC_up = tibble(Sample = character())
lax_sPTC_up = tibble(Sample = character())
full_PTC_up = tibble(Sample = character())
for (i in all_GOI) {
  assign(paste0(i,"_PTC_annotated"),
         eval(parse(text = paste0(i,"_AS_alltrans"))) %>% 
           left_join(sPTC_list, by = c("ENST.ID" = "transID")) %>% 
           left_join(all_PTC_list, by = c("ENST.ID" = "ENST-ID"))  %>%
           rename(sPTC = PTC,
                  PTC_plus = "PTC-Status") %>% 
           select("baseMean","log2FoldChange","padj","ENST.ID","Sample","ensembl_gene_id","sPTC","PTC_plus"))
  assign(paste0(i,"_sPTC_upregulated"),
         eval(parse(text = paste0(i,"_PTC_annotated"))) %>% 
           filter(sPTC == "TRUE" & baseMean > 100 & log2FoldChange > 0.58 & padj < 0.05))
  assign(paste0(i,"_sPTC_lax_up"),
         eval(parse(text = paste0(i,"_PTC_annotated"))) %>% 
           filter(sPTC == "TRUE" & baseMean > 5 & log2FoldChange > 0.58 & padj < 0.05))
  full_sPTC_up = full_sPTC_up %>% full_join(eval(parse(text = paste0(i,"_sPTC_upregulated"))))
  lax_sPTC_up = lax_sPTC_up %>% full_join(eval(parse(text = paste0(i,"_sPTC_lax_up"))))
}
full_sPTC_count = full_sPTC_up %>% filter(Sample != "UPF1" | Sample != "EIF4A3" | Sample != "MAGOH") %>% 
  group_by(ENST.ID) %>% 
  summarise(count = n(),
            med_BM = median(baseMean),
            med_FC = median(log2FoldChange),
            avg_BM = mean(baseMean),
            avg_FC = mean(log2FoldChange))
shared_sPTC_up = full_sPTC_count %>% filter(count > 5)
shared_sPTC_genes = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),
                     filters = "ensembl_transcript_id",
                     values = shared_sPTC_up$ENST.ID,
                     mart = ensembl)
shared_sPTC_up = shared_sPTC_up %>% left_join(shared_sPTC_genes, by = c("ENST.ID" = "ensembl_transcript_id")) %>% 
  rename(Number_KD = count)
write_csv(shared_sPTC_up,"sPTC_upregulated.csv")
lax_sPTC_count = lax_sPTC_up %>% filter(Sample != "UPF1" | Sample != "EIF4A3" | Sample != "MAGOH") %>% 
  group_by(ENST.ID) %>% 
  summarise(count = n(),
            med_BM = median(baseMean),
            med_FC = median(log2FoldChange),
            avg_BM = mean(baseMean),
            avg_FC = mean(log2FoldChange))
shared_lax_up = lax_sPTC_count %>% filter(count > 5)
shared_lax_genes = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),
                         filters = "ensembl_transcript_id",
                         values = shared_lax_up$ENST.ID,
                         mart = ensembl)
shared_lax_up = shared_lax_up %>% left_join(shared_lax_genes, by = c("ENST.ID" = "ensembl_transcript_id")) %>% 
  rename(Number_KD = count)
write_csv(shared_lax_up,"lax_sPTC_upregulated.csv")


#Look at the PE isoforms that are upregulated
full_PE_up = tibble(Sample = character())
lax_PE_up = tibble(Sample = character())
PE_NMD_only = PE_transcripts %>% filter(isoform == "NMD") %>% distinct(ensembl_transcript_id)
for (i in all_GOI) {
  assign(paste0(i,"_PE_upregulated"),
         eval(parse(text = paste0(i,"_AS_alltrans"))) %>% 
           filter(ENST.ID %in% PE_NMD_only$ensembl_transcript_id & baseMean > 100 & log2FoldChange > 0.58 & padj < 0.05))
  assign(paste0(i,"_PE_lax_up"),
         eval(parse(text = paste0(i,"_AS_alltrans"))) %>% 
           filter(ENST.ID %in% PE_NMD_only$ensembl_transcript_id & baseMean > 5 & log2FoldChange > 0.58 & padj < 0.05))
  full_PE_up = full_PE_up %>% full_join(eval(parse(text = paste0(i,"_PE_upregulated"))))
  lax_PE_up = lax_PE_up %>% full_join(eval(parse(text = paste0(i,"_PE_lax_up"))))
}
full_PE_count = full_PE_up %>% filter(Sample != "UPF1" | Sample != "EIF4A3" | Sample != "MAGOH") %>% 
  group_by(ENST.ID) %>% 
  summarise(count = n(),
            med_BM = median(baseMean),
            med_FC = median(log2FoldChange),
            avg_BM = mean(baseMean),
            avg_FC = mean(log2FoldChange))
shared_PE_up = full_PE_count %>% filter(count > 5)
shared_PE_genes = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),
                          filters = "ensembl_transcript_id",
                          values = shared_PE_up$ENST.ID,
                          mart = ensembl)
shared_PE_up = shared_PE_up %>% left_join(shared_PE_genes, by = c("ENST.ID" = "ensembl_transcript_id")) %>% 
  rename(Number_KD = count)
write_csv(shared_PE_up,"PE_upregulated.csv")
lax_PE_count = lax_PE_up %>% filter(Sample != "UPF1" | Sample != "EIF4A3" | Sample != "MAGOH") %>% 
  group_by(ENST.ID) %>% 
  summarise(count = n(),
            med_BM = median(baseMean),
            med_FC = median(log2FoldChange),
            avg_BM = mean(baseMean),
            avg_FC = mean(log2FoldChange))
shared_lax_PE = lax_PE_count %>% filter(count > 5)
shared_lax_PE_genes = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","external_gene_name"),
                         filters = "ensembl_transcript_id",
                         values = shared_lax_PE$ENST.ID,
                         mart = ensembl)
shared_lax_PE = shared_lax_PE %>% left_join(shared_lax_PE_genes, by = c("ENST.ID" = "ensembl_transcript_id")) %>% 
  rename(Number_KD = count)
write_csv(shared_lax_PE,"lax_PE_upregulated.csv")






#### Look at the gene level DE analysis ####
full_GL_alltrans = tibble(Sample = character())
full_GL_NMD_alltrans = tibble(Sample = character())
full_GL_NMD_res = tibble(Sample = character(), P = numeric())
for (i in all_GOI) {
  assign(paste0(i,"_GL_alltrans"),
         read_csv(paste0(i,"_geneLevel_NMD_alltrans.csv")))
  assign(paste0(i,"_GL_alltrans"),
         eval(parse(text = paste0(i,"_GL_alltrans"))) %>% 
           mutate(Sample = i))
  full_GL_alltrans = full_GL_alltrans %>% full_join(eval(parse(text = paste0(i,"_GL_alltrans"))))
  assign(paste0(i,"_GL_NMD_alltrans"),
         eval(parse(text = paste0(i,"_GL_alltrans"))) %>% 
           filter(!is.na(NMD)))
  assign(paste0(i,"_GL_NMD_res"),
         wilcox.test(log2FoldChange ~ NMD, data = eval(parse(text = paste0(i,"_GL_NMD_alltrans"))),
                     exact = FALSE, alternative = "greater")) #greater because we expect the NMD to be higher
  full_GL_NMD_res = full_GL_NMD_res %>% add_row(Sample = i,P = eval(parse(text = paste0(i,"_GL_NMD_res$p.value"))))
  full_GL_NMD_alltrans = full_GL_NMD_alltrans %>% full_join(eval(parse(text = paste0(i,"_GL_NMD_alltrans"))))
}
full_GL_NMD_summary = full_GL_NMD_alltrans %>% group_by(Sample,NMD) %>% 
  summarise(n = n(),
            med = median(log2FoldChange)) %>% 
  left_join(full_GL_NMD_res)

GL_NMD_colors = NMD_GL_colors = c("NMD" = "#F27D2E",
                                  "noNMD" = "#B27092")
GL_NMD_box = ggplot(data = full_GL_NMD_alltrans)
GL_NMD_box = GL_NMD_box + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                     y = log2FoldChange,
                                     fill = factor(NMD,levels = c("NMD","noNMD"))),
                                 position = position_dodge2(width = 0.9),
                                 width = 0.8,
                                 outlier.shape = 21,
                                 outlier.alpha = 0.5,
                                 outlier.colour = NA,
                                 linewidth = 1) +
  scale_fill_manual(values = GL_NMD_colors) +
  scale_color_manual(values = GL_NMD_colors) +
  geom_text_repel(data = full_GL_NMD_summary,
                  aes(x = factor(Sample, levels = all_GOI),
                      y = -3,
                      color = factor(NMD,levels = c("NMD","noNMD")),
                      label = n),
                  show.legend = F,
                  position = position_dodge2(width = 0.8),
                  size = 7,
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = full_GL_NMD_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = factor(NMD,levels = c("NMD","noNMD")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = full_GL_NMD_summary %>% filter(NMD == "NMD"),
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = signif(P, digits = 3)),
             show.legend = F,
             size = 7,
             color = "black") +
  labs(y = "Log2(Fold Change)",
       fill = "Gene Type") +
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()
GL_NMD_box
ggsave("Gene_level_NMD_boxplot.pdf",
       plot = GL_NMD_box,
       width = 40,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

#Gene level AS vs non-AS#
AS_genes_only = AF_AS_combined %>% select(Sample, ensembl_gene_id, AS_gene) %>% group_by(Sample) %>% 
  distinct(ensembl_gene_id,.keep_all = TRUE) %>% ungroup()
GL_AS_combined = tibble(Sample = character())
full_GL_AS_res = tibble(Sample = character(), P = numeric())
for (i in all_GOI) {
  assign(paste0(i,"_GL_AS_alltrans"),
         eval(parse(text = paste0(i,"_GL_alltrans"))) %>% 
           inner_join(AS_genes_only %>% filter(Sample == i), by = c("Sample", "ENSG.ID" = "ensembl_gene_id")))
  GL_AS_combined = GL_AS_combined %>% full_join(eval(parse(text = paste0(i,"_GL_AS_alltrans"))))
  assign(paste0(i,"_GL_AS_res"),
         wilcox.test(log2FoldChange ~ AS_gene, data = eval(parse(text = paste0(i,"_GL_AS_alltrans"))),
                     exact = FALSE, alternative = "greater")) #greater because we expect the non-spliced to be higher
  full_GL_AS_res = full_GL_AS_res %>% add_row(Sample = i,P = eval(parse(text = paste0(i,"_GL_AS_res$p.value"))))
}
full_GL_AS_summary = GL_AS_combined %>% group_by(Sample,AS_gene) %>% 
  summarise(n = n(),
            med = median(log2FoldChange)) %>% 
  left_join(full_GL_AS_res)

GL_AS_box = ggplot(data = GL_AS_combined)
GL_AS_box = GL_AS_box + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                           y = log2FoldChange,
                                           fill = AS_gene),
                                       position = position_dodge2(width = 0.9),
                                       width = 0.8,
                                       outlier.shape = 21,
                                       outlier.alpha = 0.5,
                                       outlier.colour = NA,
                                       linewidth = 1) +
  scale_fill_manual(values = AS_colors, labels = c("FALSE" = "no AS",
                                                   "TRUE" = "AS")) +
  scale_color_manual(values = AS_colors, labels = c("FALSE" = "no AS",
                                                    "TRUE" = "AS")) +
  geom_text_repel(data = full_GL_AS_summary,
                  aes(x = factor(Sample, levels = all_GOI),
                      y = -3,
                      color = AS_gene,
                      label = n),
                  show.legend = F,
                  position = position_dodge2(width = 0.8),
                  size = 8,
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = full_GL_AS_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = AS_gene,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = full_GL_AS_summary %>% filter(AS_gene == "TRUE"),
            aes(x = factor(Sample, levels = all_GOI),
                y = 3,
                label = signif(P, digits = 3)),
            show.legend = F,
            size = 8,
            color = "black") +
  labs(y = "Log2(Fold Change)",
       fill = "Gene Type",
       title = "Gene level Alternate Splicing") +
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()
GL_AS_box
ggsave("Gene_level_AS_boxplot.pdf",
       plot = GL_AS_box,
       width = 40,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)




####Do the analysis using the SMG6/SMG7 depletion PTC list####
SMG_PTC_all_isoforms <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/SMG_PTC_all_isoforms.csv")
SMG_PTC_alltrans = MF_alltrans %>% inner_join(SMG_PTC_all_isoforms, by = c("ENST.ID" = "ensembl_transcript_id"))
SMG_PTC_MANE_alltrans = SMG_PTC_alltrans %>% filter(isoform != "Other")

SMG_PTC_summary = SMG_PTC_alltrans %>% group_by(Sample,isoform) %>% 
  summarise(n = n(),
            med = median(log2FoldChange))
full_SMG_PTC_res = tibble(Sample = character(),P_Other = numeric(),P_MANE = numeric())
for (i in all_GOI) {
  assign(paste0(i,"_other_SMG_PTC_res"),
         wilcox.test(log2FoldChange ~ isoform, data = SMG_PTC_alltrans %>% filter(Sample == i & isoform != "MANE"),
                     exact = FALSE, alternative = "greater")) #Expect NMD to be greater than other
  assign(paste0(i,"_MANE_SMG_PTC_res"),
         wilcox.test(log2FoldChange ~ isoform, data = SMG_PTC_alltrans %>% filter(Sample == i & isoform != "Other"),
                     exact = FALSE, alternative = "less")) #expect MANE to be less than NMD
  full_SMG_PTC_res = full_SMG_PTC_res %>% add_row(Sample = i, P_Other = eval(parse(text = paste0(i,"_other_SMG_PTC_res$p.value"))),
                                                  P_MANE = eval(parse(text = paste0(i,"_MANE_SMG_PTC_res$p.value"))))
}
SMG_PTC_summary = SMG_PTC_summary %>% left_join(full_SMG_PTC_res)
view(SMG_PTC_summary)


all_SMG_colors = c("MANE" = "#663171",
                   "NMD" = "#EA7428",
                   "Other" = "#6FC37D")
SMG_PTC_boxplot = ggplot(data = SMG_PTC_alltrans)
SMG_PTC_boxplot = SMG_PTC_boxplot + geom_boxplot(aes(x = factor(Sample,levels = all_GOI),
                                                     y = log2FoldChange,
                                                     fill = factor(isoform, levels = c("MANE","NMD","Other"))),
                                                 position = position_dodge(width = 0.9),
                                                 width = 0.8,
                                                 outlier.shape = 21,
                                                 outlier.alpha = 0.5,
                                                 outlier.colour = NA,
                                                 linewidth = 1.5) +
  scale_fill_manual(values = all_SMG_colors) +
  geom_text_repel(data = SMG_PTC_summary %>% filter(isoform == "NMD"),
            aes(x = factor(Sample,levels = all_GOI),
                y = 4,
                label = paste0("P=",signif(P_MANE,digits = 3))),
            size = 7,
            show.legend = F,
            color = "#663171",
            direction = "y",
            segment.color = NA) +
  geom_text_repel(data = SMG_PTC_summary %>% filter(isoform == "NMD"),
            aes(x = factor(Sample,levels = all_GOI),
                y = 3,
                label = paste0("P=",signif(P_Other,digits = 3))),
            size = 7,
            show.legend = F,
            color = "#6FC37D",
            direction = "y",
            segment.color = NA) +
  geom_label(data = SMG_PTC_summary,
             aes(x = factor(Sample,levels = all_GOI),
                 y = med,
                 color = factor(isoform,levels = c("MANE","NMD","Other")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 4) +
  scale_color_manual(values = all_SMG_colors) +
  geom_text_repel(data = SMG_PTC_summary,
                  aes(x = factor(Sample,levels = all_GOI),
                      y = -3,
                      label = n,
                      color = isoform),
                  show.legend = F,
                  position = position_dodge2(width = 0.9),
                  size = 7,
                  direction = "y",
                  segment.color = NA)+
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()+
  labs(x = "Depletion",
       fill = "Isoform Type",
       caption = "SMG6/SMG7 KD sPTC list")
SMG_PTC_boxplot
ggsave("SMG_PTC_all_isoforms_boxplot.pdf",
       plot = SMG_PTC_boxplot,
       device = pdf,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

SMG_PTC_boxplot_main = ggplot(data = SMG_PTC_alltrans %>% filter(Sample %in% Pres_KD))
SMG_PTC_boxplot_main = SMG_PTC_boxplot_main + geom_boxplot(aes(x = factor(Sample,levels = all_GOI),
                                                     y = log2FoldChange,
                                                     fill = factor(isoform, levels = c("MANE","NMD","Other"))),
                                                 position = position_dodge(width = 0.9),
                                                 width = 0.8,
                                                 outlier.shape = 21,
                                                 outlier.alpha = 0.5,
                                                 outlier.colour = NA,
                                                 linewidth = 1.5) +
  scale_fill_manual(values = all_SMG_colors) +
  geom_text_repel(data = SMG_PTC_summary %>% filter(isoform == "NMD" & Sample %in% Pres_KD),
                  aes(x = factor(Sample,levels = all_GOI),
                      y = 4,
                      label = paste0("P=",signif(P_MANE,digits = 3))),
                  size = 7,
                  show.legend = F,
                  color = "#663171",
                  direction = "y",
                  segment.color = NA) +
  geom_text_repel(data = SMG_PTC_summary %>% filter(isoform == "NMD" & Sample %in% Pres_KD),
                  aes(x = factor(Sample,levels = all_GOI),
                      y = 3,
                      label = paste0("P=",signif(P_Other,digits = 3))),
                  size = 7,
                  show.legend = F,
                  color = "#6FC37D",
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = SMG_PTC_summary %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,levels = all_GOI),
                 y = med,
                 color = factor(isoform,levels = c("MANE","NMD","Other")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 4) +
  scale_color_manual(values = all_SMG_colors) +
  geom_text_repel(data = SMG_PTC_summary %>% filter(Sample %in% Pres_KD),
                  aes(x = factor(Sample,levels = all_GOI),
                      y = -3,
                      label = n,
                      color = isoform),
                  show.legend = F,
                  position = position_dodge2(width = 0.9),
                  size = 7,
                  direction = "y",
                  segment.color = NA)+
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()+
  labs(x = "Depletion",
       fill = "Isoform Type",
       caption = "SMG6/SMG7 KD sPTC list")
SMG_PTC_boxplot_main
ggsave("SMG_PTC_all_isoforms_boxplot_mainFig.pdf",
       plot = SMG_PTC_boxplot_main,
       device = pdf,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)

SMG_PTC_MANE_summary = SMG_PTC_MANE_alltrans %>% group_by(Sample,isoform) %>% 
  summarise(n = n(),
            med = median(log2FoldChange))
full_SMG_PTC_MANE_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  assign(paste0(i,"_SMG_PTC_MANE_res"),
         wilcox.test(log2FoldChange ~ isoform, data = SMG_PTC_MANE_alltrans %>% filter(Sample == i),
                     exact = FALSE, alternative = "less")) #expect MANE to be less than NMD
  full_SMG_PTC_MANE_res = full_SMG_PTC_MANE_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_SMG_PTC_MANE_res$p.value"))))
}
SMG_PTC_MANE_summary = SMG_PTC_MANE_summary %>% left_join(full_SMG_PTC_MANE_res)
view(SMG_PTC_MANE_summary)

SMG_PTC_MANE_boxplot = ggplot(data = SMG_PTC_MANE_alltrans)
SMG_PTC_MANE_boxplot = SMG_PTC_MANE_boxplot + geom_boxplot(aes(x = factor(Sample,levels = all_GOI),
                                                     y = log2FoldChange,
                                                     fill = factor(isoform, levels = c("MANE","NMD","Other"))),
                                                 position = position_dodge(width = 0.9),
                                                 width = 0.8,
                                                 outlier.shape = 21,
                                                 outlier.alpha = 0.5,
                                                 outlier.colour = NA,
                                                 linewidth = 1.5) +
  scale_fill_manual(values = all_SMG_colors) +
  geom_text_repel(data = SMG_PTC_MANE_summary %>% filter(isoform == "NMD"),
                  aes(x = factor(Sample,levels = all_GOI),
                      y = 4,
                      label = paste0("P=",signif(P,digits = 3))),
                  size = 7,
                  show.legend = F,
                  color = "black",
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = SMG_PTC_MANE_summary,
             aes(x = factor(Sample,levels = all_GOI),
                 y = med,
                 color = factor(isoform,levels = c("MANE","NMD","Other")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 4) +
  scale_color_manual(values = all_SMG_colors) +
  geom_text_repel(data = SMG_PTC_MANE_summary,
                  aes(x = factor(Sample,levels = all_GOI),
                      y = -3,
                      label = n,
                      color = isoform),
                  show.legend = F,
                  position = position_dodge2(width = 0.9),
                  size = 7,
                  direction = "y",
                  segment.color = NA)+
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()+
  labs(x = "Depletion",
       fill = "Isoform Type",
       caption = "SMG6/SMG7 KD sPTC list")
SMG_PTC_MANE_boxplot
ggsave("SMG_PTC_MANE_boxplot.pdf",
       plot = SMG_PTC_MANE_boxplot,
       device = pdf,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

SMG_PTC_MANE_boxplot_main = ggplot(data = SMG_PTC_MANE_alltrans %>% filter(Sample %in% Pres_KD))
SMG_PTC_MANE_boxplot_main = SMG_PTC_MANE_boxplot_main + geom_boxplot(aes(x = factor(Sample,levels = all_GOI),
                                                               y = log2FoldChange,
                                                               fill = factor(isoform, levels = c("MANE","NMD","Other"))),
                                                           position = position_dodge(width = 0.9),
                                                           width = 0.8,
                                                           outlier.shape = 21,
                                                           outlier.alpha = 0.5,
                                                           outlier.colour = NA,
                                                           linewidth = 1.5) +
  scale_fill_manual(values = all_SMG_colors) +
  geom_text_repel(data = SMG_PTC_MANE_summary %>% filter(isoform == "NMD" & Sample %in% Pres_KD),
                  aes(x = factor(Sample,levels = all_GOI),
                      y = 4,
                      label = paste0("P=",signif(P,digits = 3))),
                  size = 7,
                  show.legend = F,
                  color = "black",
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = SMG_PTC_MANE_summary %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,levels = all_GOI),
                 y = med,
                 color = factor(isoform,levels = c("MANE","NMD","Other")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 4) +
  scale_color_manual(values = all_SMG_colors) +
  geom_text_repel(data = SMG_PTC_MANE_summary %>% filter(Sample %in% Pres_KD),
                  aes(x = factor(Sample,levels = all_GOI),
                      y = -3,
                      label = n,
                      color = isoform),
                  show.legend = F,
                  position = position_dodge2(width = 0.9),
                  size = 7,
                  direction = "y",
                  segment.color = NA)+
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()+
  labs(x = "Depletion",
       fill = "Isoform Type",
       caption = "SMG6/SMG7 KD sPTC list")
SMG_PTC_MANE_boxplot_main
ggsave("SMG_PTC_MANE_boxplot_mainFig.pdf",
       plot = SMG_PTC_MANE_boxplot_main,
       device = pdf,
       width = 22,
       height = 10,
       units = "in",
       dpi = 300)

##Make a heatmap of the P-value of each of the KDs for the summary figure
#Rearainge the data table
hm_df =SMG_PTC_MANE_summary %>% ungroup() %>% dplyr::select(Sample,P) %>% distinct(Sample,.keep_all = T) #Filter to just the data needed for the heatmap
hm_df = hm_df %>% mutate(Sample = factor(Sample, levels = all_GOI)) %>% arrange(Sample) #Put them in the order of the other figures
hm_df = hm_df %>% column_to_rownames(var="Sample") #Put the data in the format of pheatmap
hm_small = hm_df %>% filter(P < 0.01)
hm_large = hm_df %>% filter(P >= 0.01)
#Make the breaks
breaklist_small = c(seq(min(hm_df[hm_df != min(hm_df,na.rm = T)],na.rm = T),0.01,length.out = 1000))
colors_small = colorRampPalette(c("#474ED7","#EC458D"))(999)
breaklist_large = c(seq(0.01,1,length.out = 100))
colors_large = colorRampPalette(c("#EC458D","#FFF1BF"))(99)
#Try making two heatmaps, one going from min-0.05 and one from 0.05-1

#Make heatmap
HM_PTC_small = pheatmap(hm_small,
                        color = colors_small,
                        breaks = breaklist_small,
                        na_col = "grey",
                        cluster_rows = F,
                        cluster_cols = F,
                        display_numbers = T,
                        number_format = "%.1e",
                        number_color = "white")#heatmap going from the smallest P-value to 0.01
ggsave("PTC_heatmap_smallP.pdf",
       device = pdf,
       plot = HM_PTC_small,
       width = 3,
       height = 6,
       units = "in",
       dpi = 300)
HM_PTC_large = pheatmap(hm_large,
                        color = colors_large,
                        breaks = breaklist_large,
                        na_col = "grey",
                        cluster_rows = F,
                        cluster_cols = F,
                        display_numbers = T,
                        number_format = "%.1e",
                        number_color = "black")#heatmap going from the 0.01 to 1
ggsave("PTC_heatmap_largeP.pdf",
       device = pdf,
       plot = HM_PTC_large,
       width = 3,
       height = 3,
       units = "in",
       dpi = 300)

####Look at Spliciing Inhibitor treatment####
Splicing_inhibitors = c("Risdiplam","Pladienolide","THZ531","Indisulam")
full_SI_MANE_PTC = tibble(Sample = character())
full_SI_PTC_res = tibble(Sample = character(),P_Other = numeric(),P_MANE = numeric())
full_SI_new_sPTC = tibble(Sample = character())
full_SI_newPTC_res = tibble(Sample = character(),P_Other = numeric(),P_MANE = numeric())
for (i in Splicing_inhibitors) {
  print(i)
  assign(paste0(i,"_alltrans"),
         read_csv(paste0(i,"_DBfilt_alltrans.csv")))
  assign(paste0(i,"_alltrans_annotation"),
         getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id","external_gene_name",
                              "transcript_mane_select","transcript_biotype"),
               filters = "ensembl_transcript_id",
               values = eval(parse(text = paste0(i,"_alltrans$ENST.ID"))),
               mart = ensembl))
  assign(paste0(i,"_alltrans"),
         eval(parse(text = paste0(i,"_alltrans"))) %>% 
           left_join(eval(parse(text = paste0(i,"_alltrans_annotation"))),
                     by = c("ENST.ID" = "ensembl_transcript_id")) %>% 
           mutate(Sample = i))
  assign(paste0(i,"_MANE_PTC"),
         eval(parse(text = paste0(i,"_alltrans"))) %>% 
           inner_join(all_PTCgenes, by = c("ENST.ID" = "transID")))
  full_SI_MANE_PTC = full_SI_MANE_PTC %>% full_join(eval(parse(text = paste0(i,"_MANE_PTC"))))
  assign(paste0(i,"_other_PTC_res"),
         wilcox.test(log2FoldChange ~ Type, data = eval(parse(text = (paste0(i,"_MANE_PTC")))) %>% filter(Type != "MANE"),
                     exact = FALSE, alternative = "less"))
  assign(paste0(i,"_MANE_PTC_res"),
         wilcox.test(log2FoldChange ~ Type, data = eval(parse(text = (paste0(i,"_MANE_PTC")))) %>% filter(Type != "Other"),
                     exact = FALSE, alternative = "less"))
  full_SI_PTC_res = full_SI_PTC_res %>% add_row(Sample = i, P_Other = eval(parse(text = paste0(i,"_other_PTC_res$p.value"))),
                                                P_MANE = eval(parse(text = paste0(i,"_MANE_PTC_res$p.value"))))
  assign(paste0(i,"_new_sPTC"),
         eval(parse(text = paste0(i,"_alltrans"))) %>% 
           inner_join(SMG_PTC_all_isoforms, by = c("ENST.ID" = "ensembl_transcript_id",
                                                   "ensembl_gene_id",
                                                   "transcript_biotype",
                                                   "external_gene_name")))
  full_SI_new_sPTC = full_SI_new_sPTC %>% full_join(eval(parse(text = paste0(i,"_new_sPTC"))))
  assign(paste0(i,"_new_other_PTC_res"),
         wilcox.test(log2FoldChange ~ isoform, data = eval(parse(text = (paste0(i,"_new_sPTC")))) %>% filter(isoform != "MANE"),
                     exact = FALSE, alternative = "less"))
  assign(paste0(i,"_new_MANE_PTC_res"),
         wilcox.test(log2FoldChange ~ isoform, data = eval(parse(text = (paste0(i,"_new_sPTC")))) %>% filter(isoform != "Other"),
                     exact = FALSE, alternative = "less"))
  full_SI_newPTC_res = full_SI_PTC_res %>% add_row(Sample = i, P_Other = eval(parse(text = paste0(i,"_other_PTC_res$p.value"))),
                                                P_MANE = eval(parse(text = paste0(i,"_MANE_PTC_res$p.value"))))
}

#Look at MANE vs PTC vs Other
SI_PTC_summary = full_SI_MANE_PTC %>% group_by(Sample,Type) %>%  summarise(n = n(),
                                                                           med = median(log2FoldChange)) %>% 
  left_join(full_SI_PTC_res)


SI_PTC = ggplot(data = full_SI_MANE_PTC)
SI_PTC = SI_PTC + geom_boxplot(aes(x = factor(Sample,levels = Splicing_inhibitors),
                                   y = log2FoldChange,
                                   fill = factor(Type,levels = c("MANE","PTC","Other"))),
                               position = position_dodge2(width = 0.9),
                               width = 0.8,
                               outlier.shape = 21,
                               outlier.alpha = 0.5,
                               outlier.colour = NA,
                               linewidth = 1) +
  scale_fill_manual(values = all_PTC_colors) +
  scale_color_manual(values = all_PTC_colors) +
  geom_text_repel(data = SI_PTC_summary,
                  aes(x = factor(Sample,levels = Splicing_inhibitors),
                      y = -3,
                      color = factor(Type,levels = c("MANE","PTC","Other")),
                      label = n),
                  show.legend = F,
                  position = position_dodge2(width = 0.8),
                  size = 7,
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = SI_PTC_summary,
             aes(x = factor(Sample,levels = Splicing_inhibitors),
                 y = med,
                 color = factor(Type,levels = c("MANE","PTC","Other")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 7) +
  geom_text(data = SI_PTC_summary %>% filter(Type == "PTC"),
            aes(x = factor(Sample,levels = Splicing_inhibitors),
                y = 4,
                label = paste0("P(MANE)","=",signif(P_MANE,digits = 3))),
            size = 6,
            color = "#663171") +
  geom_text(data = SI_PTC_summary %>% filter(Type == "PTC"),
            aes(x = factor(Sample,levels = Splicing_inhibitors),
                y = 3,
                label = paste0("P(Other)","=",signif(P_Other,digits = 3))),
            size = 6,
            color = "#6FC37D") +
  labs(y = "Log2(Fold Change)",
       fill = "Isoform Type",
       caption = "PTC list isoforms",
       x = "Splicing inhibitor") +
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()
SI_PTC
ggsave("Splicing_Inhibitor_sPTC_boxplot.pdf",
       plot = SI_PTC,
       width = 22,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

#Look at  new MANE vs PTC vs Other
full_SI_newPTC_res = full_SI_newPTC_res %>% distinct(Sample,.keep_all = T)
SI_newPTC_summary = full_SI_new_sPTC %>% group_by(Sample,isoform) %>%  summarise(n = n(),
                                                                           med = median(log2FoldChange)) %>% 
  left_join(full_SI_newPTC_res)


SI_newPTC = ggplot(data = full_SI_new_sPTC)
SI_newPTC = SI_newPTC + geom_boxplot(aes(x = factor(Sample,levels = Splicing_inhibitors),
                                   y = log2FoldChange,
                                   fill = factor(isoform,levels = c("MANE","NMD","Other"))),
                               position = position_dodge2(width = 0.9),
                               width = 0.8,
                               outlier.shape = 21,
                               outlier.alpha = 0.5,
                               outlier.colour = NA,
                               linewidth = 1) +
  scale_fill_manual(values = all_PTC_colors) +
  scale_color_manual(values = all_PTC_colors) +
  geom_text_repel(data = SI_newPTC_summary,
                  aes(x = factor(Sample,levels = Splicing_inhibitors),
                      y = -3,
                      color = factor(isoform,levels = c("MANE","NMD","Other")),
                      label = n),
                  show.legend = F,
                  position = position_dodge2(width = 0.8),
                  size = 7,
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = SI_newPTC_summary,
             aes(x = factor(Sample,levels = Splicing_inhibitors),
                 y = med,
                 color = factor(isoform,levels = c("MANE","NMD","Other")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 7) +
  geom_text(data = SI_newPTC_summary %>% filter(isoform == "NMD"),
            aes(x = factor(Sample,levels = Splicing_inhibitors),
                y = 4,
                label = paste0("P(MANE)","=",signif(P_MANE,digits = 3))),
            size = 6,
            color = "#663171") +
  geom_text(data = SI_newPTC_summary %>% filter(isoform == "NMD"),
            aes(x = factor(Sample,levels = Splicing_inhibitors),
                y = 3,
                label = paste0("P(Other)","=",signif(P_Other,digits = 3))),
            size = 6,
            color = "#6FC37D") +
  labs(y = "Log2(Fold Change)",
       fill = "Isoform isoform",
       caption = "PTC list isoforms",
       x = "Splicing inhibitor") +
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()
SI_newPTC
ggsave("Splicing_Inhibitor_SMG_sPTC_boxplot.pdf",
       plot = SI_newPTC,
       width = 22,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

####Look at the effect of including novel isoforms in the kallisto transcriptome ####
NK_full_alltrans = tibble(Sample = character())
for (i in Pres_KD) {
  print(i)
  assign(paste0(i,"_NK_alltrans"),
         read_csv(paste0(i,"_NK_alltrans.csv")))
  assign(paste0(i,"_NK_alltrans"),
         eval(parse(text = paste0(i,"_NK_alltrans"))) %>% mutate(Sample = i,
                                                                 Category = case_when(str_detect(ENST.ID,"MST") ~ "Novel",
                                                                                      str_detect(ENST.ID,"ENST") & str_detect(transcript_mane_select,"NM") ~ "MANE",
                                                                                      TRUE ~ "Other")))
  NK_full_alltrans = NK_full_alltrans %>% full_join(eval(parse(text = paste0(i,"_NK_alltrans"))))
}
NK_full_summary = NK_full_alltrans %>% group_by(Sample, Category) %>% summarise(n = n(),
                                                                                med = median(log2FoldChange))

NK_colors = c("MANE" = "#713E5A",
              "Other" = "#6FC37D",
              "Novel" = "#DD7373")
NK_full_boxplot = ggplot(NK_full_alltrans)
NK_full_boxplot = NK_full_boxplot + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                                     y = log2FoldChange,
                                                     fill = Category),
                                                 position = position_dodge2(width = 0.9),
                                                 width = 0.8,
                                                 outlier.shape = 21,
                                                 outlier.alpha = 0.5,
                                                 outlier.colour = NA,
                                                 linewidth = 1) +
  scale_fill_manual(values = NK_colors) +
  scale_color_manual(values = NK_colors) +
  geom_text_repel(data = NK_full_summary,
                  aes(x = factor(Sample, levels = all_GOI),
                      y = -3,
                      color = Category,
                      label = n),
                  show.legend = F,
                  position = position_dodge2(width = 0.9),
                  size = 7,
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = NK_full_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = Category,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 3) +
  labs(x = "Depletion",
       y = "Log2(Fold Change)",
       fill = "Isoform Type",
       caption = "all isoforms") +
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()
NK_full_boxplot
ggsave("novel_kallisto_all_isoforms.pdf",
       plot = NK_full_boxplot,
       width = 20,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

#Filter to only the genes on the NMD list
NK_PTC_only = NK_full_alltrans %>% inner_join(SMG_PTC_all_isoforms, by = c("ENST.ID" = "ensembl_transcript_id",
                                                                           "external_gene_name",
                                                                           "ensembl_gene_id",
                                                                           "transcript_biotype"))

NK_PTC_summary = NK_PTC_only %>% group_by(Sample, isoform) %>% summarise(n = n(),
                                                                      med = median(log2FoldChange))

NK_PTC_colors = c("MANE" = "#663171",
                  "NMD" = "#EA7428",
                  "Other" = "#6FC37D",
                  "novel NMD" = "#DD7373",
                  "novel stable" = "#0075A2")
NK_PTC_boxplot = ggplot(NK_PTC_only)
NK_PTC_boxplot = NK_PTC_boxplot + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                                   y = log2FoldChange,
                                                   fill = factor(isoform,levels = c("MANE","NMD","Other","novel NMD","novel stable"))),
                                               position = position_dodge2(width = 0.9),
                                               width = 0.8,
                                               outlier.shape = 21,
                                               outlier.alpha = 0.5,
                                               outlier.colour = NA,
                                               linewidth = 1) +
  scale_fill_manual(values = NK_PTC_colors, limits = c("MANE","NMD","Other","novel NMD","novel stable")) +
  scale_color_manual(values = NK_PTC_colors, limits = c("MANE","NMD","Other","novel NMD","novel stable")) +
  geom_text_repel(data = NK_PTC_summary,
                  aes(x = factor(Sample, levels = all_GOI),
                      y = -3,
                      color = factor(isoform,levels = c("MANE","NMD","Other","novel NMD","novel stable")),
                      label = n),
                  show.legend = F,
                  position = position_dodge2(width = 0.9),
                  size = 7,
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = NK_PTC_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = factor(isoform,levels = c("MANE","NMD","Other","novel NMD","novel stable")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 3) +
  labs(x = "Depletion",
       y = "Log2(Fold Change)",
       fill = "Isoform isoform",
       caption = "PTC list isoforms") +
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()
NK_PTC_boxplot
ggsave("novel_kallisto_PTC_isoforms.pdf",
       plot = NK_PTC_boxplot,
       width = 20,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

NK_no_novel = NK_PTC_only %>% filter(!str_detect(isoform,"novel"))
NK_no_novel_summary = NK_no_novel %>% group_by(Sample, isoform) %>% summarise(n = n(),
                                                                           med = median(log2FoldChange))
NK_no_Novel_boxplot = ggplot(NK_no_novel)
NK_no_Novel_boxplot = NK_no_Novel_boxplot + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                                             y = log2FoldChange,
                                                             fill = factor(isoform,levels = c("MANE","NMD","Other"))),
                                                         position = position_dodge2(width = 0.9),
                                                         width = 0.8,
                                                         outlier.shape = 21,
                                                         outlier.alpha = 0.5,
                                                         outlier.colour = NA,
                                                         linewidth = 1) +
  scale_fill_manual(values = NK_PTC_colors, limits = c("MANE","NMD","Other")) +
  scale_color_manual(values = NK_PTC_colors, limits = c("MANE","NMD","Other")) +
  geom_text_repel(data = NK_no_novel_summary,
                  aes(x = factor(Sample, levels = all_GOI),
                      y = -3,
                      color = factor(isoform,levels = c("MANE","NMD","Other")),
                      label = n),
                  show.legend = F,
                  position = position_dodge2(width = 0.9),
                  size = 7,
                  direction = "y",
                  segment.color = NA) +
  geom_label(data = NK_no_novel_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = factor(isoform,levels = c("MANE","NMD","Other")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 3) +
  labs(x = "Depletion",
       y = "Log2(Fold Change)",
       fill = "Isoform isoform",
       caption = "PTC list isoforms") +
  coord_cartesian(ylim = c(-4,4)) +
  theme_bw()
NK_no_Novel_boxplot
ggsave("NK_PTC_isoforms_no_novel.pdf",
       plot = NK_no_Novel_boxplot,
       width = 20,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)