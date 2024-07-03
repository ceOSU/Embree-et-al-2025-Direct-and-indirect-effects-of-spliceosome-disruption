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

all_GOI = c("UPF1","EIF4A3","MAGOH",
            "AQR","RBM22","CDC5L",
            "EFTUD2","SNRNP200","PRPF8","PRPF6",
            "SF3A1","SF3B1","SF3B3","SF3A3","U2AF1",
            "CDC40",
            "PRPF3","PRPF4",
            "GNB2L1",
            "SNRNP70","SNRPC")
Pres_KD = c("UPF1","EIF4A3","MAGOH","AQR","RBM22","EFTUD2","SNRNP200","SF3B1","SF3B3","SNRPC","SNRNP70")



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
  assign(paste0(i,"_DBfilt"), read_csv(paste0(i,"_DBfilt_alltrans.csv")))
  assign(paste0(i,"_DBfilt"), eval(parse(text = paste0(i,"_DBfilt"))) %>% mutate(Sample = paste0(i)))
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
  assign(paste0(i,"_DBfilt_sig"), eval(parse(text = paste0(i,"_DBfilt"))) %>% mutate(Sample = paste0(i)) %>% 
           filter(padj <= 0.01))
}

sig_combined = MAGOH_DBfilt_sig %>% full_join(EIF4A3_DBfilt_sig) %>% full_join(UPF1_DBfilt_sig) %>% 
  full_join(RBM22_DBfilt_sig) %>% full_join(AQR_DBfilt_sig) %>% full_join(SNRNP200_DBfilt_sig) %>% 
  full_join(EFTUD2_DBfilt_sig) %>% full_join(SF3B1_DBfilt_sig) %>% full_join(SF3B3_DBfilt_sig) %>% 
  full_join(SNRPC_DBfilt_sig) %>% full_join(SNRNP70_DBfilt_sig) %>% full_join(PRPF8_DBfilt_sig) %>%
  full_join(PRPF6_DBfilt_sig) %>% full_join(CDC5L_DBfilt_sig) %>% full_join(SF3A1_DBfilt_sig) %>% 
  full_join(SF3A3_DBfilt_sig) %>% full_join(U2AF1_DBfilt_sig) %>% full_join(CDC40_DBfilt_sig) %>% 
  full_join(PRPF3_DBfilt_sig) %>% full_join(PRPF4_DBfilt_sig) %>% full_join(GNB2L1_DBfilt_sig) %>% 
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
  geom_label(data = PTC_sum %>% filter(Sample %in% Pres_KD),
             aes(x = Sample,
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 5,
             alpha = 0.75)+
  geom_label(data = PTC_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = Pres_KD),
                 y = med,
                 color = PTC,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label(data = PTC_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = Pres_KD),
                 y = -3.5,
                 color = PTC,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Transcript")+
  coord_cartesian(y = c(-4,4))
PTC_boxplot
#Last saved and modified on 5/10/2024
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
  geom_label(data = PTC_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 5,
             alpha = 0.75)+
  geom_label(data = PTC_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = PTC,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label(data = PTC_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3.5,
                 color = PTC,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Transcript")+
  coord_cartesian(y = c(-4,4))
full_PTC_boxplot
#Last saved and modified on 5/10/2024
ggsave("PTC_boxplot_full.pdf",
       device = pdf,
       plot = full_PTC_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

##Make a heatmap of the P-value of each of the KDs for the summary figure
#Rearainge the data table
hm_df =PTC_sum %>% ungroup() %>% dplyr::select(Sample,P) %>% distinct(Sample,.keep_all = T) #Filter to just the data needed for the heatmap
hm_df = hm_df %>% mutate(Sample = factor(Sample, levels = all_GOI)) %>% arrange(Sample) #Put them in the order of the other figures
hm_df = hm_df %>% column_to_rownames(var="Sample") #Put the data in the format of pheatmap
hm_small = hm_df %>% filter(P < 0.01)
hm_large = hm_df %>% filter(P >= 0.01)
#Make the breaks
breaklist_small = c(seq(min(hm_df),0.01,length.out = 1000))
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


####Look at PTC+- vs No NMD####
#import and manipulate data (Manu's no NMD gene list"#)
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
  geom_label(data = noNMD_sum,
             aes(x = Sample,
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 5,
             alpha = 0.75)+
  geom_label(data = noNMD_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label(data = noNMD_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3.5,
                 color = isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Transcript",
       title = "No NMD biotype")+
  coord_cartesian(y = c(-4,4))
noNMD_boxplot
# Last saved and modified on 5/10/2024
ggsave("noNMD_boxplot_FC.pdf",
       device = pdf,
       plot = noNMD_boxplot,
       width = 40,
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
             size = 5,
             alpha = 0.75)+
  geom_label(data = MSnoNMD_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = MSnoNMD_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3.5,
                 color = isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
  theme_bw() + 
  labs(x = "Sample",
       y = "Log2 Fold Change",
       fill = "Transcript",
       title = "MS full list")+
  coord_cartesian(y = c(-4,4))
MSnoNMD_boxplot
# Last saved and modified on 6/7/2024
ggsave("MSnoNMD_boxplot_FC.pdf",
       device = pdf,
       plot = MSnoNMD_boxplot,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

##Look at the non-NMD transcripts from the Lykke-Anderson genes and development paper
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
total_TPM_PTC = total_TPM_PTC + geom_point(aes(x = factor(Sample, levels = all_GOI),
                                                 y = total,
                                                 color = Type),
                                             size = 3) +
  scale_color_manual(values = TPM_summary_colors) +
  labs(x = "Depletion",
       y = "Total TPM",
       title = "Annotated isoforms",
       color = "Isoform") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust=1))
total_TPM_PTC
ggsave("total_TPM_annotated.pdf",
       plot = total_TPM_PTC,
       width = 12,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

total_TPM_novel = ggplot(data = all_TPM_summary %>% filter(Treatment == "kd") %>%
                         filter(Type == "Novel_NMD" | Type == "Novel_Stable"))
total_TPM_novel = total_TPM_novel + geom_point(aes(x = factor(Sample, levels = all_GOI),
                                               y = total,
                                               color = Type),
                                           size = 3) +
  scale_color_manual(values = TPM_summary_colors) +
  labs(x = "Depletion",
       y = "Total TPM",
       title = "Novel Isoforms",
       color = "Isoform") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust=1))
total_TPM_novel
ggsave("total_TPM_novel.pdf",
       plot = total_TPM_novel,
       width = 12,
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
main_fig_total_TPM_PTC = main_fig_total_TPM_PTC + geom_point(aes(x = factor(Sample, levels = all_GOI),
                                               y = total,
                                               color = Type),
                                           size = 3) +
  scale_color_manual(values = TPM_summary_colors) +
  labs(x = "Depletion",
       y = "Total TPM",
       title = "Annotated isoforms",
       color = "Isoform") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust=1))
main_fig_total_TPM_PTC
ggsave("main_fig_total_TPM_annotated.pdf",
       plot = main_fig_total_TPM_PTC,
       width = 8,
       height = 10,
       units = "in",
       device = pdf,
       dpi = 300)

main_fig_total_TPM_novel = ggplot(data = all_TPM_summary %>% filter(Treatment == "kd" & Sample %in% Pres_KD) %>%
                           filter(Type == "Novel_NMD" | Type == "Novel_Stable"))
main_fig_total_TPM_novel = main_fig_total_TPM_novel + geom_point(aes(x = factor(Sample, levels = all_GOI),
                                                   y = total,
                                                   color = Type),
                                               size = 3) +
  scale_color_manual(values = TPM_summary_colors) +
  labs(x = "Depletion",
       y = "Total TPM",
       title = "Novel Isoforms",
       color = "Isoform") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust=1))
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

#### Look at TPM of PE NMD isoforms ####
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_full_alltrans"),
         read_csv(paste0(i,"_DBfilt_alltrans.csv")))
  assign(paste0(i,"_full_alltrans"),
         eval(parse(text = paste0(i,"_full_alltrans"))) %>% 
           mutate(Sample = i))
} 
#loads all of the log2 FC samples, filtered for TPM > 1 in both WT and KD
MF_alltrans = MAGOH_full_alltrans %>% full_join(EIF4A3_full_alltrans) %>% full_join(UPF1_full_alltrans) %>% 
  full_join(RBM22_full_alltrans) %>% full_join(AQR_full_alltrans) %>% full_join(SNRNP200_full_alltrans) %>% 
  full_join(EFTUD2_full_alltrans) %>% full_join(SF3B1_full_alltrans) %>% full_join(SF3B3_full_alltrans) %>% 
  full_join(SNRPC_full_alltrans) %>% full_join(SNRNP70_full_alltrans) %>% full_join(PRPF8_full_alltrans) %>% 
  full_join(PRPF6_full_alltrans) %>% full_join(CDC5L_full_alltrans) %>% full_join(SF3A1_full_alltrans) %>% 
  full_join(SF3A3_full_alltrans) %>% full_join(U2AF1_full_alltrans) %>% full_join(CDC40_full_alltrans) %>% 
  full_join(PRPF3_full_alltrans) %>% full_join(PRPF4_full_alltrans) %>% full_join(GNB2L1_full_alltrans)#Combines all tables
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
  geom_boxplot(aes(x = factor(Sample, levels = Pres_KD),
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
  geom_label(data = PE_pres_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = Pres_KD),
                 y = 3,
                 label = paste0("p = ",signif(P,digits = 3))),
             size = 5,
             alpha = 0.75)+
  geom_label(data = PE_pres_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = Pres_KD),
                 y = Med,
                 color = isoform,
                 label = round(Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label(data = PE_pres_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = Pres_KD),
                 y = -3,
                 color = isoform,
                 label = paste0("n =",n)),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
  geom_label(data = PE_FC_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = paste0("p = ",signif(P,digits = 3))),
             size = 5,
             alpha = 0.75)+
  geom_label(data = PE_FC_sum %>% filter(Sample %in% all_GOI),
             aes(x = factor(Sample, levels = all_GOI),
                 y = Med,
                 color = isoform,
                 label = round(Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label(data = PE_FC_sum %>% filter(Sample %in% all_GOI),
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3,
                 color = isoform,
                 label = paste0("n =",n)),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
  geom_label(data = Saltzman_alltrans_sum,
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 5,
             alpha = 0.75)+
  geom_label(data = Saltzman_alltrans_sum %>% filter(Sample %in% all_GOI),
             aes(x = factor(Sample, levels = all_GOI),
                 y = Med,
                 color = Isoform,
                 label = round(Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label(data = Saltzman_alltrans_sum %>% filter(Sample %in% all_GOI),
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
  geom_label(data = Saltzman_alltrans_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = signif(P,digits = 3)),
             size = 5,
             alpha = 0.75)+
  geom_label(data = Saltzman_alltrans_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = Med,
                 color = Isoform,
                 label = round(Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label(data = Saltzman_alltrans_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = -3,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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

#### Make TPM plot for AS genes vs not####
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
                              "transcript_biotype"),
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
             size = 8)+
  geom_label(data = AS_MANE_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = AS_gene,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = AS_MANE_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = AS_gene,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 8) +
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
             size = 8,
             alpha = 0.75)+
  geom_label(data = AS_MANE_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = AS_gene,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = AS_MANE_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = AS_gene,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 8,
             alpha = 0.75) +
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
             size = 5,
             alpha = 0.75)+
  geom_label(data = no_AS_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = no_AS_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
             size = 5,
             alpha = 0.75)+
  geom_label(data = no_AS_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = no_AS_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
             size = 5,
             alpha = 0.75)+
  geom_label(data = AS_only_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = AS_only_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
             size = 5,
             alpha = 0.75)+
  geom_label(data = AS_only_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = AS_only_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
             size = 5,
             alpha = 0.75)+
  geom_label(data = noAS_noNMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = noAS_noNMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
             size = 5,
             alpha = 0.75)+
  geom_label(data = noAS_noNMD_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = noAS_noNMD_sum %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
             size = 5,
             alpha = 0.75)+
  geom_label(data = AS_noNMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = Isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label(data = AS_noNMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = Isoform,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
PTC_list = read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/Stringent_PTC_MANE_CE.csv")
PTC_list = PTC_list %>% select(2:4)
noAS_NMD = No_AS_only %>% inner_join(PTC_list,by = c("ENST.ID" = "transID"))
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
             size = 5,
             alpha = 0.75) +
  geom_label(data = noAS_NMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = PTC,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = noAS_NMD_sum,
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = PTC,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
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
             size = 8,
             alpha = 0.75)+
  geom_label(data = Disease_PTC_sum,
             aes(x = Sample,
                 y = med,
                 color = PTC,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = Disease_PTC_sum,
             aes(x = Sample,
                 y = -3,
                 color = PTC,
                 label = n),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 8) +
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
             size = 8)+
  geom_label(data = PE_DIS_sum,
             aes(x = Sample,
                 y = Med,
                 color = isoform,
                 label = round(Med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = PE_DIS_sum,
             aes(x = Sample,
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
DIS_PE_boxplot
ggsave("Disease_PE_boxplot.pdf",
       device = pdf,
       plot = DIS_PE_boxplot,
       width = 10,
       height = 10,
       units = "in",
       dpi = 300)

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
            size = 8)+
  geom_label(data = Dis_no_NMD_sum,
             aes(x = Sample,
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = Dis_no_NMD_sum,
            aes(x = Sample,
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
            size = 8)+
  geom_label(data = disease_AS_sum,
             aes(x = Sample,
                 y = med,
                 color = factor(Gene_type,levels = c("non-AS","AS")),
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = disease_AS_sum,
            aes(x = Sample,
                y = -2.5,
                color = factor(Gene_type,levels = c("non-AS","AS")),
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 8) +
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

#### Look at long 3'UTR based NMD targets ####
long_3UTRs <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/NMD_long_3UTRs.csv")

#Make the datatables
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_3UTR"),
         eval(parse(text = paste0(i,"_full_alltrans"))) %>% 
           inner_join(long_3UTRs, by = c("ENST.ID" = "ensembl_transcript_id")))
}

combined_3UTR = MAGOH_3UTR %>% full_join(EIF4A3_3UTR) %>% full_join(UPF1_3UTR) %>% 
  full_join(RBM22_3UTR) %>% full_join(AQR_3UTR) %>% full_join(SNRNP200_3UTR) %>% 
  full_join(EFTUD2_3UTR) %>% full_join(SF3B1_3UTR) %>% full_join(SF3B3_3UTR) %>% 
  full_join(SNRPC_3UTR) %>% full_join(SNRNP70_3UTR) %>% full_join(PRPF8_3UTR) %>%
  full_join(PRPF6_3UTR) %>% full_join(CDC5L_3UTR) %>% full_join(SF3A1_3UTR) %>% 
  full_join(SF3A3_3UTR) %>% full_join(U2AF1_3UTR) %>% full_join(CDC40_3UTR) %>% 
  full_join(PRPF3_3UTR) %>% full_join(PRPF4_3UTR) %>% full_join(GNB2L1_3UTR)
summary_3UTR = combined_3UTR %>% group_by(Sample,transcript) %>% summarise(n = n(),
                                                                           med = median(log2FoldChange),
                                                                           genes = n_distinct(ensembl_gene_id))

all_UTR_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_UTR_res"),
         wilcox.test(log2FoldChange ~ transcript, data = combined_3UTR %>% filter(Sample == i),
                     exact = FALSE, alternative = "less"))
  all_UTR_res = all_UTR_res %>% add_row(Sample = i, P = eval(parse(text = paste0(i,"_UTR_res$p.value"))))
}
summary_3UTR = summary_3UTR %>% left_join(all_UTR_res)

UTR_colors = c("MANE" = "#703049",
               "NMD" = "#EB9D28")
UTR_boxplot = ggplot(data = combined_3UTR)
UTR_boxplot = UTR_boxplot + 
  geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                   y = log2FoldChange,
                   fill = transcript),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = UTR_colors, labels = c("MANE" = "MANE",
                                                   "NMD" = "long 3' UTR")) +
  scale_color_manual(values = UTR_colors) +
  geom_label(data = summary_3UTR,
             aes(x = factor(Sample, levels = all_GOI),
                 y = 3,
                 label = paste0("p = ",signif(P,digits = 3))))+
  geom_label(data = summary_3UTR,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = transcript,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8)) +
  geom_label(data = summary_3UTR,
             aes(x = factor(Sample, levels = all_GOI),
                 y = -2.5,
                 color = transcript,
                 label = paste0("n =",n)),
             show.legend = F,
             position = position_dodge2(width = 0.9)) +
  theme_bw() + 
  labs(x = "Knockdown",
       y = "Log2(Fold Change)",
       fill = "Gene Type")+
  coord_cartesian(y = c(-3,3))
UTR_boxplot
ggsave("UTR_boxplot_full.pdf",
       device = pdf,
       plot = UTR_boxplot,
       width = 35,
       height = 10,
       units = "in",
       dpi = 300)

#Limit to the main figure transcripts
UTR_main_fig = ggplot(data = combined_3UTR %>% filter(Sample %in% Pres_KD))
UTR_main_fig = UTR_main_fig + 
  geom_boxplot(aes(x = factor(Sample,
                              levels = all_GOI),
                   y = log2FoldChange,
                   fill = transcript),
               position = position_dodge2(width = 0.9),
               width = 0.8,
               outlier.shape = 21,
               outlier.alpha = 0.5,
               outlier.colour = NA,
               linewidth = 1) +
  scale_fill_manual(values = UTR_colors,
                    labels = c("MANE" = "MANE",
                               "NMD" = "long 3' UTR")) +
  scale_color_manual(values = UTR_colors) +
  geom_label(data = summary_3UTR %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = 3,
                 label = paste0("p = ",signif(P,digits = 3))),
             size = 5,
             alpha = 0.75)+
  geom_label(data = summary_3UTR %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = med,
                 color = transcript,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_label(data = summary_3UTR %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample,
                            levels = all_GOI),
                 y = -2.5,
                 color = transcript,
                 label = paste0("n =",n)),
             show.legend = F,
             position = position_dodge2(width = 0.9),
             size = 5,
             alpha = 0.75) +
  theme_bw() + 
  labs(x = "Knockdown",
       y = "Log2(Fold Change)",
       fill = "Gene Type")+
  coord_cartesian(y = c(-3,3))
UTR_main_fig
ggsave("UTR_main_fig_FC_plot.pdf",
       device = pdf,
       plot = UTR_main_fig,
       width = 22,
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
            size = 8,
            color = "black")+
  geom_label(data = noPTC_NMD_summary,
             aes(x =factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = noPTC_NMD_summary,
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
            size = 8,
            color = "black")+
  geom_label(data = noPTC_noAS_NMD_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = isoform,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = noPTC_noAS_NMD_summary,
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
            size = 8,
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
            size = 8) +
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
            size = 8,
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
            size = 8) +
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

####Compare the LA no NMD isoforms to the MANE NMD isoforms####
NMD_MANE_only = all_MANE_NMD %>% filter(NMD == "TRUE") %>% mutate(Type = "NMD")
LA_nonNMD_MANE_only = all_LAnoNMD  %>% mutate(Type = "No_NMD")
NMD_LA_MANE = NMD_MANE_only %>% full_join(LA_nonNMD_MANE_only)

NMD_LA_MANE_summary = NMD_LA_MANE %>% group_by(Sample, Type) %>% 
  summarise(n = n(),
            med = median(log2FoldChange))
NMD_LA_res = tibble(Sample = character(),P = numeric())
for (i in all_GOI) {
  print(i)
  assign(paste0(i,"_NMD_LA_res"),
         wilcox.test(log2FoldChange ~ Type, data = NMD_LA_MANE %>% filter(Sample == i),
                     exact = FALSE, alternative = "greater")) #greater because we expect the NMD to be higher
  NMD_LA_res = NMD_LA_res %>% add_row(Sample = i,P = eval(parse(text = paste0(i,"_NMD_LA_res$p.value"))))
}
NMD_LA_MANE_summary = NMD_LA_MANE_summary %>% left_join(NMD_LA_res)

NMD_LA_colors = c("NMD" = "#F27D2E",
                  "No_NMD" = "#B27092")
NMD_LA_Boxplot = ggplot(data = NMD_LA_MANE)
NMD_LA_Boxplot = NMD_LA_Boxplot + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                                   y = log2FoldChange,
                                                   fill = Type),
                                               position = position_dodge2(width = 0.9),
                                               width = 0.8,
                                               outlier.shape = 21,
                                               outlier.alpha = 0.5,
                                               outlier.colour = NA,
                                               linewidth = 1) +
  scale_fill_manual(values = NMD_LA_colors, labels = c("NMD" = "NMD",
                                                       "No_NMD" = "Non NMD")) +
  scale_color_manual(values = NMD_LA_colors, labels = c("NMD" = "NMD",
                                                        "No_NMD" = "Non NMD")) +
  geom_text(data = NMD_LA_MANE_summary %>% filter(Type == "NMD"),
            aes(x = factor(Sample, levels = all_GOI),
                y = 3,
                label = signif(P,digits = 3)),
            size = 8,
            color = "black") +
  geom_label(data = NMD_LA_MANE_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = Type,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = NMD_LA_MANE_summary,
            aes(x = factor(Sample, levels = all_GOI),
                y = -2.5,
                color = Type,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 8) +
  labs(x = "Depletion",
       y = "log2(Fold Change)",
       fill = "Transcript Type",
       title = "MANE NMD transcripts vs LA no-NMD")+
  coord_cartesian(y = c(-4,4)) +
  theme_bw()
NMD_LA_Boxplot
ggsave("NMD_MANE_vs_LA_noNMD.pdf",
       plot = NMD_LA_Boxplot,
       device = pdf,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

MainFig_NMD_LA_Boxplot = ggplot(data = NMD_LA_MANE %>% filter(Sample %in% Pres_KD))
MainFig_NMD_LA_Boxplot = MainFig_NMD_LA_Boxplot + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                                   y = log2FoldChange,
                                                   fill = Type),
                                               position = position_dodge2(width = 0.9),
                                               width = 0.8,
                                               outlier.shape = 21,
                                               outlier.alpha = 0.5,
                                               outlier.colour = NA,
                                               linewidth = 1) +
  scale_fill_manual(values = NMD_LA_colors, labels = c("NMD" = "NMD",
                                                       "No_NMD" = "Non NMD")) +
  scale_color_manual(values = NMD_LA_colors, labels = c("NMD" = "NMD",
                                                        "No_NMD" = "Non NMD")) +
  geom_text(data = NMD_LA_MANE_summary %>% filter(Type == "NMD" & Sample %in% Pres_KD),
            aes(x = factor(Sample, levels = all_GOI),
                y = 3,
                label = signif(P,digits = 3)),
            size = 8,
            color = "black") +
  geom_label(data = NMD_LA_MANE_summary %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = Type,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = NMD_LA_MANE_summary %>% filter(Sample %in% Pres_KD),
            aes(x = factor(Sample, levels = all_GOI),
                y = -2.5,
                color = Type,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 8) +
  labs(x = "Depletion",
       y = "log2(Fold Change)",
       fill = "Transcript Type",
       title = "MANE NMD transcripts vs LA no-NMD")+
  coord_cartesian(y = c(-4,4)) +
  theme_bw()
MainFig_NMD_LA_Boxplot
ggsave("MainFig_NMD_MANE_vs_LA_noNMD.pdf",
       plot = MainFig_NMD_LA_Boxplot,
       device = pdf,
       width = 20,
       height = 10,
       units = "in",
       dpi = 300)

##Look at MANE vs PTC vs LA no NMD##
PTC_only = PTC_combined %>% mutate(Type = if_else(PTC == "TRUE",
                                                  "PTC",
                                                  "MANE")) %>%
  select(1:8,10,12)
PTC_MANE_noNMD = LA_nonNMD_MANE_only %>% full_join(PTC_only)
PTC_MANE_noNMD_summary = PTC_NMDmane_noNMD %>% group_by(Sample,Type) %>% 
  summarise(n = n(),
            med = median(log2FoldChange))

NMD_PTC_LA_colors = c("MANE" = "#663171",
                      "No_NMD" = "#FE64A3",
                      "PTC" = "#EA7428")
NMD_PTC_LA_Boxplot = ggplot(data = PTC_MANE_noNMD)
NMD_PTC_LA_Boxplot = NMD_PTC_LA_Boxplot + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                                   y = log2FoldChange,
                                                   fill = Type),
                                               position = position_dodge2(width = 0.9),
                                               width = 0.8,
                                               outlier.shape = 21,
                                               outlier.alpha = 0.5,
                                               outlier.colour = NA,
                                               linewidth = 1) +
  scale_fill_manual(values = NMD_PTC_LA_colors, labels = c("NMD" = "MANE NMD",
                                                            "No_NMD" = "Non NMD",
                                                            "PTC" = "PTC NMD")) +
  scale_color_manual(values = NMD_PTC_LA_colors, labels = c("NMD" = "MANE NMD",
                                                        "No_NMD" = "Non NMD",
                                                        "PTC"="PTC NMD")) +
  geom_label(data = PTC_MANE_noNMD_summary,
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = Type,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 3) +
  geom_text(data = PTC_MANE_noNMD_summary,
            aes(x = factor(Sample, levels = all_GOI),
                y = -2.5,
                color = Type,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7) +
  labs(x = "Depletion",
       y = "log2(Fold Change)",
       fill = "Transcript Type",
       title = "MANE isoforms vs PTC isoforms vs LA no-NMD transcripts")+
  coord_cartesian(y = c(-4,4)) +
  theme_bw()
NMD_PTC_LA_Boxplot
ggsave("MANE_vs_PTC_vs_LA_noNMD.pdf",
       plot = NMD_PTC_LA_Boxplot,
       device = pdf,
       width = 40,
       height = 10,
       units = "in",
       dpi = 300)

NMD_PTC_LA_Boxplot_mainfig = ggplot(data = PTC_MANE_noNMD %>% filter(Sample %in% Pres_KD))
NMD_PTC_LA_Boxplot_mainfig = NMD_PTC_LA_Boxplot_mainfig + geom_boxplot(aes(x = factor(Sample, levels = all_GOI),
                                                           y = log2FoldChange,
                                                           fill = Type),
                                                       position = position_dodge2(width = 0.9),
                                                       width = 0.8,
                                                       outlier.shape = 21,
                                                       outlier.alpha = 0.5,
                                                       outlier.colour = NA,
                                                       linewidth = 1) +
  scale_fill_manual(values = NMD_PTC_LA_colors, labels = c("NMD" = "MANE NMD",
                                                           "No_NMD" = "Non NMD",
                                                           "PTC" = "PTC NMD")) +
  scale_color_manual(values = NMD_PTC_LA_colors, labels = c("NMD" = "MANE NMD",
                                                            "No_NMD" = "Non NMD",
                                                            "PTC"="PTC NMD")) +
  geom_label(data = PTC_MANE_noNMD_summary %>% filter(Sample %in% Pres_KD),
             aes(x = factor(Sample, levels = all_GOI),
                 y = med,
                 color = Type,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 3) +
  geom_text(data = PTC_MANE_noNMD_summary %>% filter(Sample %in% Pres_KD),
            aes(x = factor(Sample, levels = all_GOI),
                y = -2.5,
                color = Type,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 7) +
  labs(x = "Depletion",
       y = "log2(Fold Change)",
       fill = "Transcript Type",
       title = "MANE isoforms vs PTC isoforms vs LA no-NMD transcripts")+
  coord_cartesian(y = c(-4,4)) +
  theme_bw()
NMD_PTC_LA_Boxplot_mainfig
ggsave("main_fig_MANE_vs_PTC_vs_LA_noNMD.pdf",
       plot = NMD_PTC_LA_Boxplot_mainfig,
       device = pdf,
       width = 20,
       height = 10,
       units = "in",
       dpi = 300)

#### Look at MANE vs PTC vs Other genes from the same set of genes ####
PTC_genes = getBM(attributes = "ensembl_gene_id",
                  filters = "ensembl_transcript_id",
                  values = PTC_list$transID,
                  mart = ensembl)
PTC_isoforms = getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","transcript_biotype",
                                    "transcript_mane_select"),
                     filters = "ensembl_gene_id",
                     values = PTC_genes$ensembl_gene_id,
                     mart = ensembl)
all_PTCgenes = PTC_list %>% 
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
                y = -2.5,
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
                      y = -2.5,
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
                      y = -2.5,
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
                                            "transcript_biotype"),
                             filters = "ensembl_gene_id",
                             values = shared_AS$GeneID,
                             mart = ensembl)
Shared_AS_NMD_iso = shared_AS_annotation %>% filter(transcript_biotype == "nonsense_mediated_decay")
