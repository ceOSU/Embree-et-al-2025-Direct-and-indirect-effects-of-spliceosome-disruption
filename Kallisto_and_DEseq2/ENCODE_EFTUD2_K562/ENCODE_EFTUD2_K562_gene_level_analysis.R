#Written by: Caleb Embree
#Modified by: Caleb Embree
#Modified on: 9/19/2024
#Analyzing: Gene Level analysis of EFTUD2 KD

#All packages already installed#
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximport")
#BiocManager::install("rhdf5")
#BiocManager::install("biomaRt")
#BiocManager::install("DESeq2")
#BiocManager::install("pheatmap")

#Things to use find and replace on
# C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_EFTUD2_K562 --> working directory of the project
# WT1 --> Sample ID of first WT sample
# WT2 --> Sample ID of second WT sample
# KD1 --> Sample ID of first KD sample
# KD2 --> Sample ID of second KD sample
# EFTUD2 --> Ensemble name of your gene of interest

#Load all packages
library(readxl)
library(tximport)
library(biomaRt)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(MetBrewer)
library(MoMAColors)
library(pheatmap)
library(textshape)
library(janitor)
library(eulerr)
library(gghighlight)
library(tidyverse) #tidyverse should be loaded last so other packages don't mask dplyr's functions

#Load File containing sample types
getwd()
samples <- read_excel("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_EFTUD2_K562/EFTUD2_K562_samples.xlsx")
View(samples)


#### Load in the files ####
files1 <- file.path("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_EFTUD2_K562","koutput", samples$Run, "abundance.h5")
names(files1) <- paste0(c("KD1","KD2","WT1","WT2"))
files1

#load in the transcript to gene file
tx2gene <- read_csv("tx2gene.csv")

#importing output files from each library to a compiled dataframe
#txi.kallisto <- tximport(files1, type = "kallisto", txOut = TRUE)
txi.kallisto.tsv <- tximport(files1, type = "kallisto",tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
head(txi.kallisto.tsv$counts)
counts<-as.data.frame(txi.kallisto.tsv$counts)

#Filtering out transcripts with mean tpm<1
counts = counts %>% rownames_to_column(var = "ENSG_version") %>% 
  rowwise() %>% 
  mutate(WTmean = mean(c(WT1,WT2),na.rm = TRUE),
         KDmean = mean(c(KD1,KD2),na.rm = TRUE)) %>% 
  filter(WTmean >1 & KDmean >1) %>% ungroup()
counts$ENSG<- sub("\\..*", "", counts$ENSG_version)



#######Finding isoform specific counts associated with a gene in the dataset.####
ensembl <- useMart("ensembl",host = "https://feb2023.archive.ensembl.org") #use version 109 of ensembl
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
listgenes <- getBM(attributes = c("external_gene_name","ensembl_gene_id"),
                   filters = "ensembl_gene_id",
                   values = counts$ENSG,
                   mart = ensembl)

listgenes <- inner_join(counts, listgenes, by = c("ENSG" = "ensembl_gene_id"))

####Differential expression analysis#####
#sampleinfo as dataframe
colData <- data.frame(genotype = factor(c("KD","KD","WT","WT"))) #make sure these are in the right order that match your sample file

ddscounts.table <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, formula(~ genotype))
ddscounts.table = ddscounts.table[rownames(ddscounts.table) %in% counts$ENSG_version]
ddscounts <- DESeq(ddscounts.table)

#Make sure data is associated with correct genotype
colData(ddscounts)

#PCAplots
vsd <- vst(ddscounts, blind = FALSE)
rld <- rlog(ddscounts)
plotPCA(vsd, intgroup = "genotype")
#with sample labels
PCA = plotPCA(vsd, intgroup = "genotype")+
  geom_text_repel(aes(label=name))
PCA
ggsave("EFTUD2_geneLevel_PCA.pdf",
       plot = PCA,
       device = pdf,
       width = 5,
       height = 5,
       units = "in",
       dpi = 300)

#heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$genotype, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#Results
#####1st pairwise comparison
as.data.frame(colData(ddscounts))
res_counts <- results(ddscounts, contrast = c("genotype","KD","WT"))
res_counts
mcols(res_counts, use.names = TRUE)


#Histogram of p values
hist(res_counts$pvalue, col = "darkslategray1")
hist(res_counts$padj, col = "red1")

#significant results
res.sig <- res_counts[ which(res_counts$padj < 0.05), ]
hist(res.sig$pvalue, col = "green1")
head(res.sig)
plotMA(res_counts, padj = TRUE, ylim=c(-6,6), main ="MA plot: WT vs KD")

plotMA(res.sig)
sum( res.sig$log2FoldChange < 0, na.rm=TRUE )#upregulated genes=9860
sum( res.sig$log2FoldChange > 0, na.rm=TRUE )#downregulated genes=7206
head(res.sig$log2FoldChange > 0)

head( res.sig[ order( res.sig$log2FoldChange ), ] )
head( res.sig[ order( res.sig$log2FoldChange ), ],20 )
tail( res.sig[ order( res.sig$log2FoldChange ), ],20 )
summary(res.sig$log2FoldChange)
#log2fc range : -3.05793 to +4.32496

#Exporting files
write.csv(as.data.frame(res_counts), file = "KDvsWT_geneLevel_full_results.csv")
write.csv(as.data.frame(res.sig), file = "KDvsWT_geneLevel_sig_results.csv")
write.csv(counts(ddscounts, normalized = T), file = "KDvsWT_geneLevel_normalized_counts.csv")


#### Creating alltrans datatables ####
alltrans <- as.data.frame(res_counts)
alltrans = alltrans %>% filter(!is.na(padj)) %>% 
  rownames_to_column(var = "geneID_version")
alltrans = alltrans %>% mutate(ENSG.ID = str_remove(geneID_version,"\\..*")) 

alltransbiotype <- getBM(attributes = c("external_gene_name","ensembl_gene_id","gene_biotype"),
                         filters = "ensembl_gene_id",
                         values = alltrans$ENSG.ID,
                         mart = ensembl)


alltransbio <- left_join(alltrans,alltransbiotype, by=c("ENSG.ID" = "ensembl_gene_id")) #Check this for gene FC
alltrans_summary = alltransbio %>% group_by(gene_biotype) %>% summarise(n = n())

#### Look at different Classification of genes ####
noNMD <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/spliceosomeKD_analysis/noNMD_isoforms.csv")
SPTC_genes <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/spliceosomeKD_analysis/Stringent_PTC_geneID.csv")

alltrans_NMD = alltransbio %>% mutate(NMD = case_when(ENSG.ID %in% SPTC_genes$geneID ~ "NMD",
                                                      ENSG.ID %in% noNMD$ensembl_gene_id ~ "noNMD",
                                                      TRUE ~ NA))
write_csv(alltrans_NMD,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/EFTUD2_geneLevel_NMD_alltrans.csv")
alltrans__NMD_filt = alltrans_NMD %>% filter(!is.na(NMD))

NMD_summary = alltrans__NMD_filt %>% group_by(NMD) %>% summarise(n = n(),
                                                                 med = median(log2FoldChange))

NMD_res =wilcox.test(log2FoldChange ~ NMD, data = alltrans__NMD_filt,
                     exact = FALSE, alternative = "greater") #greater because we expect the NMD to be higher
NMD_summary = NMD_summary %>% mutate(P = NMD_res$p.value)


NMD_GL_colors = c("NMD" = "#F27D2E",
                  "noNMD" = "#B27092")
GL_NMD_box = ggplot(data = alltrans__NMD_filt)
GL_NMD_box = GL_NMD_box + geom_boxplot(aes(x = NMD,
                                           y = log2FoldChange,
                                           fill = NMD),
                                       position = position_dodge2(width = 0.9),
                                       width = 0.8,
                                       outlier.shape = 21,
                                       outlier.alpha = 0.5,
                                       outlier.colour = NA,
                                       linewidth = 1) +
  scale_fill_manual(values = NMD_GL_colors, labels = c("NMD" = "NMD",
                                                       "noNMD" = "Non NMD")) +
  scale_color_manual(values = NMD_GL_colors, labels = c("NMD" = "NMD",
                                                        "noNMD" = "Non NMD")) +
  geom_text(data = NMD_summary %>% filter(NMD == "NMD"),
            aes(x = NMD,
                y = 3,
                label = signif(P,digits = 3)),
            size = 8,
            color = "black") +
  geom_label(data = NMD_summary,
             aes(x = NMD,
                 y = med,
                 color = NMD,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = NMD_summary,
            aes(x = NMD,
                y = -2.5,
                color = NMD,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 8) +
  labs(y = "log2(Fold Change)",
       fill = "Contains NMD Transcript",
       title = "EFTUD2 KD Gene Level Analysis")+
  coord_cartesian(y = c(-4,4)) +
  theme_bw()
GL_NMD_box
ggsave("EFTUD2_geneLevel_NMD.pdf",
       plot = GL_NMD_box,
       device = pdf,
       width = 10,
       height = 10,
       units = "in",
       dpi = 300)
