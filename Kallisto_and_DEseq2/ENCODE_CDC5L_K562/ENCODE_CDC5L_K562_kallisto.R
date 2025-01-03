#Written by: Caleb Embree

#Modified by: Caleb Embree

#Modified on: 5/31/2023

#Analyzing: CDC5L_K562

#Kallisto and intial analysis conducted by Andreas Stephanu

#All packages already installed#

#if (!requireNamespace("BiocManager", quietly = TRUE))

#install.packages("BiocManager")

#BiocManager::install("tximport")

#BiocManager::install("rhdf5")

#BiocManager::install("biomaRt")

#BiocManager::install("DESeq2")

#BiocManager::install("pheatmap")

#Things to use find and replace on

#C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_CDC5L_K562 --> working directory of the project

#WT1 --> SRR22521386

#WT2 --> SRR22521387

#KD1 --> SRR22520919

#KD2 --> SRR22520920

#CDC5L_K562 --> CDC5L_K562

#Load all packages

library(readxl)

library(tximport)

library(rhdf5)

library(biomaRt)

library(DESeq2)

library(ggplot2)

library(ggrepel)

library(RColorBrewer)

library(MetBrewer)

library(pheatmap)

library(textshape)

library(janitor)

library(eulerr)

library(gghighlight)

library(tidyverse)

library(dplyr) #dplyr should be loaded last so other packages don't mask it's functions

setwd(dir = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_CDC5L_K562")

#Load File containing sample types

samples <- read_excel("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_CDC5L_K562/samples.xlsx")

View(samples)

#Load in the files

files1 <- file.path("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_CDC5L_K562","koutput", samples$Run, "abundance.h5")

files1

names(files1) <- paste0(c("WT1","WT2","KD1","KD2"))
files1

#importing output files from each library to a compiled dataframe

#txi.kallisto <- tximport(files1, type = "kallisto", txOut = TRUE)

txi.kallisto.tsv <- tximport(files1, type = "kallisto", txOut=TRUE, countsFromAbundance = "lengthScaledTPM")

head(txi.kallisto.tsv$counts)

counts<-as.data.frame(txi.kallisto.tsv$counts)

##split the rownames and keep ENST_ID

counts_names <- str_split(string = rownames(counts),pattern = "\\|",simplify = TRUE)

counts_names <- as.data.frame(counts_names[,1])

counts$ENST_ID <- counts_names[, 1]

rownames(counts) <- counts$ENST_ID

counts$ENST_ID <- c()

txi.kallisto.tsv$counts <- as.matrix(counts)

#Filtering out transcripts with mean tpm<1

counts$mean=rowMeans(counts[,c("WT1","WT2")], na.rm=TRUE) #Filter WT samples

counts=filter(counts, mean>=1)

counts$mean=rowMeans(counts[,c("KD1", "KD2")], na.rm=TRUE) #Fiter KD samples

counts=filter(counts, mean>=1)

counts$mean <- NULL

#######Finding isoform specific counts associated with a gene in the dataset.####

listMarts()

ensembl <- useMart("ensembl",host = "https://feb2023.archive.ensembl.org") #use version 109 of ensembl

ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

attributes <- listAttributes(ensembl)

listgenes <- getBM(attributes = c("external_gene_name","ensembl_transcript_id_version","external_transcript_name"),
                   
                   filters = "ensembl_transcript_id_version",
                   values = rownames(counts),
                   mart = ensembl)
counts$ensembl_transcript_id_version<-rownames(counts)

listgenes <- inner_join(counts, listgenes, by = "ensembl_transcript_id_version")



####Differential expression analysis#####

#sampleinfo as dataframe

colData <- data.frame(genotype = factor(c("WT","WT","KD","KD"))) #make sure these are in the right order that match your sample file

ddscounts.table <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, formula(~ genotype))

ddscounts.table = ddscounts.table[rownames(ddscounts.table) %in% rownames(counts)]

ddscounts <- DESeq(ddscounts.table)

#Make sure data is associated with correct genotype

colData(ddscounts)

#PCAplots
vsd <- vst(ddscounts, blind = FALSE)
rld <- rlog(ddscounts)
plotPCA(vsd, intgroup = "genotype")
#with sample labels
plotPCA(vsd, intgroup = "genotype")+
  geom_text_repel(aes(label=name))

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
#log2fc range : -8.31343 to +8.821

#Exporting files
write.csv(as.data.frame(res_counts), file = "KDvsWT_full_results.csv")
write.csv(as.data.frame(res.sig), file = "2KDvsWT_sig_results.csv")
write.csv(counts(ddscounts, normalized = T), file = "KDvsWT_normalized_counts.csv")

res_counts1 = read.csv("KDvsWT_full_results.csv", header = T)
res_counts1 = column_to_rownames(res_counts1, "X")
head(counts)
res_counts1=row_to_names(res_counts1, row_number=1, remove_row = TRUE, remove_rows_above = TRUE)
head(counts)
up.trans <- as.data.frame(subset(res.sig, res.sig$log2FoldChange> 0))



transcriptids=row.names(res.sig)
transcriptids1=as.vector(transcriptids)
trans.detail <- getBM(attributes = c("ensembl_transcript_id_version","ensembl_gene_id_version",
                                     "external_gene_name","transcript_biotype"),
                      filters = "ensembl_transcript_id_version",
                      values = transcriptids1,
                      mart = ensembl)
sig.trans <- as.data.frame(res.sig)
NMD.trans <- getBM(attributes = c("ensembl_transcript_id_version","ensembl_gene_id_version",
                                  "external_gene_name","transcript_biotype"),
                   filters = "transcript_biotype",
                   values = "nonsense_mediated_decay",
                   mart = ensembl)
#total NMD trans = 19738
#total significantly changing genes = 11854
sigNMD=intersect(row.names(sig.trans), NMD.trans$ensembl_transcript_id_version)
#intersect between NMDgenes and significantly changing genes in RNAseq(wt vs kd) = 950
re.alltrans<-as.data.frame(res_counts)

#### ploting venn diagrams ####
set1=intersect(NMD.trans$ensembl_transcript_id_version, rownames(re.alltrans))
set2=row.names(sig.trans)
set3=row.names(up.trans)

set.seed(1)
s <- list("NMD biotype" = set1,
          "Significant transcripts" = set2,
          "Upregulated transcripts" = set3)
plot(euler(s), quantities = TRUE,
     main="CDC5L knockout vs control")



#### Creating alltrans datatables ####
alltrans <- as.data.frame(res_counts)
alltrans = alltrans %>% filter(!is.na(padj))

alltransbiotype <- getBM(attributes = c("ensembl_transcript_id_version","external_gene_name",
                                        "transcript_biotype"),
                         filters = "ensembl_transcript_id_version",
                         values = rownames(res_counts),
                         mart = ensembl)

alltrans<-alltrans[order(rownames(alltrans)),]
alltransbiotype<-alltransbiotype[order(alltransbiotype$ensembl_transcript_id_version),]
alltrans = rownames_to_column(alltrans, var = "transcript_id")
alltransbioNMD <- inner_join(alltrans,alltransbiotype, by=c("transcript_id" = "ensembl_transcript_id_version")) #Check this for gene FC

#### CDF plot comparing NMD biotype ####
write(alltransbiotype$ensembl_transcript_id_version, "tids.txt")
alltransbioNMD = alltransbioNMD %>% filter(transcript_biotype == "protein_coding" | transcript_biotype == "nonsense_mediated_decay")

NMDBio_cdf = ggplot(alltransbioNMD, aes(log2FoldChange, colour=transcript_biotype))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "CDC5L KD vs WT", subtitle = "NMD biotype")
NMDBio_cdf

#### analysis with robert's PTC+ and PTC- list ####
ENST_PTC.EPI.TFG <- read.delim("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_CDC5L_K562/ENST_PTC-EPI-TFG.txt")
alltrans$ENST.ID <- sub("\\..*", "", alltrans$transcript_id)

alltransPTC<-inner_join(alltrans, ENST_PTC.EPI.TFG, by = "ENST.ID")

PTC_cdf = ggplot(alltransPTC, aes(log2FoldChange, colour=PTC.Status))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "CDC5L KD vs WT", subtitle = "PTC list")
PTC_cdf

resPTC <- wilcox.test(log2FoldChange ~ PTC.Status, data = alltransPTC,
                      exact = FALSE, alternative = "less")

resPTC
boxplot(log2FoldChange ~ PTC.Status, data=alltransPTC)

#### Analysis using the NMD+- transcript list ####
NMD_transcripts = read.delim("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_CDC5L_K562/NMD_transcript_list.txt")
alltransNMD = inner_join(alltrans, NMD_transcripts, by="transcript_id")

NMD_cdf = ggplot(alltransNMD, aes(log2FoldChange, colour=NMD))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "CDC5L KD vs WT", subtitle = "NMD +/- list")
NMD_cdf
NMD_res <- wilcox.test(log2FoldChange ~ NMD, data = alltransNMD,
                       exact = FALSE, alternative = "less")
NMD_res #p-vlaue = 


#### Analysis using the stringent NMD list ####
StringentNMD = read.csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_CDC5L_K562/Stringent_NMD_list_CE")
alltransSPTC = inner_join(alltrans, sPTC, by = c("ENST.ID" = "transID")) 
sPTC_res <- wilcox.test(log2FoldChange ~ PTC, data = alltransSPTC,
                        exact = FALSE, alternative = "less")
head(alltransSPTC)
count(alltransSPTC, PTC) #978 PTC false 319 PTC true
sPTC_res #p-vlaue = 2.411e-3

sPTC_cdf = ggplot(alltransSPTC, aes(log2FoldChange, colour=PTC))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "CDC5L KD vs WT", subtitle = "Stringent PTC list")
sPTC_cdf = sPTC_cdf + annotate("text", x=-2.5,y=1, label = "p-value = 2.411e-3") +
  annotate("text", x = -2.5, y = 0.95, label = "PTC+ = 319") +
  annotate("text", x = -2.5, y = 0.9, label = "PTC- = 978") +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey")
sPTC_cdf
ggsave("CDC5L_K562_sPTC_CDF.pdf", 
       plot = sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

# Make a CDF using the stringent PTC list 
sPTC = read.csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_CDC5L_K562/Stringent_PTC_CE") 
alltransSPTC = inner_join(alltrans, sPTC, by = c("ENST.ID" = "transID")) 
sPTC_res <- wilcox.test(log2FoldChange ~ PTC, data = alltransSPTC,
                        exact = FALSE, alternative = "less")
head(alltransSPTC)
count(alltransSPTC, PTC) #978 PTC false 319 PTC true
sPTC_res #p-vlaue = 2.411e-3

plot_colors = c("FALSE" = "#673272", "TRUE" = "#EA7428",
                "PE" = "#76CEB5", "MANE" = "#BCB29F")

sPTC_cdf = ggplot(alltransSPTC, aes(log2FoldChange, colour=PTC))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = plot_colors) +
  labs(y = "Cumulative Frequency")
sPTC_cdf = sPTC_cdf + annotate("text", x=-2.5,y=1, label = "CDC5L KD vs WT", size = 5, hjust = 0) +
  annotate("text", x=-2.5,y=.95, label = "Stringent PTC containing transcripts", hjust =0) +
  annotate("text", x=-2.5,y=.9, label = "p-value = 2.411e-3", hjust = 0) +
  annotate("text", x = -2.5, y = 0.85, label = "PTC+ = 319", hjust = 0) +
  annotate("text", x = -2.5, y = 0.8, label = "PTC- = 978", hjust = 0) +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey") +
  theme(legend.position = c(0.85,0.1))
sPTC_cdf
ggsave("CDC5L_K562_sPTC_CDF.pdf", 
       plot = sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)
ggsave("C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/CDC5L_K562_sPTC_CDF.pdf", 
       plot = sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

#### Make dot plot comparing the NMD plus and minus FC ####
CDC5L_x = alltransSTNMD %>% filter(NMD == "TRUE")  %>% #Makes a data table of the top NMD transcript for each gene
  select(Gene,log2FoldChange,baseMean) %>% 
  rename(plusFC = log2FoldChange) %>% 
  group_by(Gene) %>% 
  slice_max(baseMean, n = 1)

CDC5L_y  = alltransSTNMD %>% filter(NMD == "FALSE")  %>% #Makes a data table fo the top NMD- transcript for each gene
  select(Gene,log2FoldChange,baseMean) %>% 
  rename(minusFC = log2FoldChange) %>% 
  group_by(Gene) %>% 
  slice_max(baseMean, n = 1)

CDC5L_scatter = full_join(CDC5L_x, CDC5L_y, by = "Gene", suffix = c("plus","minus")) #Combine the two top tables
CDC5L_scatter = CDC5L_scatter %>% filter(!is.na(minusFC)) %>% 
  filter(!is.na(plusFC)) #Removes any NA rows


dot = ggplot(data = CDC5L_scatter)+
  geom_point(aes(x = plusFC, y = minusFC), color = "blue")+
  geom_abline(aes(intercept = 0, slope = 1))+ #Makes a graph where each dot represents the NMD+ FC (x-axis) and NMD- FC (y-axis) for a specific gene
  gghighlight(abs(plusFC) > 0.58, unhighlighted_params = list(alpha("black", 0.4)))+ #Highlights all of the genes that have more than a 1.5 log2FC in either direction in the NMD+ samples
  labs(x = "Log2FC NMD+", y = "Log2FC NMD-", title = "CDC5L gene-level FC comparison")+
  theme_minimal()+
  coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) #Use if you need to restrict the axis length
dot


#Count the number of transcripts that are significantly changed
counts = alltransSTNMD %>% filter(padj < 0.05) %>% 
  filter(abs(log2FoldChange) >0.58) #All transcripts with a fold change higher than 1.5
count(counts, NMD)
#Sig_NMD+ = 
#Sig_NMD- = 

#### Determine the effect on NMD when Alternately spliced genes are removed (Delete this section if rMATS hasn't been done) ####
# ASgenes_CDC5L <- read_delim("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/rMATS/ENCODE_CDC5L/ASgenes_CDC5L.csv", 
#                           delim = ";", escape_double = FALSE, trim_ws = TRUE)
# ASgenes_CDC5L = ASgenes_CDC5L %>% select(GeneID,geneSymbol)
# sPTC_gene <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id"),
#                    filters = "ensembl_transcript_id",
#                    values = sPTC$transID,
#                    mart = ensembl)
# sPTC_gene = sPTC_gene %>% inner_join(sPTC, by = c("ensembl_transcript_id" = "transID"))
# sPTC_AS = sPTC_gene %>% anti_join(ASgenes_CDC5L, by = c("ensembl_gene_id" = "GeneID"))
# count(sPTC_AS, PTC) #PTC false=  PTC true= 
# alltransAS = inner_join(alltrans, sPTC_AS, by = c("ENST.ID" = "ensembl_transcript_id"))
# count(alltransAS, PTC) # PTC false  PTC true
# 
# 
# asPTC_res <- wilcox.test(log2FoldChange ~ PTC, data = alltransAS,
#                          exact = FALSE, alternative = "less")
# head(alltransAS)
# asPTC_res #p-vlaue = 
# 
# asPTC_cdf = ggplot(alltransAS, aes(log2FoldChange, colour=PTC))+
#   stat_ecdf(linewidth=2)+
#   theme_minimal()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
#   coord_cartesian(xlim = c(-3,3)) +
#   scale_color_manual(values = met.brewer("Java",2)) +
#   labs(y = "Cumulative Frequency", title = "CDC5L KD vs WT", subtitle = "Non-AS transcripts")
# asPTC_cdf = asPTC_cdf + annotate("text", x=-2.5,y=1, label = "p-value = ") +
#   annotate("text", x = -2.5, y = 0.95, label = "PTC+ = ") +
#   annotate("text", x = -2.5, y = 0.9, label = "PTC- = ") +
#   geom_vline(xintercept = 0, color = "grey") + 
#   geom_hline(yintercept = 0.5, color = "grey")
# asPTC_cdf
# ggsave("CDC5L_K562_asPTC_CDF.pdf", 
#        plot = asPTC_cdf,
#        scale = 1,
#        width = 8,
#        height = 6,
#        units = "in",
#        device = "pdf",
#        dpi = 300)

#### Look at FC of NMD factors ####
alltrans_full <- as.data.frame(res_counts)
alltrans_full = alltrans_full %>% filter(!is.na(padj))
alltrans_full = rownames_to_column(alltrans_full, var = "transcript_id")
alltrans_full$ENST.ID <- sub("\\..*", "", alltrans_full$transcript_id)

alltrans_ensembl <- getBM(attributes = c("ensembl_transcript_id_version","external_gene_name",
                                         "transcript_biotype","ensembl_transcript_id","transcript_mane_select"),
                          filters = "ensembl_transcript_id",
                          values = alltrans_full$ENST.ID,
                          mart = ensembl)
alltrans_annotated = alltrans_full %>% left_join(alltrans_ensembl, by = c("ENST.ID" = "ensembl_transcript_id"))
alltrans_MANE = alltrans_annotated %>% filter(!is.na(transcript_mane_select)) %>% filter(transcript_mane_select != "")
NMD_factors <- read_excel("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/NMD_factors.xlsx")

CDC5L_NMD = alltrans_MANE %>% inner_join(NMD_factors, by = c("external_gene_name" = "NMDfactor"))
write.csv(CDC5L_NMD, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/NMD_analysis/CDC5L_NMDfactor.csv",row.names = FALSE)

####Look at just the shared stringent PTC from ISAR####
shared_sPTC <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/shared_stringent_PTC.csv")
shared_SPTC_trans = inner_join(alltrans, shared_sPTC, by = c("ENST.ID" = "transID")) 
shared_sPTC_res <- wilcox.test(log2FoldChange ~ PTC, data = shared_SPTC_trans,
                               exact = FALSE, alternative = "less")
head(shared_SPTC_trans)
count(shared_SPTC_trans, PTC)#X PTC false Y PTC true
shared_sPTC_res #p-vlaue = 

shared_sPTC_cdf = ggplot(shared_SPTC_trans, aes(log2FoldChange, colour=PTC))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = plot_colors) +
  labs(y = "Cumulative Frequency", title = "CDC5L KD vs WT", subtitle = "Shared stringent PTC transcripts")
shared_sPTC_cdf = shared_sPTC_cdf + annotate("text", x=-2.5,y=1, label = "p-value = ") +
  annotate("text", x = -2.5, y = 0.95, label = "PTC+ = ") +
  annotate("text", x = -2.5, y = 0.9, label = "PTC- = ") +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey")
shared_sPTC_cdf
ggsave("CDC5L_K562_shared_PTC_CDF.pdf", 
       plot = shared_sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

#### Look at PE from (Thomas et al 2020 Nat. Gen.) ####
PE_archive_list <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Poison_Exons/archive_PE_list.csv")
PEa_FC = alltrans %>% inner_join(PE_archive_list, by = c("ENST.ID" = "ensembl_transcript_id"))
PEa_res <- wilcox.test(log2FoldChange ~ isoform, data = PEa_FC,
                       exact = FALSE, alternative = "less")
head(PEa_FC)
count(PEa_FC, isoform)# 391 PE  224 MANE
PEa_res #p-vlaue = 0.2619

PEa_cdf = ggplot(PEa_FC, aes(log2FoldChange, colour=isoform))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = plot_colors) +
  labs(y = "Cumulative Frequency")
PEa_cdf = PEa_cdf + annotate("text", x=-2.5,y=1, label = "CDC5L KD vs WT", size = 5, hjust = 0) +
  annotate("text", x=-2.5,y=.95, label = "Transcripts with a poison exon", hjust =0) +
  annotate("text", x=-2.5,y=.9, label = "p-value = 0.2619", hjust = 0) +
  annotate("text", x = -2.5, y = 0.85, label = "PE = 391", hjust = 0) +
  annotate("text", x = -2.5, y = 0.8, label = "MANE = 224", hjust = 0) +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey") +
  theme(legend.position = c(0.85,0.1))
PEa_cdf
ggsave("CDC5L_K562_PE_archive_CDF.pdf", 
       plot = PEa_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)
ggsave("C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/CDC5L_K562_PEarchive_CDF.pdf", 
       plot = PEa_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

##Redo the sPTC analysis using the MANE list
sPTC_MANE = read.csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/Stringent_PTC_MANE_CE.csv")
alltransSPTC_MANE = inner_join(alltrans, sPTC_MANE, by = c("ENST.ID" = "transID")) 
sPTC_MANE_res <- wilcox.test(log2FoldChange ~ PTC, data = alltransSPTC_MANE,
                             exact = FALSE, alternative = "less")
count(alltransSPTC_MANE, PTC) # 261 PTC false 319 PTC true
sPTC_MANE_res #p-vlaue = 1.237e-6


MANE_sPTC_cdf = ggplot(alltransSPTC_MANE, aes(log2FoldChange, colour=PTC))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = plot_colors) +
  labs(y = "Cumulative Frequency")
MANE_sPTC_cdf = MANE_sPTC_cdf + annotate("text", x=-2.5,y=1, label = "CDC5L KD vs WT", size = 5, hjust = 0) +
  annotate("text", x=-2.5,y=.95, label = "PTC containing and MANE transcripts", hjust =0) +
  annotate("text", x=-2.5,y=.9, label = "p-value = 1.237e-6", hjust = 0) +
  annotate("text", x = -2.5, y = 0.85, label = "PTC+ = 319", hjust = 0) +
  annotate("text", x = -2.5, y = 0.8, label = "PTC- = 261", hjust = 0) +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey") +
  theme(legend.position = c(0.85,0.1))
MANE_sPTC_cdf
ggsave("CDC5L_K562_sPTC_MANE_CDF.pdf", 
       plot = MANE_sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)
ggsave("C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/CDC5L_K562_sPTC_MANE_CDF.pdf", 
       plot = MANE_sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

#### Pull the TPM of genes for the NMD TPM graph ####
#Check WT and KD are refering to the right columns in listgenes
full_tpm = listgenes %>% rowwise() %>%  mutate(WTmean = mean(c_across(1:2)), KDmean = mean(c_across(3:4))) %>% ungroup() #Check that WT and KD are referring to the right columns
full_tpm = full_tpm %>% mutate(ENSTID = str_remove(ensembl_transcript_id_version, "\\..*"))
full_tpm = full_tpm %>% select(6,8:10)
PTC_tpm = full_tpm %>% right_join(sPTC_MANE, by = c("ENSTID" = "transID"))
PTC_tpm = PTC_tpm %>% select(2:4,6,7)
write_csv(PTC_tpm, "CDC5L_PTC_MANE_TPM.csv")
write_csv(PTC_tpm, "C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/CDC5L_PTC_MANE_TPM.csv")

#### Look at the effect when only WT filtering####
WTcounts <-as.data.frame(txi.kallisto.tsv$counts)

#Filtering out transcripts with mean tpm<1
WTcounts$mean=rowMeans(WTcounts[,c("WT1","WT2")], na.rm=TRUE) #Filter WT samples
WTcounts=filter(WTcounts, mean>=1)

listMarts()
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
WTlistgenes <- getBM(attributes = c("external_gene_name","ensembl_transcript_id_version","external_transcript_name"),
                     filters = "ensembl_transcript_id_version",
                     values = rownames(WTcounts),
                     mart = ensembl)

WTcounts$ensembl_transcript_id_version<-rownames(WTcounts)
WTlistgenes <- inner_join(WTcounts, WTlistgenes, by = "ensembl_transcript_id_version")

WTddscounts.table <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, formula(~ genotype))
WTddscounts.table = WTddscounts.table[rownames(WTddscounts.table) %in% rownames(WTcounts)]
WTddscounts <- DESeq(WTddscounts.table)

#Make sure data is associated with correct genotype
colData(WTddscounts)
#PCAplots
WTvsd <- vst(WTddscounts, blind = FALSE)
WTrld <- rlog(WTddscounts)
plotPCA(WTvsd, intgroup = "genotype")
#with sample labels
plotPCA(WTvsd, intgroup = "genotype")+
  geom_text_repel(aes(label=name))

#####1st pairwise comparison
as.data.frame(colData(WTddscounts))
WTres_counts <- results(WTddscounts, contrast = c("genotype","KD","WT"))
WTres_counts
mcols(WTres_counts, use.names = TRUE)

#significant results
WTres.sig <- WTres_counts[ which(WTres_counts$padj < 0.05), ]
hist(WTres.sig$pvalue, col = "green1")
head(WTres.sig)
plotMA(WTres_counts, padj = TRUE, ylim=c(-6,6), main ="MA plot: WT vs KD")

plotMA(WTres.sig)
sum( WTres.sig$log2FoldChange < 0, na.rm=TRUE )#upregulated genes=
sum( WTres.sig$log2FoldChange > 0, na.rm=TRUE )#downregulated genes=
head(WTres.sig$log2FoldChange > 0)

head( WTres.sig[ order( WTres.sig$log2FoldChange ), ] )
head( WTres.sig[ order( WTres.sig$log2FoldChange ), ],20 )
tail( WTres.sig[ order( WTres.sig$log2FoldChange ), ],20 )
summary(WTres.sig$log2FoldChange)
#log2fc range : - to +

WTtranscriptids=row.names(WTres.sig)
WTtranscriptids1=as.vector(WTtranscriptids)
WTtrans.detail <- getBM(attributes = c("ensembl_transcript_id_version","ensembl_gene_id_version",
                                       "external_gene_name","transcript_biotype"),
                        filters = "ensembl_transcript_id_version",
                        values = WTtranscriptids1,
                        mart = ensembl)
WTsig.trans <- as.data.frame(WTres.sig)
#total significantly changing genes = 
WTsigNMD=intersect(row.names(WTsig.trans), NMD.trans$ensembl_transcript_id_version)
#intersect between NMDgenes and significantly changing genes in RNAseq(wt vs kd) = 
WTre.alltrans<-as.data.frame(WTres_counts)

WTalltrans <- as.data.frame(WTres_counts)
WTalltrans = WTalltrans %>% filter(!is.na(padj))
WTalltrans<-WTalltrans[order(rownames(WTalltrans)),]
WTalltrans = rownames_to_column(WTalltrans, var = "transcript_id")
write_csv(WTalltrans,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/CDC5L_WTfilt_alltrans.csv")
write_csv(alltrans,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/CDC5L_DBfilt_alltrans.csv")
write_csv(alltransSPTC_MANE,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/CDC5L_PTC_alltrans.csv")

####Make a list of the genes filtered to WT>5TPM ####
WT5counts <-as.data.frame(txi.kallisto.tsv$counts)

#Filtering out transcripts with mean tpm<5
WT5counts$mean=rowMeans(WT5counts[,c("WT1","WT2")], na.rm=TRUE) #Filter WT samples
WT5counts=filter(WT5counts, mean>=5)

WT5listgenes <- getBM(attributes = c("external_gene_name","ensembl_transcript_id_version","external_transcript_name"),
                      filters = "ensembl_transcript_id_version",
                      values = rownames(WT5counts),
                      mart = ensembl)

WT5counts$ensembl_transcript_id_version<-rownames(WT5counts)
WT5listgenes <- inner_join(WT5counts, WT5listgenes, by = "ensembl_transcript_id_version")

WT5ddscounts.table <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, formula(~ genotype))
WT5ddscounts.table = WT5ddscounts.table[rownames(WT5ddscounts.table) %in% rownames(WT5counts)]
WT5ddscounts <- DESeq(WT5ddscounts.table)

#Make sure data is associated with correct genotype
colData(WT5ddscounts)
#PCAplots
WT5vsd <- vst(WT5ddscounts, blind = FALSE)
WT5rld <- rlog(WT5ddscounts)
plotPCA(WT5vsd, intgroup = "genotype")
#with sample labels
plotPCA(WT5vsd, intgroup = "genotype")+
  geom_text_repel(aes(label=name))

#####1st pairwise comparison
as.data.frame(colData(WT5ddscounts))
WT5res_counts <- results(WT5ddscounts, contrast = c("genotype","KD","WT"))
WT5res_counts
mcols(WT5res_counts, use.names = TRUE)

#significant results
WT5res.sig <- WT5res_counts[ which(WT5res_counts$padj < 0.05), ]
hist(WT5res.sig$pvalue, col = "green1")
head(WT5res.sig)
plotMA(WT5res_counts, padj = TRUE, ylim=c(-6,6), main ="MA plot: WT vs KD")

plotMA(WT5res.sig)
sum( WT5res.sig$log2FoldChange < 0, na.rm=TRUE )#upregulated genes=5328
sum( WT5res.sig$log2FoldChange > 0, na.rm=TRUE )#downregulated genes=3505
head(WT5res.sig$log2FoldChange > 0)

head( WT5res.sig[ order( WT5res.sig$log2FoldChange ), ] )
head( WT5res.sig[ order( WT5res.sig$log2FoldChange ), ],20 )
tail( WT5res.sig[ order( WT5res.sig$log2FoldChange ), ],20 )
summary(WT5res.sig$log2FoldChange)
#log2fc range : --11.6032 to +7.7247

WT5transcriptids=row.names(WT5res.sig)
WT5transcriptids1=as.vector(WT5transcriptids)
WT5trans.detail <- getBM(attributes = c("ensembl_transcript_id_version","ensembl_gene_id_version",
                                        "external_gene_name","transcript_biotype"),
                         filters = "ensembl_transcript_id_version",
                         values = WT5transcriptids1,
                         mart = ensembl)
WT5sig.trans <- as.data.frame(WT5res.sig)
#total significantly changing genes = 
WT5sigNMD=intersect(row.names(WT5sig.trans), NMD.trans$ensembl_transcript_id_version)
#intersect between NMDgenes and significantly changing genes in RNAseq(WT5 vs kd) = 
WT5re.alltrans<-as.data.frame(WT5res_counts)

WT5alltrans <- as.data.frame(WT5res_counts)
WT5alltrans = WT5alltrans %>% filter(!is.na(padj))
WT5alltrans<-WT5alltrans[order(rownames(WT5alltrans)),]
WT5alltrans = rownames_to_column(WT5alltrans, var = "transcript_id")
write_csv(WT5alltrans,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/CDC5L_WT5filt_alltrans.csv")