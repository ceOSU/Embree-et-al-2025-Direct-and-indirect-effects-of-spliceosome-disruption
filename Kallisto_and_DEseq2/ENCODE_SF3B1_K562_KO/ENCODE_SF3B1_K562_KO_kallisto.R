#Written by: Caleb Embree
#Modified by: Caleb Embree
#Modified on: 3/12/2022
#Analyzing: 

#All packages already installed#
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximport")
#BiocManager::install("rhdf5")
#BiocManager::install("biomaRt")
#BiocManager::install("DESeq2")
#BiocManager::install("pheatmap")

#Things to use find and replace on
# C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_SF3B1_K562_KO --> working directory of the project
# SRR14838495 --> Sample ID of first WT sample
# SRR14838496 --> Sample ID of second WT sample
# SRR14842760 --> Sample ID of first KO sample
# SRR14842761 --> Sample ID of second KO sample
# SF3B1 --> Ensemble name of your gene of interest

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
library(tidyverse)
library(dplyr) #dplyr should be loaded last so other packages don't mask it's functions

#Load File containing sample types

samples <- read_excel("SF3B1_K562_KO_samples.xlsx")
View(samples)


#Load in the files
files1 <- file.path("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_SF3B1_K562_KO","koutput", samples$Run, "abundance.h5")
files1
names(files1) <- paste0(c("SRR14838495","SRR14838496","SRR14842760","SRR14842761"))


#importing output files from each library to a compiled dataframe
#txi.kallisto <- tximport(files1, type = "kallisto", txOut = TRUE)
txi.kallisto.tsv <- tximport(files1, type = "kallisto", txOut=TRUE, countsFromAbundance = "lengthScaledTPM")
head(txi.kallisto.tsv$counts)
counts<-as.data.frame(txi.kallisto.tsv$counts)

#Filtering out transcripts with mean tpm<1
counts$mean=rowMeans(counts[,c("SRR14838495","SRR14838496")], na.rm=TRUE) #Filter WT samples
counts=filter(counts, mean>=1)
counts$mean=rowMeans(counts[,c("SRR14842760","SRR14842761")], na.rm=TRUE) #Fiter KO samples
counts=filter(counts, mean>=1)
counts$mean <- NULL
counts$ensembl_transcript_id_version<-rownames(counts)
counts$ENSTID <- sub("\\..*", "", counts$ensembl_transcript_id_version)



#######Finding isoform specific counts associated with a gene in the dataset.####

listMarts()
ensembl <- useMart("ensembl",host = "https://feb2023.archive.ensembl.org") #use version 109 of ensembl
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
listgenes <- getBM(attributes = c("external_gene_name","ensembl_transcript_id","external_transcript_name"),
                   filters = "ensembl_transcript_id",
                   values = counts$ENSTID,
                   mart = ensembl)

listgenes <- inner_join(counts, listgenes, by = c("ENSTID" = "ensembl_transcript_id"))


#Analyzing isoform specific mean counts for a gene of interest (not essential)
SF3B1transcripts=filter(listgenes, external_gene_name == "SF3B1" )
SF3B1transcriptsWT<-SF3B1transcripts %>% dplyr::select(SRR14838495,SRR14838496, ensembl_transcript_id_version, external_gene_name, external_transcript_name)
SF3B1transcriptsWT$averagecounts = rowMeans(SF3B1transcriptsWT[,c("SRR14838495","SRR14838496")])
#put in specific library names above
SF3B1transcriptsWT$percent= SF3B1transcriptsWT$averagecounts/sum(SF3B1transcriptsWT$averagecounts)
SF3B1transcriptsWT<-SF3B1transcriptsWT %>% arrange(averagecounts)

#Repeat for the KO transcripts
SF3B1transcriptsKO<-SF3B1transcripts %>% dplyr::select(SRR14842760,SRR14842761, ensembl_transcript_id_version, external_gene_name, external_transcript_name)
SF3B1transcriptsKO$averagecounts = rowMeans(SF3B1transcriptsKO[,c("SRR14842760","SRR14842761")])
SF3B1transcriptsKO$percent= SF3B1transcriptsKO$averagecounts/sum(SF3B1transcriptsKO$averagecounts)
SF3B1transcriptsKO<-SF3B1transcriptsKO %>% arrange(averagecounts)

####Differential expression analysis#####

#sampleinfo as dataframe
colData <- data.frame(genotype = factor(c("WT","WT","KO","KO"))) #make sure these are in the right order that match your sample file


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
res_counts <- results(ddscounts, contrast = c("genotype","KO","WT"))
res_counts
mcols(res_counts, use.names = TRUE)


#Histogram of p values
hist(res_counts$pvalue, col = "darkslategray1")
hist(res_counts$padj, col = "red1")

#significant results
res.sig <- res_counts[ which(res_counts$padj < 0.05), ]
hist(res.sig$pvalue, col = "green1")
head(res.sig)
plotMA(res_counts, padj = TRUE, ylim=c(-6,6), main ="MA plot: WT vs KO")

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
write.csv(as.data.frame(res_counts), file = "KOvsWT_full_results.csv")
write.csv(as.data.frame(res.sig), file = "2KOvsWT_sig_results.csv")
write.csv(counts(ddscounts, normalized = T), file = "KOvsWT_normalized_counts.csv")

res_counts1 = read.csv("KOvsWT_full_results.csv", header = T)
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
#intersect between NMDgenes and significantly changing genes in RNAseq(wt vs KO) = 950
re.alltrans<-as.data.frame(res_counts)

##ploting venn diagrams
set1=intersect(NMD.trans$ensembl_transcript_id_version, rownames(re.alltrans))
set2=row.names(sig.trans)
set3=row.names(up.trans)

set.seed(1)
s <- list("NMD biotype" = set1,
          "Significant transcripts" = set2,
          "Upregulated transcripts" = set3)
plot(euler(s), quantities = TRUE,
     main="SF3B1 knockout vs control")




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

#CDF plot comparing NMD biotype 
write(alltransbiotype$ensembl_transcript_id_version, "tids.txt")
alltransbioNMD = alltransbioNMD %>% filter(transcript_biotype == "protein_coding" | transcript_biotype == "nonsense_mediated_decay")

NMDBio_cdf = ggplot(alltransbioNMD, aes(log2FoldChange, colour=transcript_biotype))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "SF3B1 KO vs WT", subtitle = "NMD biotype")
NMDBio_cdf

#analysis with robert's PTC+ and PTC- list
ENST_PTC.EPI.TFG <- read.delim("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_SF3B1_K562_KO/ENST_PTC-EPI-TFG.txt")
alltrans$ENST.ID <- sub("\\..*", "", alltrans$transcript_id)
test<-alltrans
test$PTC=NULL
names(test)[8]<-"ENST.ID"

alltransPTC<-inner_join(alltrans, ENST_PTC.EPI.TFG, by = "ENST.ID")

PTC_cdf = ggplot(alltransPTC, aes(log2FoldChange, colour=PTC.Status))+
  stat_ecdf(size=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "SF3B1 KO vs WT", subtitle = "PTC list")
PTC_cdf

resPTC <- wilcox.test(log2FoldChange ~ PTC.Status, data = alltransPTC,
                   exact = FALSE, alternative = "less")

resPTC
boxplot(log2FoldChange ~ PTC.Status, data=alltransPTC)

#Analysis using the NMD+- transcript list
NMD_transcripts = read.delim("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_SF3B1_K562_KO/NMD_transcript_list.txt")
alltransNMD = inner_join(alltrans, NMD_transcripts, by="transcript_id")

NMD_cdf = ggplot(alltransNMD, aes(log2FoldChange, colour=NMD))+
  stat_ecdf(size=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "SF3B1 KO vs WT", subtitle = "NMD +/- list")
NMD_cdf
NMD_res <- wilcox.test(log2FoldChange ~ NMD, data = alltransNMD,
                   exact = FALSE, alternative = "less")
NMD_res #p-vlaue = 


#Analysis using the stringent NMD list
StringentNMD = read.csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_SF3B1_K562_KO/Stringent_NMD_list_CE")
alltransSTNMD = inner_join(alltrans, StringentNMD, by = c("ENST.ID" = "transID"))

sNMD_cdf = ggplot(alltransSTNMD, aes(log2FoldChange, colour=NMD))+
  stat_ecdf(size=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "SF3B1 KO vs WT", subtitle = "Stringent NMD list")
sNMD_cdf

sNMD_res <- wilcox.test(log2FoldChange ~ NMD, data = alltransSTNMD,
                   exact = FALSE, alternative = "less")
sNMD_res #p-vlaue = 

#Make a CDF using the stringent PTC list
sPTC = read.csv("Stringent_PTC_CE") 
alltransSPTC = inner_join(alltrans, sPTC, by = c("ENST.ID" = "transID")) 
sPTC_res <- wilcox.test(log2FoldChange ~ PTC, data = alltransSPTC,
                        exact = FALSE, alternative = "less")
head(alltransSPTC)
count(alltransSPTC, PTC) #1142 PTC false 374 PTC true
sPTC_res #p-vlaue = 2.769e-6

plot_colors = c("FALSE" = "#673272", "TRUE" = "#EA7428",
                "PE" = "#76CEB5", "MANE" = "#BCB29F")

sPTC_cdf = ggplot(alltransSPTC, aes(log2FoldChange, colour=PTC))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = plot_colors) +
  labs(y = "Cumulative Frequency")
sPTC_cdf = sPTC_cdf + annotate("text", x=-2.5,y=1, label = "SF3B1 KO vs WT", size = 5, hjust = 0) +
  annotate("text", x=-2.5,y=.95, label = "Stringent PTC containing transcripts", hjust =0) +
  annotate("text", x=-2.5,y=.9, label = "p-value = 2.769e-6", hjust = 0) +
  annotate("text", x = -2.5, y = 0.85, label = "PTC+ = 374", hjust = 0) +
  annotate("text", x = -2.5, y = 0.8, label = "PTC- = 1142", hjust = 0) +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey") +
  theme(legend.position = c(0.85,0.1))
sPTC_cdf
ggsave("SF3B1_K562_sPTC_CDF.pdf", 
       plot = sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)
ggsave("C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/SF3B1_K562_sPTC_CDF.pdf", 
       plot = sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

##Make dot plot comparing the NMD plus and minus FC##
SF3B1_x = alltransSTNMD %>% filter(NMD == "TRUE")  %>% #Makes a data table of the top NMD transcript for each gene
  select(Gene,log2FoldChange,baseMean) %>% 
  rename(plusFC = log2FoldChange) %>% 
  group_by(Gene) %>% 
  slice_max(baseMean, n = 1)

SF3B1_y  = alltransSTNMD %>% filter(NMD == "FALSE")  %>% #Makes a data table fo the top NMD- transcript for each gene
  select(Gene,log2FoldChange,baseMean) %>% 
  rename(minusFC = log2FoldChange) %>% 
  group_by(Gene) %>% 
  slice_max(baseMean, n = 1)

SF3B1_scatter = full_join(SF3B1_x, SF3B1_y, by = "Gene", suffix = c("plus","minus")) #Combine the two top tables
SF3B1_scatter = SF3B1_scatter %>% filter(!is.na(minusFC)) %>% 
  filter(!is.na(plusFC)) #Removes any NA rows


dot = ggplot(data = SF3B1_scatter)+
  geom_point(aes(x = plusFC, y = minusFC), color = "blue")+
  geom_abline(aes(intercept = 0, slope = 1))+ #Makes a graph where each dot represents the NMD+ FC (x-axis) and NMD- FC (y-axis) for a specific gene
  gghighlight(abs(plusFC) > 0.58, unhighlighted_params = list(alpha("black", 0.4)))+ #Highlights all of the genes that have more than a 1.5 log2FC in either direction in the NMD+ samples
  labs(x = "Log2FC NMD+", y = "Log2FC NMD-", title = "SF3B1 gene-level FC comparison")+
  theme_minimal()+
  coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) #Use if you need to restrict the axis length
dot


#Count the number of transcripts that are significantly changed
counts = alltransSTNMD %>% filter(padj < 0.05) %>% 
  filter(abs(log2FoldChange) >0.58) #All transcripts with a fold change higher than 1.5
count(counts, NMD)
#Sig_NMD+ = 
#Sig_NMD- = 

#Determine the effect on NMD when Alternately spliced genes are removed
ASgenes_SF3B1 <- read_delim("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/rMATS/ENCODE_SF3B1/ASgenes_SF3B1.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE)
ASgenes_SF3B1 = ASgenes_SF3B1 %>% select(GeneID,geneSymbol)
sPTC_gene <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id"),
                   filters = "ensembl_transcript_id",
                   values = sPTC$transID,
                   mart = ensembl)
sPTC_gene = sPTC_gene %>% inner_join(sPTC, by = c("ensembl_transcript_id" = "transID"))
sPTC_AS = sPTC_gene %>% anti_join(ASgenes_SF3B1, by = c("ensembl_gene_id" = "GeneID"))
count(sPTC_AS, PTC) #PTC false=  170 PTC true= 39
alltransAS = inner_join(alltrans, sPTC_AS, by = c("ENST.ID" = "ensembl_transcript_id"))
count(alltransAS, PTC) # PTC false=75  PTC true=14


asPTC_res <- wilcox.test(log2FoldChange ~ PTC, data = alltransAS,
                         exact = FALSE, alternative = "less")
head(alltransAS)
asPTC_res #p-vlaue = 0.5914

asPTC_cdf = ggplot(alltransAS, aes(log2FoldChange, colour=PTC))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "SF3B1 KD vs WT", subtitle = "Non-AS transcripts")
asPTC_cdf = asPTC_cdf + annotate("text", x=-2.5,y=1, label = "p-value = 0.5914") +
  annotate("text", x = -2.5, y = 0.95, label = "PTC+ = 14") +
  annotate("text", x = -2.5, y = 0.9, label = "PTC- = 75") +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey")
asPTC_cdf
ggsave("SF3B1_K562_asPTC_CDF.pdf", 
       plot = asPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

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

SF3B1_NMD = alltrans_MANE %>% inner_join(NMD_factors, by = c("external_gene_name" = "NMDfactor"))
write.csv(SF3B1_NMD, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/NMD_analysis/SF3B1_NMDfactor.csv",row.names = FALSE)

####Look at just the shared stringent PTC from ISAR####
shared_sPTC <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/shared_stringent_PTC.csv")
shared_SPTC_trans = inner_join(alltrans, shared_sPTC, by = c("ENST.ID" = "transID")) 
shared_sPTC_res <- wilcox.test(log2FoldChange ~ PTC, data = shared_SPTC_trans,
                               exact = FALSE, alternative = "less")
head(shared_SPTC_trans)
count(shared_SPTC_trans, PTC)#253 PTC false 99 PTC true
shared_sPTC_res #p-vlaue = 4.9e-3

shared_sPTC_cdf = ggplot(shared_SPTC_trans, aes(log2FoldChange, colour=PTC))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "SF3B1 KO vs WT", subtitle = "Shared stringent PTC transcripts")
shared_sPTC_cdf = shared_sPTC_cdf + annotate("text", x=-2.5,y=1, label = "p-value = 4.9e-3") +
  annotate("text", x = -2.5, y = 0.95, label = "PTC+ = 99") +
  annotate("text", x = -2.5, y = 0.9, label = "PTC- = 253") +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey")
shared_sPTC_cdf
ggsave("SF3B1_K562_shared_PTC_CDF.pdf", 
       plot = shared_sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

#### Look at PE from (Thomas et al 2020 Nat. Gen.) ####
PE_isoform_list <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Poison_Exons/PE_isoform_list.csv")
PE_FC = alltrans %>% inner_join(PE_isoform_list, by = c("ENST.ID" = "ensembl_transcript_id"))
PE_res <- wilcox.test(log2FoldChange ~ isoform, data = PE_FC,
                      exact = FALSE, alternative = "less")
head(PE_FC)
count(PE_FC, isoform)# 83 PE 66  MANE
PE_res #p-vlaue = 0.02576

PE_cdf = ggplot(PE_FC, aes(log2FoldChange, colour=isoform))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = plot_colors) +
  labs(y = "Cumulative Frequency")
PE_cdf = PE_cdf + annotate("text", x=-2.5,y=1, label = "SF3B1 KO vs WT", size = 5, hjust = 0) +
  annotate("text", x=-2.5,y=.95, label = "Transcripts with a poison exon", hjust =0) +
  annotate("text", x=-2.5,y=.9, label = "p-value = 0.02576", hjust = 0) +
  annotate("text", x = -2.5, y = 0.85, label = "PE = 83", hjust = 0) +
  annotate("text", x = -2.5, y = 0.8, label = "MANE = 66", hjust = 0) +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey") +
  theme(legend.position = c(0.85,0.1))
PE_cdf
ggsave("SF3B1_K562_PE_CDF.pdf", 
       plot = PE_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)
ggsave("C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/SF3B1_K562_PE_CDF.pdf", 
       plot = PE_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

## Repeat with the new PE list
PE_archive_list <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Poison_Exons/archive_PE_list.csv")
PEa_FC = alltrans %>% inner_join(PE_archive_list, by = c("ENST.ID" = "ensembl_transcript_id"))
PEa_res <- wilcox.test(log2FoldChange ~ isoform, data = PEa_FC,
                       exact = FALSE, alternative = "less")
head(PEa_FC)
count(PEa_FC, isoform)# 432 PE  228 MANE
PEa_res #p-vlaue = 0.1676

PEa_cdf = ggplot(PEa_FC, aes(log2FoldChange, colour=isoform))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = plot_colors) +
  labs(y = "Cumulative Frequency")
PEa_cdf = PEa_cdf + annotate("text", x=-2.5,y=1, label = "SF3B1 KO vs WT", size = 5, hjust = 0) +
  annotate("text", x=-2.5,y=.95, label = "Transcripts with a poison exon", hjust =0) +
  annotate("text", x=-2.5,y=.9, label = "p-value = 0.1676", hjust = 0) +
  annotate("text", x = -2.5, y = 0.85, label = "PE = 432", hjust = 0) +
  annotate("text", x = -2.5, y = 0.8, label = "MANE = 228", hjust = 0) +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey") +
  theme(legend.position = c(0.85,0.1))
PEa_cdf
ggsave("SF3B1_K562_PE_archive_CDF.pdf", 
       plot = PEa_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)
ggsave("C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/SF3B1_K562_PEarchive_CDF.pdf", 
       plot = PEa_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

##Redo the sPTC analysis using the MANE list
sPTC_MANE = read.csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_SF3B1_K562_KO/Stringent_PTC_MANE_CE.csv")
alltransSPTC_MANE = inner_join(alltrans, sPTC_MANE, by = c("ENST.ID" = "transID")) 
sPTC_MANE_res <- wilcox.test(log2FoldChange ~ PTC, data = alltransSPTC_MANE,
                             exact = FALSE, alternative = "less")
count(alltransSPTC_MANE, PTC) # 263 PTC false 374 PTC true
head(alltransSPTC_MANE)
sPTC_MANE_res #p-vlaue = 2.233e-10


MANE_sPTC_cdf = ggplot(alltransSPTC_MANE, aes(log2FoldChange, colour=PTC))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = plot_colors) +
  labs(y = "Cumulative Frequency")
MANE_sPTC_cdf = MANE_sPTC_cdf + annotate("text", x=-2.5,y=1, label = "SF3B1 KD vs WT", size = 5, hjust = 0) +
  annotate("text", x=-2.5,y=.95, label = "PTC containing and MANE transcripts", hjust =0) +
  annotate("text", x=-2.5,y=.9, label = "p-value = 2.233e-10", hjust = 0) +
  annotate("text", x = -2.5, y = 0.85, label = "PTC+ = 374", hjust = 0) +
  annotate("text", x = -2.5, y = 0.8, label = "PTC- = 263", hjust = 0) +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey") +
  theme(legend.position = c(0.85,0.1))
MANE_sPTC_cdf
ggsave("SF3B1_K562_sPTC_MANE_CDF.pdf", 
       plot = MANE_sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)
ggsave("C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/SF3B1_K562_sPTC_MANE_CDF.pdf", 
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
PTC_tpm = PTC_tpm %>% select(1,3:4,6:7)
write_csv(PTC_tpm, "SF3B1_PTC_MANE_TPM.csv")
write_csv(PTC_tpm, "C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/SF3B1_PTC_MANE_TPM.csv")
write_csv(alltrans,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/SF3B1_DBfilt_alltrans.csv")
write_csv(alltransSPTC_MANE,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/SF3B1_PTC_alltrans.csv")

####Make a list of the genes filtered to WT>5TPM ####
WT5counts <-as.data.frame(txi.kallisto.tsv$counts)

#Filtering out transcripts with mean tpm<5
WT5counts$mean=rowMeans(WT5counts[,c("SRR14838495","SRR14838496")], na.rm=TRUE) #Filter WT samples
WT5counts=filter(WT5counts, mean>=5)

listMarts()
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
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
WT5res_counts <- results(WT5ddscounts, contrast = c("genotype","KO","WT"))
WT5res_counts
mcols(WT5res_counts, use.names = TRUE)

#significant results
WT5res.sig <- WT5res_counts[ which(WT5res_counts$padj < 0.05), ]
hist(WT5res.sig$pvalue, col = "green1")
head(WT5res.sig)
plotMA(WT5res_counts, padj = TRUE, ylim=c(-6,6), main ="MA plot: WT vs KO")

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
write_csv(WT5alltrans,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/SF3B1_WT5filt_alltrans.csv")
