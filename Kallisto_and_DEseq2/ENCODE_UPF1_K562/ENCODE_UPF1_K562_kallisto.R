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
# C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_UPF1_K562 --> working directory of the project
# SRR4421434 --> Sample ID of first WT sample
# SRR4421435 --> Sample ID of second WT sample
# SRR4421542 --> Sample ID of first KD sample
# SRR4421543 --> Sample ID of second KD sample
# UPF1 --> Ensemble name of your gene of interest

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

samples <- read_excel("UPF1_K562_samples.xlsx")
View(samples)


#Load in the files
files1 <- file.path("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_UPF1_K562","koutput", samples$Run, "abundance.h5")
names(files1) <- paste0(c("SRR4421434","SRR4421435","SRR4421542","SRR4421543"))
files1


#importing output files from each library to a compiled dataframe
#txi.kallisto <- tximport(files1, type = "kallisto", txOut = TRUE)
txi.kallisto.tsv <- tximport(files1, type = "kallisto", txOut=TRUE, countsFromAbundance = "lengthScaledTPM")
head(txi.kallisto.tsv$counts)
counts<-as.data.frame(txi.kallisto.tsv$counts)

#Filtering out transcripts with mean tpm<1
counts$mean=rowMeans(counts[,c("SRR4421434","SRR4421435")], na.rm=TRUE) #Filter WT samples
counts=filter(counts, mean>=1)
counts$mean=rowMeans(counts[,c("SRR4421542","SRR4421543")], na.rm=TRUE) #Fiter KD samples
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
PCA = plotPCA(vsd, intgroup = "genotype")+
  geom_text_repel(aes(label=name))
PCA
ggsave("UPF1_PCA.pdf",
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
     main="UPF1 knockout vs control")



#### Creating alltrans datatables ####
alltrans <- as.data.frame(res_counts)
alltrans = alltrans %>% filter(!is.na(padj)) %>% 
  rownames_to_column(var = "transcriptID_version")
alltrans = alltrans %>% mutate(ENST.ID = str_remove(transcriptID_version,"\\..*")) 

alltransbiotype <- getBM(attributes = c("external_gene_name","ensembl_gene_id",
                                        "ensembl_transcript_id","transcript_biotype",
                                        "transcript_mane_select"),
                         filters = "ensembl_transcript_id",
                         values = alltrans$ENST.ID,
                         mart = ensembl)


alltransbioNMD <- left_join(alltrans,alltransbiotype, by=c("ENST.ID" = "ensembl_transcript_id")) #Check this for gene FC

#### CDF plot comparing NMD biotype ####
write(alltransbiotype$ensembl_transcript_id_version, "tids.txt")
alltransbioNMD = alltransbioNMD %>% filter(transcript_biotype == "protein_coding" | transcript_biotype == "nonsense_mediated_decay")

NMDBio_cdf = ggplot(alltransbioNMD, aes(log2FoldChange, colour=transcript_biotype))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "UPF1 KD vs WT", subtitle = "NMD biotype")
NMDBio_cdf

#### analysis with robert's PTC+ and PTC- list ####
ENST_PTC.EPI.TFG <- read.delim("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_UPF1/ENST_PTC-EPI-TFG.txt")


alltransPTC<-inner_join(alltrans, ENST_PTC.EPI.TFG, by = "ENST.ID")

PTC_cdf = ggplot(alltransPTC, aes(log2FoldChange, colour=PTC.Status))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "UPF1 KD vs WT", subtitle = "PTC list")
PTC_cdf

resPTC <- wilcox.test(log2FoldChange ~ PTC.Status, data = alltransPTC,
                      exact = FALSE, alternative = "less")

resPTC
boxplot(log2FoldChange ~ PTC.Status, data=alltransPTC)


plot_colors = c("FALSE" = "#673272", "TRUE" = "#EA7428",
                "PE" = "#76CEB5", "MANE" = "#BCB29F")



#### Determine the effect on NMD when Alternately spliced genes are removed (Delete this section if rMATS hasn't been done) ####
ASgenes_UPF1 <- read_delim("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/rMATS/ENCODE_UPF1/ASgenes_UPF1.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE)
ASgenes_UPF1 = ASgenes_UPF1 %>% select(GeneID,geneSymbol)
sPTC_gene <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id"),
                   filters = "ensembl_transcript_id",
                   values = sPTC$ENST.ID,
                   mart = ensembl)
sPTC_gene = sPTC_gene %>% inner_join(sPTC, by = c("ensembl_transcript_id" = "ENST.ID"))
sPTC_AS = sPTC_gene %>% anti_join(ASgenes_UPF1, by = c("ensembl_gene_id" = "GeneID"))
count(sPTC_AS, PTC) #PTC false=  PTC true= 
alltransAS = inner_join(alltrans, sPTC_AS, by = c("ENST.ID" = "ensembl_transcript_id"))
count(alltransAS, PTC) # PTC false  PTC true


asPTC_res <- wilcox.test(log2FoldChange ~ PTC, data = alltransAS,
                         exact = FALSE, alternative = "less")
head(alltransAS)
asPTC_res #p-vlaue = 

asPTC_cdf = ggplot(alltransAS, aes(log2FoldChange, colour=PTC))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "UPF1 KD vs WT", subtitle = "Non-AS transcripts")
asPTC_cdf = asPTC_cdf + annotate("text", x=-2.5,y=1, label = "p-value = ") +
  annotate("text", x = -2.5, y = 0.95, label = "PTC+ = ") +
  annotate("text", x = -2.5, y = 0.9, label = "PTC- = ") +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey")
asPTC_cdf
ggsave("UPF1_K562_asPTC_CDF.pdf", 
       plot = asPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)

##sPTC analysis using the MANE list
sPTC_MANE = read.csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/Stringent_PTC_MANE_CE.csv")
alltransSPTC_MANE = inner_join(alltrans, sPTC_MANE, by = c("ENST.ID"= "transID")) 
sPTC_MANE_res <- wilcox.test(log2FoldChange ~ PTC, data = alltransSPTC_MANE,
                             exact = FALSE, alternative = "less")
count(alltransSPTC_MANE, PTC) # PTC false  PTC true
sPTC_MANE_res #p-vlaue = 


MANE_sPTC_cdf = ggplot(alltransSPTC_MANE, aes(log2FoldChange, colour=PTC))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = plot_colors) +
  labs(y = "Cumulative Frequency")
MANE_sPTC_cdf = MANE_sPTC_cdf + annotate("text", x=-2.5,y=1, label = "UPF1 KD vs WT", size = 5, hjust = 0) +
  annotate("text", x=-2.5,y=.95, label = "PTC containing and MANE transcripts", hjust =0) +
  annotate("text", x=-2.5,y=.9, label = "p-value = ", hjust = 0) +
  annotate("text", x = -2.5, y = 0.85, label = "PTC+ = ", hjust = 0) +
  annotate("text", x = -2.5, y = 0.8, label = "MANE = ", hjust = 0) +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey") +
  theme(legend.position = c(0.85,0.1))
MANE_sPTC_cdf
ggsave("UPF1_K562_sPTC_MANE_CDF.pdf", 
       plot = MANE_sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)
ggsave("C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/UPF1_K562_sPTC_MANE_CDF.pdf", 
       plot = MANE_sPTC_cdf,
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
                                         "transcript_biotype","ensembl_transcript_id","transcript_mane_select","ensembl_gene_id"),
                          filters = "ensembl_transcript_id",
                          values = alltrans_full$ENST.ID,
                          mart = ensembl)
alltrans_annotated = alltrans_full %>% left_join(alltrans_ensembl, by = c("ENST.ID" = "ensembl_transcript_id"))
alltrans_MANE = alltrans_annotated %>% filter(!is.na(transcript_mane_select)) %>% filter(transcript_mane_select != "")
NMD_factors <- read_excel("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/NMD_factors.xlsx")

UPF1_NMD = alltrans_MANE %>% inner_join(NMD_factors, by = c("external_gene_name" = "NMDfactor"))
write.csv(UPF1_NMD, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/NMD_analysis/UPF1_NMDfactor.csv",row.names = FALSE)


####Look at effect on the MANE transcript of spliceosome factors under study####
Splicing_factors <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/Splicing_factors.csv")

UPF1_SC_impact = alltrans_MANE %>% inner_join(Splicing_factors, by = c("external_gene_name" = "Component"))
write.csv(UPF1_SC_impact, "UPF1_splicing_impact.csv",row.names = FALSE)



#### Pull the TPM of genes for the NMD TPM graph ####
#Check WT and KD are refering to the right columns in listgenes
full_tpm = listgenes %>% rowwise() %>%  mutate(WTmean = mean(c_across(1:2)), KDmean = mean(c_across(3:4))) %>% ungroup() #Check that WT and KD are referring to the right columns
full_tpm = full_tpm %>% mutate(ENST.ID = str_remove(ensembl_transcript_id_version, "\\..*"))
full_tpm = full_tpm %>% select(9:11)
PTC_tpm = full_tpm %>% right_join(sPTC_MANE, by = c("ENST.ID" = "transID"))
PTC_tpm = PTC_tpm %>% select(1:3,5,6)
write_csv(PTC_tpm, "UPF1_PTC_MANE_TPM.csv")
write_csv(PTC_tpm, "C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/UPF1_PTC_MANE_TPM.csv")

#### Look at the effect when only WT filtering####
WTcounts <-as.data.frame(txi.kallisto.tsv$counts)

#Filtering out transcripts with mean tpm<1
WTcounts$mean=rowMeans(WTcounts[,c("SRR4421434","SRR4421435")], na.rm=TRUE) #Filter WT samples
WTcounts=filter(WTcounts, mean>=1)

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
write_csv(WTalltrans,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/UPF1_WTfilt_alltrans.csv")
write_csv(alltrans,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/UPF1_DBfilt_alltrans.csv")
write_csv(alltransSPTC_MANE,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/UPF1_PTC_alltrans.csv")

####Make a list of the genes filtered to WT>5TPM ####
WT5counts <-as.data.frame(txi.kallisto.tsv$counts)

#Filtering out transcripts with mean tpm<5
WT5counts$mean=rowMeans(WT5counts[,c("SRR4421434","SRR4421435")], na.rm=TRUE) #Filter WT samples
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
write_csv(WT5alltrans,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/UPF1_WT5filt_alltrans.csv")


##Look at the significant stress granule transcripts##
Stress_Granule_Transcriptome <- read_excel("Stress_Granule_Transcriptome.xlsx")
Stress_Granule_Transcriptome = Stress_Granule_Transcriptome %>% rename(ensembl_gene_id = test_id,
                                                                       SG_sig = significant,
                                                                       SG_local = Localization) %>% 
  select(1,6,10)
SG_alltrans_sig = alltrans_annotated %>% inner_join(Stress_Granule_Transcriptome, by = c("ensembl_gene_id")) %>% 
  filter(str_detect(transcript_mane_select,"NM") & SG_sig == "yes")
SG_summary_sig = SG_alltrans_sig %>% group_by(SG_local) %>% 
  summarise(n = n(),
            med = median(log2FoldChange))
SG_depleted_res_sig= wilcox.test(log2FoldChange ~ SG_local, data = SG_alltrans_sig %>% filter(SG_local != "Neither"),
                                 exact = FALSE, alternative = "less") #less because we expect the enriched to be higher
SG_neither_res_sig=  wilcox.test(log2FoldChange ~ SG_local, data = SG_alltrans_sig %>% filter(SG_local != "depleted"),
                                 exact = FALSE, alternative = "greater") #less because we expect the enriched to be higher

#SG_neither_res_sig$p.value

SG_Boxplot_sig = ggplot(data = SG_alltrans_sig)
SG_Boxplot_sig = SG_Boxplot_sig + geom_boxplot(aes(x = SG_local,
                                                   y = log2FoldChange,
                                                   fill = SG_local),
                                               position = position_dodge2(width = 0.9),
                                               width = 0.8,
                                               outlier.shape = 21,
                                               outlier.alpha = 0.5,
                                               outlier.colour = NA,
                                               linewidth = 1) +
  scale_color_met_d("Archambault") +
  scale_fill_met_d("Archambault") +
  geom_label(data = SG_summary_sig,
             aes(x = SG_local,
                 y = med,
                 color = SG_local,
                 label = round(med, digits = 3)),
             show.legend = F,
             position = position_dodge2(width = 0.8),
             size = 5) +
  geom_text(data = SG_summary_sig,
            aes(x = SG_local,
                y = -2.5,
                color = SG_local,
                label = n),
            show.legend = F,
            position = position_dodge2(width = 0.9),
            size = 8) +
  geom_text(aes(x = "depleted",
                y = 4,
                label = paste0("P(depleted)","=",signif(SG_depleted_res_sig$p.value,digits = 3))),
            size = 6,
            color = "#88a0dc") +
  geom_text(aes(x = "Neither",
                y = 4,
                label = paste0("P(neither)","=",signif(SG_neither_res_sig$p.value,digits = 3))),
            size = 6,
            color = "#f9d14a",) +
  labs(x = "Comparison",
       y = "log2(Fold Change)",
       fill = "Stress Granule Localization",
       title = "Stress Granule Localization",
       caption = "Sig SG genes")+
  coord_cartesian(y = c(-4,4)) +
  theme_bw()
SG_Boxplot_sig
ggsave("UPF1_KD_Stress_Granule_sig_boxplot.pdf",
       plot = SG_Boxplot_sig,
       device = pdf,
       width = 15,
       height = 10,
       units = "in",
       dpi = 300)
