#Written by: Caleb Embree
#Modified by: Caleb Embree
#Modified on: 7/10/2024
#Analyzing: AQR Kallisto results from transcriptome including novel transcripts

#All packages already installed#
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximport")
#BiocManager::install("rhdf5")
#BiocManager::install("biomaRt")
#BiocManager::install("DESeq2")
#BiocManager::install("pheatmap")

#Things to use find and replace on
# C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_AQR_K562/novel_kallisto --> working directory of the project
# SRR14846211 --> Sample ID of first WT sample
# SRR14846212 --> Sample ID of second WT sample
# SRR14844828 --> Sample ID of first KD sample
# SRR14844829 --> Sample ID of second KD sample
# AQR --> Ensemble name of your gene of interest

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
samples <- read_excel("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_AQR_K562/AQR_K562_samples.xlsx")
View(samples)


#### Load in the files ####
files1 <- file.path("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_AQR_K562/novel_kallisto","novel_koutput", samples$Run, "abundance.h5")
names(files1) <- paste0(c("SRR14844828","SRR14844829","SRR14846211","SRR14846212"))
files1

#importing output files from each library to a compiled dataframe
#txi.kallisto <- tximport(files1, type = "kallisto", txOut = TRUE)
txi.kallisto.tsv <- tximport(files1, type = "kallisto", txOut=TRUE, countsFromAbundance = "lengthScaledTPM")
head(txi.kallisto.tsv$counts)
counts<-as.data.frame(txi.kallisto.tsv$counts)

#Filtering out transcripts with mean tpm<1
counts = counts %>% mutate(wtMean = rowMeans(dplyr::select(counts,c(SRR14846211,SRR14846212)), na.rm = TRUE) ,
                           kdMean = rowMeans(dplyr::select(counts,c(SRR14844828,SRR14844829)), na.rm = TRUE)) %>% 
  filter(wtMean >= 1 | kdMean >= 1) #Filter to the transcripts with TPM >= 1 in either WT or KD
counts$ensembl_transcript_id_version<-rownames(counts)
counts = counts %>% mutate(ENSTID = if_else(str_detect(ensembl_transcript_id_version,"ENST"), #If the isoform is annotated,
                                            str_remove(ensembl_transcript_id_version,"\\..*"), #Remove the version id to get just the ENSTID
                                            ensembl_transcript_id_version)) #Otherwise use the novel isoform ID


#######Finding isoform specific counts associated with a gene in the dataset.####
listEnsemblArchives()
ensembl <- useMart("ensembl",host = "https://feb2023.archive.ensembl.org") #use version 109 of ensembl
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
listgenes <- getBM(attributes = c("external_gene_name","ensembl_transcript_id","ensembl_gene_id"),
                   filters = "ensembl_transcript_id",
                   values = counts$ENSTID,
                   mart = ensembl)

anno_genes = counts %>% filter(str_detect(ENSTID,"ENST")) %>%
  left_join(listgenes, by = c("ENSTID" = "ensembl_transcript_id"))
AQR_novel_iso_features <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_AQR_ISAR/AQR_novel_iso_features.csv")
novel_genes = counts %>% filter(str_detect(ENSTID,"MST")) %>% 
  left_join(AQR_novel_iso_features, by = c("ENSTID" = "isoform_id")) %>% 
  filter(!is.na(PTC))
full_genelist = anno_genes %>% full_join(novel_genes, by = c("ENSTID",
                                                             "ensembl_gene_id" = "gene_id",
                                                             "external_gene_name" = "gene_name",
                                                             "SRR14844828","SRR14844829","SRR14846211","SRR14846212",
                                                             "wtMean","kdMean","ensembl_transcript_id_version"))

####Differential expression analysis#####

#sampleinfo as dataframe
colData <- data.frame(genotype = factor(c("KD","KD","WT","WT"))) #make sure these are in the right order that match your sample file


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
ggsave("AQR_novel_PCA.pdf",
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
#log2fc range : -23.6439 to +24.1336

#Exporting files
write.csv(as.data.frame(res_counts), file = "KDvsWT_full_results.csv")
write.csv(as.data.frame(res.sig), file = "2KDvsWT_sig_results.csv")
write.csv(counts(ddscounts, normalized = T), file = "KDvsWT_normalized_counts.csv")

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
     main="AQR knockout vs control")



#### Creating alltrans datatables ####
alltrans <- as.data.frame(res_counts)
alltrans = alltrans %>% filter(!is.na(padj)) %>% 
  rownames_to_column(var = "transcriptID_version")
alltrans = alltrans %>% mutate(ENST.ID = if_else(str_detect(transcriptID_version,"ENST"),
                                                 str_remove(transcriptID_version,"\\..*"),
                                                 transcriptID_version)) 

alltransbiotype <- getBM(attributes = c("external_gene_name","ensembl_gene_id",
                                        "ensembl_transcript_id","transcript_biotype",
                                        "transcript_mane_select"),
                         filters = "ensembl_transcript_id",
                         values = alltrans$ENST.ID,
                         mart = ensembl)


alltrans_annotated = alltrans %>% filter(str_detect(ENST.ID,"ENST")) %>% 
  left_join(alltransbiotype, by=c("ENST.ID" = "ensembl_transcript_id")) #Add in the biomart info for the annotated isoforms
alltrans_novel = alltrans %>% filter(str_detect(ENST.ID,"MST")) %>% 
  left_join(AQR_novel_iso_features,by = c("ENST.ID" = "isoform_id"))#add in the annotation from ISAR for the novel isoforms
alltransbio = alltrans_annotated %>% full_join(alltrans_novel,by = c("ENST.ID","transcriptID_version","baseMean","log2FoldChange",
                                                                        "lfcSE","stat","pvalue","padj",
                                                                        "ensembl_gene_id" = "gene_id",
                                                                        "external_gene_name" = "gene_name")) #Combine the two tables together
write_csv(alltransbio,"C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/AQR_NK_alltrans.csv")

#### CDF plot comparing NMD biotype ####
alltransbioNMD = alltransbio %>% filter(transcript_biotype == "protein_coding" | transcript_biotype == "nonsense_mediated_decay")

NMDBio_cdf = ggplot(alltransbioNMD, aes(log2FoldChange, colour=transcript_biotype))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "AQR KD vs WT", subtitle = "NMD biotype")
NMDBio_cdf

#### analysis with robert's PTC+ and PTC- list ####
ENST_PTC.EPI.TFG <- read.delim("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/ENCODE_AQR_K562/ENST_PTC-EPI-TFG.txt")


alltransPTC<-inner_join(alltrans, ENST_PTC.EPI.TFG, by = "ENST.ID")

PTC_cdf = ggplot(alltransPTC, aes(log2FoldChange, colour=PTC.Status))+
  stat_ecdf(linewidth=2)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "grey"))+
  coord_cartesian(xlim = c(-3,3)) +
  scale_color_manual(values = met.brewer("Java",2)) +
  labs(y = "Cumulative Frequency", title = "AQR KD vs WT", subtitle = "PTC list")
PTC_cdf

resPTC <- wilcox.test(log2FoldChange ~ PTC.Status, data = alltransPTC,
                      exact = FALSE, alternative = "less")

resPTC
boxplot(log2FoldChange ~ PTC.Status, data=alltransPTC)


plot_colors = c("FALSE" = "#673272", "TRUE" = "#EA7428",
                "PE" = "#76CEB5", "MANE" = "#BCB29F")



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
MANE_sPTC_cdf = MANE_sPTC_cdf + annotate("text", x=-2.5,y=1, label = "AQR KD vs WT", size = 5, hjust = 0) +
  annotate("text", x=-2.5,y=.95, label = "PTC containing and MANE transcripts", hjust =0) +
  annotate("text", x=-2.5,y=.9, label = "p-value = ", hjust = 0) +
  annotate("text", x = -2.5, y = 0.85, label = "PTC+ = ", hjust = 0) +
  annotate("text", x = -2.5, y = 0.8, label = "MANE = ", hjust = 0) +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_hline(yintercept = 0.5, color = "grey") +
  theme(legend.position = c(0.85,0.1))
MANE_sPTC_cdf
ggsave("AQR_K562_sPTC_MANE_CDF.pdf", 
       plot = MANE_sPTC_cdf,
       scale = 1,
       width = 8,
       height = 6,
       units = "in",
       device = "pdf",
       dpi = 300)
ggsave("C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/AQR_K562_sPTC_MANE_CDF.pdf", 
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
                                         "transcript_biotype","ensembl_transcript_id","transcript_mane_select"),
                          filters = "ensembl_transcript_id",
                          values = alltrans_full$ENST.ID,
                          mart = ensembl)
alltrans_annotated = alltrans_full %>% left_join(alltrans_ensembl, by = c("ENST.ID" = "ensembl_transcript_id"))
alltrans_MANE = alltrans_annotated %>% filter(!is.na(transcript_mane_select)) %>% filter(transcript_mane_select != "")
NMD_factors <- read_excel("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/NMD_factors.xlsx")

AQR_NMD = alltrans_MANE %>% inner_join(NMD_factors, by = c("external_gene_name" = "NMDfactor"))
write.csv(AQR_NMD, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/NMD_analysis/AQR_nk_NMDfactor.csv",row.names = FALSE)


####Look at effect on the MANE transcript of spliceosome factors under study####
Splicing_factors <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/Splicing_factors.csv")

AQR_SC_impact = alltrans_MANE %>% inner_join(Splicing_factors, by = c("external_gene_name" = "Component"))
write.csv(AQR_SC_impact, "AQR_splicing_impact.csv",row.names = FALSE)



#### Pull the TPM of genes for the NMD TPM graph ####
#Check WT and KD are refering to the right columns in full_genelist
write_csv(full_genelist, "C:/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/AQR_nk_TPM.csv")
