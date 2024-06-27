# Script to analyze isoforms changes using Isoform Switch Analyzer
# Date Written: 9/6/2023
# Written by: Caleb Embree
# Modified on: 6/3/2024
# Analyzing: ENCODE KD of PRPF6 
# Go to https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html for package manual


#Processing that happened before the script (on the supercomputer):
# 1. Transcripts were mapped via HISAT2
# 2. Transcripts were assembled with Stringtie
# 3. Transcripts were merged with Stringtie and FASTA file of the merged GTF was created
# 4. Transcripts were assembled to the merged file with Stringtie


#### Copy and Paste ####
#PRPF6 = Gene of Interest
#C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_PRPF6_ISAR = Path to the working directory
#SRR4422353 = First WT sample
#SRR4422354 = Second WT sample
#SRR4422262 = First KD sample
#SRR4422263 = Second KD sample


#Setup
#if (!requireNamespace("pfamAnalyzeR", quietly = TRUE)){
#  devtools::install_github("kvittingseerup/pfamAnalyzeR")
#}
#devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR")
library(IsoformSwitchAnalyzeR)
library(rhdf5)
library(dplyr)
setwd("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_PRPF6_ISAR")

#Import datasets
PRPF6quant = importIsoformExpression(parentDir = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_PRPF6_ISAR",
                                      addIsofomIdAsColumn = TRUE, 
                                      showProgress = TRUE,
                                      readLength = 100
) #reads the t_data.ctab file in each sub directory
names(PRPF6quant)
tail(PRPF6quant$abundance)
tail(PRPF6quant$counts)

#Generate the list of isoform switches
myDesign = data.frame(sampleID = c("SRR4422353","SRR4422354","SRR4422262","SRR4422263"),
                      condition = c("WT","WT","KD","KD"))
SwitchList = importRdata(isoformCountMatrix = PRPF6quant$counts,
                         isoformRepExpression = PRPF6quant$abundance,
                         designMatrix = myDesign,
                         isoformExonAnnoation = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_PRPF6_ISAR/st_merged.gtf",
                         isoformNtFasta = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_PRPF6_ISAR/PRPF6_transcripts.fa",
                         showProgress = TRUE)
SwitchList_filt = preFilter(SwitchList,
                            geneExpressionCutoff = 1, #Cut off genes with less than 1 TPM
                            isoformExpressionCutoff = 0, #removes unused isoforms
                            removeSingleIsoformGenes = TRUE, #removes genes with a single isoform
) #The filtering removed 65968 ( 456.18% of ) transcripts. There is now 76883 isoforms left

#Perform switch test (will take a while)
AnalyzedSwitch = isoformSwitchTestDEXSeq(switchAnalyzeRlist = SwitchList_filt, #Uses DEXseq to analyze the isoform switches
                                         alpha = 0.05, #The FDR p-value has to be less than this
                                         dIFcutoff = 0.1, #THere must be at least a 10% change for the isoforms to be considered
                                         reduceToSwitchingGenes = TRUE, #This only keeps genes with one isoform significantly differently used, makes it run faster
                                         showProgress = TRUE)

#Identify ORFs
"orfAnalysis" %in% names(AnalyzedSwitch)
SwitchAndORF = addORFfromGTF(AnalyzedSwitch,
                             pathToGTF = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Genome_Files/Homo_sapiens.GRCh38.109.chr_patch_hapl_scaff.gtf",
                             ignoreAfterPeriod = TRUE)
SwitchAndORF = analyzeNovelIsoformORF(SwitchAndORF,
                                      analysisAllIsoformsWithoutORF = TRUE,
                                      orfMethod = "longest.AnnotatedWhenPossible")
"orfAnalysis" %in% names(SwitchAndORF)


#Extract Nucleotide and amino acid sequence (not needed just for NMD)
#extractSequence(AnalyzedSwitch,
#                removeLongAAseq = TRUE,
#                alsoSplitFastaFile = TRUE)
SwitchAndORF = analyzeAlternativeSplicing(switchAnalyzeRlist = SwitchAndORF,
                                          alpha = 0.05,
                                          dIFcutoff = 0.1,
                                          showProgress = TRUE)
Switch_cons = analyzeSwitchConsequences(switchAnalyzeRlist = SwitchAndORF, #Don't actually need to do any CPC2 analysis just to look at NMD status
                                        consequencesToAnalyze = c("NMD_status"))
TopSwitch = extractTopSwitches(Switch_cons,filterForConsequences = TRUE, extractGenes = FALSE) #Pulls out the top 10 switching isoforms with consequences
TopSwitch
#switchPlot(Switch_cons, gene = 'AK2') #This won't work unless all the other external datasets are added in
extractConsequenceSummary(Switch_cons)
extractConsequenceEnrichment(Switch_cons)
extractConsequenceGenomeWide(Switch_cons)
extractSplicingEnrichment(Switch_cons)

all_switch = extractTopSwitches(
  Switch_cons,
  filterForConsequences = FALSE,
  n = NA,
  extractGenes = FALSE,
  sortByQvals = FALSE)

all_switch_cons = extractTopSwitches(
  Switch_cons,
  filterForConsequences = TRUE,
  n = NA,
  extractGenes = FALSE,
  sortByQvals = FALSE)

#Determining the TPM of different classes of transcripts. 
library(dplyr)
rawTPM = PRPF6quant$abundance
annotatedTPM = rawTPM %>% mutate(kdTPM = rowMeans(dplyr::select(rawTPM,2:3), na.rm = TRUE)) %>%
  mutate(wtTPM = rowMeans(dplyr::select(rawTPM, 4:5), na.rm = TRUE))
library(biomaRt)
listMarts()
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
alltrans_ensembl <- getBM(attributes = c("ensembl_transcript_id_version","external_gene_name",
                                         "transcript_biotype","ensembl_transcript_id","transcript_mane_select","ensembl_gene_id"),
                          filters = "ensembl_transcript_id",
                          values = annotatedTPM$isoform_id,
                          mart = ensembl)
annotatedTPM = annotatedTPM %>% left_join(alltrans_ensembl, by = c("isoform_id" = "ensembl_transcript_id"))
filtTPM = annotatedTPM %>% filter(kdTPM > 1 | wtTPM >1)
MANE_TPM = filtTPM %>% filter(!is.na(transcript_mane_select)) %>% filter(transcript_mane_select != "")
library(readr)
Stringent_PTC <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/Stringent_PTC_MANE_CE.csv")
PTC_TPM = filtTPM %>% inner_join(Stringent_PTC, by = c("isoform_id" = "transID")) %>% filter(PTC == TRUE)
Novel_Switch = all_switch %>% filter(dIF < 0) %>% filter(IF2 <= 0) %>% filter(!str_detect(isoform_id,"ENST")) #Creates a list of transcripts that are switching that has NMD  consequences and 0 incluision in WT
Novel_TPM = annotatedTPM %>% inner_join(Novel_Switch, by = "isoform_id")
Novel_TPM = Novel_TPM %>% dplyr::select(1:7,9:12,23) %>% mutate(type = "Novel") %>% rename(NMD_Causing = switchConsequencesGene)
PTC_TPM = PTC_TPM %>% dplyr::select(1:7,9:12) %>% mutate(type = "PTC")
PTC_MANE_TPM = MANE_TPM %>% inner_join(Stringent_PTC, by = c("isoform_id" = "transID"))
PTC_MANE_TPM = PTC_MANE_TPM %>% dplyr::select(1:7,9:12) %>% mutate(type = "MANE")
PRPF6_TPM = PTC_TPM %>% full_join(Novel_TPM) %>% full_join(PTC_MANE_TPM)
write_csv(PRPF6_TPM, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/PRPF6_TPM.csv")

####Comparing the TPMs of different classes of novel transcript####
all_novel = all_switch %>% filter(IF2 <= 0) #X transcripts
all_novel_tpm = all_novel %>% inner_join(rawTPM, by = "isoform_id")
all_novel_tpm = all_novel_tpm %>% mutate(kdTPM = rowMeans(dplyr::select(all_novel_tpm,14:15), na.rm = TRUE))
all_novel_tpm = all_novel_tpm %>% mutate(wtTPM = rowMeans(dplyr::select(all_novel_tpm, 15:16), na.rm = TRUE))#X transcripts
all_novel_tpm_filt = all_novel %>% inner_join(filtTPM, by = "isoform_id")#X transcripts
all_novel_tpm_filt = all_novel_tpm_filt %>% dplyr::select(3:12,18:19) %>% mutate(sample = "PRPF6")
write_csv(all_novel_tpm_filt, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/PRPF6_novel_TPM.csv")

#Look at transcript length
all_switch = extractTopSwitches( #Pulls all switching genes
  Switch_cons,
  filterForConsequences = FALSE,
  n = NA,
  extractGenes = FALSE,
  sortByQvals = FALSE)
count(all_switch,switchConsequencesGene) #FALSE =  TRUE = 
ORFanalysis = SwitchAndORF$orfAnalysis
count(ORFanalysis,PTC) #False =  True =  NA = 
ORFanalysis = ORFanalysis %>% mutate(Switch = isoform_id %in% all_switch$isoform_id)
count(ORFanalysis,Switch) #False =  True= ->Same number as switching transcripts
switching_ORF = ORFanalysis %>% filter(Switch == TRUE)
count(switching_ORF,PTC) #False =  TRUE =  NA =  
switching_ORF = switching_ORF %>% left_join(all_switch, by = "isoform_id")

length_plot = ggplot(data = switching_ORF)
length_plot = length_plot + geom_jitter(aes(x = switchConsequencesGene, y = log10(orfTransciptLength),
                                            shape = switchConsequencesGene),
                                        color = "grey",
                                        show.legend = F)+
  geom_boxplot(aes(x = switchConsequencesGene, y = log10(orfTransciptLength),
                   fill = switchConsequencesGene),
               alpha = 0.5,
               outlier.shape = NA,
               show.legend = F)+
  scale_fill_manual(values = c("TRUE" = "#B0413E", "FALSE" = "#FCAA67"),
                    labels = c("TRUE" = "NMD Isoforms", "FALSE" = "No NMD")) +
  scale_color_manual(values = c("TRUE" = "#B0413E", "FALSE" = "#FCAA67"),
                     labels = c("TRUE" = "NMD Isoforms", "FALSE" = "No NMD")) +
  theme_bw() +
  labs(y = "log10(ORF length)",
       x = "Switch Consequence",
       fill = "Switch Consequence") +
  scale_x_discrete(labels = c("FALSE" = "No NMD", "TRUE" = "NMD Causing"))
length_plot
ggsave("PRPF6_ISAR_lengths.pdf",
       plot = length_plot,
       height = 10,
       width = 10,
       units = "in",
       device = pdf,
       dpi = 300)
write_csv(switching_ORF, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/PRPF6_switch_ORF_length.csv")