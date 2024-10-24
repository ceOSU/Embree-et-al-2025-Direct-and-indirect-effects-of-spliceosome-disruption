# Script to analyze isoforms changes using Isoform Switch Analyzer
# Date Written: 9/6/2023
# Written by: Caleb Embree
# Modified on: 9/8/2023
# Analyzing: U2AF1 KD from ENCODE in K562
# Go to https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html for package manual


#Processing that happened before the script (on the supercomputer):
# 1. Transcripts were mapped via HISAT2
# 2. Transcripts were assembled with Stringtie
# 3. Transcripts were merged with Stringtie and FASTA file of the merged GTF was created
# 4. Transcripts were assembled to the merged file with Stringtie


#### Copy and Paste ####
#U2AF1 = Gene of Interest
#C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/U2AF1_ISAR = Path to the working directory
#SRR3469464 = First WT sample
#SRR3469465 = Second WT sample
#SRR3469458 = First KD sample
#SRR3469459 = Second KD sample


#Setup
#if (!requireNamespace("pfamAnalyzeR", quietly = TRUE)){
#  devtools::install_github("kvittingseerup/pfamAnalyzeR")
#}
#devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR")
library(IsoformSwitchAnalyzeR)
library(rhdf5)
library(tidyverse)
library(biomaRt)
setwd("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_U2AF1_ISAR")

####Import datasets####
U2AF1quant = importIsoformExpression(parentDir = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_U2AF1_ISAR",
                                   addIsofomIdAsColumn = TRUE, 
                                   showProgress = TRUE,
                                   readLength = 100
) #reads the t_data.ctab file in each sub directory
names(U2AF1quant)
tail(U2AF1quant$abundance)
tail(U2AF1quant$counts)

#Generate the list of isoform switches
myDesign = data.frame(sampleID = c("SRR3469464","SRR3469465","SRR3469458","SRR3469459"),
                      condition = c("WT","WT","KD","KD"))
comparisons = data.frame(condition_1 = "WT",
                         condition_2 = "KD")
SwitchList = importRdata(isoformCountMatrix = U2AF1quant$counts,
                         isoformRepExpression = U2AF1quant$abundance,
                         designMatrix = myDesign,
                         isoformExonAnnoation = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_U2AF1_ISAR/st_merged.gtf",
                         isoformNtFasta = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_U2AF1_ISAR/U2AF1_transcripts.fa",
                         showProgress = TRUE,
                         comparisonsToMake = comparisons)
SwitchList_filt = preFilter(SwitchList,
                            geneExpressionCutoff = 1, #Cut off genes with less than 1 TPM
                            isoformExpressionCutoff = 0, #removes unused isoforms
                            removeSingleIsoformGenes = TRUE, #removes genes with a single isoform
) #The filtering removed 79967 ( 48.67% of ) transcripts. There is now 84332 isoforms left

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

all_switch_cons = extractTopSwitches(
  Switch_cons,
  filterForConsequences = TRUE,
  n = NA,
  extractGenes = FALSE,
  sortByQvals = FALSE)

####Determining the TPM of different classes of transcripts ####
#ID the NMD consequences of novel isoforms##
NMD_cons = Switch_cons$switchConsequence
NMD_cons = NMD_cons %>% filter(str_detect(isoformUpregulated,"MST")) %>% 
  mutate(NMD_consequence = str_remove(switchConsequence,"NMD ")) %>% 
  dplyr::select(gene_id,NMD_consequence,isoformUpregulated)

#Pull the TPM of different groups
rawTPM = U2AF1quant$abundance
annotatedTPM = rawTPM %>% mutate(kdTPM = rowMeans(dplyr::select(rawTPM,2:3), na.rm = TRUE)) %>% #Make sure these refer to the correct columns
  mutate(wtTPM = rowMeans(dplyr::select(rawTPM, 4:5), na.rm = TRUE))

listMarts()
ensembl <- useMart("ensembl",host = "https://feb2023.archive.ensembl.org") #use version 109 of ensembl
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
alltrans_ensembl <- getBM(attributes = c("ensembl_transcript_id_version","external_gene_name",
                                         "transcript_biotype","ensembl_transcript_id","transcript_mane_select","ensembl_gene_id"),
                          filters = "ensembl_transcript_id",
                          values = annotatedTPM$isoform_id,
                          mart = ensembl)
annotatedTPM = annotatedTPM %>% left_join(alltrans_ensembl, by = c("isoform_id" = "ensembl_transcript_id")) %>% 
  left_join(NMD_cons,by = c("isoform_id" = "isoformUpregulated"))
filtTPM = annotatedTPM %>% filter(kdTPM > 1 | wtTPM >1)
MANE_TPM = filtTPM %>% filter(!is.na(transcript_mane_select)) %>%
  filter(transcript_mane_select != "") %>% 
  dplyr::select(1:7,9:12)
Stringent_PTC <- read_csv("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Bioinformatics template/PTC_list_creation/Stringent_PTC_MANE_CE.csv")
PTC_TPM = filtTPM %>% inner_join(Stringent_PTC, by = c("isoform_id" = "transID")) %>%
  filter(PTC == TRUE)
Novel_TPM = annotatedTPM %>% filter(!is.na(NMD_consequence))
Novel_TPM = Novel_TPM %>% dplyr::select(1:7,9:11,13:14) %>%
  mutate(type = "Novel") %>%
  rename(ensembl_gene_id = gene_id)
PTC_TPM = PTC_TPM %>% dplyr::select(1:7,9:12,17) %>%
  mutate(type = "PTC")
PTC_MANE_TPM = MANE_TPM %>% inner_join(Stringent_PTC, by = c("isoform_id" = "transID"))
PTC_MANE_TPM = PTC_MANE_TPM %>% dplyr::select(1:11,14) %>%
  mutate(type = "MANE")
U2AF1_TPM = PTC_TPM %>% full_join(Novel_TPM) %>% full_join(PTC_MANE_TPM)
write_csv(U2AF1_TPM, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/U2AF1_TPM.csv")


####Comparing the TPMs of different classes of novel transcript####
all_novel = all_switch_cons %>% filter(IF2 <= 0) #X transcripts
all_novel_tpm = all_novel %>% inner_join(rawTPM, by = "isoform_id")
all_novel_tpm = all_novel_tpm %>% mutate(kdTPM = rowMeans(dplyr::select(all_novel_tpm,14:15), na.rm = TRUE))
all_novel_tpm = all_novel_tpm %>% mutate(wtTPM = rowMeans(dplyr::select(all_novel_tpm, 15:16), na.rm = TRUE))#X transcripts
all_novel_tpm_filt = all_novel %>% inner_join(filtTPM, by = "isoform_id")#X transcripts
all_novel_tpm_filt = all_novel_tpm_filt %>% dplyr::select(3:12,18:19) %>% mutate(sample = "U2AF1")
write_csv(all_novel_tpm_filt, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/U2AF1_novel_TPM.csv")

####Save the nt sequence of unanotated isoforms####
novel_list = subsetSwitchAnalyzeRlist(switchAnalyzeRlist = Switch_cons,
                                      subset = str_detect(Switch_cons$isoformFeatures$isoform_id,
                                                          "MST"))
head(novel_list$ntSequence)
U2AF1_novel_sequence = extractSequence(switchAnalyzeRlist = novel_list,
                                     removeORFwithStop = FALSE,
                                     extractAAseq = FALSE,
                                     removeShortAAseq = FALSE,
                                     removeLongAAseq = FALSE,
                                     outputPrefix = "U2AF1_novel")
novel_annotations = novel_list$isoformFeatures %>% dplyr::select(isoform_id,gene_id,gene_name,PTC)
write_csv(novel_annotations,"U2AF1_novel_iso_features.csv")

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
ggsave("U2AF1_ISAR_lengths.pdf",
       plot = length_plot,
       height = 10,
       width = 10,
       units = "in",
       device = pdf,
       dpi = 300)
write_csv(switching_ORF, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/U2AF1_switch_ORF_length.csv")