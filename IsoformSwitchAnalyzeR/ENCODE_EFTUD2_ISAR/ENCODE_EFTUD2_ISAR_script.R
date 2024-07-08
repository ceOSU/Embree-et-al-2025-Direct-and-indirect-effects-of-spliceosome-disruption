# Script to analyze isoforms changes using Isoform Switch Analyzer
# Date Written: 2/14/2023
# Written by: Caleb Embree
# Modified on: 
# Analyzing: KD of EFTUD2 via shRNA from ENCODE a
# Go to https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html for package manual


#Processing that happened before the script (on the supercomputer):
# 1. Transcripts were mapped via HISAT2
# 2. Transcripts were assembled with Stringtie
# 3. Transcripts were merged with Stringtie
# 4. Transcripts were assembled to the merged file with Stringtie

#Setup
#if (!requireNamespace("devtools", quietly = TRUE)){
#  install.packages("devtools")
#}
#if (!requireNamespace("pfamAnalyzeR", quietly = TRUE)){
#  devtools::install_github("kvittingseerup/pfamAnalyzeR")
#}
#devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR")
library(IsoformSwitchAnalyzeR)
library(rhdf5)
library(tidyverse)
library(biomaRt)
setwd("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_EFTUD2_ISAR")

#Import datasets
EFTUD2quant = importIsoformExpression(parentDir = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_EFTUD2_ISAR",
                                      addIsofomIdAsColumn = TRUE, 
                                      showProgress = TRUE,
                                      readLength = 100
                                      ) #reads the t_data.ctab file in each sub directory
names(EFTUD2quant)
tail(EFTUD2quant$abundance)
tail(EFTUD2quant$counts)

#Generate the list of isoform switches
myDesign = data.frame(sampleID = c("SRR4422087","SRR4422088","SRR4421357","SRR4421358"),
                      condition = c("WT","WT","KD","KD"))
comparisons = data.frame(condition_1 = "WT",
                         condition_2 = "KD")
SwitchList = importRdata(isoformCountMatrix = EFTUD2quant$counts,
                         isoformRepExpression = EFTUD2quant$abundance,
                         designMatrix = myDesign,
                         isoformExonAnnoation = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_EFTUD2_ISAR/st_merged.gtf",
                         isoformNtFasta = "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/ENCODE_EFTUD2_ISAR/EFTUD2_transcripts.fa",
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
rawTPM = EFTUD2quant$abundance
annotatedTPM = rawTPM %>% mutate(kdTPM = rowMeans(dplyr::select(rawTPM,2:3), na.rm = TRUE)) %>%
  mutate(wtTPM = rowMeans(dplyr::select(rawTPM, 4:5), na.rm = TRUE))

listMarts()
ensembl <- useMart("ensembl")
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
EFTUD2_TPM = PTC_TPM %>% full_join(Novel_TPM) %>% full_join(PTC_MANE_TPM)
write_csv(EFTUD2_TPM, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/EFTUD2_TPM.csv")

####Comparing the TPMs of different classes of novel transcript####
all_novel = all_switch %>% filter(IF2 <= 0) #829 transcripts
all_novel_tpm = all_novel %>% inner_join(rawTPM, by = "isoform_id")
all_novel_tpm = all_novel_tpm %>% mutate(kdTPM = rowMeans(select(all_novel_tpm,14:15), na.rm = TRUE)) %>% 
  mutate(wtTPM = rowMeans(select(all_novel_tpm, 15:16), na.rm = TRUE))#829 transcripts
all_novel_tpm_filt = all_novel %>% inner_join(filtTPM, by = "isoform_id")#586 transcripts
all_novel_tpm_filt = all_novel_tpm_filt %>% select(3:12,18:19) %>% mutate(sample = "EFTUD2")
write_csv(all_novel_tpm_filt, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/EFTUD2_novel_TPM.csv")

####Save the nt sequence of unanotated isoforms####
novel_list = subsetSwitchAnalyzeRlist(switchAnalyzeRlist = Switch_cons,
                                      subset = str_detect(Switch_cons$isoformFeatures$isoform_id,
                                                          "MST"))
EFTUD2_novel_sequence = extractSequence(switchAnalyzeRlist = novel_list,
                                     removeORFwithStop = FALSE,
                                     extractAAseq = FALSE,
                                     removeShortAAseq = FALSE,
                                     removeLongAAseq = FALSE,
                                     outputPrefix = "EFTUD2_novel")

#Look at transcript length
count(all_switch,switchConsequencesGene) #FALSE = 1974 TRUE = 1104
ORFanalysis = SwitchAndORF$orfAnalysis
count(ORFanalysis,PTC) #False = 11215 True = 3912 NA = 2195
ORFanalysis = ORFanalysis %>% mutate(Switch = isoform_id %in% all_switch$isoform_id)
count(ORFanalysis,Switch) #False = 14244 True=3078 ->Same number as switching transcripts
switching_ORF = ORFanalysis %>% filter(Switch == TRUE)
count(switching_ORF,PTC) #False = 2159 TRUE = 646 NA = 273 
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
ggsave("EFTUD2_ISAR_lengths.pdf",
       plot = length_plot,
       height = 10,
       width = 10,
       units = "in",
       device = pdf,
       dpi = 300)
write_csv(switching_ORF, "C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/IsoformSwitch/Combined_analysis/EFTUD2_switch_ORF_length.csv")
