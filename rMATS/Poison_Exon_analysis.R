#Created on 9/28/2023
#Created by: Caleb Embree
#Modified on: 4/26/2024
#Modified by: Caleb Embree
#Change log: 10/17/23: Added the analysis to select PE based on start coordinates being located between 
#start and end coordinates of annotated exons
#1/15/2023: Started trying to get coordinates from the Saltzman paper
#4/26/2024: Worked on doing liftover of poison exon coordinates
#4/29/2024: Pulled a list of NMD biotype isoforms from genes on Thomas PE list
#5/9/2024-5/14/2024: Rewrote the way to get PE locations from the Thomas paper


setwd("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Poison_Exons")
library(tidyverse)
library(dplyr)
library(readxl)
library(rtracklayer)
library(biomaRt)
#Import the intial dataset
Thomas_PE_list <- read_excel("C:/Users/Caleb/OneDrive - The Ohio State University/BioinfoData/Poison_Exons/Thomas_2021_PE_list.xlsx", 
                            sheet = "Table2", skip = 1)
View(Thomas_PE_list)


##Extract the Upstream exon end and downstream exon start coordinates from the list
PE_list = Thomas_PE_list %>% filter(target_type == "conservedPE") #Filter to the  poison exons
PE_wide = separate_wider_delim(PE_list,target,delim = ":", #Turns the target column into multiple columns wherever the : is
                               names = c("chr","coord_1","coord_2","strand","coord_3","coord_4","extra")) 
PE_wide = PE_wide %>% mutate(chr = str_replace(chr,"se@","chr"), #change the annotation to correctly refer to chromosome
                             strand = str_remove(strand, "\\|.*")) #Just display the strand
PE_wide = PE_wide %>% mutate(upstream_EE = if_else(strand == "+", #Pick the upstream exon end site based on if the strand is + or -
                                                   coord_1,
                                                   coord_3), #For - strand have to do the downstream ES before upstream EE because liftover doesn't work on transcripts where the second number is lower
                             downstream_ES = if_else(strand == "+", #Pick the downstream exon start site based on if the strand is + or -
                                                     coord_4,
                                                     coord_2))
PE_wide = PE_wide %>% unite("PE_range", 15:16, sep = ":", remove = FALSE)
PE_unique = PE_wide %>% group_by(gene) %>% distinct(PE_range,.keep_all = TRUE) %>% ungroup()
Thomas_sum = PE_unique %>% summarise(n = n(),
                                   nGene = n_distinct(gene))

##Convert the coordinates from hg19 to hg38 (GRch37 to GRch38)
#Have to make the coordinates into a BED file format
PE_old = PE_unique %>% mutate(chrom = chr,
                            chromStart = upstream_EE,
                            chromEnd = downstream_ES,
                             name = gene,
                            score = 0,
                            Strand = strand) %>% 
  select(18:23)


write.table(PE_old,"PE_bed_old.bed",
            row.names = F, col.names = F, quote = F)
#The BED file was provided to https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter and converted. 
PE_convert = rtracklayer::import("PE_BED_GRch38.bed",format = "bed")
newPE = as.data.frame(PE_convert)
newPE = newPE %>% select(1:4,6) %>% rename(chromosome = seqnames,
                                           intron_start = start,
                                           intron_end = end,
                                           gene_id = name) %>% 
  unite("PE_range", 2:3, sep = ":", remove = FALSE)

#Get the coordinates of exons from BioMart
listMarts()
ensembl <- useMart("ensembl",host = "https://feb2023.archive.ensembl.org") #use version 109 of ensembl
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)
exon_coordinates <- getBM(attributes = c("external_gene_name","ensembl_transcript_id","ensembl_gene_id",
                                         "ensembl_exon_id","chromosome_name","exon_chrom_start", "exon_chrom_end","strand"),
                   filters = "ensembl_gene_id",
                   values = newPE$gene_id,
                   mart = ensembl)
exon_coordinates_match = exon_coordinates %>% 
  full_join(newPE, by = c("ensembl_gene_id" = "gene_id"), relationship = "many-to-many") #Combine the coordinate table and liftover intron cordinates
exon_coordinates_summ = exon_coordinates_match %>% group_by(ensembl_exon_id) %>% 
  summarise(n = n(),
            groups = n_distinct(PE_range))
exon_coordinates_match = exon_coordinates_match %>% rowwise() %>% # Filter to exons where the start and end coordinates are in the intron range surrounding the PE
  filter(exon_chrom_start >= intron_start &
           exon_chrom_start <= intron_end &
           exon_chrom_end >= intron_start &
           exon_chrom_end <= intron_end) %>% 
  ungroup() %>% 
  distinct(ensembl_exon_id, .keep_all = TRUE)

#Now filter to exons where the annotated PE start is within the exon
PE_starts = PE_unique %>% mutate(PE_start = if_else(strand == "+", #Make a new column of start coordinates based on strand
                                                    coord_2,
                                                    coord_1))
PE_starts = PE_starts %>% mutate(end_coord = as.numeric(PE_start) + 5) %>% 
  select(chr,PE_start,end_coord,gene)#filter to just the columns we want
write.table(PE_starts,"PE_start_bed_old.bed",
            row.names = F, col.names = F, quote = F) #Write a bed file so that I can convert to GRch38
newPE_starts = rtracklayer::import("PE_start_GRch38.bed",format = "bed") #import the converted coordinates
newPE_starts = as.data.frame(newPE_starts)
newPE_starts = newPE_starts %>% select(2,6) %>% rename(PE_start = start)

exon_coordinates_match = exon_coordinates_match %>% left_join(newPE_starts,
                                                              by = c("ensembl_gene_id" = "name"),
                                                              relationship = "many-to-many") %>% #Adds in the PE_start column
  filter(PE_start >= exon_chrom_start & PE_start <= exon_chrom_end) %>% #filter to just the exons where the PE_start coordinate is within the exon coordinates
  distinct(ensembl_exon_id, .keep_all = TRUE) #remove any duplicates that may have been added from a PE having more than one start site listed
match_exon_summary = exon_coordinates_match %>% summarise(n = n(),
                                                          genes = n_distinct(ensembl_gene_id), #215 genes
                                                          exons = n_distinct(ensembl_exon_id)) #548 exons
#We have less genes than in the unique PE list, but a similar number of exons. I would assume that not all of these exons are actual poison exons
PE_coords = exon_coordinates_match %>% select(1,3:8,10,14)
write_csv(PE_coords, "Thomas_PE_Coordinates.csv")


####Get poison exons from Saltzman et al. Mol. & Cell. Bio. 2010####
SaltzmanTable <- read_excel("Saltzman_Supplemental_Table2_3Aug2010.xls", 
                            skip = 2)
SaltzmanTable %>% group_by(`PTC (for this AS event)`) %>% summarise(n())
Saltzman_PE = SaltzmanTable %>% rename(PTC = `PTC (for this AS event)`) %>% filter(PTC == "PTC upon inclusion" |
                                                                                   PTC == "PTC upon skipping") #Filters to just the AS events that result in a PTC
Saltzman_PE = Saltzman_PE %>% mutate(PTC= str_remove(PTC, "PTC upon "),
                                     Chromosome = str_remove(Chromosome,"chr")) %>%
  rename(HG18_start = A_coord_hg18start,HG18_end = A_coord_hg18end)
Saltzman_PE %>% group_by(PTC) %>% summarise(n()) #165 inclusion 331 skipping
Saltzman_PE %>% summarise(n_distinct(`LocusLink ID`)) #453 genes

##Get coordinates from biomart
Saltzman_geneID = getBM(attributes = c("external_gene_name","ensembl_gene_id"),
                   filters = "external_gene_name",
                   values = Saltzman_PE$`Gene Name`,
                   mart = ensembl)
Saltzman_geneID %>% summarise(n_distinct(external_gene_name)) #344 genes
missing_genes = Saltzman_PE %>% anti_join(Saltzman_geneID, by = c("Gene Name" = "external_gene_name")) %>% 
  select(`Gene Name`) %>% distinct()
write_csv(missing_genes, "Saltzman_missing_genes.csv")
#I manually identified the new names of the genes
Saltzman_missing_genes_filled <- read_csv("Saltzman_missing_genes_filled.csv") #110 genes
New_saltzman_geneID = getBM(attributes = c("external_gene_name","ensembl_gene_id"),
                          filters = "external_gene_name",
                          values = Saltzman_missing_genes_filled$New_name,
                          mart = ensembl)
Saltzman_geneID = Saltzman_geneID %>% full_join(New_saltzman_geneID) %>% distinct(ensembl_gene_id,.keep_all = T) #525 genes
#Get the exon IDs of all genes in the list using the HG18 build of the genome (ensembl 54)
HG18 = useMart("ensembl",host="https://may2009.archive.ensembl.org")
HG18 <- useDataset("hsapiens_gene_ensembl",mart=HG18)
attributes_HG18 <- listAttributes(HG18)
filters_HG18 = listFilters(HG18)
HG18_exons = getBM(attributes = c("external_gene_id","ensembl_gene_id","ensembl_exon_id","chromosome_name",
                                  "exon_chrom_start","exon_chrom_end"),
                   filters = "ensembl_gene_id",
                   values = Saltzman_geneID$ensembl_gene_id,
                   mart = HG18)
HG18_exons %>% summarise(n_distinct(ensembl_gene_id)) #455 genes
Saltzman_hg18_exonID = Saltzman_PE %>% inner_join(HG18_exons,by = c("HG18_start" = "exon_chrom_start"))

#Pull genes first from the locus link (AKA Enterez ID) then get exon coordinates
HG18_genes = getBM(attributes = c("external_gene_id","ensembl_gene_id","entrezgene"),
                   filters = "entrezgene",
                   values = Saltzman_PE$`LocusLink ID`,
                   mart = HG18) #459 genes
HG18_exon_coords = getBM(attributes = c(  "ensembl_exon_id","chromosome_name",
                                          "exon_chrom_start","exon_chrom_end","ensembl_gene_id"),
                         filters = "ensembl_gene_id",
                         values = HG18_genes$ensembl_gene_id,
                         mart = HG18) #459 genes
#This did not work
Saltzman_PE_coords = Saltzman_PE %>% left_join(HG18_exon_coords, by = c("HG18_start" = "exon_chrom_start"))


#### Check if each gene from the Saltzman list has an NMD biotype ####
Saltzman_inclusion = Saltzman_PE %>% filter(PTC == "inclusion")
Saltzman_inclusion %>% summarise(n = n(), genes = n_distinct(`LocusLink ID`)) #165 events, 157 genes
Inclusion_ids = getBM(attributes = c("external_gene_id","ensembl_gene_id","entrezgene"),
                      filters = "entrezgene",
                      values = Saltzman_inclusion$`LocusLink ID`,
                      mart = HG18) #157 genes
inclusion_biotype = getBM(attributes = c("external_gene_name", "ensembl_gene_id","ensembl_transcript_id",
                                        "transcript_biotype"),
                         filters = "ensembl_gene_id",
                         values = Inclusion_ids$ensembl_gene_id,
                         mart = ensembl) #157 genes
biotype_sum = inclusion_biotype %>% group_by(ensembl_gene_id) %>% 
  summarise(n = n(),
            NMD = sum(str_detect(transcript_biotype,"nonsense")))
inclusion_nmd = biotype_sum %>% filter(NMD > 0)

inclusion_NMD_transcripts = inclusion_biotype %>% filter(transcript_biotype == "nonsense_mediated_decay")
write_csv(inclusion_NMD_transcripts, file = "Saltzman_inclusion_NMD_transcripts.csv")

####Get the NMD biotype genes from Thomas list####
PE_summary = PE_list %>% summarise(n = n(), genes = n_distinct(gene)) #4099 guide RNAs, 348 genes
PE_biotype = getBM(attributes = c("external_gene_name","ensembl_gene_id","ensembl_transcript_id",
                                  "transcript_biotype","transcript_mane_select"),
                   filters = "ensembl_gene_id",
                   values = PE_list$gene,
                   mart = ensembl) #346 genes
PE_biotype_sum = PE_biotype %>% group_by(ensembl_gene_id) %>% 
  summarise(n = n(),
            NMD = sum(str_detect(transcript_biotype,"nonsense"))) %>% 
  filter(NMD > 0) #273 genes with NMD biotype transcripts
PE_NMD_transcripts = PE_biotype %>% filter(transcript_biotype == "nonsense_mediated_decay")
write_csv(PE_NMD_transcripts, file = "Thomas_NMD_transcripts.csv") #933 transcripts

####Combine the two PE gene lists to make a combined list####
Saltzman_NMD = inclusion_NMD_transcripts %>% mutate(Source = "Saltzman") %>%
  filter(!ensembl_transcript_id %in% PE_NMD_transcripts$ensembl_transcript_id) #Removes the transcripts in the other list from this annotation
Thomas_NMD = PE_NMD_transcripts %>% select(1:4) %>% 
  mutate(Source = "Thomas") %>% #Annotate where the exons come from
  filter(!ensembl_transcript_id %in% inclusion_NMD_transcripts$ensembl_transcript_id) #Removes the transcripts in the other list
Both_NMD = inclusion_NMD_transcripts %>% inner_join(PE_NMD_transcripts) %>% 
  mutate(Source = "Both") #44 transcripts
full_NMD = Saltzman_NMD %>% full_join(Thomas_NMD) %>% full_join(Both_NMD) %>% 
  mutate(isoform = "NMD")
NMD_summary = full_NMD %>% group_by(Source) %>% 
  summarise(n = n(),
            genes = n_distinct(ensembl_gene_id),
            transcripts = n_distinct(ensembl_transcript_id)) #1280 transcripts, 379 genes
NMD_mane = getBM(attributes = c("external_gene_name","ensembl_gene_id","ensembl_transcript_id",
                                "transcript_biotype","transcript_mane_select"),
                 filters = "ensembl_gene_id",
                 values = full_NMD$ensembl_gene_id,
                 mart = ensembl)
NMD_mane = NMD_mane %>% filter(str_detect(transcript_mane_select,"NM")) %>%  #filter to just the MANE isoforms
  mutate(isoform = "MANE")
full_NMD = full_NMD %>% full_join(NMD_mane)
write_csv(full_NMD, file = "PE_NMD_transcripts_both_lists.csv")


####Identify exons in Salzman list only present in NMD Biotype####
inclusion_exons = getBM(attributes = c("external_gene_name","ensembl_gene_id","ensembl_transcript_id",
                                       "ensembl_exon_id","chromosome_name","exon_chrom_start", "exon_chrom_end"),
                        filters = "ensembl_transcript_id",
                        values = inclusion_NMD_transcripts$ensembl_transcript_id,
                        mart = ensembl)
inclusion_exons = inclusion_exons %>% distinct(ensembl_exon_id,.keep_all = TRUE)
write_csv(inclusion_exons,"Saltzman_inclusion_exons.csv")
