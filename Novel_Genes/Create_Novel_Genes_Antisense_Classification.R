#!/bin/sh
# 14/08/2020: Rscript to extract the novel antisense genes from SQANTI2 classification file for mouse and human dataset

#********************** Variables and input files
library("stringr")
library("dplyr")
library("tidyr")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Input_Variables.R")

#********************** Apply function and identity novel genes 
sqanti_files()
sqanti_gtf()
antisense_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Revised_Paper/Novel_Genes/Antisense"

############# Antisense
Human_antisense <- class.files$Human[grepl("NOVEL", class.files$Human$associated_gene) & class.files$Human$structural_category == "Antisense",]
Mouse_antisense <- class.files$Mouse[grepl("NOVEL", class.files$Mouse$associated_gene) & class.files$Mouse$structural_category == "Antisense",]
write.table(Human_antisense,paste0(antisense_dir, "/Combined_Human_Novel_Antisense_classification.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(Mouse_antisense,paste0(antisense_dir, "/Mouse_Novel_Antisense_classification.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)

# gtf
gtf.files <- lapply(sqanti.gtf.names.files, function(x) read.table(x, sep = '\t') %>% mutate(pbid = word(word(.$V9, c(2), sep = fixed(";")),c(3), sep = fixed(" "))))
Human_antisense_gtf <- gtf.files$Human[gtf.files$Human$pbid %in% Human_antisense$isoform,][,c(1:9)]
Mouse_antisense_gtf <- gtf.files$Mouse[gtf.files$Mouse$pbid %in% Mouse_antisense$isoform,][,c(1:9)]
write.table(Human_antisense_gtf,paste0(antisense_dir, "/combined.NovelGenes.Antisense.gtf"), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(Mouse_antisense_gtf,paste0(antisense_dir, "/WholeIsoSeq.NovelGenes.Antisense.gtf"), quote = F, sep = "\t", col.names = F, row.names = F)

############ No Antisense
Human_noantisense <- class.files$Human[class.files$Human$structural_category != "Antisense",]
Mouse_noantisense <- class.files$Mouse[class.files$Mouse$structural_category != "Antisense",]
write.table(Human_noantisense,paste0(antisense_dir, "/Combined_Human_No_Antisense_classification.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(Mouse_noantisense ,paste0(antisense_dir, "/Mouse_Novel_No_Antisense_classification.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)

Human_noantisense_gtf <- gtf.files$Human[gtf.files$Human$pbid %in% Human_noantisense$isoform,][,c(1:9)]
Mouse_noantisense_gtf <- gtf.files$Mouse[gtf.files$Mouse$pbid %in% Mouse_noantisense$isoform,][,c(1:9)]
write.table(Human_noantisense_gtf ,paste0(antisense_dir, "/combined.No_Antisense.gtf"), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(Mouse_noantisense_gtf,paste0(antisense_dir, "/WholeIsoSeq.No_Antisense.gtf"), quote = F, sep = "\t", col.names = F, row.names = F)
