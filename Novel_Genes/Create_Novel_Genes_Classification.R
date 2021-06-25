#!/bin/sh
# 23/07/2020: Rscript to extract the novel genes from SQANTI2 classification file for mouse and human dataset 

#********************** Variables and input files
library("stringr")
library("dplyr")
library("tidyr")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Input_Variables.R")

#********************** Apply function and identity novel genes 
sqanti_files()
output_blast_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Revised_Paper/Novel_Genes/Human2Mouse"

write.table(class.files$Human[grepl("NOVEL", class.files$Human$associated_gene),],paste0(output_blast_dir, "/Combined_Human_Novel_classification.txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(class.files$Mouse[grepl("NOVEL", class.files$Mouse$associated_gene),],paste0(output_blast_dir, "/Mouse_Novel_classification.txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
