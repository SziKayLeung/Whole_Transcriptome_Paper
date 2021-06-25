# Szi Kay Leung
# Script to call functions for plots for IsoSeq Paper (2020)

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables/All_Tables_Functions.R")

### Variables (Input, Output) ########################################################################

output_plot_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Isoseq_Paper/Plots/Final"
output_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables"
input_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables/Other_Input"
output_corr_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Correlations"
reference_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/"
aaron_reference_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/ISOSEQ/reference/"
fetal_own_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Fetal/SQANTI2/own_rnaseq/"

# Input Files from Input_Variables.R
list_ccs_input_dir()        # Directory containing CCS lengths 
sqanti_files()              # SQANTI2 filtered classification files (Reference genome alignment)
junc_files()                # SQATNI2 filtered junction files
lncrna_class_files()        # SQANTI2 filtered classification files (Reference lncRNA alignment)
sqanti_gtf()                # SQANTI2 filtered gtf files (Reference genome alignment)
#rarefaction_files()         # Rarefaction curves plotted from Liz's scripts

blast_novel_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Revised_Paper/Novel_Genes/Human2Mouse/"
Mouse2HumanFlnc <- paste0(blast_novel_dir, "Mouse2HumanFlnc.blast.txt")
Human2MouseFlnc <- paste0(blast_novel_dir, "Human2MouseFlnc.blast.txt")
Human2MouseAllSqanti <- paste0(blast_novel_dir, "Human2MouseAllSqanti.blast.txt")
Mouse2HumanAllSqanti <- paste0(blast_novel_dir, "Mouse2HumanAllSqanti.blast.txt")


### Functions ########################################################################

top500Genes_enrichR(class.files$Mouse,"Mouse_Gene_Atlas","Mouse")
top500Genes_enrichR(class.files$Human,"Human_Gene_Atlas","Human")
top500Genes_enrichR(class.files$Fetal,"Human_Gene_Atlas","Fetal")
top500Genes_enrichR(class.files$Adult,"Human_Gene_Atlas","Adult")

top100Genes_enrichR(class.files$Human, "GO_Molecular_Function_2018","GOmolecularfunction", "Human")
top100Genes_enrichR(class.files$Mouse, "GO_Molecular_Function_2018","GOmolecularfunction", "Mouse")
top100Genes_enrichR(class.files$Human, "GWAS_Catalog_2019","GOgwascata", "Human")
top100Genes_enrichR(class.files$Mouse, "GWAS_Catalog_2019","GOgwascata", "Mouse")


IR_genes_enrichR(class.files$Mouse,"GO_Biological_Process_2018", "Mouse")
IR_genes_enrichR(class.files$Human,"GO_Biological_Process_2018", "Human")

### Number of annotated vs novel transcripts in annotated genes 
write.csv(annotated_genes_transcript_num()$Human, paste0(output_table_dir,"/Descriptive/Human_annotated_genes_transcript_num.csv"))
write.csv(annotated_genes_transcript_num()$Mouse, paste0(output_table_dir,"/Descriptive//Mouse_annotated_genes_transcript_num.csv"))
write.csv(annotated_genes_transcript_num()$Fetal, paste0(output_table_dir,"/Descriptive//Fetal_annotated_genes_transcript_num.csv"))
write.csv(annotated_genes_transcript_num()$Adult, paste0(output_table_dir,"/Descriptive//Adult_annotated_genes_transcript_num.csv"))

### Threshold
write.csv(threshold_gene_expression_human_mouse(), paste0(output_table_dir,"/Human_Mouse_Comparisons/threshold_gene_expression_human_mouse.csv"))
write.csv(threshold_gene_expression_human_mouse_common_threshold(), paste0(output_table_dir,"/threshold_gene_expression_human_mouse_common_threshold.csv"))


threshold <- threshold_gene_expression_human_mouse()


### Fusion
human_fusion_cg <- conjoin_prediction(paste0(input_table_dir,"/Akiva_supplementary_Table_S1.csv"),"Human")
mouse_fusion_cg <- conjoin_prediction(paste0(input_table_dir,"/Mouse_CGs.csv"),"Mouse")


write.csv(mouse_fusion_cg, paste0(output_table_dir,"/Fusion_Genes/Mouse_ConJoinedGenes_Common.csv"))
write.csv(human_fusion_cg, paste0(output_table_dir,"/Fusion_Genes/Human_ConJoinedGenes_Common.csv"))


### IR
IR_genes_only_IRtranscripts()

### NovelGenes
Mouse2HumanFlnc_hits <- novelgene_blastflnc(Mouse2HumanFlnc)
Human2MouseFlnc_hits <- novelgene_blastflnc(Human2MouseFlnc)
Human2MouseAllSqanti_hits <- blast(Human2MouseAllSqanti,"Human2Mouse")
Mouse2HumanAllSqanti_hits <- blast(Mouse2HumanAllSqanti,"Mouse2Human")

Human2MouseFlnc_hits$filteredblast %>% filter(PacBio_ID == "PB.14504.1") %>% arrange(evalue) %>% .[1,]
Human2MouseFlnc_hits$filteredblast %>% filter(PacBio_ID == "PB.21588.1") %>% arrange(evalue) %>% .[1,]

wb <- createWorkbook()
datas <- c(Human2MouseFlnc_hits,Mouse2HumanFlnc_hits)
sheetnames <- c("Human_Blast_Hits","Human_FilteredBlast_Hits","Human_Summary","Mouse_Blast_Hits","Mouse_FilteredBlast_Hits","Mouse_Summary")
sheets <- lapply(sheetnames, createSheet, wb = wb)
void <- Map(addDataFrame, datas, sheets)
saveWorkbook(wb, file = paste0(output_table_dir,"/Novel_Genes/NovelGenesAcrossSpecies.xlsx"))

Blast_Novel_Genes_Across_Genome()

