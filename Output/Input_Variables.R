# Szi Kay Leung: sl693@exeter.ac.uk
# Global files for analysis 

library("stringr")

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/SQANTI_General.R")
# SQANTI_class_preparation(input.class.file, output.class.file)
# According to SQANTI_report2.R (version 8.7) to prepare classification file in terms of ID columns etc 

#********************** Variables and input files
root <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/"
mouse_sqanti_dir <- paste0(root,"Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/WholeIsoSeq")
adult_sqanti_dir <- paste0(root,"Whole_Transcriptome/Human/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/AdultCTX")
fetal_sqanti_dir <- paste0(root,"Whole_Transcriptome/Human/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/FetalCTX")
human_sqanti_dir <- paste0(root,"Whole_Transcriptome/Human/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/HumanCTX")
humanHIP_sqanti_dir <- paste0(root,"Whole_Transcriptome/Human/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/FetalHIP")
humanSTR_sqanti_dir <- paste0(root,"Whole_Transcriptome/Human/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/FetalSTR")

mouse_lnc_sqanti_dir <- paste0(root,"Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/LNCRNA/WholeIsoSeq")
human_lnc_sqanti_dir <- paste0(root,"Whole_Transcriptome/Human/Post_IsoSeq/SQANTI_TAMA_FILTER/LNCRNA/HumanCTX")

# RNASeq Stringtie Transcriptome
rnaseq_transcriptome <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/SQANTI2/rnaseq_stringtie_merged_final_classification.filtered_lite_classification.txt"
cuffrefmap_input <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/SQANTI2/long_short_transcripts.rnaseq_stringtie_merged_final_classification.filtered_lite.gtf.refmap"
cufftmap_input <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/SQANTI2/long_short_transcripts.rnaseq_stringtie_merged_final_classification.filtered_lite.gtf.tmap"


class_suffix <- "_sqantitamafiltered.classification.txt"
gtf_suffix <- "_sqantitamafiltered.classification.gtf"
junc_suffix <- "_sqantitamafiltered.junction.txt" 

# Mouse ERCC
ERCC_sqanti_dir <- paste0(root,"Whole_Transcriptome/All_Tg4510/ERCC/WholeIsoSeq")

# INPUT: CCS Directory 
# Directory containing ccs read lengths (From CCS_Length.sh)
list_ccs_input_dir <- function(){
  ccs_input_dir <- list(paste0(root,"Whole_Transcriptome/All_Tg4510/IsoSeq/CCS/Lengths"),
                        paste0(root, "Whole_Transcriptome/Human/IsoSeq/CCS/Lengths/FetalCTX"),
                        paste0(root, "Whole_Transcriptome/Human/IsoSeq/CCS/Lengths/AdultCTX"))
  
  names(ccs_input_dir) <- c("Mouse","Fetal","Adult")
  
  ccs_input_dir <<- ccs_input_dir
}

# INPUT: SQANTI FILTERED GTF 
sqanti_gtf <- function(){
  
  sqanti.gtf.names.files <- list(
    paste0(mouse_sqanti_dir, gtf_suffix),
    paste0(adult_sqanti_dir , gtf_suffix),
    paste0(fetal_sqanti_dir , gtf_suffix),
    paste0(human_sqanti_dir ,gtf_suffix)
  )
  names(sqanti.gtf.names.files) <- c("Mouse","Adult","Fetal","Human")
  sqanti.gtf.names.files <<- sqanti.gtf.names.files
  
}

ERCC_sqanti_files <- function(){
  Ercc.class.file <- read.table(paste0(ERCC_sqanti_dir,class_suffix), header=T, as.is=T, sep="\t")
  Ercc.class.file <<- Ercc.class.file
}

# INPUT: SQANTI FILTERED DATA!!!
sqanti_files <- function(){
  
  mouse.class.file <- paste0(mouse_sqanti_dir,class_suffix)
  adult.class.file <- paste0(adult_sqanti_dir,class_suffix)
  fetal.class.file <- paste0(fetal_sqanti_dir,class_suffix)
  human.class.file <- paste0(human_sqanti_dir,class_suffix)
  
  class.names.files <- list(mouse.class.file, adult.class.file, fetal.class.file, human.class.file)
  names(class.names.files) <- c("Mouse","Adult","Fetal","Human")
  
  class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"standard"))
  
  class.files$Mouse$associated_gene <- toupper(class.files$Mouse$associated_gene)
  class.files$Adult$associated_gene <- toupper(class.files$Adult$associated_gene)
  class.files$Fetal$associated_gene <- toupper(class.files$Fetal$associated_gene)
  class.files$Human$associated_gene <- toupper(class.files$Huma$associated_gene)
  
  class.files$Mouse$Sample <- "Mouse"
  class.files$Adult$Sample <- "Human (Adult)"
  class.files$Fetal$Sample <- "Human (Fetal)"
  class.files$Human$Sample <- "Human"
  
  #assign(class.names.files, class.names.files, envir=.GlobalEnv)
  #assign(class.files, class.files, envir = .GlobalEnv, inherits = T)
  class.names.files <<- class.names.files
  class.files <<- class.files
  
}

# fetal hippocampus and fetal striatum SQANTI2 files
hipstr_sqanti_files <- function(){
  
  fetalhip.class.file <- human.class.file <- paste0(humanHIP_sqanti_dir,class_suffix)
  fetalstr.class.file <- human.class.file <- paste0(humanSTR_sqanti_dir,class_suffix)
  
  class.names.files <- list(fetalhip.class.file, fetalstr.class.file)
  names(class.names.files) <- c("FetalHip","FetalStr")
  
  class.files <- lapply(class.names.files, function(x) SQANTI_class_preparation(x,"standard"))
  
  class.files$FetalHip$associated_gene <- toupper(class.files$FetalHip$associated_gene)
  class.files$FetalStr$associated_gene <- toupper(class.files$FetalStr$associated_gene)
  
  hipstr.class.files <<- class.files
  
}

rnaseq_sqanti_files <- function(){
  rnaseq.class.files <- SQANTI_class_preparation(rnaseq_transcriptome,"rnaseq") %>% mutate(associated_gene = toupper(.$associated_gene))
  rnaseq.class.files$Sample <- "RNA-Seq"
  rnaseq.class.files <<- rnaseq.class.files
}

# INPUT: JUNC Files
junc_files <- function(){
  mouse.junc.file <- paste0(mouse_sqanti_dir,junc_suffix)
  adult.junc.file <- paste0(adult_sqanti_dir,junc_suffix)
  fetal.junc.file <-  paste0(fetal_sqanti_dir,junc_suffix)
  human.junc.file <- paste0(human_sqanti_dir,junc_suffix)
  
  junc.names.files <- list(mouse.junc.file,adult.junc.file,fetal.junc.file,human.junc.file)
  names(junc.names.files) <- c("Mouse","Adult","Fetal","Human")
  
  # Process Junction file 
  junc.files <- lapply(junc.names.files, function(x) SQANTI_junction_preparation(x))
  junc.files <<- junc.files
  
}

# INPUT: lncRNA Classification files 
lncrna_class_files <- function(){
  
  mouse.lnc.file <- paste0(mouse_lnc_sqanti_dir, class_suffix)
  adult.lnc.file <- paste0(human_lnc_sqanti_dir, class_suffix)
  fetal.lnc.file <- paste0(human_lnc_sqanti_dir, class_suffix)
  human.lnc.file <- paste0(human_lnc_sqanti_dir, class_suffix)
  
  lnc.names.files <- list(mouse.lnc.file,adult.lnc.file,fetal.lnc.file,human.lnc.file)
  names(lnc.names.files) <- c("Mouse","Adult","Fetal","Human")
  
  # Process Junction file 
  lnc.files <- lapply(lnc.names.files, function(x) SQANTI_class_preparation(x,"standard"))
  
  lnc.files$Mouse$associated_gene <- toupper(lnc.files$Mouse$associated_gene)
  lnc.files$Adult$associated_gene <- toupper(lnc.files$Adult$associated_gene)
  lnc.files$Fetal$associated_gene <- toupper(lnc.files$Fetal$associated_gene)
  lnc.files$Human$associated_gene <- toupper(lnc.files$Huma$associated_gene)
  
  lnc.files <<- lnc.files
}

mytheme <- function(){
  mytheme <- theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   text=element_text(size=20,  family="ArialMT", colour = "black"),
                   axis.title.x = element_text(vjust=-0.5),
                   axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 15, b = 0, l = 0)), 
                   legend.position = "bottom") 
  
  mytheme <<- mytheme 
}



# INPUT: REFERENCE 
# Gencode gtf file was extracted using Aaron's code: 
#cat gencode.vM20.primary_assembly.annotation.gtf | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"transcript") print a[2]"\t"a[4]"\t"$1":"$4"-"$5"\t"a[3]"\t"$7}' |sed 's/transcript_id "//' | sed 's/transcript_id "//' | sed 's/gene_type "//'| sed 's/gene_name "//' | sed 's/"//g' | awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | sed "1i\Transcriptid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength"  > gencode.vM20_gene_annotation_table.txt
reference_files <- function(){
  
  mouse_reference <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/ISOSEQ/reference/gencode.vM22_gene_annotation_table.txt", header = TRUE)
  human_reference <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/ISOSEQ/reference/gencode.v31_gene_annotation_table.txt", header = TRUE)
  
  reference_files <- list(mouse_reference, human_reference)
  names(reference_files) <- c("Mouse","Human")
  
  # label list of dataframes with sample name for distinguishing in later aggregation
  reference_files$Mouse$Sample <- "Mouse"
  reference_files$Human$Sample <- "Human"
  
  # Capitalise Mouse reference files 
  reference_files$Mouse$GeneSymbol <- toupper(reference_files$Mouse$GeneSymbol)
  
  reference_files <<- reference_files
}

# RAREFACTION at gene level and isoform level 
rarefaction_files <- function(){
  human_root_rarefaction <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/Post_IsoSeq/RAREFACTION/"
  # list input directory of rarefaction files (prepared)
  Mouse_rarefaction_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/RAREFACTION/WholeIsoSeq"
  Adult_rarefaction_dir <- paste0(human_root_rarefaction,"AdultCTX")
  Fetal_rarefaction_dir <- paste0(human_root_rarefaction,"FetalCTX")
  Human_rarefaction_dir <- paste0(human_root_rarefaction,"HumanCTX")
  rarefaction_dir <- list(Mouse_rarefaction_dir, Adult_rarefaction_dir, Fetal_rarefaction_dir, Human_rarefaction_dir)
  names(rarefaction_dir) <- c("Mouse","Human (Adult)","Human (Fetal)","Human")
  
  # Read input files pertaining to genes across all datasets
  all_rarefaction_genes <- list()
  for(i in 1:length(rarefaction_dir)){
    all_rarefaction_genes[[i]] <- read.table(paste0(rarefaction_dir[i],'.rarefaction.by_refgene.min_fl_2.txt'),sep=' ',header=T,skip=1)
    all_rarefaction_genes[[i]]$Sample <- names(rarefaction_dir)[i]
    all_rarefaction_genes[[i]]$type <- "Genes"
  }
  names(all_rarefaction_genes) <- c("Mouse","Human (Adult)","Human (Fetal)","Human")
  
  # Read input files pertaining to isoforms across all datasets
  all_rarefaction_isoforms <- list()
  for(i in 1:length(rarefaction_dir)){
    all_rarefaction_isoforms [[i]] <- read.table(paste0(rarefaction_dir[i],'.rarefaction.by_refisoform.min_fl_2.txt'),sep=' ',header=T,skip=1)
    all_rarefaction_isoforms [[i]]$Sample <- names(rarefaction_dir)[i]
    all_rarefaction_isoforms[[i]]$type <- "Isoforms"
  }
  names(all_rarefaction_isoforms) <- c("Mouse","Human (Adult)","Human (Fetal)","Human")
  
  # Read input files pertaining to isoforms by category across all datasets
  # Differentiate whether the type of category of transcripts is Annotated, Novel, or others for plotting 
  all_rarefaction_isoforms_category <- list()
  for(i in 1:length(rarefaction_dir)){
    all_rarefaction_isoforms_category[[i]] <- read.table(paste0(rarefaction_dir[i],'.rarefaction.by_refisoform.min_fl_2.by_category.txt'),sep=' ',header=T,skip=1)
    all_rarefaction_isoforms_category[[i]]$Sample <- names(rarefaction_dir)[i]
    all_rarefaction_isoforms_category[[i]]$category <- factor(all_rarefaction_isoforms_category[[i]]$category, 
                                                              levels = c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", 
                                                                         "novel_not_in_catalog", "fusion", "antisense", "genic", "intergenic"),
                                                              labels = c("FSM", "ISM", "NIC", "NNC", "Fusion", "Antisense", "Genic", "Intergenic"))
    
    for(j in 1:nrow(all_rarefaction_isoforms_category[[i]])){
      all_rarefaction_isoforms_category[[i]]$type[j] <- if(all_rarefaction_isoforms_category[[i]]$category[j] %in% c("FSM", "ISM")){
        "Annotated" 
      } else if (all_rarefaction_isoforms_category[[i]]$category[j] %in% c("NIC", "NNC")){
        "Novel"
      } else {
        "Others"
      }
    }
  }
  
  names(all_rarefaction_isoforms_category) <- c("Mouse","Adult","Fetal","Human")
  
  all_rarefaction_genes <<- all_rarefaction_genes
  all_rarefaction_isoforms <<- all_rarefaction_isoforms
  all_rarefaction_isoforms_category <<- all_rarefaction_isoforms_category 
  
}


# output file from TAMA merge 
TAMA_transfile <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/Mouse/Whole_vs_Targeted/merged_whole_targeted_trans_report.txt", header = T) %>% mutate(gene_id = word(transcript_id, c(1), sep = fixed("."))) %>% mutate(gene_name = word(word(all_source_trans, c(3), sep = fixed("_")),c(1), sep = ","))

######################################################################
# Merge all class.output for plot

# pre_and_post 
# prerequisite to run function: pre_sqanti_files and sqanti_files 
pre_and_post <- function(){
  
  all_class.output <- bind_rows(class.files[1:3],pre.class.files[1:3])
  all_class.output$SQANTI <- factor(as.character(all_class.output$SQANTI),
                                    levels = c("pre","post"))
  
  all_class.output <<- all_class.output
}

# prerequisite to run function: sqanti_files 
# work with all structural cateogory minus transcripts that are monoexonic
multiexonic <- function(){
  multiexonic <- lapply(class.files, function(x) x %>% filter(subcategory != "mono-exon"))
  multiexonic <<- multiexonic
}

FSM <- function(){
  FSM <- lapply(class.files, function(x) x %>% filter(structural_category == "FSM"))
  FSM <<- FSM
}


##### FeatureCounts 
# FeatureCounts output from RNA2IsoSeq (aligned RNASeq to Iso-Seq defined transcriptome)
#IsoSeq_Def <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsovsRnaseq/IsoSeq_Def/RNA2IsoSeq.transcript_id.tsv", header = T)
# FeatureCounts output from RNA2RNASeq (aligned RNASeq to RNA-Seq defined transcriptome)
#RNASeq_Def <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsovsRnaseq/RNASeq_Def/RNA2RNASeq.transcript_id.tsv", header = T)
##### Kallisto
RNASeq_Def <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq/KALLISTO/WholeIsoSeq.abundance.tsv", header = T) 
IsoSeq_Def <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME/KALLISTO/WholeIsoSeq.abundance.tsv", header = T) %>% mutate(target_id = word(target_id, c(1), sep = fixed("<")))
