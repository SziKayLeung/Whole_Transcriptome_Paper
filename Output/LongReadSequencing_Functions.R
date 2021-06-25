## Functions for LongReadSequencing

suppressMessages(library(DT))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


#source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/human_mouse_isoseq/DownStream/QC_Figures/All_Plots_Functions.R")


# tabulating_sqanti_num: tabulating sqanti number of genes and isoforms
# input: list of classification file of all datasets (already read.table) 
# output: aggregated table of stats across all datasets
tabulating_sqanti_num <- function(type_class_file){
  dat <- data.frame()
  count=1
  for(i in type_class_file){
    # Total Unique Isoforms: tabulated by number of rows
    isoforms <- dim(i)[1]
    # Total Unique Genes: Remove novel genes and count 
    annotated_genes <- i[!grepl("NOVELGENE",i$associated_gene),] %>% count(associated_gene) %>% nrow(.)
    # Total Novel Genes: 
    novel_genes <- i[grepl("NOVELGENE",i$associated_gene),] %>% count(associated_gene) %>% nrow(.)
    # Total Number of genes 
    total_genes <- annotated_genes + novel_genes
    # FSM % of isoforms 
    FSM <- round(nrow(i[i$structural_category == "FSM",])/dim(i)[1]*100,2)
    # Number of isoforms (Range per Gene)
    isoform_count <- i %>% count(associated_gene) 
    min_isoform_count <- min(isoform_count$n)
    max_isoform_count <- max(isoform_count$n)
    # Number of Genes with >1 isoform
    num_isoform_morethanone <- nrow(isoform_count[isoform_count$n > 1,])
    num_isoform_morethanten <- nrow(isoform_count[isoform_count$n > 10,]) 
    perc_isoform_morethanone <- round((isoform_count[isoform_count$n > 1,] %>% nrow())/total_genes * 100, 2)
    perc_isoform_morethanten <- round((isoform_count[isoform_count$n > 10,] %>% nrow())/total_genes * 100, 2)
    
    # percentage of annotated and novel genes 
    annotated_genes <- paste0(annotated_genes, " (", round(annotated_genes/total_genes * 100,2), "%)")
    novel_genes <- paste0(novel_genes, " (", round(novel_genes/total_genes * 100,2), "%)")
    
    # Number of Annotated Isoforms (FSM, ISM) 
    annotated_isoforms <- paste0(nrow(i[i$associated_transcript != "novel",])," (",
                                 round(nrow(i[i$associated_transcript != "novel",])/isoforms * 100,2),"%)")
    novel_isoforms <- paste0(nrow(i[i$associated_transcript == "novel",])," (",
                             round(nrow(i[i$associated_transcript == "novel",])/isoforms * 100,2),"%)")
    
    # Number of Novel Isoforms
    
    # 9 levels of structural cateogory 
    struct <- vector("numeric", 9)
    for(num in 1:length(levels(i$structural_category))){
      cate <- nrow(i[i$structural_category == levels(i$structural_category)[num],])
      struct[num] <- paste0(cate, " (",
                            round(cate/dim(i)[1]*100,2),"%)") 
    }
    
    dat[1:10,count] <- rbind(total_genes,annotated_genes, novel_genes, isoforms, FSM, 
                             paste(min_isoform_count, max_isoform_count, sep = "-"),
                             paste0(num_isoform_morethanone, " (", perc_isoform_morethanone,"%)"),
                             paste0(num_isoform_morethanten, " (", perc_isoform_morethanten, "%)"),
                             annotated_isoforms, novel_isoforms)
    dat[11:19,count] <- struct
    colnames(dat)[count] <- names(type_class_file)[count]
    count = count + 1
  }
  
  row.names(dat) <- append(c("Total Number of Genes", "Annotated Genes","Novel Genes","Total Number of Isoforms", "Percentage of FSM", "Number of Isoforms per Gene", "Genes with > 1 Isoform", "Genes with > 10 Isoforms","Annotated Isoforms", "Novel Isoforms"),
                           levels(i$structural_category))
  dat
}

prepare_merge_numbers <- function(){
  # REFERENCE 
  # Gencode gtf file was extracted using Aaron's code: 
  #cat gencode.vM20.primary_assembly.annotation.gtf | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"transcript") print a[2]"\t"a[4]"\t"$1":"$4"-"$5"\t"a[3]"\t"$7}' |sed 's/transcript_id "//' | sed 's/transcript_id "//' | sed 's/gene_type "//'| sed 's/gene_name "//' | sed 's/"//g' | awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | sed "1i\Transcriptid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength"  > gencode.vM20_gene_annotation_table.txt
  
  FSM <- lapply(class.files, function(x) x%>% filter(structural_category == "FSM"))
  mouse_reference <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/ISOSEQ/reference/gencode.vM22_gene_annotation_table.txt", header = TRUE)
  human_reference <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/ISOSEQ/reference/gencode.v31_gene_annotation_table.txt", header = TRUE)
  
  reference_files <- list(mouse_reference, human_reference)
  names(reference_files) <- c("Mouse","Human")
  
  # label list of dataframes with sample name for distinguishing in later aggregation
  reference_files$Mouse$Sample <- "Mouse"
  reference_files$Human$Sample <- "Human"
  
  # Capitalise Mouse reference files 
  reference_files$Mouse$GeneSymbol <- toupper(reference_files$Mouse$GeneSymbol)
  
  ## Tabulate sum of transcripts per gene in reference files and isoseq files 
  gencode_transcripts <- function(df) df %>% group_by(GeneSymbol, Sample) %>% tally()
  isoseq_transcripts  <- function(df) df %>% group_by(associated_gene, Sample) %>% tally()
  
  reference_freq <- lapply(reference_files, function(x) gencode_transcripts(x))
  class_freq <- lapply(FSM, function(x) isoseq_transcripts(x))
  
  reference_freq <<- reference_freq
  class_freq <<- class_freq
}

## Merge
Merge_Numbers <- function(species1.df.classfreq, species1.df.referencefreq, species2.df.classfreq, species2.df.referencefreq, output_file){
  # Merge species 1 (i.e Mouse) numbers from reference and isoseq files by gene names
  # Normalise by dividing detected/known
  species1.df <- merge(species1.df.classfreq, species1.df.referencefreq, by.x = "associated_gene", by.y = "GeneSymbol", all = TRUE)
  species1.df$normalised <- species1.df$n.x/species1.df$n.y
  
  # Merge species 2 (i.e Human) numbers from reference and isoseq files by gene names
  # Normalise by dividing detected/known
  species2.df <- merge(species2.df.classfreq, species2.df.referencefreq, by.x = "associated_gene", by.y = "GeneSymbol", all = TRUE)
  species2.df$normalised <- species2.df$n.x/species2.df$n.y
  
  # Merge species 1 and species 2 by gene names for comparison
  merge_df <- merge(species1.df, species2.df, by = "associated_gene", all = TRUE)
  
  species1 <- species1.df.classfreq$Sample[1]
  species2 <- species2.df.classfreq$Sample[1]
  #print(species1)
  #print(species2)
  
  merge_df$normalised.x[is.na(merge_df$normalised.x)] <- 0
  merge_df$normalised.y[is.na(merge_df$normalised.y)] <- 0
  
  for(row in 1:nrow(merge_df)){
    merge_df$Detected[row] <- 
      if(merge_df[row,"normalised.x"] == 0){
        paste(species2)
      }else{
        if(merge_df[row,"normalised.y"] == 0){
          paste(species1)
        }else{
          paste("Both")
        }
      }
  }
  
  
  # Trajectory, histogram of Normalised counts 
  df1 <- data.frame(table(merge_df$normalised.x))
  df1$Var1 <-as.numeric(as.character(df1$Var1))
  df1$Sample <- paste(species1)
  
  df2 <- data.frame(table(merge_df$normalised.y))
  df2$Var1 <-as.numeric(as.character(df2$Var1))
  df2$Sample <- paste(species2)
  
  df3 <- rbind(df1,df2)
  
  # remove repeated Sample name columns
  merge_df <- merge_df[,-c(4,9)]
  
  colnames(merge_df) <- c("associated_gene", species1, paste0(species1,"_Detected_Num_Isoseq"), paste0(species1,"_Known_Num_Gencode"), paste0(species1,"_Normalised"), species2, paste0(species2,"_Detected_Num_Isoseq"), paste0(species2,"_Known_Num_Gencode"), paste0(species2,"_Normalised"), "Detected")
  
  
  assign(output_file, merge_df, envir=.GlobalEnv)
  assign(paste0(output_file,"_hist"), df3, envir=.GlobalEnv)
}



