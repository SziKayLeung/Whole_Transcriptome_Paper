# Szi Kay Leung
# 26/07/2019: readapted Aaron's script from Isoseq_Vs_RNASeq_Abundance.R 
# 08/08/2019: modified scirpt to allow processing directly from linux with only input directory/sample_name
# 08/08/2019: converted SQANTI2 FL input to TPM for correlation (and code for normalising to length) and removed use of TOFU
# 26/11/2019: Merge_Transcript_Kallisto_Input and Merge_Gene_Kallisto_Input function

suppressMessages(library(dplyr))# dataframe splitting
#suppressMessages(library(ggpubr)) # correlation plot
suppressMessages(library(ggplot2))
suppressMessages(library(grid)) # ggplot add corr value
suppressMessages(library(gridExtra)) # ggplot add corr value
#suppressMessages(library(ggthemes)) # ggplot minimum theme
suppressMessages(library(tidyr)) # ggplot minimum theme
suppressMessages(library(extrafont))

#********************** Function for input and plot
### Merge_Transcript_Input: Function for correlating RNASeq expression matrix from featurecount against IsoSeq from SQANTI 
# Input Sqanti classification file, tofu's abundance filtered file, and TPM calculated using script processfeaturecounts.R
Merge_Transcript_Input <- function(sqanti_dir, featurecount_dir,prefix_postisoseq3_name,prefix_ftrcnt_name,sample_name,rnaseq_threshold){
  sqanti <- read.table(print(paste0(sqanti_dir,"/",prefix_postisoseq3_name,".collapsed.filtered.rep_classification.filtered_lite_classification.txt")),header=T, stringsAsFactors=F, sep="\t")
  All_TPM <- read.table(print(paste0(featurecount_dir,"/",prefix_ftrcnt_name,"_tpm.txt")),header=T, stringsAsFactors=F, sep="\t")
  
  sample_TPM <- data.frame(All_TPM[, grep(sample_name, names(All_TPM))]) 
  TPM <- cbind(All_TPM$Geneid,sample_TPM)
  colnames(TPM) <- c("Geneid","RNASeq_TPM")

  # convert SQANTI FL to TPM (based on E.Tseng's SQANTI2.report https://github.com/Magdoll/SQANTI2)
  total_fl <- sum(sqanti$FL, na.rm=T)
  sqanti$ISOSEQ_TPM <- sqanti$FL*(10**6)/total_fl
  ## Alternatively to normalise SQANTI FL and to lenghts
  #tpm <- function(counts, lengths) {
  #rate <- counts / lengths
  #rate / sum(rate) * 1e6
  #}
  #sqanti$ISOSEQ_TPM_Normalised <- tpm(sqanti$FL,sqanti$length)

  # merge and format
  FSM <- subset(sqanti, sqanti$structural_category=="full-splice_match")
  FSM_TPM <- merge(FSM,TPM,by.x="associated_transcript", by.y="Geneid")
  
  # to avoid log output of "-inf" for 0 reads, set 0 reads in RNASeq expression data as NA 
  FSM_TPM[FSM_TPM$RNASeq_TPM <= rnaseq_threshold, "RNASeq_TPM"] <- NA
  FSM_TPM$RNASeq_TPM <- as.numeric(as.character(FSM_TPM$RNASeq_TPM ))
  
  # log counts
  FSM_TPM$log_ISOSEQ_TPM <- log10(FSM_TPM$ISOSEQ_TPM)
  FSM_TPM$log_RNASEQ_TPM <- log10(FSM_TPM$RNASeq_TPM)
  FSM_TPM <<- FSM_TPM
  
  # FSM, only multi-exonic
  FSM_TPM_multi <- FSM_TPM[FSM_TPM$subcategory=="multi-exon",]
  FSM_TPM_multi <<- FSM_TPM_multi  

  # FSM, only specific monoexonic within cage peak
  FSM_TPM_specific_mono <- FSM_TPM[!(FSM_TPM$subcategory =="mono-exon" & FSM_TPM$within_cage_peak == "False"),]
  FSM_TPM_specific_mono <<- FSM_TPM_specific_mono
}

# Merge_Transcript_Kallisto_Input: Function for correlating RNASeq expression and Isoseq expression in SQANTI output file
Merge_Transcript_Kallisto_Input <- function(input_sqanti, rnaseq_threshold){
  
  # sqanti <- read.table(input_sqanti_file, header=T, as.is=T, sep="\t")
  # convert SQANTI FL to TPM (based on E.Tseng's SQANTI2.report https://github.com/Magdoll/SQANTI2)
  # recalculate by removing all monoexonic counts even within total_fl 
  sqanti <- input_sqanti %>% filter(subcategory != "mono-exon")
  total_fl <- sum(sqanti$FL, na.rm=T)
  sqanti$ISOSEQ_TPM <- sqanti$FL*(10**6)/total_fl
  
  # log counts
  sqanti$log_ISOSEQ_TPM <- log10(sqanti$ISOSEQ_TPM)
  # to avoid log output of "-inf" for 0 reads, set 0 reads in RNASeq expression data as NA 
  sqanti[sqanti$iso_exp <= rnaseq_threshold, "iso_exp"] <- NA
  sqanti$iso_exp <- as.numeric(as.character(sqanti$iso_exp))
  sqanti$log_RNASEQ_TPM <- log10(sqanti$iso_exp)
  
  # remove rows with NA as should only matching genes with expression in both RNASeq and IsoSeq
  sqanti <- sqanti[!is.na(sqanti$iso_exp),]
  cat("Number of Transcripts used in correlation:", nrow(sqanti))
  
  return(sqanti)
}

# Merge_Gene_Kallisto_Input: gene expression correlation
Merge_Gene_Kallisto_Input <- function(sqanti, rnaseq_threshold, isoseq_cateogory){
  #sqanti <- read.table(input_sqanti_file, header=T, as.is=T, sep="\t")
  
  # convert SQANTI FL to TPM (based on E.Tseng's SQANTI2.report https://github.com/Magdoll/SQANTI2)
  # recalculate by removing all monoexonic counts even within total_fl 
  sqanti_multiexonic <- sqanti %>% filter(subcategory != "mono-exon")
  total_fl <- sum(sqanti_multiexonic$FL, na.rm=T)
  
  # Gene Count from IsoSeq
  # sum the transcript FL/gene and calculate TPM
  if(isoseq_cateogory == "FSM"){
    iso_counts <- sqanti %>%
      filter(structural_category == "full-splice_match") %>%
      group_by(associated_gene) %>%
      summarise(ISOSEQ_GENE_FL = sum(FL),) %>%
      mutate(ISOSEQ_GENE_TPM = ISOSEQ_GENE_FL*(10**6)/total_fl)
  }else{
    iso_counts <- sqanti %>%
      filter(subcategory != "mono-exon") %>%
      group_by(associated_gene) %>%
      summarise(ISOSEQ_GENE_FL = sum(FL),) %>%
      mutate(ISOSEQ_GENE_TPM = ISOSEQ_GENE_FL*(10**6)/total_fl)
  }
  
  
  # Gene Count from RNASeq
  # Note: Kallisto input, gene expression is the same for all the transcripts from the same gene
  rna_counts <- unique(sqanti_multiexonic[,c("associated_gene","gene_exp")])
  
  # Merge Counts at Gene level
  gene_counts <- merge(iso_counts, rna_counts, by = "associated_gene", all = TRUE) 
  # to avoid log output of "-inf" for 0 reads, set 0 reads in RNASeq expression data as NA, or set threshold 
  gene_counts[gene_counts$gene_exp <= rnaseq_threshold, "gene_exp"] <- NA
  gene_counts$gene_exp <- as.numeric(as.character(gene_counts$gene_exp))
  
  # remove rows with NA as should only matching genes with expression in both RNASeq and IsoSeq
  gene_counts <- gene_counts[complete.cases(gene_counts), ]
  
  # log counts
  gene_counts$log_ISOSEQ_GENE_TPM <- log10(gene_counts$ISOSEQ_GENE_TPM)
  gene_counts$log_RNASEQ_GENE_TPM <- log10(gene_counts$gene_exp)
  
  cat("Number of Transcripts used in correlation:", nrow(gene_counts))
  return(gene_counts)
}

multiple_transcript_corr_plots <- function(dat){
  
  p.1 <- density_plot(dat$ALL,"log_ISOSEQ_TPM","log_RNASEQ_TPM","Log(Isoseq TPM)","Log(RNASeq TPM)","Transcript: All_noMonoExonic" )
  p.2 <- density_plot(dat$FSM,"log_ISOSEQ_TPM","log_RNASEQ_TPM","Log(Isoseq TPM)","Log(RNASeq TPM)","Transcript: FSM" )
  p.3 <- density_plot(dat$FSM_Multiexonic,"log_ISOSEQ_TPM","log_RNASEQ_TPM","Log(Isoseq TPM)","Log(RNASeq TPM)","Transcript: FSM,Multi-exonic only" )
  p.4 <- density_plot(dat$FSM_Multiexonic_CAGE,"log_ISOSEQ_TPM","log_RNASEQ_TPM","Log(Isoseq TPM)","Log(RNASeq TPM)","Transcript: FSM,Multi-exonic, CAGE only" )
  
  print(p.1)
  print(p.2)
  print(p.3)
  print(p.4)
}

