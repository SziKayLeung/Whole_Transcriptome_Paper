# Szi Kay Leung
# Script to functions for tables for IsoSeq Paper (2020)


suppressMessages(library(enrichR))
suppressMessages(library(dplyr))
suppressMessages(library(xlsx))

### Source additional scripts #############################################################################
# Read in input files 
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Input_Variables.R")


output_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables/"


### General ##############################################################
# summary_info <read.clasification.file> 
# Aim: summarise the number of isoforms, minimum number of exons and max etc per gene 
summary_info <- function(dat){
  
  total_fl <- sum(dat$FL, na.rm=T)
  
  info <- list(
    # Number of isoforms 
    dat %>% count(associated_gene, name = "Num_of_Isoforms"),
    # Min Exons 
    dat %>% group_by(associated_gene) %>% summarise(Min_exon = min(exons)),
    # Max Exons 
    dat %>% group_by(associated_gene) %>% summarise(Max_exon = max(exons)),
    # Min Length 
    dat %>% group_by(associated_gene) %>% summarise(Min_length = min(length)),
    # Max Length
    dat %>% group_by(associated_gene) %>% summarise(Max_length = max(length)), 
    # Num of FL reads 
    dat %>% group_by(associated_gene) %>% summarise(Total_Reads = sum(FL))
  )
  
  final <- Reduce(function(...) merge(..., by='associated_gene', all.x=TRUE), info) %>%
    filter() 
  
  # TPM by gene (sum of FL reads)
  final <- final %>% mutate("FL_TPM" = round(Total_Reads*(10**6)/total_fl)) %>%
    mutate("Log10_FL_TPM" = log10(FL_TPM))
  
  return(final)
}

# obtain genome coordinates for isoforms in sqanti classification file 
gtf_obtain <- function(){
  gtf.class.files <- list()
  for(i in 1:length(class.files)){
    gtf.class.files[[i]] <- read_merge_gtf(sqanti.gtf.names.files[[i]], class.files[[i]])
  }
  names(gtf.class.files) <- names(class.files)
}

### GO ##############################################################
## Top 500 genes ranked by total reads
top500Genes_enrichR <- function(dataset, enrichr_category,species){
  
  # run the summary per gene on dataset
  info <- summary_info(dataset)
  # order by reads 
  info_exp <- info[order(info$Total_Reads, decreasing = T),]
  top500 <- info_exp$associated_gene[1:500]
  # run enrichR
  enriched <- enrichr(top500, enrichr_category)
  write.csv(enriched[1], paste0(output_table_dir,"/GO/",species,"_500topGenes.csv"))
}

# Top 100 Genes ranked by isoform diversity
top100Genes_enrichR <- function(dataset, enrichr_category,species,title){
  
  # run the summary per gene on dataset
  info <- summary_info(dataset)
  # order by reads 
  info_num <- info[order(info$Num_of_Isoforms, decreasing = T),]
  top100 <- info_num$associated_gene[1:100]
  # run enrichR
  enriched <- enrichr(top100, enrichr_category)
  write.csv(enriched[1], paste0(output_table_dir,"/GO/",title,"_",species,"_100topGenes.csv"))
}


IR_genes_enrichR <- function(dat,enrichr_category,species){
  # remove first row with X
  IR <- dat[dat$subcategory == "intron_retention",]
  IR_gene <- unique(IR$associated_gene)
  # run enrichR
  enriched <- enrichr(IR_gene, enrichr_category)
  write.csv(enriched[1], paste0(output_table_dir,"/GO/",species,"_IRGenes.csv"))
}

### Numbers and Expression Analysis #################################################################

# annotated_genes_transcript_num
# Aim: Extract the number of novel and known transcripts from annotated genes in each dataset 
# Input: classification files (SQANTI2) filtered 
# Output: txt. file 
# caveat: removal of monoexonic transcripts
annotated_genes_transcript_num <- function(){
  # subset classification files (filtered) by annotated genes, and novel transcripts 
  annotated.class.files <- lapply(class.files, function(x) x[!grepl("NOVELGENE",x$associated_gene),])
  annotated.class.files.novel.transcripts <- lapply(annotated.class.files, function(x) x[grepl("novel",x$associated_transcript),] %>% filter(subcategory != "mono-exon"))
  annotated.class.files.annotated.transcripts <- lapply(annotated.class.files, function(x) x[!grepl("novel",x$associated_transcript),] %>% filter(subcategory != "mono-exon"))
  
  # function to tally the number of transcripts per gene 
  transcripts_tally <- function(dat){
    dat1 <- dat %>% bind_rows() %>% 
      group_by(associated_gene, Sample) %>% 
      tally() %>% 
      mutate(ID = paste(associated_gene, Sample, sep = "_"))
    
    return(dat1)
  }
  
  # apply function to novel transcripts and annotated transcripts subsetted from known genes
  novel_tally <- transcripts_tally(annotated.class.files.novel.transcripts)
  known_tally <- transcripts_tally(annotated.class.files.annotated.transcripts)
  
  # Merge and datawrangle 
  full <- merge(novel_tally,known_tally, by = "ID", all = TRUE) %>% 
    select(ID, associated_gene.x, Sample.x, novel_transcripts = n.x, associated_gene.y, Sample.y, annotated_transcripts = n.y) %>%
    mutate(Sample = ifelse(!is.na(Sample.x),Sample.x,Sample.y)) %>%
    mutate(associated_gene = ifelse(!is.na(associated_gene.x),associated_gene.x,associated_gene.y)) %>% 
    rowwise() %>% mutate(total = sum(novel_transcripts, annotated_transcripts, na.rm = TRUE)) %>%
    .[,c("Sample","associated_gene","annotated_transcripts","novel_transcripts","total")]
  
  Mouse <- full[full$Sample == "Mouse",]
  Human <- full[full$Sample == "Human",]
  Adult <- full[full$Sample == "Human (Adult)",]
  Fetal <- full[full$Sample == "Human (Fetal)",]
  
  output <- list(Mouse,Human,Adult,Fetal)
  names(output) <- c("Mouse","Human","Adult","Fetal")
  return(output)
}

# threshold_gene_expression_human_mouse_common_threshold
# Prerequisite: files generated from Human_Mouse_Comparison.Rmd
# both genes in human and mouse expression have to be above the threshold
threshold_gene_expression_human_mouse_common_threshold <- function(){
  # read in gene expression (output from Human_Mouse_Comparison.Rmd)
  gene_expression <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/human_mouse_multiexonic_TPM_homologs_fullymatched.csv")) 
  num_isoform <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/human_mouse_multiexonic_NUM_homology_fullymatched.csv"))
  
  # determine threshold of IsoSeq TPM Gene Expression (Log10TPM)
  threshold <- seq(0,4,0.5) 
  
  # loop through the threshold value for genes with Log10 TPM gene expression > threshold value 
  dat <- data.frame()
  for(i in 1:length(threshold)){
    # filter genes with human and mouse threshold values > 10 TPM
    cutoff <- gene_expression %>% filter(Human_Detected_LogTPM_Isoseq > threshold[i] & 
                                           Mouse_Detected_LogTPM_Isoseq > threshold[i]) 
    # only include genes passing through threshold
    cut_off_num_isoform <- num_isoform[num_isoform$associated_gene %in% cutoff$associated_gene,] 
    corr_value <- ifelse(length(cut_off_num_isoform$Human_Detected_Num_Isoseq) > 5 & 
                           length(cut_off_num_isoform$Mouse_Detected_Num_Isoseq) > 5, 
                         cor(cut_off_num_isoform$Human_Detected_Num_Isoseq, cut_off_num_isoform$Mouse_Detected_Num_Isoseq, 
                             use = "pairwise.complete.obs"), 
                         "NA")
    
    # p value only determined if more than 5 values
    p.value <- ifelse(length(cut_off_num_isoform$Human_Detected_Num_Isoseq) > 5 & 
                        length(cut_off_num_isoform$Mouse_Detected_Num_Isoseq) > 5, 
                      cor.test(cut_off_num_isoform$Human_Detected_Num_Isoseq, cut_off_num_isoform$Mouse_Detected_Num_Isoseq, 
                               use = "pairwise.complete.obs")$p.value, 
                      "NA")
    
    dat[i,1] <- threshold[i]
    dat[i,2] <- round(corr_value,3)
    dat[i,3] <- p.value
    dat[i,4] <- length(cut_off_num_isoform$Human_Detected_Num_Isoseq) # Note matched number of mouse genes
    
    colnames(dat) <- c("Threshold_IsoSeq_GeneExpression_LogTPM","Correlation_Value", "p_value", "Number of common genes")
  }
  return(dat)
}

# threhsold_gene_expression_human_mouse
# Prerequisite: files generated from Human_Mouse_Comparison.Rmd
# either genes in human and mouse expression have to be above the threshold
threshold_gene_expression_human_mouse <- function(){
  # read in gene expression (output from Human_Mouse_Comparison.Rmd)
  gene_expression <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/human_mouse_multiexonic_TPM_homologs_fullymatched.csv")) 
  num_isoform <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/human_mouse_multiexonic_NUM_homology_fullymatched.csv"))
  
  # determine threshold of IsoSeq TPM Gene Expression (Log10TPM)
  threshold <- seq(0,4,0.5) 
  
  # loop through the threshold value for genes with Log10 TPM gene expression > threshold value 
  dat <- data.frame()
  for(i in 1:length(threshold)){
    # filter genes with human and mouse threshold values > 10 TPM
    human_cutoff <- gene_expression %>% filter(Human_Detected_LogTPM_Isoseq > threshold[i]) 
    mouse_cutoff <- gene_expression %>% filter(Mouse_Detected_LogTPM_Isoseq > threshold[i]) 
    
    # only include genes passing through threshold
    cut_off_num_isoform <- bind_rows(num_isoform[num_isoform$associated_gene %in% human_cutoff$associated_gene,],
                                     num_isoform[num_isoform$associated_gene %in% mouse_cutoff$associated_gene,])
    corr_value <- cor(cut_off_num_isoform$Human_Detected_Num_Isoseq, cut_off_num_isoform$Mouse_Detected_Num_Isoseq, 
                      use = "pairwise.complete.obs")
    
    p.value <- ifelse(length(cut_off_num_isoform$Human_Detected_Num_Isoseq) > 5 & 
                        length(cut_off_num_isoform$Mouse_Detected_Num_Isoseq) > 5, 
                      cor.test(cut_off_num_isoform$Human_Detected_Num_Isoseq, 
                               cut_off_num_isoform$Mouse_Detected_Num_Isoseq, 
                               use = "pairwise.complete.obs")$p.value, 
                      "NA")
    
    dat[i,1] <- threshold[i]
    dat[i,2] <- round(corr_value,3)
    dat[i,3] <- p.value
    dat[i,4] <- length(cut_off_num_isoform$Human_Detected_Num_Isoseq) # Note matched number of mouse genes
    
    colnames(dat) <- c("Threshold_IsoSeq_GeneExpression_LogTPM","Correlation_Value", "p_value", "Number of common genes")
  }
  return(dat)
}


### Fusion Genes #################################################################
fusion_genes_findmorethan2 <- function(){
  
  perdataset <- function(pre_or_post_sqanti){
    fusion_genes <- bind_rows(pre_or_post_sqanti) %>% 
      filter(structural_category == "Fusion") %>% 
      mutate(num_fusion_genes = count.fields(textConnection(associated_gene), sep = "_")) 
    
    stats <- fusion_genes %>%
      group_by(Sample,num_fusion_genes) %>% 
      tally() %>% 
      arrange(num_fusion_genes)
    
    print(stats)
    
    print(fusion_genes[fusion_genes$num_fusion_genes == "3",])
  }
  
  print("Working with SQANTI filtered files")
  perdataset(class.files)
  
  print("Working with SQANTI prefiltered files")
  perdataset(pre.class.files)
}

conjoin_prediction <- function(input_table, sample){
  
  # subset SQANTI2 classification file 
  Fusion_class <- class.files[[sample]][class.files[[sample]]$structural_category == "Fusion",] %>% 
    mutate(Fusion_Gene_1 = word(associated_gene,c(1),  sep = fixed ('_'))) %>% 
    mutate(Fusion_Gene_2 = word(associated_gene,c(2),  sep = fixed ('_')))
  
  
  
  # read input table and fuse Gene 1 and Gene, both ways to check against fusion genes in SQANTI2 classification file
  if(sample == "Human"){
    fusion_conjoin <- read.csv(input_table) %>% 
      mutate(Fused_Gene_1 = paste0(Gene.symbol..Gene.1.,"_",Gene.symbol..Gene.2.)) %>% 
      mutate(Fused_Gene_2 = paste0(Gene.symbol..Gene.2.,"_",Gene.symbol..Gene.1.))
    
    # Double check not misisng any fusion genes 
    print("grep successful of genes in fusion and from supplementary table")
    for(i in fusion_conjoin$Gene.symbol..Gene.1.){
      if(nrow(Fusion_class[grepl(i,Fusion_class$associated_gene) & Fusion_class$structural_category == "Fusion",]) > 0){print(i)}
    }
    # Note:MAD is from SMAD4_ELAC1 which is missed on the supplementary table but not in the database
    for(i in fusion_conjoin$Gene.symbol..Gene.2.){
      if(nrow(Fusion_class[grepl(i,Fusion_class$associated_gene) & Fusion_class$structural_category == "Fusion",]) > 0){print(i)}
    }
    
  }else if(sample == "Mouse"){
    fusion_conjoin <- read.csv(input_table) %>% 
      mutate(Fused_Gene_1 = paste0(X5..Gene.Symbol,"_",X3..Gene.Symbol)) %>% 
      mutate(Fused_Gene_2 = paste0(X3..Gene.Symbol,"_",X5..Gene.Symbol)) 
    fusion_conjoin$Fused_Gene_1 <- toupper(fusion_conjoin$Fused_Gene_1)
    fusion_conjoin$Fused_Gene_2 <- toupper(fusion_conjoin$Fused_Gene_2)
    
    # Double check not misisng any fusion genes 
    print("grep successful of genes in fusion and from supplementary table")
    for(i in toupper(fusion_conjoin$X5..Gene.Symbol)){
      if(nrow(Fusion_class[grepl(i,Fusion_class$associated_gene) & Fusion_class$structural_category == "Fusion",]) > 0){print(i)}
    }
    # Note:number of genes that are detected in database as fusion but not necessarily the same 2 genes as fusion 
    # SERPINA3G, NXF1, RIOK2, KCNK4, 1600002K03RIK, RHBDL1, GP1BB
    for(i in toupper(fusion_conjoin$X3..Gene.Symbol)){
      if(nrow(Fusion_class[grepl(i,Fusion_class$associated_gene) & Fusion_class$structural_category == "Fusion",]) > 0){print(i)}
    }
  }
  
  
  Fusion_class$fusion_conjoin <- ifelse(Fusion_class$associated_gene %in% fusion_conjoin$Fused_Gene_1 | Fusion_class$associated_gene %in% fusion_conjoin$Fused_Gene_2,
                                        "Yes","No")
  
  
  return(Fusion_class[Fusion_class$fusion_conjoin == "Yes",])
  
  
}

### IR #################################################################
IR_genes_only_IRtranscripts <- function(){
  IR_genes_only_IRtranscripts_dataset <- function(input.class.files){
    Total <- input.class.files %>% group_by(associated_gene) %>% tally() %>% dplyr::rename(., total_num_transcripts = n)
    IR <- input.class.files %>% filter(subcategory == "intron_retention") %>% group_by(associated_gene) %>% tally() %>% 
      dplyr::rename(.,IR_Transcripts = n) 
    
    perc_IR_Transcripts <- merge(IR, Total, by = "associated_gene", all.x = T) %>% 
      mutate(perc_IR_trascripts = round(IR_Transcripts/ total_num_transcripts * 100,2)) 
    
    num_genes_only_IR <- perc_IR_Transcripts %>% filter(perc_IR_trascripts == 100) 
    return(list(nrow(Total), nrow(IR), nrow(num_genes_only_IR), as.data.frame(num_genes_only_IR)))
  }
  
  dat <- data.frame()
  for(i in 1:length(class.files)){
    dat[1,i] <- IR_genes_only_IRtranscripts_dataset(class.files[[i]])[[1]]
    dat[2,i] <- IR_genes_only_IRtranscripts_dataset(class.files[[i]])[[2]]
    dat[3,i] <- IR_genes_only_IRtranscripts_dataset(class.files[[i]])[[3]]
    dat[4,i] <- round(dat[3,i]/dat[1,i] * 100,2)
    dat[5,i] <- round(dat[3,i]/dat[2,i] * 100,2)
    
    colnames(dat)[i] <- names(class.files)[i] 
    # write.csv 
    IR_genes_only_IRtranscripts_output <- IR_genes_only_IRtranscripts_dataset(class.files[[i]])[[4]]
    write.table(IR_genes_only_IRtranscripts_output[-1,1], paste0(output_table_dir, "AS_IR/", names(class.files)[[i]],
                                                                 "_IR_genes_only_IR_GOList.txt"), row.names=FALSE,quote = FALSE)
    write.csv(IR_genes_only_IRtranscripts_output, paste0(output_table_dir, "AS_IR/", names(class.files)[[i]],
                                                         "_IR_genes_only_IR_Transcripts.csv"))
  }
  rownames(dat) <- c("Number of Total Genes","Number of Genes with IR transcripts", "Number of Genes with only IR transcripts", "Percentage of Genes with only IR transcripts out of all Total Genes", "Percentage of Genes with only IR transcripts out of IR Genes")
  
  return(dat)
}


### NovelGenes #################################################################
novelgene_blastflnc <- function(blast_input){
  dat <- read.table(blast_input, as.is=T, sep="\t", header = F)
  colnames(dat) <- c("CCS","PacBio_ID","Perc_Identity","alignment_length","mismatches", "gap","q.start","q.end","s.start","s.end","evalue","bit_score")
  
  ## Filtering requirements
  # > 200bp length 
  # > 90% blast identity 
  
  # filter blast hits > 200bp 
  filtered <- dat[dat$alignment_length > 200 & dat$Perc_Identity > 90, ]
  novelgene_hit <- filtered %>% group_by(PacBio_ID) %>% tally()
  
  output <- list(dat,filtered,novelgene_hit)
  names(output) <- c("blast","filteredblast","novelgene")
  print(output$novelgene)
  return(output)
}

# prepare_blast
# Aim: to read in output from blast and relabel headings, and merge output with description from sqanti classification file 
# Important to get gtf.class.filesfirst 
# Assumption: # Query ID => Human, Subject ID => Mouse, as blasting Human to mouse novel genes 
# Input (and Output): Blast output table 
blast <- function(blast_input, type){
  gtf_obtain()
  dat <- read.table(blast_input, as.is=T, sep="\t", header = F)
  
  ## Filtering requirements
  # > 200bp length 
  # > 90% blast identity 
  
  if(type == "Human2Mouse"){
    colnames(dat) <- c("Human_PacBio_ID","Mouse_PacBio_ID","Perc_Identity","alignment_length","mismatches", "gap","q.start","q.end","s.start","s.end","evalue","bit_score")
    unfiltered_hitcounts <- dat %>% group_by(Human_PacBio_ID) %>% tally()
    dat <- dat[dat$alignment_length > 200 & dat$Perc_Identity > 90, ]
    hitcounts <- dat %>% group_by(Human_PacBio_ID) %>% tally()
  } else {
    colnames(dat) <- c("Mouse_PacBio_ID","Human_PacBio_ID","Perc_Identity","alignment_length","mismatches", "gap","q.start","q.end","s.start","s.end","evalue","bit_score")
    unfiltered_hitcounts <- dat %>% group_by(Mouse_PacBio_ID) %>% tally()
    dat <- dat[dat$alignment_length > 200 & dat$Perc_Identity > 90, ]
    hitcounts <- dat %>% group_by(Mouse_PacBio_ID) %>% tally()
  }
  
  # blast output with description from sqanti clasification file (human)
  final_human <- merge(dat, gtf.class.files$Human, 
                       by.x = "Human_PacBio_ID", by.y = "isoform", all.x = T) %>% 
    select(-"Mouse_PacBio_ID")
  
  # blast output with description from sqanti clasification file (mouse)
  final_mouse <- merge(dat, gtf.class.files$Mouse, 
                       by.x = "Mouse_PacBio_ID", by.y = "isoform", all.x = T) %>% 
    select(-"Human_PacBio_ID")
  
  print(unfiltered_hitcounts)
  print(hitcounts)
  # return as list
  output <- list(dat, unfiltered_hitcounts,hitcounts, final_human, final_mouse)
  names(output) <- c("Novel_Genes_Blast_output","unfiltered_hitcounts","Hits","Human", "Mouse")
  return(output)
}



# Prerequisite: Files generated from Characterising_Novel_Genes.sh
Blast_Novel_Genes_Across_Genome <- function(){
  # dir 
  AcrossGenome <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Revised_Paper/Novel_Genes/AcrossGenome/"
  
  Human_Novel_Human <- read.table(paste0(AcrossGenome,"Human_Novel_Human_Genome_blast.txt"))
  Human_Novel_gtf <- read.csv(paste0(output_table_dir,"/Novel_Genes/Novel_Genes.csv"))[,-1] %>% filter(Sample == "human") 
  Mouse_Novel_Mouse <- read.table(paste0(AcrossGenome,"Mouse_Novel_Mouse_Genome_blast.txt"))
  Mouse_Novel_gtf <- read.csv(paste0(output_table_dir,"/Novel_Genes/Novel_Genes.csv"))[,-1] %>% filter(Sample == "mouse") 
  
  #QC: Check that the blast number of genes is equivalent to the number of novel genes from sqanti
  #setdiff(unique(Human_Novel_gtf$isoform), unique(Human_Novel_Human$V1))
  #length(unique(Human_Novel_Human$V1))
  #setdiff(unique(Mouse_Novel_gtf$isoform), unique(Mouse_Novel_Mouse$V1))
  #length(unique(Mouse_Novel_Mouse$V1))
  
  curated_list <- function(input_blast_dat, input_gtf){
    
    ## Filtering requirements
    # > 500bp length 
    # > 90% blast identity 
    
    # filter blast hits > 500bp 
    input_blast_dat_500 <- input_blast_dat[input_blast_dat$V4 > 500, ] 
    # restrict to length > 500 bp; different chromosomes
    #hist(input_blast_dat$V4)
    #hist(input_blast_dat$V3)
    
    # merge with input gtf based on isoform to get original gtf coordinates 
    input_blast_dat_500 <- merge(input_blast_dat_500, input_gtf[,c("isoform","gtf_coordinates")], by.x = "V1", by.y = "isoform")
    
    # note coordinates are sometimes in the wrong order for start and end 
    blast_directionality <- function(blast){
      # created new column of START and END 
      blast$START <- ifelse(blast$V9 < blast$V10, blast$V9, blast$V10) # for START column: if START coordinates < END coordinates, keep START, otherwise switch END
      blast$END <- ifelse(blast$V9 < blast$V10, blast$V10, blast$V9) # for END column: if START coordinates < END coordinates, keep END, otherwise switch START
      
      return(blast)
    }
    
    corrected_input_blast_500 <- blast_directionality(input_blast_dat_500)
    # Further filtering to only extract blast hits that are not of original sequence using gtf coordinates 
    corrected_input_blast_500  <- corrected_input_blast_500 %>%
      mutate(chr_gtf = word(corrected_input_blast_500$gtf_coordinates,c(1),  sep = fixed (' '))) %>%
      mutate(start_gtf = word(corrected_input_blast_500$gtf_coordinates,c(3),  sep = fixed (' '))) %>%
      mutate(end_gtf = word(corrected_input_blast_500$gtf_coordinates,c(5),  sep = fixed (' '))) %>% 
      # If blast start coordinate is bigger than start coordinate and less than end coordinate, put "Yes" 
      # If blast end coordinate is bigger than start coordinate and less than end coordinate, put "Yes"
      # as this means that the sequence of interest is within the same coordinate range as the original sequence
      # START is start sequence, END is end sequence 
      mutate(blat_start_classifer = ifelse(.$START >= .$start_gtf & .$START <=.$end_gtf, "Yes","No")) %>% 
      mutate(blat_end_classifer = ifelse(.$END >= .$start_gtf & .$END <=.$end_gtf, "Yes","No")) %>%
      mutate(blat_start_end_classifier = ifelse(.$blat_start_classifer == .$blat_end_classifer, "Identical","Different")) %>% 
      mutate(blat_gtf_coordinate = paste0(V2,":",START,"-",END)) %>% 
      mutate(blat_identified = ifelse(.$blat_gtf_coordinate == .$gtf_coordinates, "Identified","Not_Identified")) %>%
      # filter for sequences with coordinates outside of the same coordinate range as the original sequence 
      filter(blat_end_classifer == "No" & blat_start_classifer == "No") %>% 
      mutate(blat_coordinates = paste(.$V2,":",.$START, "-", .$END)) %>% 
      # further filter for blast sequence identity > 90% 
      filter(V3 > 90)  
    # filter for different chromosomes 
    #mutate(chr_classifier = ifelse(.$V2 == .$chr_gtf, "Identical","Different")) %>% 
    #filter(chr_classifier == "Different")
    
    # curated with the tally of the number of blast hits for that isoform 
    curated <- corrected_input_blast_500 %>% group_by(V1) %>% tally() %>% 
      `colnames<-`(c("isoform", "num_blast_hits")) %>%
      arrange(desc(num_blast_hits)) %>%
      left_join(., input_gtf[,c("isoform","gtf_coordinates", "associated_gene","cagepeak")], by = "isoform")
    
    colnames(corrected_input_blast_500)[1:12] <- c("isoform","chr","%_Identity","alignment_length","mismatches","gap_openings","q.start","q.end","s.start","s.end","e_value","bit_score")
    output <- list(corrected_input_blast_500,curated)
    names(output) <- c("input_blast_dat_500","tally_curated_list")
    return(output)
  }
  
  # In Mouse: any hits with less than 5 is from just spurious hits that actually refer to other transcripts of the same novel gene 
  # however visually checked human and some of the hits are of different novel gene
  # NOTE: Inflation of blast hits
  
  
  Human <- curated_list(Human_Novel_Human, Human_Novel_gtf)
  Mouse <- curated_list(Mouse_Novel_Mouse, Mouse_Novel_gtf)
  
  write.csv(Human$input_blast_dat_500, paste0(output_table_dir,"/Novel_Genes/Human_filtered_novelblastgenome.csv"), quote = F, row.names = F)
  write.csv(Mouse$input_blast_dat_500, paste0(output_table_dir,"/Novel_Genes/Mouse_filtered_novelblastgenome.csv"), quote = F, row.names = F)
  write.csv(Human$tally_curated_list, paste0(output_table_dir,"/Novel_Genes/Human_filtered_novelblastgenome_number.csv"), quote = F, row.names = F)
  write.csv(Mouse$tally_curated_list, paste0(output_table_dir,"/Novel_Genes/Mouse_filtered_novelblastgenome_number.csv"), quote = F, row.names = F)
  
}



