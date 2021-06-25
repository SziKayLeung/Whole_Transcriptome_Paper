# global parameters
# need disease_list.csv for AD genes, SZ genes 
disease_list <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables/Disease/DiseaseList.csv")
fusion_genes <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables/Fusion_Genes/Fusion_Genes.csv")

# original files for disease 
# copied from Aaron's folder
#PGC <- read.table(paste0(diseasegenelists,"/SZ_PGC2.txt"), sep="\t", header=T, stringsAsFactors = F)
#CLOZ <- read.table(paste0(diseasegenelists,"/SZ_CLOZUK_GENE.txt"), sep="\t", header=T, stringsAsFactors = F)

# updated autism SFARI list 
diseasegenelists <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables/Disease/disease_gene_lists"
SFARI <- read.csv(paste0(diseasegenelists,"/SFARI-Gene_genes_07-29-2020release_08-07-2020export.csv"),header=T)
SFARI_CLASS_1_2_S <- subset(SFARI$gene.symbol, SFARI$gene.score == 1 | SFARI$gene.score == 2 |SFARI$syndromic == 1)
SFARI_CLASS_1_2 <- subset(SFARI$gene.symbol, SFARI$gene.score == 1 | SFARI$gene.score == 2)

# 1. Datawrangle the disease list  
# aim: remove white space or empty records from each column of the DiseaseList.csv and save genes as new list 
# output: rehash_disease_list  
rehash_disease_list <- list(
  AD = disease_list$AD[disease_list$AD != ""], 
  SZ = unique(as.factor(c(as.character(disease_list$CLOZ[disease_list$CLOZ != ""]),
                          as.character(disease_list$PGC_SZ[disease_list$PGC_SZ != ""])))),
  Autism = SFARI_CLASS_1_2
)

# 2. Number of transcripts and genes in disease
disease_all <- function(input.class.files){
  disease <- unique(unlist(rehash_disease_list))
  alldisease_class.files <- input.class.files[input.class.files$associated_gene %in% disease,]
  
  totaltranscripts <- nrow(alldisease_class.files)
  totalgenes <- length(unique(alldisease_class.files$associated_gene))
  cat("Number of transcripts to disease genes:",totaltranscripts,"\n")
  cat("Number of genes to disease genes:",totalgenes,"\n")
  
  # more than one isoform
  morethan1 <- alldisease_class.files %>% group_by(associated_gene) %>% tally() %>% filter(n > 1) %>% nrow()
  cat("Number of disease genes with more than isoform:",morethan1,"(",round(morethan1/totalgenes*100,2),"%)","\n")
  
  # novel transcripts 
  novel <- alldisease_class.files %>% filter(associated_transcript == "novel") %>% nrow()
  cat("Number of novel transcripts of disease genes:",novel,"(",round(novel/totaltranscripts*100,2),"%)","\n")
}



# 3. Subset classification file of just the genes interested for each disease 
# aim: Loop through the rehash_disease_list and extract from classification file the rows mapped to those genes 
# input: rehash_disease_list 
# output: disease_class.files 
subset_disease <- function(input.class.files){
  # use count to loop and save each diseaes gene list into new list entry 
  disease_class.files <- list()
  count=1
  for(i in rehash_disease_list){
    disease_class.files[[count]] <- input.class.files[input.class.files$associated_gene %in% i,]
    disease_class.files[[count]]$disease <- i[1]
    count = count + 1
  }
  names(disease_class.files) <- names(rehash_disease_list)
  return(disease_class.files)
}



# 4. General Basic Stats about genes associated with disease 
# aim: Extract the number of associated isoforms, novel and known and the differen SQANTI2 categories 
# input: disease_class.files 
# output: output basic stats
tabulating_sqanti_num <- function(type_class_file){
  dat <- data.frame()
  count=1
  for(i in type_class_file){
    # Total Unique Isoforms: tabulated by number of rows
    isoforms <- dim(i)[1]
    # Total Unique Genes: Remove novel genes and count 
    annotated_genes <- i[!grepl("NOVELGENE",i$associated_gene),] %>% count(associated_gene) %>% nrow(.)
    # Number of Annotated Isoforms (FSM, ISM) 
    annotated_isoforms <- paste0(nrow(i[i$associated_transcript != "novel",])," (",
                                 round(nrow(i[i$associated_transcript != "novel",])/isoforms * 100,2),"%)")
    novel_isoforms <- paste0(nrow(i[i$associated_transcript == "novel",])," (",
                             round(nrow(i[i$associated_transcript == "novel",])/isoforms * 100,2),"%)")
    
    # 9 levels of structural cateogory 
    struct <- vector("numeric", 9)
    for(num in 1:length(levels(i$structural_category))){
      struct[num] <- nrow(i[i$structural_category == levels(i$structural_category)[num],]) 
    }
    
    dat[1:4,count] <- rbind(annotated_genes, isoforms, annotated_isoforms, novel_isoforms)
    dat[5:13,count] <- struct
    colnames(dat)[count] <- names(type_class_file)[count]
    count = count + 1
  }
  row.names(dat) <- append(c("Total Number of Detected Genes", "Total Number of Isoforms", 
                             "Number and % of Annotated Isoforms", "Number and % of Novel Isoforms"),
                           levels(i$structural_category))
  
  dat
}



# 5. IR vs NMD vs NMD_IR
# aim: Tabulating for each gene in the disease list the number of IR, NMD and Fusion transcripts 
# input: rehash_disease_list 
# output: number of associated IR vs NMD transcripts for each gene: disease_num
num_disease <- function(disease_list, input.class.files, dataset){
  df <- input.class.files
  df_fusion <- fusion_genes[fusion_genes[[dataset]] != 0,]
  total_num_isoforms <- data.frame()
  total_num_IR <- data.frame()
  total_num_NMD <- data.frame()
  total_num_NMD_IR <- data.frame()
  total_num_fusion <- data.frame()
  for(i in 1:length(disease_list)){
    total_num_isoforms[i,1] <- disease_list[i]
    total_num_isoforms[i,2] <- nrow(df[df$associated_gene == disease_list[i],])
    
    # intron retention
    total_num_IR[i,1] <- disease_list[i]
    total_num_IR[i,2] <-  nrow(df[df$associated_gene == disease_list[i] & df$subcategory == "intron_retention",]) 
    
    # NMD 
    total_num_NMD[i,1] <- disease_list[i]
    total_num_NMD[i,2] <- nrow(df %>% filter(associated_gene == disease_list[i] & predicted_NMD == "TRUE"))
    
    
    # NMD and IR 
    total_num_NMD_IR[i,1] <- disease_list[i]
    total_num_NMD_IR[i,2] <- nrow(df %>% filter(associated_gene == disease_list[i] & predicted_NMD == "TRUE" & subcategory == "intron_retention"))
    
    # Fusion Genes 
    total_num_fusion[i,1] <- disease_list[i]
    total_num_fusion[i,2] <- ifelse(nrow(df_fusion[grep(disease_list[i],df_fusion$associated_gene),]) > 0, "YAY","NAY")
    
  }
  
  colnames(total_num_isoforms) <- c("associated_gene", "Total_Num_Isoforms")
  colnames(total_num_IR) <- c("associated_gene","Total_Num_IRTranscripts")
  colnames(total_num_NMD) <- c("associated_gene","Total_Num_NMDTranscripts")
  colnames(total_num_NMD_IR) <- c("associated_gene","Total_Num_NMDIRTranscripts")
  colnames(total_num_fusion) <- c("associated_gene","Total_Num_Fusiontranscripts")
  
  full_list <- list(total_num_isoforms, total_num_IR, total_num_NMD, total_num_NMD_IR, total_num_fusion)
  
  disease_num <- Reduce(function(...) merge(..., by='associated_gene', all.x=TRUE), full_list) 
  disease_num$rate_IR_NMD <- paste0(round(disease_num$Total_Num_NMDIRTranscripts/disease_num$Total_Num_IRTranscripts * 100, 2),"%")
  return(disease_num)
}

# 6. Tabulate IR vs NMD vs NMD_IR 
# aim: Tabulate the general stats about IR, NMD and fusion 
# input: output from function(4): disease_num 
# output: stats 
num_disease_stats <- function(input.class.files,num_disease_output_all){
  stats <- data.frame()
  for(i in 1:length(num_disease_output_all)){
    # remove duplicated rows 
    df <- num_disease_output_all[[i]] %>% distinct()
    # Total Number of genes 
    # Total Number of detected genes 
    
    
    total_num_IRgenes <- length(unique(input.class.files[input.class.files$subcategory == "intron_retention", "associated_gene"]))
    
    stats[1,i] <- nrow(df)
    stats[2,i] <- paste0(nrow(df[df$Total_Num_Isoforms !=0,]), "(",round(nrow(df[df$Total_Num_Isoforms !=0,])/nrow(df) * 100,2), "%)")
    stats[3,i] <- paste0(nrow(df[df$Total_Num_IRTranscripts != 0,]), "(", round(nrow(df[df$Total_Num_IRTranscripts != 0,])/total_num_IRgenes * 100,2), "%)")
    stats[4,i] <- paste0(nrow(df[df$Total_Num_IRTranscripts != 0,]), "(", round(nrow(df[df$Total_Num_IRTranscripts != 0,])/
                                                                                  nrow(df[df$Total_Num_Isoforms !=0,])* 100,2), "%)")
    stats[5,i] <- paste0(nrow(df[df$Total_Num_NMDTranscripts !=0,]), "(", 
                         round(nrow(df[df$Total_Num_NMDTranscripts !=0,])/ nrow(df[df$Total_Num_Isoforms !=0,]) * 100,2), "%)")
    stats[6,i] <- paste0(nrow(df[df$Total_Num_NMDIRTranscripts !=0,]), "(",
                         round(nrow(df[df$Total_Num_NMDIRTranscripts !=0,])/ nrow(df[df$Total_Num_Isoforms !=0,]) * 100,2), "%)")
    stats[7,i] <- paste0(nrow(df[df$Total_Num_Fusiontranscripts == "YAY",]))
    stats[8,i] <- paste0(nrow(df[df$Total_Num_Isoforms > 1,]), "(",
                         round(nrow(df[df$Total_Num_Isoforms > 1,])/nrow(df[df$Total_Num_Isoforms !=0,]) * 100,2),"%)")
    stats[9,i] <- paste0(min(df[df$Total_Num_Isoforms != 0,"Total_Num_Isoforms"]), "-", max(df$Total_Num_Isoforms))
    colnames(stats)[i] <- names(num_disease_output_all)[[i]]
    
  }
  rownames(stats) <- c("Total Genes", "Detected Genes", 
                       "IR Genes (% from all IR Genes)", 
                       "IR Genes (% from Detected Genes)",
                       "NMD Genes (% from Detected Genes)",
                       "IR and NMD genes (% from Detected Genes)",
                       "Fusion Genes",
                       "Number of Detected Genes with more than one isoform",
                       "Range of isoform diversity")
  print(stats)
  return(stats)
}

### Apply all functions into one 
run_disease_stats_all <- function(input.class.files, dataset){
  
  cat("Working with:", dataset,"\n")
  # 1. Datawrangle the disease list
  
  # 2. All disease
  disease_all(input.class.files)
  
  # 3. Subset classification file of just the genes interested for each disease 
  disease_class.files <- subset_disease(input.class.files)
  
  # 4. General Basic Stats about genes associated with disease 
  disease_stats <- tabulating_sqanti_num(disease_class.files)
  print(disease_stats)
  
  # 5. IR vs NMD vs NMD_IR
  num_disease_output_all <- lapply(rehash_disease_list, function(x) num_disease(x, input.class.files, dataset))
  num_disease_output_all <<- num_disease_output_all
  
  # 5. Tabulate IR vs NMD vs NMD_IR 
  num_disease_stat <- num_disease_stats(input.class.files, num_disease_output_all)
  
  all_stats <- bind_rows(num_disease_stat,disease_stats)
  
  # write.csv output from function 3 of the subset classification file 
  #write.csv(do.call(bind_rows, disease_class.files), paste0(output_table_dir, "/Disease/",dataset,"_disease_class.csv"))
  
  # write.csv output from function 5
  #write.csv(do.call(rbind, num_disease_output_all) %>% mutate(disease = word(rownames(.),c(1),  sep = fixed ('.'))) %>% distinct(), 
  #          paste0(output_table_dir, "/Disease/",dataset,"_disease_numbers.csv"))
  
  return(all_stats)
}

MouseStats <- run_disease_stats_all(class.files$Mouse,"Mouse")
HumanStats <- run_disease_stats_all(class.files$Human,"Human")
run_disease_stats_all(class.files$Fetal,"Fetal")
run_disease_stats_all(class.files$Fetal,"Adult")
