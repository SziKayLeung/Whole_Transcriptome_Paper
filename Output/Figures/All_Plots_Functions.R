# Szi Kay Leung (sl693@exeter.ac.uk)
# Functions Script for plots for IsoSeq Paper (2020)

### Packages #########################################################################################

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(reshape))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(stringr)) 
suppressMessages(library(viridis)) 
suppressMessages(library(wesanderson)) 
suppressMessages(library(extrafont))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(tibble))
suppressMessages(library(VennDiagram))
suppressMessages(library(directlabels))
suppressMessages(library(cowplot))
suppressMessages(library(wesanderson))


### Source additional scripts #############################################################################
# Read in input files 
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Input_Variables.R")
# RNASeq 
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/RNAseqvsISOseq/RNASeqvsIsoSeq_Expression.R")

### Generic Functions and Variables for Plots ##############################################################

# theme for plots
mymaintheme <- theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 text=element_text(size=20,  family="ArialMT"),
                 axis.title.x = element_text(vjust=-0.5, colour = "black"),
                 axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                 legend.position = c(.90, 0.95),
                 #legend.justification = c(1,1),
                 legend.box.just = "right",
                 legend.margin = margin(6, 6, 6, 6), 
                 legend.text = element_text(size = 20,family="ArialMT"),
                 axis.text.x= element_text(size=15,  family="ArialMT"),
                 axis.text.y= element_text(size=15,  family="ArialMT"))

mytheme <- theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   text=element_text(size=14,  family="ArialMT"),
                   axis.title.x = element_text(vjust=-0.5, colour = "black"),
                   axis.title.y = element_text(vjust=0.5, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                   legend.position = c(.90, 0.95),
                   #legend.justification = c(1,1),
                   legend.box.just = "right",
                   legend.margin = margin(6, 6, 6, 6), 
                   legend.text = element_text(size = 14,family="ArialMT"),
                   axis.text.x= element_text(size=12,  family="ArialMT"),
                   axis.text.y= element_text(size=12,  family="ArialMT"))


# Colour scheme: #440154FF (Human Adult) , #287D8EFF (Human Fetal) , #20A387FF (Human), #FDE725FF (Mouse)
# https://cran.r-project.org/web/packages/viridis/viridis.pdf
label_colour <- function(dataset){
  if (dataset == "Mouse" || dataset == "Mouse Cortex"){"#FDE725FF"} 
  else if (dataset == "Human (Fetal)" || dataset == "Fetal" || dataset == "Human (Fetal) Cortex"){"#287D8EFF"} 
  else if (dataset == "Human (Adult)" || dataset == "Adult" || dataset == "Human (Adult) Cortex"){"#440154FF"} 
  else if (dataset == "Human" || dataset == "Human Cortex"){"#20A387FF"}
  else if (dataset == "targeted"){wes_palette("Darjeeling1")[2]}
  else if (dataset == "whole"){wes_palette("Darjeeling1")[1]}
  else if (dataset == "whole+targeted"){wes_palette("Darjeeling2")[1]}else{
    print("NA")}}

# To scale axis into 1000s
ks <- function(x){ format(x/1000, big.mark=",")} 


density_plot <- function(dat,x.var,y.var, x_lab, y_lab,title){
  
  
  print(paste0(title))
  print(paste0("Correlation between", x.var, "and", y.var))
  
  print(cor.test(dat[[x.var]],dat[[y.var]]))
  cor(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")
  
  corr.value <- cor(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")
  p.value <- cor.test(dat[[x.var]],dat[[y.var]], use = "pairwise.complete.obs")$p.value 
  
  
  # corr.value <- cor(FSM_TPM$ISOSEQ_TPM_Normalised,FSM_TPM$RNASeq_TPM) # normalised ISOSEQ FL counts to length
  corr <- grobTree(textGrob(paste("r = ", round(corr.value, 2)), 
                            x = 0.05, y = 0.97, hjust = 0, 
                            gp = gpar(col = "black", fontsize = 20, fontface = "italic")))
  
  x.var <- rlang::sym(quo_name(enquo(x.var)))
  y.var <- rlang::sym(quo_name(enquo(y.var)))
  
  print(paste0("corr.value", corr.value))
  print(paste0("p.value", p.value))
  
  p <- ggplot(dat, aes(x = !! x.var, y = !! y.var)) +
    annotation_custom(corr) +
    stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
    geom_point(size = 0.4, alpha = 0.25) +
    scale_fill_distiller(palette=4, direction=1, name = "Density") +
    theme_bw() +
    labs(x = x_lab, y = y_lab, title = paste(title,"\n\n")) + 
    geom_smooth(method=lm, colour = "black") + 
    mytheme + 
    theme(legend.position = "none")
  
  return(p)
}

geom_plot <- function(dat, x_var, y_var, group_var, xlabel, ylabel){
  x_var <- rlang::sym(quo_name(enquo(x_var)))
  y_var <- rlang::sym(quo_name(enquo(y_var)))
  group_var <- rlang::sym(quo_name(enquo(group_var)))
  
  p <- ggplot(dat, aes(x = !!x_var, y = !!y_var, fill = !! group_var)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    mytheme +
    labs(x = xlabel, y = ylabel) + 
    scale_fill_viridis(discrete = TRUE, 
                       name="",
                       labels=c("Human(Adult)", "Human(Fetal)", "Mouse")) + 
    scale_y_continuous(labels = percent, limits = c(0,1))
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_y_continuous(labels = comma) +
  
  
  return(p)
}

venn_diagram_plot_twocircles <- function(set1, set2, label_set1, label_set2){
  
  p <- venn.diagram(
    x = list(set1, set2),
    category.names = c(paste(label_set1,"Cortex"),paste(label_set2,"Cortex")),
    filename = NULL,
    output=TRUE,
    
    # Circles
    lwd = 0.2,
    lty = 'blank',
    fill = c(label_colour(label_set1),label_colour(label_set2)),
    
    # Numbers
    cex = 3,
    fontface = "bold",
    fontfamily = "ArialMT",
    
    # Set names
    cat.cex = 2,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "ArialMT",
    #rotation = 1, 
    main = "\n\n\n\n",
    
    print.mode = "raw"
  )
  
  return(p)
  
}


venn_diagram_plot_twocircles_basic <- function(set1, set2, label_set1, label_set2){
  
  p <- venn.diagram(
    x = list(set1, set2),
    category.names = c(label_set1,label_set2),
    filename = NULL,
    output=TRUE,
    
    # Circles
    lwd = 0.2,
    lty = 'blank',
    fill = c(label_colour(label_set1),label_colour(label_set2)),
    
    # Numbers
    cex = 3,
    fontface = "bold",
    fontfamily = "ArialMT",
    
    # Set names
    cat.cex = 2,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "ArialMT",
    #rotation = 1, 
    main = "\n\n\n\n",
    
    print.mode = "raw"
  )
  
  return(p)
  
}

# summary_info <read.clasification.file> 
# Aim: summarise the number of isoforms, minimum number of exons and max etc per gene 
summary_info <- function(dat){
  
  total_fl <- sum(dat$FL, na.rm=T)
  
  info <- list(
    # Number of isoforms 
    dat %>% count(associated_gene) %>% select(associated_gene, "Num_of_Isoforms" = "n"),
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

# homologs_genes 
# Aim: Remove genes from the homologous gene list that are present in both human and mouse genome, eliminating number overestimation or underestimation 
# Input: homologs gene list, mm10 and hg38 genome list 
# Output: same format as homologs gene list but removed common genes 
homologs_genes_filter <- function(){
  # read in list of homologous genes between mouse and human; and capitalise for consistency
  homologs_genes <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/mousehumangeneconversion.csv")
  homologs_genes$mouse <- toupper(homologs_genes$mouse)
  homologs_genes$human <- toupper(homologs_genes$human)
  cat("Number of homologous genes listed between human and mouse:", nrow(homologs_genes), "\n") 
  
  # list only the unique genes with different names across human and mouse (from read in homologous list)
  # if homologous human and does not match with the homologous gene, take note of the row number and append to vector list 
  # use vector list of row numbers to extract the rows with mismatching names 
  unique_naming_row <- c()
  for(i in 1:nrow(homologs_genes)){
    mismatching_row <- if(homologs_genes$human[i] != homologs_genes$mouse[i]){i}
    unique_naming_row <- append(mismatching_row, unique_naming_row)
  }
  unique_homologs_genes <- homologs_genes[unique_naming_row,]
  cat("Number of mismatching gene names between human and mouse:", nrow(unique_homologs_genes), "\n") 
  
  
  ## Remove any genes from the unique genes list that contain genes that are present across both species
  # i.e SMS gene is present in both mouse and human genome; in unique gene list: 
  # SMS-PS (mouse) --> SMS (human); 
  # if convert SMS-PS to SMS, there may be two separate SMS records (one from converted SMS-PS and SMS original gene) 
  # result in number overestimation
  # i.e RPL29 gene is present in both mouse and human genome; in uniquge gene list: 
  # RPL29 (mouse) --> LINC01621 (human); if convert RPL29 to LINC01621
  # missed RPL29 when matching human and mouse genome even though RPL29 is present in human -----> number underestimation
  # duplicate_mouse <- genes present in unique human list for conversion also present in mouse genome
  # duplicate_human <- genes present in unique mouse list for conversion also present in human genome 
  duplicate_mouse <- intersect(unique(mm10$Gene), unique_homologs_genes$human) 
  duplicate_human <- intersect(unique(hg38$GeneSymbol), unique_homologs_genes$mouse)
  cat("Number of genes from human homologous list present in mouse genome, thus removed for conversion", length(duplicate_mouse), "\n") 
  cat("Number of genes from mouse homologous list present in human genome, thus removed for conversion", length(duplicate_human), "\n") 
  # filter list of genes from unique homologous genes so no conversion
  homologs_genes_filtered <- homologs_genes %>% 
    filter(!human %in% duplicate_mouse) %>% 
    filter(!mouse %in% duplicate_human)
  
  return(homologs_genes_filtered)
}



### QC Plots ########################################################################################

# ccs_length_distribution <ccs_input_dir> <label_for_plot, i.e. "Mouse"> <type = "individual", "Merge> 
# Aim: Plot the distribution of CCS read length per dataset (individual) or merged dataset
# Input: path directory of fasta.seqlength generated from CCS_Length_Stats.sh
# Output of Length Distribution of CCS reads based for each dataset (based on ccs_input_dir)
ccs_length_distribution <- function(ccs_input_dir_path, plot_label, type){
  
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  fetal_filenames <- list.files(path = ccs_input_dir$Fetal, pattern = "fasta.seqlengths.txt$", 
                                full.names = TRUE)
  adult_filenames <- list.files(path = ccs_input_dir$Adult, pattern = "fasta.seqlengths.txt$", 
                                full.names = TRUE)
  mouse_filenames <- list.files(path = ccs_input_dir$Mouse, pattern = "fasta.seqlengths.txt$", 
                                full.names = TRUE)
  
  # plot of histogram distribution for all samples
  plot_distribution <- function(dat, dataset1, dataset2){
    p <- ggplot(dat, aes(V2, colour = Sample)) +
      geom_freqpoly(binwidth = 500,size = 1.5 ) + 
      labs(x = "CCS Read Length (kb)", y = "Number of Reads (Thousand)") + 
      scale_fill_manual(values = cbbPalette) + 
      theme_bw() + mytheme +
      scale_y_continuous(labels = ks) + 
      scale_x_continuous(labels = ks) +
      theme(legend.title = element_blank(), legend.position = c(0.7, 0.6)) +
      scale_color_manual(values=c(label_colour(dataset1), label_colour(dataset2)))
    return(p)
  }
  
  # distinguish data from the "Sample" column and rbind to one dataframe as "all"
  if(type == "Individual"){
    # read in all ccs files  as list 
    filenames <- list.files(path = ccs_input_dir_path, pattern = "fasta.seqlengths.txt$", full.names = TRUE)
    ccs_files <- lapply(filenames, read.table)
    cat("Working with the following files", filenames)
    
    for (i in 1:length(ccs_files)){
      ccs_files[[i]]$Sample <- print(paste(plot_label, i))
    }
    
    # all = dataframe of the individual CCS (V1), CCS read length (V2), Sample (V3)
    all <- bind_rows(ccs_files)
    cat("Average CCS Read Length across all individuals", mean(all$V2))
    
    # plot of histogram distribution for all samples in dataset
    output_plot <- ggplot(all, aes(V2, colour = Sample, linetype = Sample)) + 
      geom_freqpoly(binwidth = 500,size = 1.5 ) + 
      labs(x = "CCS Read Length (kb)", y = "Number of Reads (Thousand)") + 
      scale_fill_manual(values = cbbPalette) + 
      theme_bw() + mytheme +
      scale_y_continuous(labels = ks) + 
      scale_x_continuous(labels = ks) +
      theme(legend.title = element_blank(), legend.position = c(0.85,0.85))
    
    # average length across human 
    #ccs_input_dir_path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/IsoSeq/CCS/Lengths/AdultCTX"
    #filenames <- list.files(path = ccs_input_dir_path, pattern = "fasta.seqlengths.txt$", full.names = TRUE)
    #ccs_files_adult <- lapply(filenames, read.table)
    #for (i in 1:length(ccs_files_adult)){ccs_files_adult[[i]]$Sample <- print(paste("Adult", i))}
    
    #ccs_input_dir_path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/IsoSeq/CCS/Lengths/FetalCTX"
    #filenames <- list.files(path = ccs_input_dir_path, pattern = "fasta.seqlengths.txt$", full.names = TRUE)
    #ccs_files_fetal <- lapply(filenames, read.table)
    #for (i in 1:length(ccs_files_fetal)){ccs_files_fetal[[i]]$Sample <- print(paste("Fetal", i))}
    #all <- bind_rows(ccs_files_adult,ccs_files_fetal)
    #cat("Average CCS Read Length across all human samples", mean(all$V2))
    
    
  }else if (type == "Human_Mouse_Merge"){
    fetal_ccs_files <- lapply(fetal_filenames, function(x) read.table(x) %>% mutate(Sample = "Human Cortex"))
    adult_ccs_files <- lapply(adult_filenames, function(x) read.table(x) %>% mutate(Sample = "Human Cortex"))
    mouse_ccs_files <- lapply(mouse_filenames, function(x) read.table(x) %>% mutate(Sample = "Mouse Cortex"))
    all <- bind_rows(fetal_ccs_files, adult_ccs_files, mouse_ccs_files)
    output_plot <- plot_distribution(all, "Human Cortex","Mouse Cortex")
  }else if (type == "Fetal_Adult_Merge"){
    fetal_ccs_files <- lapply(fetal_filenames, function(x) read.table(x) %>% mutate(Sample = "Human (Fetal) Cortex"))
    adult_ccs_files <- lapply(adult_filenames, function(x) read.table(x) %>% mutate(Sample = "Human (Adult) Cortex"))
    all <- bind_rows(fetal_ccs_files, adult_ccs_files)
    output_plot <- plot_distribution(all, "Human (Fetal)","Human (Adult)")
  }else{
    paste0("Nothing")
  }
  
  return(output_plot)
}

# rarefaction_distribution
# Aim: Plot the rarefaction curves for gene and isoform level for each dataset and merged
# Input: path directory of rarefaction.by_refgene.min_fl_2.txt 
rarefaction_distribution <- function(){
  
  # Merge input files pertainng to genes and isoforms 
  genes <- bind_rows(all_rarefaction_genes)
  isoforms <- bind_rows(all_rarefaction_isoforms)
  all_rarefaction_levels <- bind_rows(genes,isoforms)
  
  ## Plots 
  # p1 <- isoform and gene level for human and mouse 
  # p2 <- isoform and gene level for adult and fetal 
  # p3 <- isoform per category for individual datasets
  plot_distribution <- function(dataset1, dataset2){
    p <- all_rarefaction_levels %>% 
      filter(Sample %in% c(dataset1, dataset2)) %>% 
      ggplot(., aes(x = size, y = mean, color = Sample, linetype = type)) + 
      geom_line(size = 1.5) + 
      labs(x ="Number of Subsampled Reads (Thousand)", y = "Number of Genes/Isoforms (Thousand)") + 
      theme_bw() + mytheme + 
      scale_y_continuous(labels = ks) + scale_x_continuous(labels = ks) + 
      scale_color_manual(values=c(label_colour(dataset1), label_colour(dataset2)), 
                         labels=c(paste(dataset1,"Cortex"), paste(dataset2, "Cortex")), name = "") + 
      scale_linetype_manual(values=c("longdash","solid")) +
      theme(legend.position = c(0.8, 0.6), legend.spacing.y = unit(-0.1, "cm"),legend.title = element_blank())
    
    return(p)
  }
  
  p1 <- plot_distribution("Human","Mouse")
  p2 <- plot_distribution("Human (Adult)","Human (Fetal)")
  
  plot_category <- list()
  count = 1 
  for(i in all_rarefaction_isoforms_category){
    # plot order: Mouse, Adult, Fetal, Human
    p <- ggplot(i, aes(x=size, y=mean, color=category)) + geom_line(aes(linetype = type), size = 1.5) + 
      labs(x = "Number of Subsampled Reads (Thousand)", 
           y = "Number of Isoforms (Thousand)", title = paste0("\n\n")) +
      mytheme + 
      scale_y_continuous(labels = ks, limits = c(0, 20000)) + 
      scale_x_continuous(labels = ks, expand = c(0.15, 0)) + 
      theme(legend.position = "none") +
      scale_linetype_manual(values=c("solid", "longdash","dotted")) +
      geom_dl(aes(label = category),  method = list(dl.trans(x = x + 0.5), "last.bumpup", cex = 1.3, hjust = .1))
    
    plot_category[[count]] <- p 
    count = count + 1
  }
  
  names(plot_category) <- c("Mouse","Adult","Fetal","Human")
  
  output <- list(p1, p2, plot_category)
  return(output)
}

# Plot results from GO analysis of most abundant 500 Genes
most_abundant500Genes <- function(){
  # Read in GO results from Top 500 most abundant genes from human
  human_gene_atlas <- read.csv(paste0(output_table_dir, "/GO/Human_500topGenes.csv")) %>% arrange(Human_Gene_Atlas.Adjusted.P.value) %>% .[1:5,] %>% mutate(Sample = "Human Cortex")
  
  # Read in GO results from Top 500 most abundant genes from mouse
  mouse_gene_atlas <- read.csv(paste0(output_table_dir, "/GO/Mouse_500topGenes.csv")) %>% arrange(Mouse_Gene_Atlas.Adjusted.P.value) %>% .[1:5,] %>% mutate(Sample = "Mouse Cortex") 
  
  # PLOT
  # rename after checking the terms
  mouse_gene_atlas$Term <- c("Cerebral Cortex", "Cerebral Cortex Prefontal", "Spinal Chord", "Nucleus Accumbens", "Olfactory Bulb")
  human_gene_atlas$Term <- c("Prefrontal Cortex", "Amygdala", "Whole Brain", "Fetal Brain", "Spinal Chord")
  human_gene_atlas <- human_gene_atlas %>% select(Adjusted.P.value = Human_Gene_Atlas.Adjusted.P.value, Sample,Term)
  mouse_gene_atlas <- mouse_gene_atlas %>% select(Adjusted.P.value = Mouse_Gene_Atlas.Adjusted.P.value, Sample,Term)
  
  p <- rbind(human_gene_atlas, mouse_gene_atlas) %>% 
    ggplot(., aes(x = reorder(Term, Adjusted.P.value), y = -log(Adjusted.P.value), fill = Sample)) + 
    geom_bar(stat = "identity") + 
    mytheme + 
    theme(axis.text.x = element_text(angle = 45,hjust=1), 
          legend.position = "none") + 
    facet_grid(~Sample,scales="free_x", space = "free") + 
    labs(x = "", y = "-log(adj P value)") + 
    scale_fill_manual(values=c(label_colour("Human (Fetal)"), label_colour("Mouse")), name = "")  +
    theme(legend.position="none", strip.background = element_rect(colour="white", fill="white"),
          strip.text.x = element_text(size=20,  family="ArialMT"))
  
  return(p)
  
}

iso_length <- function(class){
  class <- class %>% filter(subcategory != "mono-exon")
  p <- ggplot(class, aes(x = length)) + geom_histogram(bins = 15, fill="gray", col="black") + 
    labs(x = "Transcript Length (kb)", y = "Number of Isoforms (Thousand)") + mytheme +
    scale_x_continuous(labels = ks) + 
    scale_y_continuous(labels = ks) 
  
  two <- class[which(class$length >= 2000 & class$length <= 4000),] %>% nrow()
  print(paste0("Number of isoforms 2-4kb:", two, "(",round(two/nrow(class),2) *100,"%)"))
  
  return(p)
}

### RNASeq ################################################################################

# rnaseq_isoseq_correlation <read_sqanti_file>
# Aim: Correlate the gene expression and transcript expression of RNASeq vs IsoSeq 
# RNASeq expression is from Kallisto alignment of merged reads to IsoSeq transcriptome (in SQANTI2)
# IsoSeq expression is FL reads 
# Input: Functions from RNASeqvsIsoSeq_Transcript.R
# Note: there are isoforms and genes that are not detected by RNASeq and thus removed in correlation
rnaseq_isoseq_correlation <- function(sqanti_file){
  
  # 0.005 for better plot
  Gene_corr_cut <- Merge_Gene_Kallisto_Input(sqanti_file, 0, "All")
  Gene_Expression_Cut_off <- Gene_corr_cut  %>% filter(log_ISOSEQ_GENE_TPM > 2.5)
  Transcript_corr_cut <- Merge_Transcript_Kallisto_Input(sqanti_file, 0.01)
  
  ## Plots 
  # p1 <- Gene Level Correlation
  # p2 <- Gene Level Correlation with cutoff expression (Log Isoseq Gene TPM > 2.5)
  # p3 <- Isoform Level Correlation
  
  p1 <- density_plot(Gene_corr_cut,"log_ISOSEQ_GENE_TPM","log_RNASEQ_GENE_TPM",
                     "Iso-Seq Expression (Log10 TPM)",
                     "RNA-Seq Expression (Log10 TPM)","") +
    theme(legend.position = "none") 
  
  p2 <- density_plot(Gene_Expression_Cut_off,"log_ISOSEQ_GENE_TPM","log_RNASEQ_GENE_TPM",
                     "Iso-Seq Expression (Log10 TPM)",
                     "RNA-Seq Expression (Log10 TPM)","" ) +
    theme(legend.position = "none") 
  
  p3 <- density_plot(Transcript_corr_cut,"log_ISOSEQ_TPM","log_RNASEQ_TPM",
                     "Iso-Seq Expression (Log10 TPM)",
                     "RNA-Seq Expression (Log10 TPM)","" ) +
    theme(legend.position = "none") 
  
  return(list(p1,p2,p3))
}

# rnaseq_junction_support 
# Aim: Plot the number of unique junctions supported by RNASeq (>1 read)
# Output: 2 plots 
# Plot1: Junctions supported by RNASeq coverage 
# Plot2: Junctions in general and whether supported by RNASeq coverage
rnaseq_junction_support <- function(){
  junc.files <- list(junc.files$Mouse, junc.files$Fetal)
  names(junc.files) <- c("Mouse","Fetal")
  
  # function for counts 
  calculate_count_junction_coverage <- function(dat){
    uniqJuncCov <- unique(dat[,c("junctionLabel","SJ_type", "total_coverage", "splice_site")])
    
    # total frequency of type of splice junctions regardless of coverage in IsoSeq dataset
    e <- data.frame(table(uniqJuncCov$SJ_type))
    f <- data.frame(table(uniqJuncCov[which(uniqJuncCov$total_coverage>0),"SJ_type"]))
    
    df.SJcov <- merge(e, f, by="Var1")
    # Junction suppported and unsupported
    df.juncSupport <- data.frame(type=e$Var1, count=e$Freq-f$Freq, name='Not Supported')
    df.juncSupport <- rbind(df.juncSupport, data.frame(type=f$Var1, count=f$Freq, name='Supported'))
    
    print("Further subset of novel canonical junctions that are not supported (1st - Mouse, 2nd Fetal)")
    novelcanonical <- uniqJuncCov %>% filter(SJ_type == "Novel\nCanonical " & total_coverage == 0) %>% group_by(splice_site) %>% tally()
    print(novelcanonical)
    
    return(df.juncSupport)
  }
  
  # function for percentage 
  calculate_perc_junction_coverage <- function(dat){
    # unique junctions 
    # examples of junctions with more than one transcript having that junction i.e chr1_-_10025064_10027102
    uniqJuncCov <- unique(dat[,c("junctionLabel","SJ_type", "total_coverage","splice_site")])
    
    # total frequency of type of splice junctions regardless of coverage in IsoSeq dataset
    e <- data.frame(table(uniqJuncCov$SJ_type))
    f <- data.frame(table(uniqJuncCov[which(uniqJuncCov$total_coverage>0),"SJ_type"]))
    
    df.SJcov <- merge(e, f, by="Var1")
    # calculate the percentage of junctions that have more short read junction coverage
    df.SJcov$perc <- df.SJcov$Freq.y / df.SJcov$Freq.x ;
    df.SJcov[is.na(df.SJcov$perc), "perc"] <- 0
    
    return(df.SJcov)
  }
  
  # return table of stats if necessary
  table1 <- lapply(junc.files, function(x) calculate_count_junction_coverage(x)) %>%
    do.call("rbind", . ) %>%
    rownames_to_column(., var = "Sample") %>%
    mutate(Sample = word(Sample, c(1), sep = fixed ('.'))) %>% 
    group_by(Sample) %>% mutate(perc = count/sum(count) * 100)
  table1$name <-  factor(table1$name, levels = c("Supported","Not Supported"))
  
  #sum(table1[table1$Sample == "Mouse","count"]) #sum of total junctions 
  #sum(table1[table1$Sample == "Mouse" & table1$name == "Supported","count"]) # sum of supported junctions
  #sum(table1[table1$Sample == "Fetal","count"]) #sum of total junctions 
  #sum(table1[table1$Sample == "Fetal" & table1$name == "Supported","count"]) # sum of supported junctions
  print(table1)
  
  # plotting percentage support by short-reads 
  table1$Sample[table1$Sample == "Fetal"] <- "Human (Cortex)"
  table1$Sample[table1$Sample == "Mouse"] <- "Mouse (Cortex)"
  p1 <- lapply(junc.files, function(x) calculate_perc_junction_coverage(x)) %>%
    do.call("rbind", . ) %>%
    rownames_to_column(., var = "Sample") %>%
    mutate(Sample = word(Sample, c(1), sep = fixed ('.'))) %>%
    ggplot(., aes(x = Var1, y = perc, fill = Sample)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "", y = "Percentage of Unique Junctions with RNASeq coverage") +
    scale_fill_manual(values=c(label_colour("Human"), label_colour("Mouse")), 
                      labels=c("Human", "Mouse"), name = "") +
    scale_y_continuous(labels = percent, limits = c(0,1)) + 
    mytheme + theme(legend.position="bottom")
  
  
  # plotting percentage support and not support 
  p2 <- ggplot(table1, aes(x = name, y = perc, fill = type)) + 
    geom_bar(stat = "identity") +
    facet_grid(~Sample) + 
    theme_bw() + 
    mytheme + 
    labs(x = "", y = "Unique Junctions (%)") + 
    scale_fill_manual(values = wes_palette("Darjeeling2")) +
    theme(legend.position="bottom",legend.title = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          strip.text.x = element_text(size=20,  family="ArialMT"))
  
  return(list(p1,p2))
}

junction_support_subset <- function(genome_type){
  # data.junction files contain only multiexonic transcripts (as monoexonic transcripts do not contain junctions)
  # Create new column of whether RNASeq supported of the junction based on total coverage
  if(genome_type == "intropolis"){
    input.junc.files <- lapply(junc.files, function(x) x %>% mutate(RNASeq = ifelse(.$total_coverage > 0, "Supported", "Not_Supported")))
    
    # Filter only multiexonic transcripts 
    multiexonic_genes <- lapply(class.files, function(x) x[x$subcategory != "mono-exon",])
    
    # Extract isoform id for novel and annotated genes
    # isoform id of novel genes
    novel_genes_isoform_id <- lapply(multiexonic_genes, function(x) x[grepl("NOVEL", x$associated_gene),"isoform"])
    # isoform id of annotated genes (known and novel transcripts)
    annotated_genes_isoform_id <- lapply(multiexonic_genes, function(x) x[!grepl("NOVEL", x$associated_gene),"isoform"])
    annotated_genes_novel_transcripts_isoform_id <- lapply(multiexonic_genes, function(x) x[!grepl("NOVEL", x$associated_gene) & 
                                                                                              x$associated_transcript == "novel","isoform"])
    
  } else {
    
    input.junc.files <- list(junc.files$Mouse,junc.files$Fetal)
    input.junc.files <- lapply(input.junc.files, function(x) x %>% mutate(RNASeq = ifelse(.$total_coverage > 0, "Supported", "Not_Supported")))
    
    # Filter only multiexonic transcripts 
    multiexonic_genes <- lapply(list(class.files$Mouse, class.files$Fetal), function(x) x[x$subcategory != "mono-exon",])
    names(multiexonic_genes) <- c("Mouse","Fetal")
    
    # Extract isoform id for novel and annotated genes
    # isoform id of novel genes
    novel_genes_isoform_id <- lapply(multiexonic_genes, function(x) x[grepl("NOVEL", x$associated_gene),"isoform"])
    # isoform id of annotated genes (known and novel transcripts)
    annotated_genes_isoform_id <- lapply(multiexonic_genes, function(x) x[!grepl("NOVEL", x$associated_gene),"isoform"])
    annotated_genes_novel_transcripts_isoform_id <- lapply(multiexonic_genes, function(x) x[!grepl("NOVEL", x$associated_gene) & 
                                                                                              x$associated_transcript == "novel","isoform"])
    
  }
  
  
  # Internal function: 
  # subset <junction> file of isoforms (based on <isoform_id>)
  # know the matching gene to isoform id using <multiexonic_class_file>
  # output: table of each isoform and the number of junctions supported or not supported by RNASeq 
  count_and_support <- function(multiexonic_class_file, input_junction, isoform_id){
    subset_input_junction <- input_junction[input_junction$isoform %in% isoform_id,]
    subset_input_junction$Sample <- multiexonic_class_file$Sample[1]
    
    subset_input_junction_support <- 
      subset_input_junction  %>% group_by(isoform, RNASeq, Sample) %>% tally() %>%
      # merge to match the isoform id for the transcript from the associated gene 
      left_join(., multiexonic_class_file[,c("isoform","associated_gene")], by = "isoform") %>%
      spread(RNASeq, n) %>%   
      replace(is.na(.), 0)
    
    subset_input_junction_support$total_num_input_junctions <- rowSums(subset_input_junction_support[4:5])
    subset_input_junction_support$support_coverage <- subset_input_junction_support$Supported/subset_input_junction_support$total_num_input_junctions
    
    return(subset_input_junction_support)
  }
  
  
  # Plots 
  # Plots the cumulative junction coverage of the <input_support_file> generated from count_and_support
  # Transcript type refers to the title of the plot
  plot_junction_support <- function(input_support_file, transcript_type){
    
    input_support_file <- bind_rows(input_support_file) %>%
      group_by(support_coverage, Sample) %>%
      count() 
    
    p1 <- input_support_file %>%
      group_by(Sample) %>% arrange(support_coverage) %>% mutate(csum = cumsum(n)) %>%
      ggplot(., aes(x = support_coverage, y = csum, color = Sample)) + 
      geom_line() + 
      labs(y = paste0("Number of ", transcript_type), x = "Junction Coverage by RNA-Seq") + 
      theme_bw() + 
      mytheme + 
      theme(legend.position = c(.20, 0.90), legend.title = element_blank()) 
    
    total <- aggregate(input_support_file$n, by=list(Sample=input_support_file$Sample), FUN=sum)
    
    p2 <- input_support_file %>%
      left_join(total, by = "Sample") %>%
      mutate(perc = n/x * 100) %>% 
      group_by(Sample) %>% arrange(support_coverage) %>% mutate(cperc = cumsum(perc)) %>%
      ggplot(., aes(x = support_coverage, y = cperc, color = Sample)) + 
      geom_line() + 
      labs(y = paste0(transcript_type, " (%)"), x = "Junction Coverage by RNA-Seq") + 
      theme_bw() + 
      mytheme + 
      theme(legend.position = c(.20, 0.90), legend.title = element_blank()) 
    
    # note the percentages do not sum to 100 as also missing the other transcripts with in between RNASeq support coverage 
    stats <- input_support_file %>%
      left_join(total, by = "Sample") %>%
      mutate(perc = n/x * 100, ) %>% 
      group_by(Sample) %>% arrange(support_coverage) %>% 
      filter(support_coverage %in% c(0, 1))%>% 
      `colnames<-`(c("RNASeq support coverage of junctions", "Sample", "number of transcripts", "number of total transcripts", "percentage")) %>%
      as.data.frame()
    
    if(genome_type == "intropolis"){
      p1 <- p1 +  scale_color_manual(labels = c("Human(Adult)","Human(Fetal)","Human","Mouse"), 
                                     values = c("#440154FF","#287D8EFF","#29AF7FFF", "#FDE725FF"), 
                                     name = "")
      p2 <- p2 + scale_color_manual(labels = c("Human(Adult)","Human(Fetal)","Human","Mouse"), 
                                    values = c("#440154FF","#287D8EFF","#29AF7FFF", "#FDE725FF"), 
                                    name = "")
    }else{
      p1 <- p1 +  scale_color_manual(labels = c("Human(Fetal) Cortex ","Mouse Cortex"), values = c("#29AF7FFF", "#FDE725FF"), 
                                     name = "")
      p2 <- p2 +   scale_color_manual(labels = c("Human(Fetal) Cortex","Mouse Cortex"), values = c("#29AF7FFF", "#FDE725FF"), 
                                      name = "")
      
    }
    
    print(stats)
    return(list(p1,p2,stats))
  }
  
  
  annotated_genes_novel_transcripts_support <- list()
  novel_genes_all_transcripts_support <- list()
  for(i in 1:length(names(multiexonic_genes))){
    # RNASeq support of annotated genes, novel transcripts
    annotated_genes_novel_transcripts_support[[i]] <- 
      count_and_support(multiexonic_genes[[i]], input.junc.files[[i]], annotated_genes_novel_transcripts_isoform_id[[i]])
    
    # RNASeq support of novel genes
    novel_genes_all_transcripts_support[[i]] <- 
      count_and_support(multiexonic_genes[[i]], input.junc.files[[i]], novel_genes_isoform_id[[i]])
  }
  names(annotated_genes_novel_transcripts_support) <- names(multiexonic_genes)
  
  cat("Annotated Genes, Novel Transcripts")
  output1 <- plot_junction_support(annotated_genes_novel_transcripts_support, "Annotated Genes, Novel Transcripts")
  print(output1)
  
  cat("Novel Genes")
  output2 <- plot_junction_support(novel_genes_all_transcripts_support, "Novel Genes")
  print(output2)
  
  return(list(output1,output2))
  
}

fetal_mouse_rnaseq_distribution <- function(){
  # mouse_annotated_genes_novel_transcript_isoforms
  mouse <- class.files$Mouse[!grepl("NOVEL", class.files$Mouse$associated_gene) & class.files$Mouse$associated_transcript == "novel",]
  own_fetal_class_file <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Fetal/SQANTI2/own_rnaseq/fetalFL.sample.collapsed.filtered.rep_classification.filtered_lite_classification.txt", as.is = T, sep = "\t", header = T)
  
  # own_fetal_annotated_genes_novel_transcript_isoforms 
  fetal <- own_fetal_class_file[!grepl("novel", own_fetal_class_file$associated_gene) & own_fetal_class_file$associated_transcript == "novel",]
  
  
  p1 <- rbind(data.frame("Num" = fetal$iso_exp, "Sample" = "human"), 
              data.frame("Num" = mouse$iso_exp, "Sample" = "mouse")) %>%
    ggplot(., aes(x = log10(Num), fill = Sample)) + geom_density(alpha = 0.5) + 
    theme_bw() + 
    labs(x = "RNASeq Isoform Expression (Log10TPM)", y = "Density") +
    mytheme + 
    theme(legend.position = c(.20, 0.90), legend.title = element_blank()) +
    scale_fill_manual(labels = c("Human(Fetal)","Mouse"), values = c("#287D8EFF", "#FDE725FF"), name = "")
  
  stats <- data.frame()
  stats[1,1:2] <- c(nrow(mouse[mouse$iso_exp == 0,]),nrow(mouse))
  stats[2,1:2] <- c(nrow(fetal[fetal$iso_exp == 0,]),nrow(fetal))
  colnames(stats) <- c("Num_transcripts_no_RNASeq coverage","Total_Number_of_novel_transcripts")
  stats$perc <- round(stats$`Num_transcripts_no_RNASeq coverage`/stats$Total_Number_of_novel_transcripts * 100,2)
  rownames(stats) <- c("Mouse","Fetal")
  
  return(p1)
}

rnaseq_isoseq_transcriptome <- function(cuffrefmap_input,cufftmap_input){
  #### READ Files from Gff compare 
  # note only include isoforms from Iso-Seq that are partially or fully matched to RNA-Seq, not all isoforms from Iso-Seq dataset 
  cuffrefmap <- read.table(cuffrefmap_input, header = T)
  # no duplicates of isoform - QC checked
  # cuffrefmap[cuffrefmap$class_code == "=",] %>% group_by(ref_id) %>% tally() %>% filter(n > 1)
  
  cufftmap <- read.table(cufftmap_input, header = T)
  # Note not all isoforms are detected by RNA-Seq therefore not listed in cufftmap
  # novel_gene <- class.files$Mouse[grepl("NOVE", class.files$Mouse$associated_gene),"isoform"]
  # for(i in novel_gene){print(cufftmap[cufftmap $ref_id == i,])}
  
  # classification of rnaseq reads to isoseq 
  # cufftmap %>% group_by(class_code) %>% tally() %>% mutate(perc = n/sum(n) * 100) %>% ggplot(., aes(x = reorder(class_code, -perc), y = perc)) + geom_bar(stat = "identity")
  # cufftmap %>% filter(class_code %in% c("u")) %>% ggplot(.,aes(x = num_exons)) + geom_bar(aes(y = (..count..)/sum(..count..)))
  
  ## 
  # Replace id of those rnaseq isoforms that fully match with isoseq isoforms with pbid for venn diagram
  rnaseq <- c(as.character(cufftmap[cufftmap$class_code != "=","qry_id"]),as.character(cuffrefmap[cuffrefmap$class_code == "=","ref_id"])) 
  p1 <- venn.diagram(
    x = list(rnaseq, class.files$Mouse$isoform), category.names = c("RNA-Seq","Iso-Seq"), filename = NULL, output=TRUE,
    lwd = 0.2,lty = 'blank', fill = c("#B3E2CD", "#FDCDAC"), main = "\n",
    cex = 1,fontface = "bold",fontfamily = "ArialMT",
    cat.cex = 1,  cat.default.pos = "outer",  cat.pos = c(-27, 27),  cat.dist = c(0.055, 0.055),  cat.fontfamily = "ArialMT",  #rotation = 1,   main = "\n\n\n\n",
    print.mode = "raw"
  )
  
  
  ## Number of Isoforms per dataset 
  num_iso <- list(class.files$Mouse %>% mutate(Sample = "Iso-Seq") %>%  .[,c("isoform", "associated_gene", "novelGene","FSM_class","gene_exp","Sample")],
                  rnaseq.class.files %>%  .[,c("isoform","associated_gene", "novelGene","FSM_class","gene_exp","Sample")])
  isoPerGene <- lapply(num_iso, function(x) SQANTI_gene_preparation(x)) %>% bind_rows()
  
  # Total Number of Genes per Type 
  Total_Num <- isoPerGene %>% group_by(Sample) %>% count(Sample)
  # Number of isoform cateogories per Sample
  p2 <- isoPerGene %>% group_by(Sample) %>% count(nIsoCat) %>% full_join(Total_Num,., by = "Sample") %>% mutate(Perc = n.y/n.x * 100) %>% 
    ggplot(., aes(x=nIsoCat, fill=Sample)) +
    geom_bar(stat="identity", aes(y= Perc, group = as.factor(Sample)), color="black", size=0.3, width=0.7, 
             position="dodge") + 
    labs(x ="Number of Isoforms", y = "Genes (%)", fill = "", title = "\n") +
    #scale_fill_manual(values=c(label_colour(dataset1), label_colour(dataset2)), 
    #                 labels=c(paste(dataset1,"Cortex"), paste(dataset2,"Cortex")))  + 
    mytheme + 
    theme(legend.position = c(0.75,0.95))
  
  # Length
  p3 <- bind_rows(class.files$Mouse %>% mutate(Sample = "Iso-Seq") %>% .[,c("length","Sample")], rnaseq.class.files %>% .[,c("length","Sample")]) %>% 
    ggplot(., aes(x = Sample, y = log10(length/1000))) + geom_boxplot() + mytheme +
    labs(x = "", y = "Isoform Length (Log10 kb)", title = "\n")
  
  # Exons
  p4 <- bind_rows(class.files$Mouse %>% mutate(Sample = "Iso-Seq") %>% .[,c("exons","Sample")], rnaseq.class.files %>% .[,c("exons","Sample")]) %>% 
    ggplot(., aes(x = Sample, y = log10(exons))) + geom_boxplot() + mytheme +
    labs(x = "", y = "Number of Exons (Log10)", title = "\n")
  
  # cage peak
  p5 <- data.frame(dataset = c("Iso-Seq","RNA-Seq"),
                   within_50 = c(nrow(class.files$Mouse %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(class.files$Mouse) * 100,
                                 nrow(rnaseq.class.files %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(rnaseq.class.files) * 100),
                   # include NAs
                   without_50 = c(100 - nrow(class.files$Mouse %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(class.files$Mouse) * 100,
                                  100 - nrow(rnaseq.class.files %>% filter(abs(dist_to_cage_peak) <= 50))/nrow(rnaseq.class.files) * 100)) %>%
    melt() %>%
    mutate(variable = factor(variable, levels = c("without_50","within_50"))) %>%
    ggplot(., aes(x = dataset, y = value, fill = variable)) + geom_bar(stat = "identity") + mytheme +
    labs(y ="Isoforms within 50bp CAGE (%)", x = "", fill = "", title = "\n") + theme(legend.position = "right") +
    scale_fill_discrete(labels = c("No","Yes"))
  
  # structural_category
  p6 <- bind_rows(class.files$Mouse %>% mutate(Sample = "Iso-Seq") %>% .[,c("structural_category","Sample")] %>% 
                    group_by(structural_category,Sample) %>% tally() %>% group_by(Sample) %>% mutate(perc = n/sum(n)*100),
                  rnaseq.class.files %>% .[,c("structural_category","Sample")] %>% group_by(structural_category,Sample) %>% tally() %>% group_by(Sample) %>% mutate(perc = n/sum(n)*100)) %>% 
    ggplot(., aes(x = Sample, y = perc, fill = structural_category)) + geom_bar(stat = "identity") + mytheme +
    labs(x = "", y = "Isoforms (%)", fill = "Structural Category", title = "\n") + theme(legend.position = "bottom")
  
  output <- plot_grid(grobTree(p1),p2,p3,p4,p5,p6, labels = "auto", label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
  
  return(output)
}

# rnaseq_isoseq_counts
# correlate the number of reads from RNASeq alignment to Iso-Seq and RNA-Seq defined transcriptome using only the 23,761 isoforms considered matching from gffcompare (cuffrefmap file from gffcompare output)
rnaseq_isoseq_counts <- function(class_file){
  # only consider the isoforms that are "matching" between Iso-Seq and RNA-Seq defined transcriptome (output from gff compare)
  matching <- cuffrefmap[cuffrefmap$class_code == "=",] %>% mutate(qry_id = word(qry_id_list,c(2),sep = fixed("|")))
  
  # number of isoforms matching = 23,761
  # nrow(cuffrefmap[cuffrefmap$class_code == "=",]) 
  
  # ref_id = pacbio id from isoseq-defined transcriptome gtf
  matching_pbisoform <- cuffrefmap[cuffrefmap$class_code == "=","ref_id"]
  
  # qry_id_list = id from rnaseq-defined transcriptome gtf
  matching_rnaseqisoform <- unique(word(cuffrefmap[cuffrefmap$class_code == "=","qry_id_list"],c(2),sep = fixed("|")))
  
  # QC no repeated refid (PBid from isoforms)
  #n_occur <- data.frame(table(cuffrefmap[cuffrefmap$class_code == "=","ref_id"]))
  #n_occur[n_occur$Freq > 1,]
  
  # repeated qry_id due to imperfect match from cuffrefmap i.e. two PB isoforms refer to the same RNA-Seq isoform if 5' and 3' end different but internal junction the same
  #n_occur <- data.frame(table(word(cuffrefmap[cuffrefmap$class_code == "=","qry_id_list"],c(2),sep = fixed("|"))))
  #n_occur[n_occur$Freq > 1,]
  
  ### FeatureCounts
  # subset Iso-Seq Defined and RNA-Seq defined transcriptome with corresponding matching isoform id for counts  
  #IsoSeq_FL <- class_file[class_file$isoform %in% matching_pbisoform,c("isoform","FL","ISOSEQ_TPM")]
  #IsoSeq_Defmatching <- IsoSeq_Def %>% filter(Geneid %in% matching_pbisoform) %>% select(Geneid, RNA2IsoSeq.sorted.bam) %>% `colnames<-`(c("PBID", "RNA2Isoseq_Counts"))
  #RNASeq_Defmatching <- RNASeq_Def %>% filter(Geneid %in% matching_rnaseqisoform) %>% select(Geneid, RNA2RNASeq.sorted.bam) %>% `colnames<-`(c("RNASeqID", "RNA2RNAseq_Counts"))
  
  ### Kallisto
  IsoSeq_FL <- class_file[class_file$isoform %in% matching_pbisoform,c("isoform","FL","ISOSEQ_TPM")]
  IsoSeq_Defmatching <- IsoSeq_Def %>% filter(target_id %in% matching_pbisoform) %>% select(target_id, est_counts,tpm) %>% `colnames<-`(c("PBID", "RNA2Isoseq_Counts","RNA2Isoseq_TPM"))
  RNASeq_Defmatching <- RNASeq_Def %>% filter(target_id %in% matching_rnaseqisoform) %>% select(target_id, est_counts,tpm) %>% `colnames<-`(c("RNASeqID", "RNA2RNAseq_Counts","RNA2RNAseq_TPM"))
  
  
  # Merge all the counts with the isoforms that are considered matching 
  final <- merge(matching,IsoSeq_Defmatching ,by.x = "ref_id", by.y = "PBID", all = T)
  final <- merge(final,RNASeq_Defmatching ,by.x = "qry_id", by.y = "RNASeqID", all = T)
  final <- merge(final,IsoSeq_FL ,by.x = "ref_id", by.y = "isoform", all = T)
  
  # Correlation tests 
  ## Featurecounts
  #cor.test(final$RNA2Isoseq_Counts,final$RNA2RNAseq_Counts)
  #cor.test(final$RNA2Isoseq_Counts,final$FL)
  #cor.test(final$RNA2RNAseq_Counts,final$FL)
  
  cor.test(final$RNA2Isoseq_TPM,final$RNA2RNAseq_TPM, method = "pearson")
  cor.test(final$RNA2Isoseq_TPM,final$ISOSEQ_TPM, method = "pearson")
  cor.test(final$RNA2RNAseq_TPM,final$ISOSEQ_TPM,method = "pearson")
  
  
  p1 <- density_plot(final,"RNA2Isoseq_TPM","ISOSEQ_TPM", "RNA-Seq Transcript Expression from \n Iso-Seq-defined transcriptome (TPM)", "Iso-Seq Transcript Expression \n from Iso-Seq transcriptome (TPM)","")
  p2 <- density_plot(final,"RNA2Isoseq_TPM","RNA2RNAseq_TPM", "RNA-Seq Expression from \n Iso-Seq-defined transcriptome (TPM)", "RNA-Seq expression from RNA-Seq transcriptome (TPM)","")
  p3 <- density_plot(final,"RNA2RNAseq_TPM","ISOSEQ_TPM", "RNA-Seq expression from RNA-Seq transcriptome (TPM)", "Iso-Seq Transcript Expression \n from Iso-Seq transcriptome (TPM)","")
  
  return(list(p1,p2,p3))
}


### Isoforms Descriptive Plots #####################################################################

# no_of_isoforms_group
# Aim: Plot the range of number of isoforms for human and mouse, adult and fetal 
no_of_isoforms_group <- function(){
  isoPerGene <- lapply(class.files, function(x) SQANTI_gene_preparation(x)) %>% bind_rows()
  
  # Total Number of Genes per Sample 
  Total_Num <- isoPerGene %>% group_by(Sample) %>% count(Sample)
  # Number of isoform cateogories per Sample
  Counts <- isoPerGene %>% group_by(Sample) %>% count(nIsoCat) %>% 
    full_join(Total_Num,., by = "Sample")
  Counts$Perc <- Counts$n.y/Counts$n.x * 100
  
  # plot distribution for human vs mouse, adult vs fetal
  plot_distribution <- function(dataset1, dataset2){
    Counts$Sample <- factor(Counts$Sample, levels=c("Human (Adult)", "Human (Fetal)", "Human", "Mouse"))
    
    p <- Counts %>% 
      filter(Sample %in% c(dataset1, dataset2)) %>% 
      ggplot(., aes(x=nIsoCat, fill=Sample)) +
      geom_bar(stat="identity", aes(y= Perc, group = as.factor(Sample)), color="black", size=0.3, width=0.7, 
               position="dodge") + 
      labs(x ="Number of Isoforms", y = "Genes (%)", fill = "") +
      scale_fill_manual(values=c(label_colour(dataset1), label_colour(dataset2)), 
                        labels=c(paste(dataset1,"Cortex"), paste(dataset2,"Cortex")))  + 
      mytheme + 
      theme(legend.position = c(0.85,0.85))
    
    return(p)
  }
  
  p1 <- plot_distribution("Human","Mouse")
  p2 <- plot_distribution("Human (Adult)","Human (Fetal)")
  
  return(list(p1,p2))
  
  #write output of number of isoforms (cateogorised)
  #write.csv(Counts, paste0(output_table_dir,"/Num_Isoform_categories.csv"), quote = F)
  
}  

# no_of_isoforms_gene
# plot the number of isoforms per unique gene (not filtered by monoexonic) in each dataset
# input: the list of class.files 
no_of_isoforms_gene <- function(){
  
  # bind all datasets as one dataframe 
  Counts <- bind_rows(class.files) %>% 
    group_by(Sample,associated_gene) %>% tally() %>%
    mutate(Sample_name = paste(Sample, "Cortex"))
  
  # box-plot the number of isoforms per unique gene, all data, not filtered to monoexonic isoforms
  plot_distribution <- function(dataset1, dataset2){
    p <- Counts %>% 
      filter(Sample %in% c(dataset1, dataset2)) %>% 
      ggplot(., aes(x = Sample_name, y = log10(n), fill = Sample_name)) + 
      geom_boxplot() + theme_bw() + mytheme + 
      labs(x = "", y= "Number of Isoforms per Gene (Log10)", title = paste("\n\n\n")) + 
      theme(legend.position = "none") + 
      scale_fill_manual(values=c(alpha(label_colour(dataset1),0.3), alpha(label_colour(dataset2), 0.3))) 
    
    return(p)
  }
  
  p1 <- plot_distribution("Human","Mouse")
  p2 <- plot_distribution("Human (Adult)","Human (Fetal)")
  
  return(list(p1, p2))
}

# structural_categories_distribution <individual/merged>
# Aim: Plot percentage of isoforms classified by SQANTI2 per dataset --> sample = "individual
# Plot percentage of isoforms as merged comparison --> sample = "merged" 
structural_categories_distribution <- function(sample){
  
  if(sample == "Individual"){
    
    ## Individual Plot per Sample 
    for(i in c("Mouse", "Adult","Fetal", "Human")){
      p <- ggplot(class.files[[i]], aes(x=structural_category)) +
        geom_bar(aes(y = (..count..)/sum(..count..) * 100, alpha=coding, fill=structural_category), 
                 color="black", size=0.3, width=0.7) +
        scale_x_discrete(drop=FALSE) + scale_alpha_manual(values=c(1,0.3), 
                                                          name = "", labels = c("Protein-coding", "Non-coding"))+ 
        labs(x = "", y = "Isoforms (%)") +
        mytheme + 
        geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
        theme(axis.text.x = element_text(angle = 45)) +
        theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=margin(22,0,0,0), size=15, colour = "black")) +
        theme(legend.justification=c(1,1), legend.position=c(1,1)) + 
        guides(fill=FALSE) 
      
      print(p)
    }
  }else{
    if(sample == "Merged"){
      ## Across species
      structural_cateogory_labels <- c("FSM", "ISM", "NIC", "NNC")
      all <- bind_rows(class.files)
      df <- all[all$structural_category %in% c(structural_cateogory_labels),]
      df$structural_category <- factor(df$structural_category,
                                       labels = structural_cateogory_labels,
                                       levels = structural_cateogory_labels,
                                       ordered=TRUE)
      
      # Total Number of Transcripts
      Total_num <- df %>% group_by(Sample) %>% tally()
      
      # Count of Transcripts
      Counts <- df %>% group_by(structural_category, coding, Sample) %>%  tally()
      
      proportions <- merge(Counts,Total_num, by = "Sample")
      proportions$perc_covered <- proportions$n.x/proportions$n.y * 100
      
      # plot distribution for human vs mouse, adult vs fetal
      plot_distribution <- function(dataset1, dataset2){
        p <- proportions %>%
          filter(Sample %in% c(dataset1, dataset2)) %>%
          ggplot(., aes(x=structural_category, y = perc_covered, fill = structural_category, group = Sample)) +
          geom_bar(stat = "identity", aes(alpha=coding)) +
          facet_wrap( ~ paste(Sample, "Cortex"), nrow = 1) +    
          scale_x_discrete(drop=FALSE) +
          scale_alpha_manual(values=c(1,0.3), 
                             name = "", labels = c("Protein-coding", "Non-coding")) +
          labs(x = "", y = "Isoforms (%)") +
          mytheme +
          theme(legend.position="bottom") + 
          guides(fill=FALSE) +
          scale_fill_manual(values = wes_palette("Darjeeling2")) +
          theme(legend.position="bottom",legend.title = element_blank(),
                strip.background = element_rect(colour="white", fill="white"),
                strip.text.x = element_text(size=20,  family="ArialMT"))
        
        return(p)
      }
      
      output <- list(plot_distribution("Human","Mouse"),plot_distribution("Human (Adult)","Human (Fetal)"))
      return(output)
    }
  }
}

# exon_length_isoform_correlation <read_classification_file>
# Aim: Correlate the number of exons/length with number of isoforms (only multiexonic only)
# Mulitple plots and .txt file for the output of correlation of stats
exon_length_isoform_correlation <- function(input_class_files){
  
  # remove mono-exonic transcripts 
  dat1 <- input_class_files %>% filter(subcategory != "mono-exon")
  
  # summary_info = function to summarise the number of isoforms, etc
  # note TPM calculated with the removal of monoexonic transcripts 
  dat2 <- data.frame(summary_info(dat1))
  
  # Gene expression cut off threshold at Log10_FL_TPM at >2.5 
  dat3 <- dat2 %>% filter(Log10_FL_TPM > 2.5)
  
  sample <- input_class_files$Sample[1]
  
  # Plots 
  # P1: transcript length vs number of exons for all transcripts 
  # P2: transcrpt length vs number of exons (using only representative transcript per gene)
  # P3: same plot as P2 but filtered at expression level 
  # P4: transcript length vs number of isoforms (using only representative length) --> SUPP
  # P5: same plot as P4 but filtered at expression level --> SUPP
  # P6: exon vs number of isoforms (exon defined by representative transcript) --> SUPP
  # P7: same plot as P7 but filtered at expression level --> SUPP
  p1 <- density_plot(dat1, "length","exons","Transcript length (bp)","Number of exons", 
                     paste(sample, "Transcript length\nvs Number of Exons (all transcripts)"))
  
  p2 <- density_plot(dat2, "Max_length","Max_exon", "Transcript length (bases)", "Number of exons",
                     paste(sample, "Transcript length \nvs Number of Exons (representative transcript per gene)"))
  
  p3 <- density_plot(dat3, "Max_length","Max_exon", "Transcript length (bases)", "Number of exons",
                     paste(sample, "Transcript length\nvs Number of Exons (representative transcript per gene):\nFiltered at 25TPM"))
  
  sink(paste0(output_corr_dir,"/",sample,"_exon_length_isoform_correlation.txt"))
  # Transcript length (max representative)\nvs Number of Isoforms
  p4 <- density_plot(dat2, "Max_length","Num_of_Isoforms", 
                     "Gene Length (kb)", "Number of Isoforms", "Transcript length (max representative) vs Number of Isoforms") + 
    scale_x_continuous(labels = ks) + labs(title="\n\n\n")
  
  # Transcript length (max representative)\nvs Number of Isoforms:\nFiltered for high gene expression
  p5 <- density_plot(dat3, "Max_length","Num_of_Isoforms", 
                     "Gene Length (kb)", "Number of Isoforms",
                     " Transcript length (max representative) vs Number of Isoforms:Filtered for high gene expression") + 
    scale_x_continuous(labels = ks) + labs(title="\n\n\n")
  
  # Number of Exons (max representative)\nvs Number of Isoforms
  p6 <- density_plot(dat2, "Max_exon","Num_of_Isoforms", "Number of Exons", "Number of Isoforms","Number of Exons (max representative) vs Number of Isoforms") + labs(title="\n\n\n")
  
  # Number of Exons (max representative)\nvs Number of Isoforms:\nFiltered for high gene expression
  p7 <- density_plot(dat3, "Max_exon","Num_of_Isoforms", "Number of Exons", "Number of Isoforms","Number of Exons (max representative) vs Number of Isoforms:Filtered for high gene expression") + labs(title="\n\n\n")
  sink()
  
  
  # linear regression 
  # summary(lm(Num_of_Isoforms ~ Max_length + Max_exon, dat2))
  return(list(p1,p2,p3,p4,p5,p6,p7))
}

# gene_expression 
# Aim: plot histogram of TPM of genes (aggregating by sum)
gene_expression <- function(class.files){
  
  dat <- aggregate(class.files$ISOSEQ_TPM, by=list(Category=class.files$associated_gene), FUN=sum)
  
  # plot histogram of all genes with no threshold 
  p1 <- ggplot(dat, aes(x)) + 
    geom_histogram(binwidth = 50) + 
    labs(x = "Gene Expression (TPM)", y = "Counts") + 
    theme_bw()
  
  # plot histogram of genes with expression < 250 TPM 
  p2 <- dat %>% 
    mutate(x_new = ifelse(x > 250, 250, x)) %>% 
    ggplot(aes(x_new)) +
    geom_histogram(binwidth = .1, col = "black", fill = "blue") + 
    labs(x = "Gene Expression (TPM)", y = "Counts") + 
    theme_bw()
  
  return(list(p1, p2))
}

all_hist_peaks <- function(type){
  
  prepare_hist_breaks <- function(dat, feature){
    dat <- dat %>% .[!is.na(.[[feature]]),] # remove NAs for that specific feature 
    
    # settting the breaks for the histogram
    # threshold is for replacing >200, whic differs for each feature and dataset
    diff_max <- max(max(abs(dat[[feature]])), max(abs(dat[[feature]])))
    diff_breaks <- c(-(diff_max+1), seq(-200, 200, by = 20), diff_max+1)
    dat$diff <- cut(-(dat[[feature]]), breaks = diff_breaks) 
    threshold <- paste(formatC(diff_breaks[1], format = "g", digits = 3))
    
    # formatting of the x-axis to grab the first number, between the brackets
    dat$diff2 <- gsub("]", "", paste(dat$diff))
    dat$diff2 <- word(dat$diff2,c(2),  sep = fixed ('('))
    dat$diff2 <- word(dat$diff2,c(1),  sep = fixed (','))
    dat$diff2 <- as.factor(dat$diff2)
    
    levels(dat$diff2)[match(threshold,levels(dat$diff2))] <- ">-220"
    dat$diff2 <- factor(dat$diff2,
                        levels = c(">-220","-200","-180","-160","-140","-120",
                                   "-100","-80","-60","-40","-20","0",
                                   "20","40", "60", "80","100",
                                   "120","140","160","180","200" ))
    return(dat)
  }
  
  
  
  hist_TSS_Cage_TTS_peak <- function(dat, sample, feature,presentation){
    
    # plot 
    x_label <- 
      if(feature == "diff_to_TSS"){paste("Distance to Annotated Transcription Start Site (bp)")}else{
        if(feature == "dist_to_cage_peak"){paste("Distance to Annotated CAGE Peak (bp)")}else{
          if(feature == "diff_to_TTS"){paste("Distance to Annotated Transcription Termination Site (bp)")}else{
            if(feature == "diff_to_gene_TTS"){paste("Distance to any Annotated Transcription Termination Site (bp)")}else{
              if(feature == "diff_to_gene_TSS"){paste("Distance to any Annotated Transcription Start Site (bp)")}else{
                paste("NA")
              }
            }
          }
        }
      }
    
    
    if(presentation == "individual"){
      dat <- prepare_hist_breaks(dat, feature)
      p <- ggplot(dat, aes(x=diff2)) +
        geom_bar(aes(y = (..count..)/sum(..count..)*100), fill=label_colour(sample), color="black", size=0.3, 
                 position = position_nudge(x = 0.5)) +
        mytheme +
        labs(x = x_label, 
             y = "Isoforms (%)", 
             title = paste("\n\n")) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
      
      return(p)
    } else if (presentation == "merged"){
      human <- prepare_hist_breaks(class.files$Human, feature)
      mouse <- prepare_hist_breaks(class.files$Mouse, feature)
      
      merge <- rbind(mouse[,c("diff2","Sample")], human[,c("diff2","Sample")]) 
      p <- ggplot(merge, aes(x=diff2, fill = Sample)) +
        geom_bar(aes(y = (..count..)/sum(..count..)*100), color="black", size=0.3, 
                 position = position_dodge())+
        mytheme +
        scale_fill_manual(values=c(label_colour("Human"), label_colour("Mouse")), 
                          labels=c("Human Cortex", "Mouse Cortex"), name = "") +
        labs(x = x_label, y = "Isoforms (%)") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank(), 
              legend.position = c(0.85,0.85))
      
      return(p)
    }
  }
  
  cage_plots <- list()
  tss_plots <- list()
  tts_plots <- list()
  
  # individual plots
  count = 1
  for(i in c("Fetal","Mouse","Adult","Human")){
    
    # Subset based on input defined classification
    if(type == "All_Transcripts"){
      input_dat <- class.files[[i]]
      TSS_variable <- "diff_to_TSS"
      TTS_variable <- "diff_to_TTS"
    }else{
      if(type == "Novel_Transcripts_Annotated_Genes"){
        input_dat <- class.files[[i]] %>% .[!grepl("NOVELGENE",.$associated_gene),] %>% .[.$associated_transcript == "novel",]
        TSS_variable <- "diff_to_gene_TSS"
        TTS_variable <- "diff_to_gene_TTS"
      }else{paste("NA")}}
    
    cage_plots[[count]] <- hist_TSS_Cage_TTS_peak(input_dat,i, "dist_to_cage_peak","individual")
    tss_plots[[count]] <- hist_TSS_Cage_TTS_peak(input_dat,i, TSS_variable,"individual")
    tts_plots[[count]] <- hist_TSS_Cage_TTS_peak(input_dat,i, TTS_variable ,"individual")
    
    cage_plots_merged <- hist_TSS_Cage_TTS_peak(input_dat,"NA", "dist_to_cage_peak","merged")
    tss_plots_merged <- hist_TSS_Cage_TTS_peak(input_dat,"NA", TSS_variable,"merged")
    tts_plots_merged <- hist_TSS_Cage_TTS_peak(input_dat,"NA", TTS_variable ,"merged")
    
    count = count + 1
  }

  return(list(cage_plots,tss_plots,tts_plots, cage_plots_merged, tss_plots_merged, tts_plots_merged))
}

overlap_noveltranscripts_annotated_genes <- function(){
  annotated.class.files <- lapply(class.files, function(x) x[!grepl("NOVELGENE",x$associated_gene),])
  annotated.class.files.novel.transcripts <- lapply(annotated.class.files, function(x) x[grepl("novel",x$associated_transcript),])
  fetal.annotated.class.files.novel.transcripts <- annotated.class.files.novel.transcripts$Fetal %>% group_by(associated_gene) %>% tally()
  adult.annotated.class.files.novel.transcripts <- annotated.class.files.novel.transcripts$Adult %>% group_by(associated_gene) %>% tally()
  
  p1 <- grid.draw(venn_diagram_plot_twocircles(fetal.annotated.class.files.novel.transcripts$associated_gene, adult.annotated.class.files.novel.transcripts$associated_gene, 
                                               "Fetal(Human)", "Adult(Human)"))
  
  return(p1)
}

no_of_isoforms_sample <- function(class){
  # number of samples with detected expression of isoform
  # prerequisite: demultiplexing with numver of counts per sample
  # across each row (i.e isoform, count the number of occurences where reads are != 0)
  dat <- class %>% select(starts_with("FL.")) %>% 
    mutate(median_FL = apply(.,1, function(x) median(x)), num_samples = apply(.,1, function(x) length(x[which(x != "0")]))) 
  table(dat$num_samples)
  table(dat$num_samples)/sum(table(dat$num_samples))
  p1 <- ggplot(dat, aes(x = as.factor(num_samples))) + geom_bar(aes(y = (..count..)/sum(..count..))) + 
     scale_y_continuous(labels=scales::percent) + mytheme + labs(x = "Number of Samples", y = "Isoforms (%)")
  
  p2 <- ggplot(dat, aes(x = as.factor(num_samples), y = log(median_FL))) + geom_boxplot() + 
    mytheme + labs(x = "Number of Samples", y = "Median FL Read Count(Log10)")
  
  return(list(p1,p2))
}

### Annotated Genes Plots ##################################################################


novel_transcripts_annotated_genes_cate <- function(){
  # subset classification files (filtered) by annotated genes, and novel transcripts 
  annotated.class.files <- lapply(class.files, function(x) x[!grepl("NOVELGENE",x$associated_gene),])
  annotated.class.files.novel.transcripts <- lapply(annotated.class.files, function(x) x[grepl("novel",x$associated_transcript),])
  
  # sum number of novel transcripts from annotated genes (for downstream percentage plotting)
  sum_novel_cat <- lapply(annotated.class.files.novel.transcripts, function(x) x %>% group_by(Sample) %>% tally()) %>% bind_rows()
  
  # merge and % of categories of novel transcripts 
  novel <- annotated.class.files.novel.transcripts  %>% 
    bind_rows() %>%
    group_by(structural_category, Sample) %>%
    tally() %>%
    full_join(., sum_novel_cat, by = "Sample") %>%
    mutate(perc = n.x/n.y * 100) %>%
    rename(num_novel_transcript_type = n.x) %>% 
    rename(total_num_novel_transcript = n.y)
  
  # plot the percentage of categories from novel transcripts 
  plot_distribution <- function(dataset1, dataset2){
    p <- novel %>% 
      filter(Sample %in% c(dataset1, dataset2)) %>% 
      ggplot(., aes(x = structural_category, y = perc, fill = Sample)) + 
      geom_bar(stat = "identity", position = position_dodge()) + 
      theme_bw() + mytheme + 
      labs(x = "", y = "Novel Transcripts of Annotated Genes (%)", title = "\n\n") + 
      theme(legend.justification=c(1,1), legend.position=c(0.95,0.95), legend.title = element_blank()) + 
      scale_fill_manual(values=c(label_colour(dataset1), label_colour(dataset2)), 
                        labels=c(paste(dataset1,"Cortex"), paste(dataset2, "Cortex")), name = "") 
    
    return(p)
  }
  # output percentages for values in paper 
  #sink(paste0(output_plot_dir, "/num_novel_transcripts_annotated_genes_cate.txt"))
  #print(novel)
  #sink()
  
  p1 <- plot_distribution("Human","Mouse")
  p2 <- plot_distribution("Human (Adult)","Human (Fetal)")
  return(list(p1,p2))
}

subset_feature <- function(col_name_feature, ylabel, category){
  
  # subset classification files (filtered) by annotated genes, and novel transcripts 
  annotated.class.files <- lapply(class.files, function(x) x[!grepl("NOVELGENE",x$associated_gene),])
  annotated.class.files.novel.transcripts <- lapply(annotated.class.files, 
                                                    function(x) x[grepl("novel",x$associated_transcript),] %>% 
                                                      filter(subcategory != "mono-exon"))
  annotated.class.files.annotated.transcripts <- lapply(annotated.class.files, 
                                                        function(x) x[!grepl("novel",x$associated_transcript),] %>% 
                                                          filter(subcategory != "mono-exon"))
  
  subset_transcripts <- function(dataset, type){
    if(dataset == "Novel"){
      dat <- lapply(annotated.class.files.novel.transcripts, function(x) x %>% filter(structural_category == type))
    } else {
      dat <- lapply(annotated.class.files.annotated.transcripts, function(x) x %>% filter(structural_category == type))
    }
    return(dat)
  }
  
  #annotated.class.files.NIC.transcripts <- lapply(annotated.class.files.novel.transcripts, function(x) x %>% filter(structural_category == "NIC"))
  #annotated.class.files.NNC.transcripts <- lapply(annotated.class.files.novel.transcripts, function(x) x %>% filter(structural_category == "NNC"))
  # annotated.class.files.Fusion.transcripts <- lapply(annotated.class.files.novel.transcripts, function(x) x %>% filter(structural_category == "Fusion"))
  
  
  # Subset length from classification file (bind_rows as input is list of dataframes)
  extract_feature <- function(dat, dat_name){
    dat1 <- dat %>% bind_rows() %>% 
      mutate(Transcripts = dat_name) %>% .[,c(col_name_feature, "Transcripts", "Sample")]
  }
  
  if(category == "all_novel") { 
    set1 <- extract_feature(annotated.class.files.annotated.transcripts, "Known")
    set2 <- extract_feature(annotated.class.files.novel.transcripts, "Novel") 
    merge <- bind_rows(list(set1, set2))
  } else { 
    set1 <- extract_feature(subset_transcripts("Known", "FSM"), "FSM (Known)")
    set2 <- extract_feature(subset_transcripts("Known", "ISM"), "ISM (Known)")
    set3 <- extract_feature(subset_transcripts("Novel", "NIC"), "NIC (Novel)")
    set4 <- extract_feature(subset_transcripts("Novel", "NNC"), "NNC (Novel)")
    set5 <- extract_feature(subset_transcripts("Novel", "Fusion"), "Fusion (Novel)")
    merge <- bind_rows(list(set1, set2, set3, set4, set5)) 
    merge$Transcripts <- factor(merge$Transcripts, levels = c("FSM (Known)", "ISM (Known)", "NIC (Novel)", "NNC (Novel)","Fusion (Novel)"))
  }
  
  plot_distribution <- function(dataset1, dataset2){
    y.var <- rlang::sym(quo_name(enquo(col_name_feature)))
    p <- merge %>% filter(Sample %in% c(dataset1, dataset2)) %>%  
      ggplot(., aes(x = Sample, y = !! y.var, fill= Transcripts)) + 
      geom_boxplot() + 
      theme_bw() + 
      mytheme + 
      labs(x = "", y= ylabel, fill = "Isoforms", title = "\n\n") + 
      theme(legend.justification=c(1,1), legend.position=c(0.95,0.95),legend.title = element_blank(),
            legend.direction = "horizontal") + 
      scale_x_discrete(labels = c(paste(dataset1,"Cortex"),paste(dataset2,"Cortex"))) + 
      guides(fill=guide_legend(nrow=3,byrow=TRUE))
    
    return(p)
  }
  
  p1 <- plot_distribution("Human","Mouse")
  p2 <- plot_distribution("Human (Adult)","Human (Fetal)")
  
  output <- list(p1,p2, merge)
  names(output) <- c("human_mouse_plot","adult_fetal_plot","merge")
  return(output)
}

novel_annotated_transcript_expression <- function(corr_output_name){
  output <- subset_feature("Log_ISOSEQ_TPM", "Iso-Seq Expression of Annotated Genes (Log10 TPM)", "all_novel")
  split_output <- subset_feature("Log_ISOSEQ_TPM", "Iso-Seq Expression of Annotated Genes (Log10 TPM)", "all_split")
  
  # mann whitney test
  sink(paste0(output_corr_dir, "/", corr_output_name,".txt"))
  for(i in unique(factor(output$merge$Sample))){
    cat("Mann Whitney for Expression of Novel vs Annotated in Annotated Genes of dataset:", i,"\n")
    dat <- output$merge[output$merge$Sample == paste(i),]
    test <- wilcox.test(Log_ISOSEQ_TPM ~ Transcripts,data=dat)
    print(test)
    cat("Exact p value:", test$p.value,"\n")
  }
  sink()
  
  return(list(output$human_mouse_plot, output$adult_fetal_plot, split_output$human_mouse_plot, split_output$adult_fetal_plot))
}

novel_annotated_transcript_length <- function(corr_output_name){
  output <- subset_feature("length", "Transcript Length (kB)", "all_novel")
  split_output <- subset_feature("length", "Transcript Length (kB)", "all_split")
  
  sink(paste0(output_corr_dir, "/", corr_output_name))
  for(i in unique(factor(output$merge$Sample))){
    cat("Mann Whitney for length of Novel vs Annotated in Annotated Genes of dataset:", i,"\n")
    dat <- output$merge[output$merge$Sample == paste(i),]
    test <- wilcox.test(length ~ Transcripts,data=dat)
    print(test)
    cat("Exact p value:", test$p.value,"\n")
  }
  sink()
  
  # plot 
  return(list(output$human_mouse_plot + scale_y_continuous(labels = ks), 
              output$adult_fetal_plot + scale_y_continuous(labels = ks), 
              split_output$human_mouse_plot + scale_y_continuous(labels = ks), 
              split_output$adult_fetal_plot + scale_y_continuous(labels = ks)))
}

novel_annotated_transcript_exon <- function(corr_output_name){
  output <- subset_feature("exons", "Number of Exons", "all_novel")
  split_output <- subset_feature("exons", "Number of Exons", "all_split")
  
  sink(paste0(output_corr_dir, "/", corr_output_name))
  for(i in unique(factor(output$merge$Sample))){
    cat("Mann Whitney for Number of Exons of Novel vs Annotated in Annotated Genes of dataset:", i,"\n")
    dat <- output$merge[output$merge$Sample == paste(i),]
    test <- wilcox.test(exons ~ Transcripts,data=dat)
    print(test)
    cat("Exact p value:", test$p.value,"\n")
  }
  sink()
  
  return(list(output$human_mouse_plot, output$adult_fetal_plot, split_output$human_mouse_plot, split_output$adult_fetal_plot))
}

# splice_site_frequency
# Aim: plot the number of splice sites of annotated genes 
# Output: plot and stats (saved in list)
splice_site_frequency <- function(){
  subset_junctions <- function(type_class_file,type_junc_file, dataset){
    # subset by annotated genes only 
    annotated_genes <- type_class_file[!grepl("NOVEL", type_class_file$associated_gene),]
    
    # from junction file, subset only with PBID of annotated genes
    output_junc_file <- type_junc_file[type_junc_file$isoform %in% annotated_genes$isoform,] 
    output_junc_file$Sample <- dataset
    return(output_junc_file)
  }
  
  # QC: the number of junctions extracted are of those novel transcripts (annotated genes only), exact same number of transcripts
  #length(unique(junc.files$Mouse[junc.files$Mouse$isoform %in% annotated_genes_novel_transcripts$isoform,"isoform"]))
  #annotated_genes_novel_transcripts %>% nrow()
  
  # QC the number of junctions extracted are of those multiexonic transcripts (annotated genes only)
  #length(unique(junc.files$Mouse[junc.files$Mouse$isoform %in% annotated_genes$isoform,"isoform"]))
  #annotated_genes %>% filter(subcategory != "mono-exon") %>% nrow()
  
  
  count <- bind_rows(subset_junctions(class.files$Mouse, junc.files$Mouse, "Mouse"), 
                     subset_junctions(class.files$Fetal, junc.files$Fetal, "Fetal"))
  
  
  count$canonical[count$canonical == "canonical"] <- "Canonical"
  count$canonical[count$canonical == "non_canonical"] <- "Non-canonical"
  
  p1 <- count %>% 
    group_by(splice_site, canonical, Sample) %>% 
    count() %>% 
    ggplot(., aes(x = reorder(splice_site, -n), y = log10(n), group = canonical, colour = Sample)) + 
    geom_point(stat = "identity") + 
    facet_grid(~ canonical,scales="free_x",space = "free") + 
    theme_bw() + 
    labs(x = "", y = "Number of Splice Sites (Log10)", title = "\n\n") + 
    mytheme +  
    scale_color_manual(values=c(label_colour("Human"), label_colour("Mouse")),
                       name="",
                       labels=c( "Human Cortex", "Mouse Cortex")) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = c(.90, 0.85),
          strip.background = element_rect(colour="white", fill="white"),
          strip.text.x = element_text(size=20,  family="ArialMT"))
  
  
  # summarise the percentage of each splice site across all splice sites per dataset
  stats <- count %>% 
    group_by(splice_site, canonical, Sample) %>% 
    count()
  
  stats[stats$splice_site %in% c("GTAG","GCAG","ATAC"),]
  
  total_stats <- aggregate(stats$n, by=list(Sample=stats$Sample), FUN=sum)
  stats <- stats %>% left_join(., total_stats, by = "Sample") %>%
    mutate(perc = round(n/x * 100,3))
  
  stats_canonical <- count %>% 
    group_by(canonical, Sample) %>% 
    count() %>%
    left_join(., total_stats, by = "Sample") %>%
    mutate(perc = round(n/x * 100,3))
  
  # Support of noncanonical splice sites by total RNASeq coverage
  p2 <- count %>% 
    .[.$canonical == "Non-canonical",c("splice_site", "total_coverage", "Sample")] %>% 
    aggregate(total_coverage ~ splice_site + Sample, ., sum) %>%
    ggplot(., aes(x = reorder(splice_site, - total_coverage), y = log10(total_coverage), color = Sample)) +
    geom_jitter() + 
    mytheme + 
    labs(x = "Non-canonical Splice Sites", y = "RNA-Seq Coverage (Log10)", title = "\n\n") + 
    theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = c(.90, 0.85)) + 
    scale_color_manual(values=c(label_colour("Human"), label_colour("Mouse")),
                       name="",
                       labels=c("Human Cortex", "Mouse Cortex"))
  
  count_supported <- count %>% 
    .[.$canonical == "Non-canonical",c("splice_site", "total_coverage", "Sample")] %>% 
    aggregate(total_coverage ~ splice_site + Sample, ., sum) 
  
  cat("Number of non canonical splice sites supported by RNA-Seq:", length(unique(count_supported[count_supported$total_coverage != 0,"splice_site"])),"\n")
  cat("Number of non canonical splice sites in Mouse supported by RNA-Seq:", 
      length(unique(count_supported[count_supported$Sample == "Mouse" & count_supported$total_coverage != 0,"splice_site"])),"\n")
  cat("Number of non canonical splice sites in Fetal supported by RNA-Seq:", 
      length(unique(count_supported[count_supported$Sample == "Fetal" & count_supported$total_coverage != 0,"splice_site"])),"\n")
  cat("Number of common non canonical splice sites in Mouse and Fetal supported by RNA-Seq:", 
      length(intersect(unique(count_supported[count_supported$Sample == "Mouse","splice_site"]),
                       unique(count_supported[count_supported$Sample == "Fetal","splice_site"]))),"\n")
  
  
  
  
  # tally the number of samples and occurence of that for non_canonical splice site
  dat <- count %>% 
    .[.$canonical == "Non-canonical",c("splice_site", "sample_with_cov", "Sample")] %>% 
    group_by(splice_site, sample_with_cov, Sample) %>% 
    tally()
  
 
  # lowly expressed non_canonical junctions
  # dat[dat$splice_site %in% c("TCCT", "GGGC", "GTGC", "TGGC"),]
  
  
  return(list(p1,p2,stats))
}

### Intron Retention ##################################################################
# venn_diagram_IR
# Input: homologues genes between human and mouse (already read),
# List of genes with intron retention generated from Intron_Retention.Rmd (pre-requisite)
venn_diagram_IR <- function(){
  homologs_genes_filtered <- homologs_genes_filter()
  
  # Intron retention transcripts
  ALL_IR <- read.csv(paste0(output_table_dir,"/AS_IR/Intron_Retention.csv"))
  human <- ALL_IR %>% select(associated_gene, Human) %>% filter(!is.na(Human))
  mouse <- ALL_IR %>% select(associated_gene, Mouse) %>% filter(!is.na(Mouse))
  
  human_homologs <- human[human$associated_gene %in% homologs_genes_filtered$human,]
  mouse_homologs <- mouse[mouse$associated_gene %in% homologs_genes_filtered$mouse,]
  
  ################ human only (no need for considering homologues)
  fetal <- ALL_IR %>% select(associated_gene, Fetal) %>% filter(!is.na(Fetal))
  adult <- ALL_IR %>% select(associated_gene, Adult) %>% filter(!is.na(Adult))
  
  dat <- data.frame()
  count <- 1
  for(i in list(human_homologs, mouse_homologs, fetal, adult)){
    dat[1, count] <- sum(i[[2]])
    dat[2, count] <- length(unique(i$associated_gene))
    count <- count + 1
  }
  rownames(dat) <- c("Number of IR transcripts", "Number of Genes associated with IR transcripts")
  colnames(dat) <- c("human_homologs", "mouse_homologs", "fetal", "adult")
  print(dat)
  
  #sink(paste0(output_plot_dir,"/venndiagram_IRtests"), append = TRUE)
  cat("Statstical Tests for number of IR transcripts in adult vs fetal:\n")
  dat <- data.frame("Fetal" = c(class.files$Fetal %>% filter(subcategory == "intron_retention") %>% nrow(),
                                class.files$Fetal %>% filter(subcategory != "intron_retention") %>% nrow()),
                    "Adult" = c(class.files$Adult %>% filter(subcategory == "intron_retention") %>% nrow(),
                                class.files$Adult %>% filter(subcategory != "intron_retention") %>% nrow()))
  rownames(dat) <- c("IR","Non-IR")                  
  
  f.test <- fisher.test(dat)
  print(f.test)
  print(f.test$p.value)
  
  
  output <- list(human, mouse, human_homologs, mouse_homologs, fetal, adult)
  names(output) <- c("human","mouse","human_homologs","mouse_homologs","fetal","adult")
  return(output)
}

IR_rnaseq_support <- function(){
  # using fetal classification file with input of fetal RNASeq rather than intropolis 
  # note therfore different number of filtered transcripts (fetal RNASeq: 3077 IR transcripts, intropolis: 3053 IR transcripts)
  own_fetal_class <- read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Fetal/SQANTI2/own_rnaseq/fetalFL.sample.collapsed.filtered.rep_classification.filtered_lite_classification.txt", as.is = T, sep = "\t", header = T)
  own_fetal_class$Sample <- "Human(Fetal)"
  
  # select for IR transcripts
  df <- bind_rows(class.files$Mouse, own_fetal_class)
  df_IR <- df[df$subcategory == "intron_retention",]
  
  # Coverage is defined by > 0 TPM
  df_IR <- df_IR %>% mutate(IR_iso_esp = ifelse(.$iso_exp == 0., "No RNASeq Coverage", "RNASeq Coverage"))
  
  # stats 
  stats <- df_IR %>% group_by(Sample, IR_iso_esp) %>% tally() %>% 
    group_by(Sample) %>% mutate(perc = n/sum(n))
  print(stats)
  
  # plot supported by RNASeq Isoform Expression 
  p <- ggplot(stats, aes(x = Sample, y = perc, fill = IR_iso_esp)) + 
    geom_bar(stat = "identity") + 
    mytheme + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
    labs(x = "", y = "RNASeq Isoform Expression of IR-transcripts") + 
    theme(legend.position = "bottom", legend.title = element_blank()) 
  
  return(p)
}

# IR_rate 
# Aim: to plot the number of IR transcripts associated with highly expressed vs lowly expressed genes 
# Output: 
# Plot1: IR-transcript isoseq expression across each dataset
# Plot2: Proportion of IR transcripts across gene expression
# Plot2b: Proportion of genes out of total genes with IR transcripts
# Plot3: Proportion of IR transcripts using gene expression cut off (previously estabilshed)
# Plot4: Density of gene expression of genes with IR transcripts
# Plot5: IsoSeq Gene expression with IR transcripts between adult vs fetal
IR_rate <- function(x){ 
  # Filter each dataset with only transcipts classified as intron retention, and bind to one big dataframe 
  IR_dataset <- lapply(class.files, function(x) x[x$subcategory == "intron_retention",])
  IR <- bind_rows(IR_dataset)
  
  ###
  # Proportion of IR transcripts across gene expression 
  # Hypothesis: highly expressed genes contain fewer IR transcripts
  ###
  # Determine gene expression for each gene (with removal of monoexonic transcripts)
  # note TPM calculated with the removal of monoexonic transcripts 
  multiexonic <- lapply(class.files, function(x) x %>% filter(subcategory != "mono-exon"))
  summary_info_dataset <- lapply(multiexonic, function(x) data.frame(summary_info(x)))
  
  # Create function
  # Merge for each gene, the gene expression and number of IR transcripts associated 
  merge_expression_num_IR <- function(IR_dataset_input, summary_info_dataset_input){
    # in the IR classification file, count the number of IR transcripts for each gene
    # merge with expression for that gene
    output <- IR_dataset_input %>% group_by(associated_gene) %>% count() %>% 
      left_join(.,summary_info_dataset_input[,c("associated_gene","Log10_FL_TPM")], by = "associated_gene")
    
    # classify the gene expression into three classes for plot
    for(i in 1:nrow(output)){
      output$FL_TPM_cate[i] <- if ( output$Log10_FL_TPM[i] >= 0 & output$Log10_FL_TPM[i] <= 2) {
        "< 2"
      } else if ( output$Log10_FL_TPM[i] > 2 & output$Log10_FL_TPM[i] <= 3) {
        "2 - 3"
      } else {
        "3 <"
      }
    }
    
    output$Sample <- IR_dataset_input$Sample[1]
    return(output)
  }
  
  
  # Apply function to each dataset to merge the number of associated transcripts and expression per gene
  IR_gene_expression <- list(
    merge_expression_num_IR(IR_dataset$Adult, summary_info_dataset$Adult),
    merge_expression_num_IR(IR_dataset$Fetal, summary_info_dataset$Fetal),
    merge_expression_num_IR(IR_dataset$Human, summary_info_dataset$Human),
    merge_expression_num_IR(IR_dataset$Mouse, summary_info_dataset$Mouse)
  )
  names(IR_gene_expression) <- c("Adult","Fetal","Human","Mouse")
  IR_gene_expression_bind <- bind_rows(IR_gene_expression)
  
  # QC: check the number of genes in IR gene expression = the number of genes with IR transcripts 
  # length(unique(class.files$Adult[class.files$Adult$subcategory == "intron_retention","associated_gene"]))
  # IR_gene_expression_bind[IR_gene_expression_bind$Sample == "Adult",] %>% nrow()
  
  # Count the number of IR transcripts per dataset for percentage for plot
  total <- aggregate(IR_gene_expression_bind$n, by=list(Sample=IR_gene_expression_bind$Sample), FUN=sum)
  colnames(total)[2] <- "total_num"
  
  class.files.bind <- bind_rows(class.files)
  total_genes <- lapply(class.files, function(x) length(unique(x$associated_gene))) %>% bind_rows() 
  total_genes <- as.data.frame(gather(total_genes, Sample, total_num_genes, Mouse:Human, factor_key=TRUE))
  total_genes$Sample <- as.character(total_genes$Sample)
  
  
  ###
  # Proportion of IR transcripts using gene expression cut off (previously estabilshed)
  # Hypothesis: highly expressed genes contain fewer IR transcripts
  ###
  # tally the number of genes with IR transcripts in gene expression with Log10_FL_TPM > 2.5 and < 2.5
  lowly_expressed <- lapply(IR_gene_expression, function(x) x[x$Log10_FL_TPM < 2.5,"associated_gene"] %>% nrow()) %>% bind_rows() %>% gather(., Sample, total_num_genes, Adult:Mouse, factor_key=TRUE)
  
  highly_expressed <- lapply(IR_gene_expression, function(x) x[x$Log10_FL_TPM > 2.5,"associated_gene"] %>% nrow()) %>% bind_rows() %>% gather(., Sample, total_num_genes, Adult:Mouse, factor_key=TRUE)
  
  
  # Plot
  # merge total_num.genes.x == highly expressed, total_num_genes.y = lowly expressed genes
  # table for stats before plotting
  merge <- merge(highly_expressed,lowly_expressed, by = "Sample") %>% melt(.)
  colnames(merge)[2] <- "gene_expression"
  colnames(merge)[3] <- "num_genes_IR_transcript"
  print(merge)
  
  # Plot Gene expression categories and the number of IR transcripts for that category  
  # merge IR-transcript expression and total per dataset to calcualte the percentage 
  plot_distribution_A <- function(dataset1, dataset2){
    
    pA <- aggregate(.~FL_TPM_cate+Sample, IR_gene_expression_bind[,c("n","FL_TPM_cate","Sample")], sum) %>% 
      left_join(., total, by = "Sample") %>% 
      filter(Sample %in% c(dataset1, dataset2)) %>% 
      mutate(FL_TPM_cate = factor(FL_TPM_cate,levels = c("< 2", "2 - 3", "3 <"))) %>%
      mutate(perc = n/total_num * 100) %>% 
      ggplot(., aes(x = FL_TPM_cate, y =perc, colour = Sample, group = Sample)) + geom_point() + 
      geom_line() + theme_bw() +  mytheme + 
      labs(x = "Iso-Seq Gene Expression (Log10 TPM)", y = "IR-Transcripts (%)", title = "\n\n") + 
      scale_color_manual(values=c(label_colour(dataset1), label_colour(dataset2)),
                         labels=c(paste(dataset1,"Cortex"), paste(dataset2, "Cortex")), name = "") +
      theme(legend.title = element_blank(), legend.position = c(0.85,0.85))
    
    
    return(pA)
  }
  
  plot_distribution_B <- function(dataset1, dataset2){
    pB <- merge %>% 
      filter(Sample %in% c(dataset1, dataset2)) %>% 
      left_join(., total_genes, by = "Sample") %>% 
      mutate(perc = num_genes_IR_transcript/total_num_genes * 100) %>% 
      ggplot(., aes(x = Sample, y = perc, fill = gene_expression)) + 
      geom_bar(stat = "identity") + 
      labs(x = "", y = "Genes with IR-Transcripts (%)", title = "\n\n") + 
      theme_bw() + 
      mytheme +  
      theme(legend.position = "top") + 
      scale_fill_manual(values = c("#FC4E07","#E7B800"),
                        name="Iso-Seq Gene Expression",
                        labels=c("High (>2.5 Log10 TPM)", "Low (<2.5 Log10 TPM)")) 
    
    return(pB)
  }
  
  p1 <- plot_distribution_A("Human","Mouse")
  p2 <- plot_distribution_A("Human (Adult)","Human (Fetal)")
  p3 <- plot_distribution_B("Human","Mouse") + scale_x_discrete(labels = c('Human Cortex','Mouse Cortex'))
  p4 <- plot_distribution_B("Adult","Fetal") + scale_x_discrete(labels = c('Human (Adult) Cortex','Human (Fetal) Cortex'))
  # Statistical analysis for intron retention expression (at gene level) between adult vs fetal
  #sink(paste0(output_plot_dir, "/Mann_Whitney_adult_fetal_IRexpression"))
  cat("\n\nMann Whitney for intron rentetion rate in gene expression between adult vs fetal:")
  test <- wilcox.test(IR_gene_expression$Adult$Log10_FL_TPM,IR_gene_expression$Fetal$Log10_FL_TPM)
  print(test)
  cat("Exact p value:", test$p.value)
  #sink()
  
  # QC to check the % of genes with IR transcripts refer to that draw in plot p3.b
  #length(unique(class.files$Mouse[class.files$Mouse$subcategory == "intron_retention","associated_gene"]))/length(unique(class.files$Mouse$associated_gene))
  
  ###
  # Density of expression of genes with IR transcripts
  # Hypothesis: highly expressed genes contain fewer IR transcripts
  ###
  
  p5 <- rbind(data.frame("Num" = IR_gene_expression$Adult$Log10_FL_TPM, "Sample" = "Adult"),
              data.frame("Num" = IR_gene_expression$Fetal$Log10_FL_TPM, "Sample" = "Fetal")) %>% 
    ggplot(., aes(y = Num, x = Sample)) + geom_boxplot() + 
    theme_bw() + 
    labs(y = "Iso-Seq Gene Expression (Log10 TPM)", x = "") +
    mytheme + 
    theme(legend.position = c(.90, 0.90), legend.title = element_blank()) #+
  #scale_fill_manual(values=c(label_colour("Adult"),label_colour("Fetal")),
  # name="",
  # labels=c("Human (Adult) Cortex", "Human (Fetal) Cortex"))
  
  
  ###
  # IsoSeq gene expression with IR transcripts between Fetal and Adult 
  ###
  subset_expression_rename <- function(dat,sample){dat1 <- dat[,c("associated_gene","Log10_FL_TPM")];colnames(dat1)[2] <- sample; return(dat1)}
  IR_gene_expression_Fetal_Adult <- merge(subset_expression_rename(IR_gene_expression$Adult, "Adult"), 
                                          subset_expression_rename(IR_gene_expression$Fetal, "Fetal"), all = TRUE) %>% 
    mutate(Adult_TPM = 10^.$Adult) %>% 
    mutate(Fetal_TPM = 10^.$Fetal) %>%
    mutate(diff = abs(.$Adult_TPM - .$Fetal_TPM)) 
  
  # only common gene with IR in fetal and adult 
  p6 <- IR_gene_expression_Fetal_Adult [complete.cases(IR_gene_expression_Fetal_Adult ), ] %>% 
    density_plot(., "Fetal", "Adult", "Iso-Seq Gene Expression with IR-transcripts (Log10 TPM) - Fetal","Iso-Seq Gene Expression with IR-transcripts (Log10 TPM) - Fetal","Fetal vs Adult")
  
  # write.csv absolute differences of gene expression with IR between fetal and adult 
  write.csv(IR_gene_expression_Fetal_Adult, paste0(output_table_dir,"/IR_gene_expression_Fetal_Adult.csv"), quote = F)
  
  return(list(p1,p2,p3,p4,p5,p6))
}



### Alternative Splicing events ##################################################################
# human_mouse_suppa2
# Aim: plots (venn diagram and bar plot) of common genes with spliced events between human and mouse
# Prerequisite: input homologues genes from Aaron between human and mouse 
# input: type == "homology" 
human_mouse_suppa2 <- function(type){
  
  # homologs_genes
  homologs_genes_filtered <- homologs_genes_filter()
  
  # Input SUPPA2 output from SUPPA2_Output.Rmd of a table of list of genes with splicing events 
  SUPPA2_Genes <- read.csv(paste0(output_table_dir, "/AS_IR/ALL_SUPPA2_Genes_Output_Updated.csv"))
  
  # OLD 
  #SUPPA2_Genes <- read.table(paste0(output_table_dir, "/ALL_SUPPA2_Genes_Output.txt"))
  
  # Subset on dataset, of genes in each dataset with each event (therefore note replicated gene name) 
  Human_SUPPA2 <- SUPPA2_Genes[SUPPA2_Genes$Sample == "Human",] %>% .[,c(-1)]
  Mouse_SUPPA2 <- SUPPA2_Genes[SUPPA2_Genes$Sample == "Mouse",] %>% .[,c(-1)]
  Mouse_SUPPA2$Gene <- toupper(Mouse_SUPPA2$Gene)
  
  ########################## Plots by removing non-homologous genes (no conversion of names )
  #### Plot of genes and splicing events separated by two datasets
  merge_human_mouse <- bind_rows(list(Human_SUPPA2,Mouse_SUPPA2)) %>% spread(., Sample, No_of_Events)
  cat("Number of unique associated genes with splicing events in Human:", 
      length(unique(merge_human_mouse[!is.na(merge_human_mouse$Human),"Gene"])), "\n")
  cat("Number of unique associated genes with splicing events in Mouse:", 
      length(unique(merge_human_mouse[!is.na(merge_human_mouse$Mouse),"Gene"])), "\n")
  
  # Detection based on splice events on two datasets
  for(row in 1:nrow(merge_human_mouse)){
    merge_human_mouse$Detected[row] <- 
      if(is.na(merge_human_mouse[row,"Human"])){paste("Mouse Only")
      }else{
        if(is.na(merge_human_mouse[row,"Mouse"])){paste("Human Only")
        }else{paste("Common")}}}
  
  
  ##### Venn Diagram of unique genes with splicing events
  p1 <- venn_diagram_plot_twocircles(Human_SUPPA2$Gene,Mouse_SUPPA2$Gene,"Human","Mouse")
  
  p2 <- merge_human_mouse %>% group_by(Event, Detected) %>% tally() %>%
    ggplot(., aes(x = Event, fill = Detected, y = n)) + 
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"), 
                      labels=c("Common", "Human Only", "Mouse Only")) + 
    theme_bw() + 
    mytheme + 
    labs(x = "Splicing Events", y = "Number of Genes") + 
    theme(legend.position = c(.15, 0.85))
  
  
  ########################## Plots by converting mouse ogous gene names to human 
  cat("Working with homologous genes; Converting mouse genes to human equivalent genes", "\n")
  homologs_genes_filtered <- homologs_genes_filter()
  
  # Convert mouse gene to human homology gene before matching  
  for(i in 1:nrow(Mouse_SUPPA2)){
    Mouse_SUPPA2$human_homolog[i] <- 
      ifelse(Mouse_SUPPA2$Gene[i] %in% homologs_genes_filtered$mouse,
             paste(as.character(homologs_genes_filtered$human[which(homologs_genes_filtered$mouse == 
                                                                      Mouse_SUPPA2$Gene[i])])), 
             Mouse_SUPPA2$Gene[i] )
  }
  
  
  # Mouse_SUPPA2[Mouse_SUPPA2$Gene == "GPI1",] # to test accurate conversion 
  # datawrangle to dataframe from conversion for rowbinding with human genes 
  Mouse_SUPPA2_converted <- Mouse_SUPPA2 %>% select(human_homolog, No_of_Events, Event, Sample)
  colnames(Mouse_SUPPA2_converted)[1] <- "Gene" 
  
  merge_human_mouse_converted <- bind_rows(list(Human_SUPPA2, Mouse_SUPPA2_converted)) %>% 
    group_by_at(vars(-No_of_Events)) %>% # https://github.com/tidyverse/tidyr/issues/426 (issue otherwise)
    mutate(row_id=1:n()) %>% ungroup() %>% 
    spread(., Sample, No_of_Events)
  
  # Detection based on splice events on two datasets
  merge_human_mouse_converted$Detected <- "NA"
  for(row in 1:nrow(merge_human_mouse_converted)){
    merge_human_mouse_converted$Detected[row] <- 
      if(is.na(merge_human_mouse_converted[row,"Human"])){paste("Mouse Only")
      }else{
        if(is.na(merge_human_mouse_converted[row,"Mouse"])){paste("Human Only")
        }else{paste("Common")}}}
  
  ##### Venn Diagram of unique genes with splicing events
  Human_SUPPA2_homologues <- merge_human_mouse_converted[merge_human_mouse_converted$Detected %in% c("Human Only", "Common"), ] 
  Mouse_SUPPA2_homologues <- merge_human_mouse_converted[merge_human_mouse_converted$Detected %in% c("Mouse Only", "Common"), ] 
  
  p4 <- venn_diagram_plot_twocircles(Human_SUPPA2_homologues$Gene, 
                                     Mouse_SUPPA2_homologues$Gene, "Human", "Mouse")
  
  p3 <- merge_human_mouse_converted %>% 
    group_by(Event, Detected) %>% tally() %>%
    ggplot(., aes(x = Event, fill = Detected, y = n)) + 
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values=c("#999999", label_colour("Human"), label_colour("Mouse")), 
                      labels=c("Common", "Human Only", "Mouse Only")) + 
    theme_bw() + 
    mytheme + 
    labs(x = "Splicing Events", y = "Number of Genes", title = "\n\n") + 
    theme(legend.position = c(.90, 0.90),legend.title = element_blank())
  
  # unique(Human_SUPPA2$Gene[Human_SUPPA2$Gene%in% Mouse_SUPPA2_homologues_match$Gene])
  if(type == "not_homology"){return(list(p1,p2))}
  if(type == "homology"){return(list(p4,p3))}
  
  
  number_differences <- function(){
    
    p5 <- venn.diagram(
      x = list(Human_SUPPA2_homologues$Gene, Mouse_SUPPA2_converted_homologues$Gene, Mouse_SUPPA2_homologues$Gene),
      category.names = c("Human_SUPPA2_homologues" , "Mouse_SUPPA2_converted_homologues", "Mouse_SUPPA2_homologues"),
      filename = NULL,
      output=TRUE)
    
    ncommon_human_mouse <- unique(Mouse_SUPPA2_converted_homologues[!Mouse_SUPPA2_converted_homologues$Gene %in% Human_SUPPA2_homologues$Gene, "Gene"])
    ncommon_mouse <- unique(Mouse_SUPPA2_converted_homologues[!Mouse_SUPPA2_converted_homologues$Gene %in% Mouse_SUPPA2_homologues$Gene, "Gene"])
    n3376 <- unique(Mouse_SUPPA2_homologues[Mouse_SUPPA2_homologues$Gene %in% Human_SUPPA2_homologues$Gene, "Gene"])
    n2422 <- setdiff(ncommon_human_mouse, ncommon_mouse)
    n104 <- setdiff(ncommon_mouse, ncommon_human_mouse)
    combined_common <- c(n2422, n104, n3376)
    n183 <- unique(Mouse_SUPPA2_converted_homologues[!Mouse_SUPPA2_converted_homologues$Gene %in% combined_common,"Gene"])
    
    hg38[hg38$Gene %in% n183,]
    mm10[mm10$Gene %in% n183,]
  }
}

# fetal_adult_suppa2 
# Aim: plots (venn diagram and bar plot) of common genes with spliced events between fetal and adult  
fetal_adult_suppa2 <- function(){
  
  # Input SUPPA2 output from SUPPA2_Output.Rmd of a table of list of genes with splicing events 
  SUPPA2_Genes <- read.csv(paste0(output_table_dir, "/AS_IR/ALL_SUPPA2_Genes_Output_Updated.csv"))
  
  # Subset on dataset 
  Adult_SUPPA2 <- SUPPA2_Genes[SUPPA2_Genes$Sample == "Human (Adult)",] %>% .[,c(2:5)]
  Fetal_SUPPA2 <- SUPPA2_Genes[SUPPA2_Genes$Sample == "Human (Fetal)",] %>% .[,c(2:5)]
  
  ##### Venn Diagram of unique genes with splicing events (3919)
  p1 <- venn_diagram_plot_twocircles(Adult_SUPPA2$Gene, Fetal_SUPPA2$Gene,"Human (Adult)","Human (Fetal)")
  Adult_Unique <- setdiff(Adult_SUPPA2$Gene, Fetal_SUPPA2$Gene)
  Fetal_Unique <- setdiff(Fetal_SUPPA2$Gene,Adult_SUPPA2$Gene)
  write.table(Adult_Unique, paste0(output_table_dir, "/Adult_Only_AS.txt"), col.names=FALSE, row.names=FALSE,quote = FALSE)
  write.table(Fetal_Unique, paste0(output_table_dir, "/Fetal_Only_AS.txt"), col.names=FALSE, row.names=FALSE,quote = FALSE)
  
  cat("Number of unique AS genes in Fetal", length(unique(Fetal_SUPPA2$Gene)))
  cat("Number of unique AS genes in Adult", length(unique(Adult_SUPPA2$Gene)))
  
  #### Plot of genes and splicing events separated by two datasets
  
  merge_human <- bind_rows(Fetal_SUPPA2, Adult_SUPPA2) %>% spread(., Sample, No_of_Events)
  
  # Detection based on splice events on two datasets 
  for(row in 1:nrow(merge_human)){
    merge_human$Detected[row] <- 
      if(is.na(merge_human[row,"Human (Adult)"])){paste("Fetal Only")
      }else{
        if(is.na(merge_human[row,"Human (Fetal)"])){paste("Adult Only")
        }else{paste("Common")}}}
  
  merge_human$Detected <- factor(merge_human$Detected, levels=c("Common", "Adult Only", "Fetal Only"))
  
  # Note there will be duplicates of genes in the bars 
  p2 <- merge_human %>% group_by(Event, Detected) %>% tally() %>% 
    ggplot(., aes(x = Event, fill = Detected, y = n)) + 
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values=c("#999999", label_colour("Adult"), label_colour("Fetal"), 
                               labels=c("Common", "Human (Adult) Cortex Only", "Human (Fetal) Cortex Only"))) + 
    theme_bw() + 
    labs(x = "Splicing Events", y = "Number of AS Genes", title = "\n\n") + 
    mytheme + 
    theme(legend.position = c(.90, 0.90), legend.title = element_blank())
  
  return(list(p1,p2))
  
}

AS_genes_events <- function(dataset1, dataset2){
  cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  splicing_events <- read.csv(paste0(output_table_dir,"/AS_IR/ALL_SUPPA2_Genes_Output_Updated2.csv"))
  dataset_tally_events <- splicing_events %>% group_by(Sample) %>% tally(n)  
  
  # Number of splicing events
  splicing_events_tally <- splicing_events %>% group_by(Event, Sample) %>% tally(n) %>% left_join(dataset_tally_events, by = "Sample") %>% mutate(perc = n.x/n.y * 100)
  cat("Number of splicing events: \n")
  splicing_events_tally_present <- splicing_events_tally %>% 
    mutate(Num_Perc = paste0(n.x," (",round(perc,2),"%)")) %>% select(Event, Sample, Num_Perc) %>% spread(., Sample, Num_Perc)
  print(splicing_events_tally_present)
  splicing_events_tally_present <<- splicing_events_tally_present
  
  p1 <- splicing_events_tally %>%
    filter(Sample %in% c(dataset1,dataset2)) %>%
    mutate(Sample = paste(Sample, "Cortex")) %>%
    ggplot(., aes(x = Sample, y = perc, fill = reorder(Event, -perc))) + geom_bar(stat = "identity") + 
    theme_bw() + labs(y = "Splicing Events (%)", x = "", title = "\n\n") + mytheme + 
    theme(legend.position = "bottom",legend.title = element_blank()) +
    guides(fill = guide_legend(nrow = 1)) 
  
  # Number of genes
  total_gene <- class.files %>% bind_rows() %>% group_by(Sample, associated_gene) %>% tally() %>% group_by(Sample) %>% tally()  
  splicing_events_genes <- splicing_events %>% group_by(Event, Sample) %>% tally() %>% left_join(total_gene, by = "Sample") %>% mutate(perc = n.x/n.y * 100)
  cat("Number of splicing genes: \n")
  splicing_events_gene_present <- splicing_events_genes %>% 
    mutate(Num_Perc = paste0(n.x," (",round(perc,2),"%)")) %>% select(Event, Sample, Num_Perc) %>% spread(., Sample, Num_Perc)
  print(splicing_events_gene_present)
  splicing_events_gene_present <<- splicing_events_gene_present
  
  p2 <- splicing_events_genes %>%
    filter(Sample %in% c(dataset1,dataset2)) %>%
    mutate(Sample = paste(Sample, "Cortex")) %>%
    ggplot(., aes(x = Sample, y = perc, fill = reorder(Event, -perc),)) + geom_bar(stat = "identity", position = position_dodge()) + 
    theme_bw() + labs(y = "Genes (%)", x = "", title = "\n\n") + mytheme + 
    theme(legend.position = "bottom",legend.title = element_blank()) +
    guides(fill = guide_legend(nrow = 1)) 
  
  dataset_tally_gene <- splicing_events %>% group_by(Sample) %>% count(associated_gene) %>% tally()
  #nrow((unique(splicing_events[splicing_events$Sample == "Mouse","associated_gene"])))
  splicing_events_number <- splicing_events %>% group_by(Sample, associated_gene) %>% tally() %>% group_by(n, Sample) %>% tally() %>% left_join(dataset_tally_gene, by = "Sample") %>%
    mutate(perc = nn/n.y * 100)  %>% `colnames<-`(c("Number_of_Splicing_Events", "Sample", "Genes", "Total_Genes","perc")) %>% as.data.frame(.) 
  print(splicing_events_number)
  p3 <- splicing_events_number %>%
    filter(Sample %in% c(dataset1,dataset2)) %>%
    mutate(Sample = paste(Sample, "Cortex")) %>%
    ggplot(., aes(x = Number_of_Splicing_Events, y = perc, fill = Sample)) + 
    geom_bar(stat = "identity", position = position_dodge()) + 
    theme_bw() + labs(y = "AS Genes (%)", x = "Number of Splicing Events", title = "\n\n") + mytheme + 
    theme(legend.position = c(0.85,0.85), legend.title = element_blank()) + 
    scale_x_continuous(breaks = 1:7) +
    scale_fill_manual(values=c(label_colour(dataset1), label_colour(dataset2)))
  
  return(list(p1,p2,p3))
}

SUPPA2_NumAS_PerGene <- function(){
  # file generated from Alternative Splicing.Rmd
  num_events_merge <- read.csv(paste0(output_table_dir, "/AS_IR/SUPPA2_NumEvents_PerGene.csv"))
  
  total_num_genes <- num_events_merge %>% group_by(dataset) %>% tally() 
  plot <- num_events_merge %>% group_by(Number.of.Splicing.Events, dataset) %>% tally() %>% 
    left_join(., total_num_genes, by = "dataset") %>%
    mutate(perc = n.x/n.y) %>%
    mutate(dataset = factor(dataset, levels = c("AdultCTX", "FetalCTX", "HumanCTX","WholeIsoSeq"))) 
  
  plot_fun <- function(plot, d1,d2){
    
    p <- plot %>% filter(dataset %in% c(d1,d2)) %>%
    ggplot(., aes(x = Number.of.Splicing.Events, y = perc, fill = dataset)) + 
    geom_bar(stat = "identity", position = position_dodge()) + 
    theme_bw() + 
    mytheme + 
    labs(x = "Number of Splicing Events", y = "AS Genes (%)") + 
    scale_x_continuous(breaks = 1:7) + theme(legend.position =c(0.85,0.85))
    
    if(d1 == "AdultCTX"){
      p <- p + scale_fill_manual(labels = c("Human (Adult) Cortex","Human (Fetal) Cortex"), values = c("#440154FF","#287D8EFF"), name = "")
    } else {
      p <- p + scale_fill_manual(labels = c("Human Cortex","Mouse Cortex"), values = c("#29AF7FFF", "#FDE725FF"), name = "")
    }
    
  }
  
  p1 <- plot_fun(plot,"AdultCTX", "FetalCTX")
  p2 <- plot_fun(plot,"HumanCTX", "WholeIsoSeq")
  
  stats <- num_events_merge %>% group_by(Number.of.Splicing.Events, dataset) %>% tally() %>% 
    left_join(., total_num_genes, by = "dataset") %>%
    mutate(perc = n.x/n.y * 100) 
  
  colnames(stats)[2] <- "Num_of_AS_Genes"
  colnames(stats)[3] <- "Total_AS_Genes"
  
  return(list(p1,p2,stats))
}

SUPPA2_events_genes_plot <- function(){
  # file generated from Alternative Splicing.Rmd
  Sum_No_Genes <- read.csv(paste0(output_table_dir, "/AS_IR/SUPPA2_Sum_No_Genes.csv"))
  Sum_No_Events <- read.csv(paste0(output_table_dir, "/AS_IR/SUPPA2_Sum_No_Events.csv")) %>% rename(Percentage = Perc)
  
  plot_distribution <- function(input_df, dataset1, dataset2, y_label){
    p <- input_df %>% 
      mutate(dataset = ifelse(as.character(dataset) == "WholeIsoSeq", "Mouse", as.character(dataset))) %>% 
      mutate(dataset = ifelse(as.character(dataset) == "AdultCTX", "Human(Adult)", as.character(dataset))) %>%
      mutate(dataset = ifelse(as.character(dataset) == "FetalCTX", "Human(Fetal)", as.character(dataset))) %>%
      mutate(dataset = ifelse(as.character(dataset) == "HumanCTX", "Human", as.character(dataset))) %>%
      mutate(percentage_plot = Percentage/100) %>%
      filter(dataset %in% c(dataset1, dataset2)) %>%
      ggplot(., aes(x = dataset, y = `percentage_plot`, fill = Events)) + 
      geom_bar(stat = "identity") + 
      theme_bw() +
      mytheme +
      scale_fill_viridis(discrete=TRUE) + 
      scale_y_continuous(labels = scales::percent) +
      labs(x = "", y = y_label) 
    
    return(p)
  }
  
  
  p1 <- plot_distribution(Sum_No_Genes,"Human","Mouse", "Genes")
  p2 <- plot_distribution(Sum_No_Genes,"Human(Adult)","Human(Fetal)", "Genes")
  p3 <- plot_distribution(Sum_No_Events,"Human","Mouse", "Splicing Events")
  p4 <- plot_distribution(Sum_No_Events,"Human(Adult)","Human(Fetal)", "Splicing Events")
  
  return(list(p1,p2,p3,p4))
}


### Numbers comparison#################################################################

# human_mouse_comparison_num
# Aim: plot the number of isoforms and gene expression between mouse and adult 
# Input: Number of matched isoforms and gene expression 
# Prerequisite: Run Human_Mouse_Comparison.Rmd for input 
# note the number of homologous genes between huamn and mouse: 9,823 (nrow(num_isoform))
human_mouse_comparison_num <- function(homology_human_mouse){
  
  # read in the number of matched isoforms and gene expression (output from Human_Mouse_Comparison.Rmd)
  if(homology_human_mouse == "homology_yes"){
    cat("Working with homologous genes; Converting mouse genes to human equivalent genes", "\n")
    num_isoform <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/human_mouse_multiexonic_NUM_homology_fullymatched.csv"))
    num_isoform_all <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/human_mouse_multiexonic_NUM_homology.csv"))
    gene_expression <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/human_mouse_multiexonic_TPM_homologs_fullymatched.csv"))
    num_isoform_all_includemonoexonic <-  read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/human_mouse_all_NUM_homology.csv"))
  } else {
    num_isoform <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/human_mouse_multiexonic_NUM_fullymatched.csv"))
    gene_expression <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/human_mouse_multiexonic_TPM_fullymatched.csv"))
  }
  
  cat("Number of genes between human and mouse", nrow(num_isoform))
  
  # Correlation of matched isoseq numbers between human and mouse 
  cat("Correlation of Number of detected isoforms per Gene (Human) and Number of detected Isoforms per Gene (Mouse)")
  p1 <- density_plot(num_isoform, "Human_Detected_Num_Isoseq", "Mouse_Detected_Num_Isoseq", 
                     "Number of Isoforms - Human Cortex", "Number of Isoforms - Mouse Cortex", "") +
    ylim(0,80) + xlim(0,80)
  
  # Correlation of gene expression between human and mouse (only of matched genes)
  cat("Correlation of IsoSeq Gene Expression (Human) and IsoSeq Gene Expression (Mouse)")
  p2 <- density_plot(gene_expression, "Human_Detected_LogTPM_Isoseq","Mouse_Detected_LogTPM_Isoseq", 
                     "Iso-Seq Gene Expression (Log10TPM) - Human Cortex", "Iso-Seq Gene Expression (Log10TPM) - Mouse Cortex", "")
  
  # Histogram of gene expression (for determining threshold) (only of matched genes)
  p3 <- rbind(data.frame("Num" = gene_expression$Human_Detected_LogTPM_Isoseq, "Sample" = "human"), 
              data.frame("Num" = gene_expression$Mouse_Detected_LogTPM_Isoseq, "Sample" = "mouse")) %>% 
    ggplot(., aes(Num, fill = Sample)) + geom_density(alpha = 0.5) + 
    theme_bw() + 
    labs(x = "IsoSeq Gene Expression (Log10TPM)", y = "Density") +
    mytheme + 
    theme(legend.position = c(.90, 0.90), legend.title = element_blank()) +
    scale_fill_manual(labels = c("Human","Mouse"), values = c("#29AF7FFF", "#FDE725FF"), name = "")
  
  # Correlation of gene expression between human and mouse with IsoSeq Gene Expression < Log 10 TPM removed 
  #human_cutoff <- gene_expression %>% filter(Human_Detected_LogTPM_Isoseq > 2.5) 
  #mouse_cutoff <- gene_expression %>% filter(Mouse_Detected_LogTPM_Isoseq > 2.5) 
  both_cutoff <- gene_expression %>% filter(Human_Detected_LogTPM_Isoseq > 2.5 & Mouse_Detected_LogTPM_Isoseq > 2.5)
  
  # only include genes passing through threshold
  #cut_off_num_isoform <- bind_rows(num_isoform[num_isoform$associated_gene %in% human_cutoff$associated_gene,],
  #num_isoform[num_isoform$associated_gene %in% mouse_cutoff$associated_gene,])
  cut_off_num_isoform <- num_isoform[num_isoform$associated_gene %in% both_cutoff$associated_gene,]
  
  cat("Correlation of IsoSeq Number of isoforms (Human) and IsoSeq Number of isoforms(Mouse) with Gene Expression > 2.5 (Log10 TPM)")
  p4 <- density_plot(cut_off_num_isoform, "Human_Detected_Num_Isoseq", "Mouse_Detected_Num_Isoseq", 
                     "Number of Isoforms - Human Cortex", "Number of Isoforms - Mouse Cortex", "") +
    ylim(0,80) + xlim(0,80)
  
  
  # venn diagram of genes between human and mouse using the num_isoform_all datasheet
  # remove novel genes 
  human_both_gene <- num_isoform_all_includemonoexonic %>% filter(Detected %in% c("Human","Both")) %>% .[!grepl("NOVEL", .$associated_gene),] %>% select(associated_gene)
  mouse_both_gene <- num_isoform_all_includemonoexonic %>% filter(Detected %in% c("Mouse","Both")) %>% .[!grepl("NOVEL", .$associated_gene),] %>% select(associated_gene)
  
  # note non-overlapping gene could be homologues that were just not detected in other dataset
  # num_isoform_all[num_isoform_all$Detected == "Mouse" & num_isoform_all$associated_gene %in% homologs_genes$human,]
  # num_isoform_all[num_isoform_all$Detected == "Human" %in% homologs_genes$human,]
  return(list(p1,p2,p3,p4,human_both_gene,mouse_both_gene))
}

overlap_annotated_genes <- function(){
  annotated.class.files <- lapply(class.files, function(x) x[!grepl("NOVELGENE",x$associated_gene),])
  p1 <- grid.draw(venn_diagram_plot_twocircles(annotated.class.files$Fetal$associated_gene, annotated.class.files$Adult$associated_gene, 
                                               "Fetal(Human)", "Adult(Human)"))
  return(p1)
}

adult_fetal_comparison_num <- function(){
  
  # read in the number of matched isoforms and gene expression (output from Human_Mouse_Comparison.Rmd)
  num_isoform <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/adult_fetal_multiexonic_NUM_fullymatched.csv"))
  num_isoform_all <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/adult_fetal_multiexonic_NUM.csv"))
  gene_expression <- read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/adult_fetal_multiexonic_TPM_fullymatched.csv"))
  num_isoform_all_includemonoexonic <-  read.csv(paste0(output_table_dir,"/Human_Mouse_Comparisons/adult_fetal_all_NUM.csv"))
  
  # Correlation of matched isoseq numbers between adult and fetal 
  cat("Correlation of Number of detected isoforms per Gene (Adult) and Number of detected Isoforms per Gene (fetal)")
  p1 <- density_plot(num_isoform, "Human..Adult._Detected_Num_Isoseq", "Human..Fetal._Detected_Num_Isoseq", 
                     "Number of Isoforms - Human (Adult) Cortex", "Number of Isoforms - Human (Fetal) Cortex","") +
    ylim(0,60) + xlim(0,60)
  
  # Correlation of gene expression between human and mouse (only of matched genes)
  cat("Correlation of IsoSeq Gene Expression (Human) and IsoSeq Gene Expression (Mouse)")
  p2 <- density_plot(gene_expression, "Human..Adult._Detected_LogTPM_Isoseq","Human..Fetal._Detected_LogTPM_Isoseq", 
                     "Iso-Seq Gene Expression (Log10TPM) - Human (Adult) Cortex", "Iso-Seq Gene Expression (Log10TPM) - Human (Fetal) Cortex", "")
  
  # Histogram of gene expression (for determining threshold) (only of matched genes)
  p3 <- rbind(data.frame("Num" = gene_expression$Human..Adult._Detected_LogTPM_Isoseq, "Sample" = "Human (Adult) Cortex"), 
              data.frame("Num" = gene_expression$Human..Fetal._Detected_LogTPM_Isoseq, "Sample" = "Human (Fetal) Cortex")) %>% 
    ggplot(., aes(Num, fill = Sample)) + geom_density(alpha = 0.5) + 
    theme_bw() + 
    labs(x = "Iso-Seq Gene Expression (Log10TPM)", y = "Density") +
    mytheme + 
    theme(legend.position = c(.80, 0.90), legend.title = element_blank()) +
    scale_fill_manual(labels = c("Human (Adult) Cortex","Human (Fetal) Cortex"), values = c(label_colour("Adult"),label_colour(("Fetal"))), name = "")
  
  # Correlation of gene expression between human and mouse with IsoSeq Gene Expression < Log 10 TPM removed 
  #adult_cutoff <- gene_expression %>% filter(Adult_Detected_LogTPM_Isoseq > 2.5) 
  #fetal_cutoff <- gene_expression %>% filter(Fetal_Detected_LogTPM_Isoseq > 2.5) 
  both_cutoff <- gene_expression %>% filter(Human..Fetal._Detected_LogTPM_Isoseq > 2.5 & Human..Adult._Detected_LogTPM_Isoseq > 2.5) 
  
  # only include genes passing through threshold
  #cut_off_num_isoform <- bind_rows(num_isoform[num_isoform$associated_gene %in% adult_cutoff$associated_gene,],
  #num_isoform[num_isoform$associated_gene %in% fetal_cutoff$associated_gene,])
  cut_off_num_isoform <- num_isoform[num_isoform$associated_gene %in% both_cutoff$associated_gene,]
  
  cat("Correlation of IsoSeq Gene Expression (Adult) and IsoSeq Gene Expression (Fetal) with Gene Expression > 2.5 (Log10 TPM)")
  p4 <- density_plot(cut_off_num_isoform, "Human..Adult._Detected_Num_Isoseq", "Human..Fetal._Detected_Num_Isoseq", 
                     "Number of Isoforms - Human (Adult) Cortex", "Number of Isoforms - Human (Fetal) Cortex", "") +
    ylim(0,60) + xlim(0,60)
  
  
  # venn diagram of genes between human and mouse using the num_isoform_all datasheet
  # remove novel genes 
  adult_both_gene <- num_isoform_all_includemonoexonic %>% filter(Detected %in% c("Human (Adult)","Both")) %>% .[!grepl("NOVEL", .$associated_gene),] %>% select(associated_gene)
  fetal_both_gene <- num_isoform_all_includemonoexonic %>% filter(Detected %in% c("Human (Fetal)","Both")) %>% .[!grepl("NOVEL", .$associated_gene),] %>% select(associated_gene)
  
  return(list(p1,p2,p3,p4,adult_both_gene,fetal_both_gene))
}




### IR and NMD#################################################################
IR_NMD_run <- function(){
  
  library(RColorBrewer)
  myCol <- brewer.pal(3, "Set2")
  
  dat <- data.frame()
  list_IR_NMD_genes <- list()
  plot_IR_NMD <- list()
  for(i in 1:length(class.files)){
    df <- class.files[[i]]
    annotated_genes <- df[!grepl("NOVEL", df$associated_gene), ]
    IR_NMD <- df[!grepl("NOVEL", df$associated_gene), ] %>% filter(subcategory == "intron_retention" & predicted_NMD == "TRUE")
    IR_Not_NMD <- df[!grepl("NOVEL", df$associated_gene), ] %>% filter(subcategory == "intron_retention" & predicted_NMD == "FALSE")
    NMD <- df[!grepl("NOVEL", df$associated_gene), ] %>% filter(predicted_NMD == "TRUE")
    NMD_Not_IR <- df[!grepl("NOVEL", df$associated_gene), ] %>% filter(predicted_NMD == "TRUE") %>% filter(subcategory != "intron_retention")
    IR_coding <- df[!grepl("NOVEL", df$associated_gene), ] %>% filter(subcategory == "intron_retention" & coding == "coding")
    # number of IR_transcripts identified as predicted NMD (coding)
    dat[1,i] <- nrow(IR_NMD)
    # number of IR_transcripts identified as predicted NMD (coding)
    dat[2,i] <- nrow(IR_Not_NMD)
    # number of those IR-NMD transcripts classified as annotated 
    dat[3,i] <- nrow(IR_NMD[IR_NMD$associated_transcript != "novel", ])
    # number of those IR-NMD transcripts classified as novel
    dat[4,i] <- nrow(IR_NMD[IR_NMD$associated_transcript == "novel", ])
    # IR_transcripts (coding)
    dat[5,i] <- nrow(IR_coding) 
    # NMD transcripts
    dat[6,i] <- paste0(nrow(NMD)," (",round(nrow(NMD)/nrow(annotated_genes)*100,2),"%)")
    # Genes associated with IR_NMD transcripts 
    dat[7,i] <- paste0(length(unique(IR_NMD$associated_gene))," (", round(length(unique(IR_NMD$associated_gene))/
                                                                            length(unique(annotated_genes$associated_gene))*100,2), "%)")
    # genes associated with IR_transcripts (coding)
    dat[8,i] <- length(unique(IR_coding$associated_gene))
    # genes of NMD transcripts
    dat[9,i] <- paste0(length(unique(NMD$associated_gene))," (",round(length(unique(NMD$associated_gene))/
                                                                        length(unique(annotated_genes$associated_gene))*100,2),"%)")
    colnames(dat)[i] <- df$Sample[1]
    # % of IR_NMD out of NMD transcripts 
    dat[10,i] <- round(nrow(IR_NMD)/nrow(NMD)* 100,2)
    # % of IR_NMD out of IR transcripts 
    dat[11,i] <- round(nrow(IR_NMD)/nrow(IR_coding)* 100,2)
    # NMD transcripts but not intron retained 
    dat[12,i] <- nrow(NMD_Not_IR)
    # Number of genes with IR and NMD mutually exclusive (i.e. transcipt with IR but not NMD and vice versa)
    dat[13,i] <- paste0(length(setdiff(intersect(NMD$associated_gene,IR_coding$associated_gene),IR_NMD$associated_gene))," (", 
                        round(length(setdiff(intersect(NMD$associated_gene,IR_coding$associated_gene),IR_NMD$associated_gene))/
                                length(unique(annotated_genes$associated_gene))*100,2), "%)")  
    
    #plot[[i]] <- venn_diagram_plot_twocircles(unique(IR_coding$associated_gene), unique(NMD$associated_gene), "Genes with Intron Retention", "Genes with NMD")
    plot_IR_NMD[[i]] <- venn.diagram(
      x = list(unique(IR_coding$associated_gene), unique(NMD$associated_gene), unique(IR_NMD$associated_gene)),
      category.names = c("IR", "NMD", "IR&NMD"),
      filename = NULL,
      output=TRUE, 
      print.mode = "raw",
      fill = myCol, 
      cex = 2,
      fontface = "bold",
      fontfamily = "ArialMT",
      
      # Set names
      cat.cex = 2,
      cat.fontfamily = "ArialMT",
      
      main = paste0("\n\n"))
    
    list_IR_NMD_genes[[i]] <- data.frame(IR_NMD %>% group_by(associated_gene) %>% tally())
    colnames(list_IR_NMD_genes[[i]]) <- c("associated_gene", names(class.files)[[i]])
  }
  
  
  rownames(dat) <- c("number of IR_transcripts identified as predicted NMD (coding)",
                     "number of IR_transcripts not identified as predicted NMD (coding)",
                     "number of those IR-NMD transcripts classified as annotated", 
                     "number of those IR-NMD transcripts classified as novel",
                     "IR_transcripts (coding)",
                     "NMD transcripts",
                     "Genes associated with IR_NMD transcripts",
                     "Genes associated with IR_transcripts (coding)",
                     "Genes of NMD transcripts", 
                     "% of IR_NMD transcripts out of total NMD transcripts",
                     "% of IR_NMD transcripts out of total IR transcripts (coding)",
                     "Number of NMD transcripts not intron retained",
                     "Genes associated with mutually exclusive transcripts")
  names(plot_IR_NMD) <- names(class.files)
  print(dat)
  
  # List of IR genes per dataset 
  #write.csv(list_IR_NMD_genes %>% reduce(full_join, by = "associated_gene", all = T), paste0(output_table_dir, "/IR_NMD_Genes.csv"), row.names = F, quote = F)
  
  return(plot_IR_NMD)
}


NMD_vs_NonNMD <- function(type_class_file){
  count=1
  plots <- list()
  for(i in type_class_file){
    sink(paste0(output_corr_dir, "/Mann_Whitney_NMD.txt"), append = T)
    cat("\n#############################################################\n")
    cat("Statstical Tests for dataset:", names(class.files)[[count]], "\n")
    cat("#############################################################\n\n")
    
    # subset transcripts whether from annotated or novel genes
    annotated_genes <- i[!grepl("NOVEL", i$associated_gene),]
    novel_genes <- i[grepl("NOVEL", i$associated_gene),]
    
    # Annotated Genes, NMD transcripts
    NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "TRUE")
    non_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "FALSE")
    
    IR_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "TRUE" & subcategory == "intron_retention")
    Non_IR_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "TRUE" & subcategory != "intron_retention")
    IR_Non_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "FALSE" & subcategory == "intron_retention")
    Non_IR_Non_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "FALSE" & subcategory != "intron_retention")
    
    
    # Novel Genes, NMD transcripts 
    NMD_transcripts_novel_genes <- novel_genes %>% filter(predicted_NMD == "TRUE")
    
    cat("Test on Expression of NMD transcripts vs Non-NMD transcripts in Annotated genes\n")
    cat("Expression of NMD Transcripts in Annotated genes (IsoSeq TPM)\n", 
        "mean = ", mean(NMD_transcripts_annotated_genes$ISOSEQ_TPM), "s.d. =",  sd(NMD_transcripts_annotated_genes$ISOSEQ_TPM),"\n")
    cat("Expression of Non-NMD Transcripts in Annotated genes (IsoSeq TPM) \n", 
        "mean = ", mean(non_NMD_transcripts_annotated_genes$ISOSEQ_TPM), "s.d. =",  sd(non_NMD_transcripts_annotated_genes$ISOSEQ_TPM))
    test <- wilcox.test(NMD_transcripts_annotated_genes$ISOSEQ_TPM, non_NMD_transcripts_annotated_genes$ISOSEQ_TPM)
    print(test)
    print(test$p.value)
    
    cat("Test on Expression of IR-NMD transcripts vs IR, but no NMD transcripts in Annotated genes or NMD, but no IR\n")
    test1 <- wilcox.test(IR_NMD_transcripts_annotated_genes$ISOSEQ_TPM, c(IR_Non_NMD_transcripts_annotated_genes$ISOSEQ_TPM,
                                                                          Non_IR_NMD_transcripts_annotated_genes$ISOSEQ_TPM))
    print(test1)
    print(test1$p.value)
    
    p1 <- rbind(data.frame("Num" = IR_NMD_transcripts_annotated_genes$Log_ISOSEQ_TPM, "Type" = "IR-NMD"),
                data.frame("Num" = IR_Non_NMD_transcripts_annotated_genes$Log_ISOSEQ_TPM, "Type" = "IR"),
                data.frame("Num" = Non_IR_NMD_transcripts_annotated_genes$Log_ISOSEQ_TPM,"Type" = "NMD"),
                data.frame("Num" = Non_IR_Non_NMD_transcripts_annotated_genes$Log_ISOSEQ_TPM, "Type" = "Non IR-NMD")) %>%
      ggplot(., aes(y = Num, x = Type)) + geom_boxplot() +
      theme_bw() +
      labs(x = "", y = "Iso-Seq Transcript Expression (Log10 TPM)", title = paste0("\n\n")) + mytheme +
      theme(legend.position = c(.90, 0.90), legend.title = element_blank()) 
    
    plots[[count]] <- p1
    count = count + 1
  }
  names(plots) <- names(class.files) 
  sink()
  
  
  # Plot of Expression 
  # Further filtering
  #annotated_genes <- lapply(class.files, function(x) x[!grepl("NOVEL", x$associated_gene),]) %>% bind_rows()
  #IR_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "TRUE" & subcategory == "intron_retention")
  #Non_IR_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "TRUE" & subcategory != "intron_retention")
  #IR_Non_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "FALSE" & subcategory == "intron_retention")
  #Non_IR_Non_NMD_transcripts_annotated_genes <- annotated_genes %>% filter(predicted_NMD == "FALSE" & subcategory != "intron_retention")
  
  #p1 <- rbind(data.frame("Num" = IR_NMD_transcripts_annotated_genes$ISOSEQ_TPM, "Sample" = IR_NMD_transcripts_annotated_genes$Sample, "Type" = "IR and NMD"),
  #data.frame("Num" = Non_IR_NMD_transcripts_annotated_genes$ISOSEQ_TPM,"Sample" = Non_IR_NMD_transcripts_annotated_genes$Sample,"Type" = "NMD Only, No IR"),
  #data.frame("Num" = IR_Non_NMD_transcripts_annotated_genes$ISOSEQ_TPM, "Sample" = IR_Non_NMD_transcripts_annotated_genes$Sample,"Type" = "IR Only, No NMD"),
  #data.frame("Num" = Non_IR_Non_NMD_transcripts_annotated_genes$ISOSEQ_TPM, "Sample" = Non_IR_Non_NMD_transcripts_annotated_genes$Sample, "Type" = "No IR,No NMD")) %>%
  #filter(Sample == c("Human","Mouse")) %>%
  #ggplot(., aes(y = Num, x = Sample, fill = Type)) + geom_boxplot() +
  #theme_bw() +
  #labs(x = "", y = "Transcript IsoSeq Expression (Log10 TPM)") +
  # mytheme +
  #theme(legend.position = c(.90, 0.90), legend.title = element_blank()) 
  #print(p1)
  
  return(plots)
}


### lncRNA #################################################################
lncRNA <- function(){
  
  # Internal Function
  # statstical test between lncRNA and non-lncRNA
  lncRNA_tests <- function(all_dat, lnc_dat, title, sample){
    
    cat(title,"\n")
    test <- wilcox.test(all_dat, lnc_dat)
    
    if(title == "Transcript Expression (Annotated Genes only)"){
      p <- rbind(data.frame("Num" = all_dat, "Sample" = "Non-lncRNA"),
                 data.frame("Num" = lnc_dat, "Sample" = "lncRNA")) %>%
        ggplot(., aes(y = Num, x = Sample)) + geom_boxplot() +
        theme_bw() +
        labs(x = "", y = "Iso-Seq Transcript Expression (Log10 TPM)", title = paste0("\n\n")) +
        mytheme + scale_fill_manual(values = wes_palette("Darjeeling2")) +
        theme(legend.position = c(.85, 0.85), legend.title = element_blank())
    }else if(title == "Number of Exons (Annotated Genes only)"){
      p <- rbind(data.frame("Num" = all_dat, "Sample" = "Non-lncRNA"),
                 data.frame("Num" = lnc_dat, "Sample" = "lncRNA")) %>%
        ggplot(., aes(Num, fill = Sample)) + geom_density(alpha = 0.2) +
        theme_bw() +
        labs(x = "Number of Exons", y = "Density", title = paste0("\n\n")) +
        theme_bw() + 
        mytheme + 
        theme(legend.position = c(.85, 0.85), legend.title = element_blank())
    }else if(title == "ORF Length (Annotated Genes only)"){
      p <- rbind(data.frame("Num" = all_dat, "Sample" = "Non-lncRNA"),
                 data.frame("Num" = lnc_dat, "Sample" = "lncRNA")) %>%
        ggplot(., aes(Num, fill = Sample)) + geom_density(alpha = 0.2) +
        theme_bw() +
        labs(x = "ORF Length (kb)", y = "Density", title = paste0("\n\n")) +
        theme_bw() + 
        mytheme + 
        theme(legend.position = c(.85, 0.85), legend.title = element_blank())
    }else if(title == "Number of Isoforms per Gene"){
      p <- rbind(data.frame("Num" = all_dat, "Sample" = "Non-lncRNA"),
                 data.frame("Num" = lnc_dat, "Sample" = "lncRNA")) %>%
        ggplot(., aes(y = Num, x = Sample)) + geom_boxplot() +
        theme_bw() +
        labs(x = "", y = "Number of Isoforms per Gene (Log10)", title = paste0("\n\n")) +
        theme_bw() + 
        mytheme + 
        theme(legend.position = c(.85, 0.85), legend.title = element_blank())
    }
    else{
      p <- rbind(data.frame("Num" = all_dat, "Sample" = "Non-lncRNA"),
                 data.frame("Num" = lnc_dat, "Sample" = "lncRNA")) %>%
        ggplot(., aes(Num, fill = Sample)) + geom_density(alpha = 0.2) +
        theme_bw() +
        labs(x = "Transcript Length (kb)", y = "Density", title = paste0("\n\n")) +
        theme_bw() + 
        mytheme + 
        theme(legend.position = c(.85, 0.85), legend.title = element_blank())
    }
    
    print(test)
    cat("Exact p value:", test$p.value,"\n")
    return(p)
  }
  
  
  dat <- data.frame()
  plots <- list()
  for(i in 1:length(class.files)){
    # input file = annotated genes only 
    input_file <- class.files[[i]][!grepl("NOVEL", class.files[[i]]$associated_gene),]
    # lncfiles = classifiction files aligned to lncRNA reference annotation
    input_lnc_file <- lnc.files[[i]]
    # transcrips annotated genes = known transcripts of annotated genes 
    transcripts_annotated_genes <- input_file[!grepl("NOVEL", input_file$associated_gene),]
    # genes_annotated_genes = tally of the known transcripts of annotated genes
    genes_annotated_genes <- transcripts_annotated_genes %>% group_by(associated_gene) %>% tally()
    # coding_transcripts_annotated_genes = coding, known transcripts of annotated genes 
    coding_transcripts_annotated_genes <- input_file[!grepl("NOVEL", input_file$associated_gene) & input_file$coding == "coding",]
    coding_genes_annotated_genes <- coding_transcripts_annotated_genes %>% group_by(associated_gene) %>% tally()
    
    # annotated lncRNA transcripts and genes from lnc.files 
    lncRNA_transcripts <- input_lnc_file[!grepl("novel",input_lnc_file$associated_gene),]
    lncrna_genes <- unique(input_lnc_file[!grepl("novel",input_lnc_file$associated_gene),"associated_gene"])
    
    # subset of lncRNA transcripts from original classification file using genes extracted from lnc.files 
    class_non_lncRNA_transcripts <- input_file[!input_file$associated_gene %in% lncrna_genes,]
    class_lncRNA_transcripts <- input_file[input_file$associated_gene %in% lncrna_genes,]
    
    # number of isoforms (n) of lncRNA genes vs nonlncRNA genes for each gene
    isoform_diversity <- bind_rows(class_non_lncRNA_transcripts %>% group_by(associated_gene) %>% tally() %>% mutate(gene_type = "non_lncRNA"),
                                   class_lncRNA_transcripts %>% group_by(associated_gene) %>% tally() %>% mutate(gene_type = "lncRNA")) 
    isoform_diversity_lncRNA <- isoform_diversity[isoform_diversity$gene_type == "lncRNA",]
    isoform_diversity_nonlncRNA <- isoform_diversity[isoform_diversity$gene_type == "non_lncRNA",]
    
    dat[1,i] <- nrow(transcripts_annotated_genes)
    dat[2,i] <- nrow(genes_annotated_genes)
    dat[3,i] <- paste0(nrow(coding_transcripts_annotated_genes), "(",
                       round(nrow(coding_transcripts_annotated_genes)/nrow(transcripts_annotated_genes) * 100, 2), "%)")
    dat[4,i] <- nrow(coding_genes_annotated_genes)
    dat[5,i] <- nrow(class_lncRNA_transcripts)
    dat[6,i] <- nrow(class_lncRNA_transcripts %>% group_by(associated_gene) %>% tally())
    dat[7,i] <- paste0(nrow(class_lncRNA_transcripts %>% filter(coding == "coding")),"(",
                       round(nrow(class_lncRNA_transcripts  %>% filter(coding == "coding"))/nrow(class_lncRNA_transcripts ) * 100,2),"%)")
    dat[8,i] <- nrow(class_lncRNA_transcripts  %>% filter(coding == "non_coding"))
    dat[9,i] <- paste0(round(mean(class_lncRNA_transcripts$ORF_length, na.rm = TRUE),1),"/",
                       round(sd(class_lncRNA_transcripts$ORF_length, na.rm = TRUE),1))  # NA for non-coding
    dat[10,i] <- paste0(round(mean(class_non_lncRNA_transcripts$ORF_length,  na.rm = TRUE),1),"/",
                        round(sd(class_non_lncRNA_transcripts$ORF_length, na.rm = TRUE),1))# NA for non-coding
    # known lncRNA genes, known lncRNA transcripts
    dat[11,i] <- nrow(class_lncRNA_transcripts[!grepl("novel", class_lncRNA_transcripts $associated_transcript),]) 
    # known lncRNA genes, novel lncRNA transcripts
    dat[12,i] <- nrow(class_lncRNA_transcripts [grepl("novel", class_lncRNA_transcripts $associated_transcript),]) 
    dat[13,i] <- paste0(nrow(class_lncRNA_transcripts %>% filter(subcategory == "mono-exon")),"/",
                        round(nrow(class_lncRNA_transcripts %>% filter(subcategory == "mono-exon"))/nrow(class_lncRNA_transcripts) * 100,2))
    dat[14,i] <- paste0(nrow(class_non_lncRNA_transcripts %>% filter(subcategory == "mono-exon")),"/",
                        round(nrow(class_non_lncRNA_transcripts %>% filter(subcategory == "mono-exon"))/nrow(class_non_lncRNA_transcripts) * 100, 2))
    dat[15,i] <- nrow(class_non_lncRNA_transcripts)
    dat[16,i] <- paste0(round(mean(class_lncRNA_transcripts$length),1),"(s.d = ",
                        round(sd(class_lncRNA_transcripts$length),1),", range = ",
                        round(min(class_lncRNA_transcripts$length),1),"-",
                        round(max(class_lncRNA_transcripts$length),1),")")
    dat[17,i] <- paste0(round(mean(class_non_lncRNA_transcripts$length),1)," (s.d = ",
                        round(sd(class_non_lncRNA_transcripts$length),1),", range = ",
                        round(min(class_non_lncRNA_transcripts$length),1),"-",
                        round(max(class_non_lncRNA_transcripts$length),1),")")
    dat[18,i] <- paste0("mean = ", round(mean(isoform_diversity_nonlncRNA$n),2),
                        ",min = ", min(isoform_diversity_nonlncRNA[isoform_diversity_nonlncRNA$gene_type == "non_lncRNA","n"]),
                        ",max = ",max(isoform_diversity_nonlncRNA[isoform_diversity_nonlncRNA$gene_type == "non_lncRNA","n"]))
    dat[19,i] <- paste0("mean = ", round(mean(isoform_diversity_lncRNA$n),2),
                        ",min = ", min(isoform_diversity_lncRNA[isoform_diversity_lncRNA$gene_type == "lncRNA","n"]),
                        ",max = ",max(isoform_diversity_lncRNA[isoform_diversity_lncRNA$gene_type == "lncRNA","n"]))
    
    colnames(dat)[i] <- names(class.files)[i]
    
    # Run Statstical Test per dataset 
    sink(paste0(output_corr_dir,"/Mann_Whitney_LncRNA.txt"), append = TRUE)
    cat("#############################################################\n")
    cat("Statstical Tests for dataset:", names(class.files)[[i]], "\n")
    cat("#############################################################\n\n")
    
    # Ensure non-lncRNA dataset first
    # p1: Transcripts Length (Annotated Genes only)
    # p2: Number of Exons (Annotated Genes only)
    # p3: Non-Coding Transcripts Length (Annotated Genes only)
    # p4: Coding Transcripts Length (Annotated Genes only)
    # p5: Transcript Expression (Annotated Genes only)
    # p6: ORF Length (Annotated Genes only)
    # p7: Number of Isoforms per Gene
    
    
    # comparing transcript length of lncRNA vs non-lncRNA (irrespective of coding or not)
    p1 <- lncRNA_tests(class_non_lncRNA_transcripts$length, class_lncRNA_transcripts$length, 
                       "Transcripts Length (Annotated Genes only)",names(class.files)[[i]]) + 
      scale_x_continuous(labels = ks) 
    
    p2 <- lncRNA_tests(class_non_lncRNA_transcripts$exons, class_lncRNA_transcripts$exons, 
                       "Number of Exons (Annotated Genes only)",
                       names(class.files)[[i]])
    
    p3 <- lncRNA_tests(class_non_lncRNA_transcripts[class_non_lncRNA_transcripts$coding == "non_coding","length"],
                       class_lncRNA_transcripts[class_lncRNA_transcripts$coding == "non_coding","length"], 
                       "Non-Coding Transcripts Length (Annotated Genes only)", names(class.files)[[i]]) +
      scale_x_continuous(labels = ks) 
    
    p4 <- lncRNA_tests(class_non_lncRNA_transcripts[class_non_lncRNA_transcripts$coding == "coding","length"],
                       class_lncRNA_transcripts[class_lncRNA_transcripts$coding == "coding","length"], 
                       "Coding Transcripts Length (Annotated Genes only)",names(class.files)[[i]]) +
      scale_x_continuous(labels = ks) 
    
    
    p5 <- lncRNA_tests(class_non_lncRNA_transcripts$Log_ISOSEQ_TPM, class_lncRNA_transcripts$Log_ISOSEQ_TPM, 
                       "Transcript Expression (Annotated Genes only)",
                       names(class.files)[[i]])
    
    p6 <- lncRNA_tests(class_non_lncRNA_transcripts$ORF_length, class_lncRNA_transcripts$ORF_length, 
                       "ORF Length (Annotated Genes only)",names(class.files)[[i]]) +
      scale_x_continuous(labels = ks) 
    
    # Log number for plotting, makes no difference to statstical test
    p7 <-lncRNA_tests(log10(isoform_diversity_nonlncRNA$n), log10(isoform_diversity_lncRNA$n), 
                      "Number of Isoforms per Gene",names(class.files)[[i]])
    
    
    plots[[i]] <- list(p1,p2,p3,p4,p5,p6,p7)
    names(plots)[[i]] <- names(class.files)[i]
    sink()
  }
  
  rownames(dat) <- c("Number of total transcripts of annotated genes", 
                     "Number of genes of annotated genes",
                     "Number and % of coding transcripts of annotated genes", 
                     "Number of coding genes of annnotated genes", 
                     "Number of lncRNA transcripts",
                     "Number of lncRNA genes", 
                     "Number of lncRNA coding transcripts", 
                     "Number of lncRNA non_coding transcripts", 
                     "Mean/standard deviation ORF Length of lncRNA transcripts",
                     "Mean/standard deviation ORF Length of non-lncRNA transcripts", 
                     "Number of known lncRNA transcripts from known lncRNA genes", 
                     "Number of novel lncRNA transcripts from known lncRNA genes", 
                     "Number and % of monoexonic lncRNA transcripts", 
                     "Number and % of monoexonic non-lncRNA transcripts",
                     "Number of non lncRNA transcripts",
                     "Transcript Length of lncRNA transcripts",
                     "Transcript Length of non lncRNA transcripts", 
                     "Isoform diversity per nonlncRNA gene",
                     "Isoform diversity per lncRNA gene")
  
  print(dat)
  dat <<- dat
  return(plots)
}



### Novel Genes#################################################################
# subset classification files (filtered) by annotated genes, and novel transcripts 
expression_novel_annotated_genes <- function(){
  
  # internal functions
  # extract specific columns for plotting after merging all datasets and creating classifier 
  extract_feature <- function(dat, dat_name, col_name_feature){
    dat1 <- dat %>% bind_rows() %>% 
      mutate(Transcripts = dat_name) %>% .[,c(col_name_feature, "Transcripts", "Sample")] 
  }
  
  # summary_info_TPM 
  # abridged function of summary_info, but specific to calculating the expression of novel genes by aggregation with PB.ID instead of "novelgene category"
  # recalculate TPM by removing monoexonic transcripts from total FL reads 
  # use PB.ID rather than novel gene category as overestimation of numbers --> Liz renames the same novel gene differently (despite same first part of PB.ID)
  summary_info_TPM <- function(dat){
    total_fl <- sum(dat$FL, na.rm=T)
    dat$PB.ID <- word(dat$isoform,c(2),  sep = fixed ('.'))
    novel_genes <- dat[grepl("NOVELGENE",dat$associated_gene),] 
    
    final <- novel_genes %>% group_by(PB.ID) %>% summarise(Total_Reads = sum(FL)) %>% 
      mutate("FL_TPM" = round(Total_Reads*(10**6)/total_fl)) %>%
      mutate("Log10_FL_TPM" = log10(FL_TPM))
    
    return(final)
  }
  
  
  # boxplot 
  # plot for expression of novel transcripts (of novel genes) vs all transcripts (of annotated genes) at transcript level 
  # or plot for expression of novel genes vs annotated genes at gene level 
  boxplot <- function(annotated_set1, novel_set2, col_name_feature, ylabel, cate_name,dataset1,dataset2){
    set1 <- extract_feature(annotated_set1, "Annotated Genes", col_name_feature)
    set2 <- extract_feature(novel_set2, "Novel Genes", col_name_feature) 
    merge <- bind_rows(list(set1, set2))
    
    y.var <- rlang::sym(quo_name(enquo(col_name_feature)))
    p <- merge %>% 
      filter(Sample %in% c(dataset1, dataset2)) %>%
      ggplot(., aes(x = Sample, y = !! y.var, fill= Transcripts)) + 
      geom_boxplot() + 
      theme_bw() + 
      mytheme + 
      labs(x = "", y= ylabel, title = "\n\n") + 
      theme(legend.justification=c(1,1), legend.position=c(1,1)) + 
      scale_fill_discrete(name = cate_name)
    
    return(list(p, merge))
  }
  
  
  # subset classification files at transcript level to capture expression
  multiexonic <- lapply(class.files, function(x) x %>% filter(subcategory != "mono-exon"))
  annotated.gene.class.files <- lapply(multiexonic, function(x) x[!grepl("NOVELGENE",x$associated_gene),])
  novel.gene.class.files <- lapply(multiexonic, function(x) x[grepl("NOVELGENE",x$associated_gene),])
  
  # subset classification file at gene level, using summary_info (for annotated_genes) and summary_info_TPM (novel_genes)
  # able to use summary_info for annotated genes as overestimation from double renaming isn't likely 
  annotated.transcript.class.files <- lapply(multiexonic, function(x) data.frame(summary_info(x)) %>% .[!grepl("NOVELGENE",.$associated_gene),])
  novel.transcript.class.files <- lapply(multiexonic, function(x) summary_info_TPM(x))
  annotated.transcript.class.files$Mouse$Sample <- "Mouse Cortex"
  annotated.transcript.class.files$Adult$Sample <- "Human (Adult) Cortex"
  annotated.transcript.class.files$Fetal$Sample <- "Human (Fetal) Cortex"
  annotated.transcript.class.files$Human$Sample <- "Human Cortex"
  novel.transcript.class.files$Mouse$Sample <- "Mouse Cortex"
  novel.transcript.class.files$Adult$Sample <- "Human (Adult) Cortex"
  novel.transcript.class.files$Fetal$Sample <- "Human (Fetal) Cortex"
  novel.transcript.class.files$Human$Sample <- "Human Cortex"
  
  # plot 
  transcript <- boxplot(annotated.gene.class.files,novel.gene.class.files,"Log_ISOSEQ_TPM",
                        "Iso-Seq Transcript Expression (Log10 TPM)", "Transcripts from","Human Cortex","Mouse Cortex")
  gene <- boxplot(annotated.transcript.class.files,novel.transcript.class.files,
                  "Log10_FL_TPM","Iso-Seq Gene Expression (Log10 TPM)", "", "Human Cortex","Mouse Cortex")
  
  output <- list(transcript[[1]], gene[[1]], gene[[2]])
  names(output) <- c("Transcript_Plot", "Gene_Plot","Gene_Stats")
  
  # mann whitney test
  sink(paste0(output_corr_dir, "/Mann_Whitney_Transcript_Expression_Novel_genes"))
  for(i in unique(factor(output$Gene_Stats$Sample))){
    cat("\n\nMann Whitney for Expression of Novel Genes vs Annotated Genes of dataset:", i)
    dat <- output$Gene_Stats[output$Gene_Stats$Sample == paste(i),]
    test <- wilcox.test(Log10_FL_TPM ~ Transcripts,data=dat)
    print(test)
    cat("Exact p value:", test$p.value)
  }
  sink()
  
  return(output)
  
  # rationale for using PBId for aggregating number of novel genes rahter than novelgene category
  # check <- summary_info_TPM(multiexonic$Mouse)
  #check[check$PB.ID == "1291",]
  #class.files$Mouse[class.files$Mouse$isoform %in% c("PB.1291.1","PB.1291.2"),"FL"] 
}

length_novel_annotated_genes <- function(){
  
  # internal functions
  # extract specific columns for plotting after merging all datasets and creating classifier 
  extract_feature <- function(dat, dat_name, col_name_feature){
    dat1 <- dat %>% bind_rows() %>% 
      mutate(Transcripts = dat_name) %>% .[,c(col_name_feature, "Transcripts", "Sample")] 
    return(dat1)
  }
  
  boxplot <- function(annotated_set1, novel_set2, col_name_feature, ylabel, cate_name){
    set1 <- extract_feature(annotated_set1, "Annotated Genes", col_name_feature)
    set2 <- extract_feature(novel_set2, "Novel Genes", col_name_feature) 
    merge <- bind_rows(list(set1, set2))
    
    y.var <- rlang::sym(quo_name(enquo(col_name_feature)))
    p <- merge %>% mutate(Sample = paste(Sample, "Cortex")) %>%
      ggplot(., aes(x = Sample, y = !! y.var, fill= Transcripts)) + 
      geom_boxplot() + 
      theme_bw() + 
      mytheme + 
      labs(x = "", y= ylabel, title = "\n\n") + 
      theme(legend.justification=c(1,1), legend.position=c(1,1)) +
      scale_fill_discrete(name = cate_name) + 
      scale_y_continuous(labels = ks)
    
    output <- list(p, merge)
    names(output) <- c("plot","stats")
    return(output)
  }
  
  
  # subset classification files at transcript level to capture expression
  multiexonic <- lapply(class.files, function(x) x %>% filter(subcategory != "mono-exon"))
  annotated.gene.class.files <- lapply(multiexonic, function(x) x[!grepl("NOVELGENE",x$associated_gene),])
  novel.gene.class.files <- lapply(multiexonic, function(x) x[grepl("NOVELGENE",x$associated_gene),])
  
  
  # plot 
  human_mouse <- boxplot(list(annotated.gene.class.files$Human, annotated.gene.class.files$Mouse),
                         list(novel.gene.class.files$Human, novel.gene.class.files$Mouse),"length",
                         "Iso-Seq Transcript Length (kb)", "")[["plot"]]
  
  fetal_adult <- boxplot(list(annotated.gene.class.files$Adult, annotated.gene.class.files$Fetal),
                         list(novel.gene.class.files$Adult, novel.gene.class.files$Fetal),"length",
                         "Iso-Seq Transcript Length (kb)", "") [["plot"]]
  
  all <-  boxplot(annotated.gene.class.files, novel.gene.class.files,"length",
                  "Iso-Seq Transcript Length (kb)", "") 
  
  sink(paste0(output_corr_dir, "/Mann_Whitney_Transcript_Length_Novel_Genes.txt"))
  for(i in unique(factor(all$stats$Sample))){
    cat("\n\nMann Whitney for Length of Novel Genes vs Annotated Genes of dataset:", i)
    dat <- all$stats[all$stats$Sample == paste(i),]
    test <- wilcox.test(length ~ Transcripts,data=dat)
    print(test)
    cat("Exact p value:", test$p.value)
  }
  sink()
  
  return(list(human_mouse,fetal_adult))
}


### ONT_valdation#################################################################

ONT_validation <- function(){
  ONT_validated_df <- data.frame(categories = c("Total_Isoforms","Novel_Isoforms"), 
                                 PacBio_Total = c(42645, 13931),
                                 ONT_validated = c(27715, 7081))
  
  
  p1 <- ONT_validated_df %>% 
    # number of PacBio transcripts not supported by ONT
    mutate(PacBio_Only = PacBio_Total - ONT_validated) %>% 
    # calculate percentages and select columns for plots
    mutate(perc_supported = ONT_validated / PacBio_Total * 100) %>% 
    mutate(perc_not_supported = PacBio_Only/ PacBio_Total * 100) %>%
    select(categories, perc_not_supported, perc_supported) %>%
    melt(., id = "categories") %>% 
    # plot datawrangle 
    mutate(categories = factor(.$categories,levels = c("Total_Isoforms", "Novel_Isoforms"), label = c("All Transcripts", "Novel Transcripts, Annotated Genes"))) %>%
    ggplot(., aes(x = categories, y = value, fill = variable)) + geom_bar(stat = "identity") + 
    labs(x = "", y = "Iso-Seq Transcripts (%)") + 
    mytheme + 
    scale_fill_manual(labels=c("Not supported by ONT", "Supported by ONT"),
                      values = c(alpha(wes_palette("Royal1")[1],0.4), wes_palette("Darjeeling1")[2]), 
                      name = "", 
                      guide = guide_legend(reverse=TRUE)) +
    theme(legend.position = "bottom")
  
  return(p1)
}


### ERCC #################################################################
run_ERCC <- function(class){
  
  ERCC_conc <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/ERCC/ERCC_calculations.csv", header = T)[-1,]
  
  cat("Total unique ERCCs:", length(unique(class$chrom)), paste0("(", round(length(unique(class$chrom))/92 * 100,2), "%)"))
  
  # redundant 
  redundant <- class[,c("chrom")] %>% table() %>% melt 
  colnames(redundant) <- c("ERCC","num_isoforms")
  
  p1 <- redundant %>% 
    group_by(num_isoforms) %>% tally() %>% ggplot(., aes(x = as.factor(num_isoforms), y = n)) + geom_bar(stat = "identity") + 
    mytheme + labs(x = "ERCC", y = "Number of Isoforms")
  
  # isoform vs concentration
  isoform_conc <- merge(redundant, ERCC_conc, by.x = "ERCC", by.y = "ERCC_ID", all = TRUE) %>%
    mutate(log2_amount_of_ERCC = log2(amount_of_ERCC)) %>%
    replace_na(list(num_isoforms = 0)) %>% 
    mutate(num_isoforms = as.factor(num_isoforms))
  
  p2 <- ggplot(isoform_conc, aes(x = num_isoforms, y = log2_amount_of_ERCC, colour = num_isoforms)) + geom_jitter(width = 0.2) + theme(legend.position = "none") + 
    mytheme + labs(x = "Number of Isoforms", y = "Amount of ERCC (log2)") + theme(legend.position = "none")
  
  # correlation 
  ERCC_corr <- merge(class,ERCC_conc, by.x = "chrom", by.y = "ERCC_ID") %>% 
    mutate(log2_amount_of_ERCC = log2(amount_of_ERCC)) %>%
    mutate(log2_FL_reads = log2(FL)) 
  
  p3 <- density_plot(ERCC_corr,"log2_amount_of_ERCC","log2_FL_reads", "Amount of ERCC (Log2)", "Number of Full-Length Reads (Log2)","")
  p3
  
  return(list(p1,p2,p3))
}

### vs Targeted #################################################################

whole_vs_targeted_plots <- function(){
  
  TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")
  # from TAMA merge transcript file, filter transcripts associated with target gene 
  ADtrans <- TAMA_transfile %>% filter(toupper(gene_name) %in% TargetGene)
  # from TAMA merge transcript file, filter transcripts not associated with target gene 
  nonADtrans <- TAMA_transfile %>% filter(!toupper(gene_name) %in% TargetGene)
  
  #### plots ###
  #p1 = number of isoforms per target gene identified in whole transcriptome, targeted transcriptome, and in both
  #p2 = FL read counts in targeted transcriptome of isoforms annotated to AD genes (target) that are detected uniquely in targeted transcriptome and also in whole
  #p3 = FL read counts in whole transcriptome of isoforms not annotated to AD genes (off target) that are detected uniquely in whole and also in targeted
  
  ##### p1
  # calculate the total number of merged transcripts for plot in reorder 
  ADtrans_ttgene <- TAMA_transfile %>% group_by(gene_name) %>% tally()
  
  # Tally the merged transcripts by where it is from and also by gene name 
  # source from tama merge transcript file would either be "Targeted" or "Targeted,Whole" or "Whole"
  p1 <- ADtrans %>% group_by(sources,gene_name) %>% tally() %>%
    mutate(sources = recode(sources,"Targeted,Whole" = "Both")) %>%
    mutate(sources = factor(sources, levels = c("Whole","Both","Targeted"))) %>%
    left_join(.,ADtrans_ttgene, by = "gene_name") %>%
    ggplot(.,aes(x = reorder(gene_name,-n.y), y = n.x, fill = sources)) + geom_bar(stat = "identity") + mytheme + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(x = "", y = "Number of Isoforms") +
    scale_fill_manual(name = "", values = c(label_colour("whole"),label_colour("whole+targeted"),label_colour("targeted")))
  
  return(p1)
}
