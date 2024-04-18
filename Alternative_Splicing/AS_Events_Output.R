#!/usr/bin/env Rscript
# Szi Kay Leung: sl693@exeter.ac.uk
# Aim: Classification of AS events in Iso-Seq data (IR, SE, MX, A5', A3', AF, AL) 
	# IR(Intron Retention) - identified from SQANTI2 
	# MX(Mutually Exclusive), SE (Skipped Exons) - identified from SUPPA2 
	# A5'(Alternative 5'), A3'(Alternative 3'), AF(Alternative First), AL(Alternative Last) - identified from custom scripts (AS_Events.R)

library("dplyr")
library("tidyr")
library("ggplot2")
mouse_suppa_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SUPPA"
human_suppa_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/Post_IsoSeq/SUPPA"
output_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables"

# theme for plots
mytheme <- theme(axis.line = element_line(colour = "black"),
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

### SQANTI2 input for Intron Retention ###################################################################
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Figures/All_Plots_Functions.R")
# Read in input files 
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Input_Variables.R")
sqanti_files()
sqanti_gtf()

### SUPPA2 input for Mutually Exclusive Events and Skipped Exon ###################################################################

# Suppa2_input 
# Input: path of SUPPA2 output files, and prefix name 
# Output: output table for MX and SE 
Suppa2_input <- function(SUPPA2_input_dir, dataset){
  SUPPA2_MX_output_files <- list.files(path = SUPPA2_input_dir, pattern = paste0('^',dataset,'_MX_strict.ioe.txt'), full.names = TRUE)
  SUPPA2_SE_output_files <- list.files(path = SUPPA2_input_dir, pattern = paste0('^',dataset,'_SE_strict.ioe.txt'), full.names = TRUE)
  
  extract_gene <- function(output_file){
    dat <- read.table(output_file, header = TRUE)
    # Assumption: Each splicing event is only involving the one same gene, and no overlap between genes/fusion genes
    transcript <- word(dat$total_transcripts, c(1),  sep = fixed (','))
    gene <- word(transcript, c(1), sep = fixed('_'))
    dat$associated_gene <- gene 
    
    # Splicing eent name taken from first entry of each row of event_id
    event <- word(dat$event_id, c(1),  sep = fixed (':'))
    event <- word(event, c(2),  sep = fixed (';'))
    dat$Event <- event
    
    return(dat)
  }
  
  MX_SE <- rbind(extract_gene(SUPPA2_MX_output_files), extract_gene(SUPPA2_SE_output_files))
  return(MX_SE)
}

# Read in suppa2 output files (modified from python script)
suppa2_output <- list(Suppa2_input(mouse_suppa_dir, "WholeIsoSeq"),
                      Suppa2_input(human_suppa_dir, "AdultCTX"),
                      Suppa2_input(human_suppa_dir, "FetalCTX"),
                      Suppa2_input(human_suppa_dir, "HumanCTX"),
                      Suppa2_input(human_suppa_dir, "FetalHIP"),
                      Suppa2_input(human_suppa_dir, "FetalSTR"))

names(suppa2_output) <- c("WholeIsoSeq","AdultCTX","FetalCTX","HumanCTX","FetalHIP","FetalSTR") # also prefix of files



### Other AS events (personal script) ###################################################################

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Alternative_Splicing/AS_Events.R")
sqanti_gtf <- lapply(sqanti.gtf.names.files, 
                     function(x) read.table(x, as.is = T) %>% 
                       .[,c("V1","V3","V4","V5","V12")] %>% 
                       `colnames<-`(c("chromosome","source","start","end","isoform")) %>% 
                       mutate(gtf_coordinates = paste0(chromosome,":",start,"-",end)) %>% 
                       mutate(isoform = word(isoform,c(1), sep = ";"))
)
# names(sqanti.gtf.names.files): [1] "Mouse" "Adult" "Fetal" "Human"
names(sqanti_gtf) <- c("Mouse","Human (Adult)","Human (Fetal)","Human")

# Filter monoexonic transcripts as no alternative splicing events with just one exon, thus merge all.x = T
class.files <- lapply(class.files, function(x) x %>% filter(subcategory != "mono-exon") 
                      %>% .[,c("isoform","associated_gene","structural_category","strand","subcategory","Sample")])
names(class.files) <- c("Mouse","Human (Adult)","Human (Fetal)","Human")

annotated_sqanti_gtf <- lapply(names(sqanti_gtf), function(x) merge(class.files[[x]], sqanti_gtf[[x]], by = "isoform", all.x = T))
names(annotated_sqanti_gtf) <- c("Mouse","Human (Adult)","Human (Fetal)","Human")
names(suppa2_output)[1:4] <- c("Mouse","Human (Adult)","Human (Fetal)","Human")

dataset_IR <- list()
dataset_AF_AL <- list()
dataset_A5_A3 <- list()
dataset_tally <- list()
count = 1 
for(i in c("Mouse","Human (Adult)","Human (Fetal)","Human")){
  dataset_IR[[count]] <- class.files[[i]][class.files[[i]]$subcategory == "intron_retention",] %>% mutate(Event = "IR") %>% mutate(Sample = paste(i))
  dataset_AF_AL[[count]] <- AF_AL(annotated_sqanti_gtf[[i]])
  dataset_A5_A3[[count]] <- A5_A3(annotated_sqanti_gtf[[i]],250)
  dataset_tally[[count]] <- lapply(list(suppa2_output[[i]],dataset_AF_AL[[count]],dataset_A5_A3[[count]],dataset_IR[[count]]), 
                                   function(x) x %>% group_by(associated_gene, Event) %>% 
                                     tally()) %>% bind_rows() %>% mutate(Sample = paste(i))
  count = count + 1 
}
names(dataset_AF_AL) <- c("Mouse","Human (Adult)","Human (Fetal)","Human")
names(dataset_A5_A3) <- c("Mouse","Human (Adult)","Human (Fetal)","Human")
names(dataset_IR) <- c("Mouse","Human (Adult)","Human (Fetal)","Human")

plot_distribution <- function(dataset1, dataset2){
  cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  splicing_events <- dataset_tally %>% bind_rows()  
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

# Write output table of splicing events across all datasets
output <- dataset_tally %>% bind_rows() %>% `colnames<-`(c("Gene", "Event", "No_of_Events", "Sample")) %>% as.data.frame()
write.csv(output, paste0(output_table_dir,"/AS_IR/ALL_SUPPA2_Genes_Output_Updated.csv"))

#output2 <- dataset_tally %>% bind_rows() %>% as.data.frame()
#write.csv(output2, paste0(output_table_dir,"/AS_IR/ALL_SUPPA2_Genes_Output_Updated2.csv"))

#plot_distribution("Human","Mouse")
#plot_distribution("Human (Adult)","Human (Fetal)")

#write.csv(splicing_events_gene_present, paste0(output_table_dir,"/AS_IR/Splicing_Events_Genes.csv"))
#write.csv(splicing_events_tally_present, paste0(output_table_dir,"/AS_IR/Splicing_Events_Number.csv"))

#suppA <- plot_distribution("Human","Mouse")[[1]]
#suppB <- plot_distribution("Human (Adult)","Human (Fetal)")[[1]]
#plot_grid(suppA, suppB, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", ncol = 2)

### Annotated vs Novel - Mouse only #############################################
subset_isoform_events <- function(type){
  # Subset suppa events into novel and known 
  if(type == "known"){suppa2 <- suppa2_output$Mouse[!grepl("novel",suppa2_output$Mouse$alternative_transcripts),]}else 
      if(type == "novel"){suppa2 <-suppa2_output$Mouse[grepl("novel",suppa2_output$Mouse$alternative_transcripts),]}else{print("NA")}
  
  if(type == "known"){cate = c("FSM","ISM")}else if(type == "novel"){cate = c("NNC","NIC","Fusion","Genic\nGenomic")}else{print("NA")}
  
  print(cate)
  AS_events <- lapply(list(dataset_A5_A3$Mouse, dataset_AF_AL$Mouse, dataset_IR$Mouse), 
                      function(x) x %>% filter(as.factor(structural_category) %in% cate))
  tally <- bind_rows(AS_events) %>% group_by(Event) %>% tally() %>% 
    bind_rows(suppa2 %>% group_by(Event) %>% tally()) %>% mutate(class_isoform = type)
  
  return(tally)
}


anno_novel_events <- function(event, stats){
  cat("##################### Processing:",event,"\n\n")
  stats$known <- as.numeric(as.character(stats$known))
  stats$novel <- as.numeric(as.character(stats$novel))
  
  dat <- bind_rows(stats[stats$Event == event,], data.frame(Event = c("non-event"),
                                                            novel=sum(stats[stats$Event != event,"novel"]),
                                                            known = sum(stats[stats$Event != event,"known"]))) %>% 
    remove_rownames %>% column_to_rownames(var="Event")
  
  dat <- dat[,c("novel","known")]
  print(dat)
  #print(fisher.test(dat, alternative = "greater"))
  print(fisher.test(dat, alternative = "less"))
}

anno_novel_AS <- bind_rows(subset_isoform_events("novel"),subset_isoform_events("known")) %>% spread(.,class_isoform,n) 
for(e in c("A3","A5","AF","AL","IR","MX","SE")){anno_novel_events(e,anno_novel_AS)}
write.csv(anno_novel_AS, paste0(output_table_dir,"/AS_IR/Mouse_AS_NovelAnnoIsoforms.csv"))


knowngenes_AS <- dataset_tally[[1]] %>% filter(!grepl("NOVEL",associated_gene)) %>% group_by(Event) %>% tally(n) %>% mutate(perc = n/ sum(n) * 100)
write.csv(knowngenes_AS, paste0(output_table_dir,"/AS_IR/Mouse_AS_KnownGenes.csv"))


