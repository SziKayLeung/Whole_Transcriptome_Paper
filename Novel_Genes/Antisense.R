
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Figures/All_Plots_Functions.R")

antisense <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Revised_Paper/Novel_Genes/Antisense"
sqanti_files()              # SQANTI2 filtered classification files (Reference genome alignment)

# Note one novel gene in human dataframe that is not antisense but classified as antisense (chr1 : 144716079 - 144758287)
class.files$Mouse[class.files$Mouse$associated_gene == "NOVELGENE_RPS2_AS","associated_gene"] <- "SNHG9"
class.files$Mouse[class.files$Mouse$associated_gene == "SNHG9","structural_category"] <- "FSM"

dat <- data.frame()
count=1
for(i in class.files){
  novel_genes <- i[grepl("NOVEL", i$associated_gene),]
  antisense_genes <- i[i$structural_category == "Antisense",]
  novel_antisense_genes <- i[grepl("NOVEL", i$associated_gene) & i$structural_category == "Antisense",]
  
  dat[1:4,count] <- rbind(nrow(novel_genes), nrow(antisense_genes), nrow(novel_antisense_genes), 
                          paste0(round(nrow(novel_antisense_genes)/ nrow(novel_genes) * 100, 2),"%"))
  dat[5,count] <- length(unique(novel_antisense_genes$associated_gene))
  
  colnames(dat)[count] <- names(class.files)[count]
  count = count + 1
}
rownames(dat) <- c("Number of Novel Transcripts from Novel Genes","Number of Antisense Transcripts", "Number of Antisense Transcripts from Genes",
                   "% of Antisense transcripts from Novel Genes","Number of Antisense, Novel Genes ")
dat


antisense_sense_pairs <- function(path_inputdf_novelgenes_antisense, path_inputdf_annotated_genes_overlap,path_inputdf_novelgenes_antisense_reference, dataset){
  
  antisense_transcripts <- class.files[[dataset]] %>% filter(structural_category == "Antisense")
  antisense_genes <- unique(antisense_transcripts$associated_gene)
  
  # Extract from bedtools the overlapping exonic regions
  bedtools_datawrangle <- function(path_inputdf, type_label){
    if(type_label == "Novel Antisense Genes"){
      input_df <- read.table(path_inputdf) %>% .[,c(1:4,6,8,14,ncol(.))] 
    }else{
      input_df <- read.table(path_inputdf, as.is = T, sep = "\t")  %>% .[,c(1:4,6,8,10,ncol(.))] 
    }
    
    colnames(input_df) <- c("chr","start","end","associated_gene","strand","exonic_structure","isoform","no_exon_overlap")
    output_df <- input_df %>%
      filter(no_exon_overlap > 0) %>% 
      filter(exonic_structure != "transcript") %>%  
      mutate(gtf_coordinates = paste0(chr,":",start,"-",end)) %>% mutate(type = type_label)
    
    return(output_df)
  }
  
  novel_genes_antisense_overlap <- bedtools_datawrangle(path_inputdf_novelgenes_antisense,"Novel Antisense Genes")
  annotated_genes_overlap <- bedtools_datawrangle(path_inputdf_annotated_genes_overlap,"Annotated Genes")
  novel_genes_antisense_genome_overlap <- bedtools_datawrangle(path_inputdf_novelgenes_antisense_reference,"Novel Antisense Genes Reference") 
  
  # Human
  # Unique genes that are not identified as overlap of annotated to novel and vice versa 
  # manually checked "CDKAL1","AC069288.1","CR769767.2"; and "CR769767.2" refers to "BX664615.2"
  # therefore the annotated_genes_overlap more accurate and picked up all
  # false positive from novel_genes_antisense_overlap = AC069288.1 & CDKAL1 (but are within gene)
  setdiff(unique(annotated_genes_overlap$associated_gene),word(unique(novel_genes_antisense_overlap$associated_gene),c(2), sep = "_"))
  setdiff(word(unique(novel_genes_antisense_overlap$associated_gene),c(2), sep = "_"),unique(annotated_genes_overlap$associated_gene))

  
  # Mouse
  # in annotated genes, but not novel genes: "Rfx4","Lgalsl_Gm12037","Trim39","Idua","Sh2b2"  
  # in novel genes, but not annotated genes: "Camk2b","Gm12037","Cog5","Prkcd","H2-M10.2","Atrnl1","Gpsm1","Mbnl1","Igsf21","Fgfrl1","Gm5050","Il16"     "Zbtb16" 
  # Manual Curation 
  # Rfx4 = already detected as Ric8b in novel, but also fusion antisense novel gene so picked up as unique to annotated
  # Lgalsl_Gm12037 = GM12037 = within gene but not overlapping exon
  # Trim39 = H2-M10.2
  # Idua = Fgfrl1
  # Sh2b2 = Gm5050
  # NOVELGENE_TM9SF4_AS overlap two genes
  # "NOVELGENE_RPS2_AS" not a novel gene PB.6773.1
  
  
  if (dataset == "Human"){
    within_gene_not_exon_overlap <- c("novelGene_CDKAL1_AS","novelGene_AC069288.1_AS","novelGene_LRRC27_AS",
                                      "novelGene_RN7SL786P_AS","novelGene_DGKB_AS","novelGene_BRD9P2_AS")
    novel_genes_antisense_overlap <- novel_genes_antisense_overlap %>% filter(!associated_gene %in% within_gene_not_exon_overlap )
    
  } else if (dataset == "Mouse"){
    within_gene_not_exon_overlap <-c("NOVELGENE_DUSP22_AS","NOVELGENE_NRROS_AS","NOVELGENE_ENTPD1_AS","NOVELGENE_ATRNL1_AS","NOVELGENE_ETL4_AS",
                                     "NOVELGENE_GPSM1_AS","NOVELGENE_LHX2_AS","NOVELGENE_DCAF12_AS","NOVELGENE_IGSF21_AS","NOVELGENE_RNF216_AS",
                                     "NOVELGENE_AKR1B10_AS","NOVELGENE_PARD3_AS","NOVELGENE_HIST3H2BA_AS")
    
    novel_genes_antisense_overlap <- novel_genes_antisense_overlap %>% filter(!toupper(associated_gene) %in% within_gene_not_exon_overlap)

  }else {
    print("dataset must be as human or mouse")
  }
  
  # list of genes detected from bedtools
  isoseq_detected <- rbind(as.data.frame(toupper(within_gene_not_exon_overlap)) %>% mutate(type = "isoseq_exon_not_overlap") %>% `colnames<-`(c("Gene", "type")),
                           as.data.frame(unique(toupper(novel_genes_antisense_overlap$associated_gene))) %>% 
                             mutate(type = "isoseq_exon_overlap") %>% `colnames<-`(c("Gene", "type")))
  
  # manually checked all genes
  all_novel <- rbind(isoseq_detected, as.data.frame(setdiff(antisense_genes,isoseq_detected$Gene)) %>% 
                       mutate(type = "NA") %>% `colnames<-`(c("Gene", "type")))
  
  # Reference
  # Detected list of genes (overlapping exon and within genes, but not overlapping)
  Isoseq_detected <- c(toupper(unique(novel_genes_antisense_overlap$associated_gene)),paste0("NOVELGENE_",toupper(within_gene_not_exon_overlap),"_AS"))
  # list of genes that are classified as antisense novel genes, but not detected as antisense, overlap or within gene from dataset
  # However note these novel genes could still be antisense to known genes, but just not that detected in our dataset 
  undetected_external_gene <- setdiff(antisense_genes, Isoseq_detected)
  # Genes that are not detected in PacBio Isoseq dataset as overlapping/within and also not detected in genome dataset 
  genome_undetected_external_gene <- setdiff(undetected_external_gene, unique(toupper(novel_genes_antisense_genome_overlap$associated_gene))) 
  # Genes that are not detected in PacBio Isoseq dataset, but detected as overlapping/within genome reference 
  # Manually check the genomes that are detected in genome reference whether overlaping or within gene only 
  genome_overlapping_within <- intersect(undetected_external_gene, unique(toupper(novel_genes_antisense_genome_overlap$associated_gene)))
  
  
  if (dataset == "Human"){
   
    novel_genes_antisense_genome_overlap_exonconfirmed <- c("NOVELGENE_RPEL1_AS","NOVELGENE_SLC16A11_AS",
                                                            "NOVELGENE_S1PR4_AS","NOVELGENE_AC093866.1_AS")
    
    novel_genes_antisense_genome_overlap_withingeneconfirmed <- c("NOVELGENE_LINC02291_AS","NOVELGENE_ATP2C2_AS",
                                                                  "NOVELGENE_AC018709.1_AS", "NOVELGENE_AC113167.1_AS",
                                                                  "NOVELGENE_E2F3_AS")
    
    pseudogene <- c("NOVELGENE_RASA4DP_AS","NOVELGENE_DUX4L50_AS")
    
    notwithin_gene <- c("NOVELGENE_KIF28P_AS","NOVELGENE_PIGZ_AS","NOVELGENE_AC092047.1_AS")
    
  } else if (dataset == "Mouse"){
    novel_genes_antisense_genome_overlap_exonconfirmed <- c("NOVELGENE_GM28577_AS", 
                                                            "NOVELGENE_CFAP44_AS", "NOVELGENE_SLIT3_AS","NOVELGENE_ASIC2_AS","NOVELGENE_GM49200_AS","NOVELGENE_LBHD1_AS",
                                                            "NOVELGENE_ADRA2A_AS",
                                                            "NOVELGENE_CHRNA1_AS","NOVELGENE_CFAP100_AS","NOVELGENE_CYP2B19_AS","NOVELGENE_LIPT2_AS","NOVELGENE_GM3985_AS",
                                                            "NOVELGENE_PPP2R3D_AS")
    
    novel_genes_antisense_genome_overlap_withingeneconfirmed <- c("NOVELGENE_GM48845_AS", "NOVELGENE_RGS6_AS", "NOVELGENE_CFAP299_AS", 
                                                                  "NOVELGENE_CPVL_AS", "NOVELGENE_FILIP1_AS","NOVELGENE_CDH23_AS","NOVELGENE_SDK2_AS",
                                                                  "NOVELGENE_LNCPPARA_AS","NOVELGENE_LSAMP_AS","NOVELGENE_GM14620_AS")
    
    notwithin_gene <- c("NOVELGENE_GM13192_AS","NOVELGENE_GM43833_AS")
    
    pseudogene <- numeric(0)
    
  }else {
    print("dataset must be as human or mouse")
  }
  
  antisense_genes <<- antisense_genes
  within_gene_not_exon_overlap  <<- within_gene_not_exon_overlap
  novel_genes_antisense_overlap <<- novel_genes_antisense_overlap
  annotated_genes_overlap <<- annotated_genes_overlap
  genome_overlapping_within <<- genome_overlapping_within
  novel_genes_antisense_genome_overlap <<- novel_genes_antisense_genome_overlap
  
  stats <- function(){
    
    dat <- data.frame()
    dat[1:7,1] <- c(length(antisense_genes), length(within_gene_not_exon_overlap), length(unique(novel_genes_antisense_overlap$associated_gene)),
                    length(novel_genes_antisense_genome_overlap_withingeneconfirmed), length(novel_genes_antisense_genome_overlap_exonconfirmed), length(pseudogene),
                    length(notwithin_gene))
    
    row.names(dat) <- c("Antisense Genes Detected",
                        "Antisense Genes within Another Detected gene, but not exonic overlap",
                        "Antisense Genes within Another Detected gene, and exonic overlap",
                        "Antisense genes within known gene (but not detected IsoSeq), but not exonic overlap",
                        "Antisense Genes within known gene (but not detected in Isoseq, and exonic overlap",
                        "Antisense Genes within pseudogene (but not detected in Isoseq",
                        "Antisense Genes not within or overlapping another gene")
    colnames(dat) <- c("Number")
    print(dat)
  }
  
  stats()
  
  p1 <- rbind(novel_genes_antisense_overlap[,c("exonic_structure","no_exon_overlap","type")], 
        annotated_genes_overlap[,c("exonic_structure","no_exon_overlap","type")]) %>% 
    ggplot(., aes(x = type, y = no_exon_overlap, fill = exonic_structure)) + geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(color = exonic_structure), position=position_jitterdodge(jitter.height = 0.1)) + 
    theme_bw() + 
    mytheme 
  

  
  final <- rbind(isoseq_detected,
                 as.data.frame(novel_genes_antisense_genome_overlap_withingeneconfirmed) %>% mutate(type = "genome_exon_not_overlap") %>% `colnames<-`(c("Gene", "type")),
                 as.data.frame(novel_genes_antisense_genome_overlap_exonconfirmed) %>% mutate(type = "genome_exon_overlap") %>% `colnames<-`(c("Gene", "type")),
                 as.data.frame(pseudogene) %>% mutate(type = "pseudogene_overlap") %>% `colnames<-`(c("Gene", "type")),
                 as.data.frame(notwithin_gene) %>% mutate(type = "no_gene_overlap") %>% `colnames<-`(c("Gene", "type")))
  
  
  #output <- list(novel_genes_antisense_overlap,annotated_genes_overlap)
  #names(output) <- c("novel_genes_antisense_overlap","annotated_genes_overlap")
  return(final)
}

human_antisense <- antisense_sense_pairs(paste0(antisense, "/combined.NovelGenes.Antisense.output"),
                      paste0(antisense, "/combined.NovelGenes.Antisense.output_collapsed"),
                      paste0(antisense, "/combined.NovelGenes.Antisense.genome.output"),"Human")

mouse_antisense <- antisense_sense_pairs(paste0(antisense, "/WholeIsoSeq.NovelGenes.Antisense.output"),
                      paste0(antisense, "/WholeIsoSeq.NovelGenes.Antisense.output_collapsed"),
                      paste0(antisense, "/WholeIsoSeq.NovelGenes.Antisense.genome.output"),"Mouse")


novelgenes <- read.csv(paste0(output_table_dir,"/Novel_Genes.csv"))
