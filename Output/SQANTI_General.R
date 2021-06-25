# Szi Kay Leung: sl693@exeter.ac.uk
# Functions to data wrangle output files from SQANTI2 (classification file, junction file, gtf, reason file)

# SQANTI_class_preparation <classification_file> 
# Aim: Read in classification file generated from SQANTI2
SQANTI_class_preparation <- function(class.file,standard) {
  data.class = read.table(class.file, header=T, as.is=T, sep="\t")
  rownames(data.class) <- data.class$isoform
  
  # (Liz) not sorting by expression
  #if (!all(is.na(data.class$iso_exp))){
  #  sorted <- data.class[order(data.class$iso_exp, decreasing = T),]
  #  FSMhighestExpIsoPerGene <- sorted[(!duplicated(sorted$associated_gene) & sorted$structural_category=="full-splice_match"),"isoform"]
  #  data.class[which(data.class$isoform%in%FSMhighestExpIsoPerGene),"RTS_stage"] <- FALSE
  #  write.table(data.class, file=class.file, row.names=FALSE, quote=F, sep="\t")
  #}
  
  xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
  xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
  
  legendLabelF1 <- levels(as.factor(data.class$coding));
  
  data.class$structural_category = factor(data.class$structural_category,
                                          labels = xaxislabelsF1, 
                                          levels = xaxislevelsF1,
                                          ordered=TRUE)
  
  data.FSMISM <- subset(data.class, structural_category %in% c('FSM', 'ISM'))
  data.FSM <- subset(data.class, (structural_category=="FSM" & exons>1))
  data.ISM <- subset(data.class, (structural_category=="ISM" & exons>1))
  
  
  
  # Label Empty blanks in associated_gene column as "Novel Genes_PB_<isoform_ID>"
  data.class[data.class$associated_gene == "",]
  data.class$associated_gene[data.class$associated_gene == ""] <- paste0("novelGene_PB.",
                                                                         word(data.class$isoform[data.class$associated_gene == ""],c(2), 
                                                                              sep = fixed ('.')
                                                                         ))
  
  # Create a new attribute called "novelGene"
  
  data.class$novelGene <- "Annotated Genes"
  data.class[grep("novelGene", data.class$associated_gene), "novelGene"] <- "Novel Genes"
  data.class$novelGene = factor(data.class$novelGene,
                                levels = c("Novel Genes","Annotated Genes"),
                                ordered=TRUE)
  
  # Create a new attribute called "exonCat"
  
  data.class[which(data.class$exons>1), "exonCat"] <- "Multi-Exon"
  data.class[which(data.class$exons==1), "exonCat"] <- "Mono-Exon"
  data.class$exonCat = factor(data.class$exonCat,
                              levels = c("Multi-Exon","Mono-Exon"),
                              ordered=TRUE)
  
  data.class$all_canonical = factor(data.class$all_canonical,
                                    levels = c("canonical","non_canonical"),
                                    ordered=TRUE)
  
  data.class$within_cage_peak = factor(data.class$within_cage_peak)
  data.class$within_cage_peak <- factor(data.class$within_cage_peak, c("True","False"))
  
  if(standard == "standard"){
    # relabel the sum of the samples FL due to demultiplexing 
    # starts with FL refer to columns with samples
    # can append total directly to FL as do not change the order of the rows
    dat <- data.class %>% dplyr::select(starts_with("FL.")) %>% mutate(total =  rowSums(.[1:ncol(.)]))
    data.class$FL <- dat$total
    
  
    # convert SQANTI FL to TPM (based on E.Tseng's SQANTI2.report https://github.com/Magdoll/SQANTI2)
    total_fl <- sum(data.class$FL, na.rm=T)
    #print(paste0("Total FL counts:", total_fl))
    data.class$ISOSEQ_TPM <- data.class$FL*(10**6)/total_fl
    data.class$Log_ISOSEQ_TPM <- log10(data.class$ISOSEQ_TPM)
  }
  
  print(paste0("Loading classification file:",class.file))
  #assign(data_class_output_file, data.class, envir=.GlobalEnv)
  return(data.class)
}

TPM_Calculation <- function(dat){
  # recalculate Isoseq TPM based on removal of monoexons 
  total_fl <- sum(dat$FL, na.rm=T)
  dat$ISOSEQ_TPM <- dat$FL*(10**6)/total_fl
  dat$Log_ISOSEQ_TPM <- log10(dat$ISOSEQ_TPM)
  return(dat)
}


# read_merge_gtf
# Aim: Extract coordinates with output sqanti2 filtered gtf file and merge with sqanti2 classification.txt using PBID
# Input: sqanti2.classification.gtf and sqanti2_classification file (already read) 
# Output: dataframe with 3 columns: V9 from input gtf file, genome coordinates, isoform (pacbio id)
read_merge_gtf <- function(gtf_input, sqanti2_class){
  # read in gff
  input <- read.delim(gtf_input, header=F, comment.char="#") %>%
    # filter only transcripts in column 3
    filter(V3 == "transcript") %>%
    # take only chromosome (V1), start coordinates (V4), end coordinates (v5)
    mutate(gtf_coordinates = paste(V1,":", V4,"-", V5)) %>%
    select(-(V1:V8))
  
  # create separate column for merging 
  input$isoform <- word(input$V9,2, sep = ";")
  input$isoform <- word(input$isoform,3, sep = " ")
  input$isoform <- gsub('"', '', input$isoform)
  
  final_merge <- merge(input, sqanti2_class, by = "isoform", all.y = TRUE)
  return(final_merge)
}

SQANTI_gene_preparation <- function(data_class_output_file){
  
  # ----------------------------------------------------------
  # Make "isoPerGene" which is aggregated information by gene
  #  $associatedGene - either the ref gene name or novelGene_<index>
  #  $novelGene      - either "Novel Genes" or "Annotated Genes"
  #  $FSM_class      - "A", "B", or "C"
  #  $geneExp        - gene expression info
  #  $nIso           - number of isoforms associated with this gene
  #  $nIsoCat        - splicing complexity based on number of isoforms
  # ----------------------------------------------------------
  
  if (!all(is.na(data_class_output_file$gene_exp))){
    isoPerGene = aggregate(data_class_output_file$isoform,
                           by = list("associatedGene" = data_class_output_file$associated_gene,
                                     "novelGene" = data_class_output_file$novelGene,
                                     "FSM_class" = data_class_output_file$FSM_class,
                                     "geneExp"=data_class_output_file$gene_exp),
                           length)
  } else {
    isoPerGene = aggregate(data_class_output_file$isoform,
                           by = list("associatedGene" = data_class_output_file$associated_gene,
                                     "novelGene" = data_class_output_file$novelGene,
                                     "FSM_class" = data_class_output_file$FSM_class),
                           length)
  }
  # assign the last column with the colname "nIso" (number of isoforms)
  colnames(isoPerGene)[ncol(isoPerGene)] <- "nIso"
  
  
  isoPerGene$FSM_class2 = factor(isoPerGene$FSM_class, 
                                 levels = c("A", "B", "C"), 
                                 labels = c("MonoIsoform Gene", "MultiIsoform Genes\nwithout expression\nof a FSM", "MultiIsoform Genes\nexpressing at least\none FSM"), 
                                 ordered=TRUE)
  
  isoPerGene$novelGene = factor(isoPerGene$novelGene, 
                                levels = c("Annotated Genes", "Novel Genes"), 
                                ordered=TRUE)
  
  # Altered to extend out the graphs
  isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3,5,7,9, max(isoPerGene$nIso)+1), labels = c("1", "2-3", "4-5", "6-7","8-9",">=10"))
  Sample_Type <- data_class_output_file$Sample[1]
  isoPerGene$Sample <- Sample_Type
  
  return(isoPerGene)
  #assign(isoPerGene_output_file, isoPerGene, envir=.GlobalEnv)
}

# SQANTI_reason_preparation <reason_file> 
# Aim: Read in reason file generated from SQANTI2 of rationale for filtering isoforms in SQANTI2 filter
SQANTI_reason_preparation <- function(reason.file){
  reason.class = read.table(reason.file, header=T, as.is=T, sep=",")
  return(reason.class)
}

# SQANTI_junction_preparation <junction_file> 
# Aim: Read in junction file generated from SQANTI2 
SQANTI_junction_preparation <- function(junc.file){
  data.junction <- read.table(junc.file, header=T, as.is=T, sep="\t")
  
  # (Liz) don't sort junction file by expression
  #if (!all(is.na(data.class$iso_exp))){
  #  data.junction[which(data.junction$isoform%in%FSMhighestExpIsoPerGene),"RTS_junction"] <- FALSE
  #  write.table(data.junction, file=junc.file, row.names=FALSE, quote=F, sep="\t")
  #}
  
  # make a unique identifier using chrom_strand_start_end
  data.junction$junctionLabel = with(data.junction, paste(chrom, strand, genomic_start_coord, genomic_end_coord, sep="_"))
  
  data.junction$SJ_type <- with(data.junction, paste(junction_category,canonical,"SJ", sep="_"))
  data.junction$SJ_type <- factor(data.junction$SJ_type, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                  labels=c("Known\nCanonical ", "Known\nNon-canonical ", "Novel\nCanonical ", "Novel\nNon-canonical "), order=T)
  
  #data.junction$structural_category = data.class[data.junction$isoform, "structural_category"]
  
  uniqJunc <- unique(data.junction[,c("junctionLabel", "SJ_type", "total_coverage")]);
  uniqJuncRTS <- unique(data.junction[,c("junctionLabel","SJ_type", "RTS_junction")])
  
  return(data.junction)
}
