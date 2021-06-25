# Szi Kay Leung
# Script to call functions for plots for IsoSeq Paper (2020)


source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Figures/All_Plots_Functions.R")

### Variables (Input, Output) ########################################################################

output_plot_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Figures"
output_table_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Tables"
output_corr_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Correlations"
reference_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/"
aaron_reference_dir <- "/gpfs/mrc0/projects/Research_Project-MRC148213/ISOSEQ/reference/"

# Input Files from Input_Variables.R
list_ccs_input_dir()        # Directory containing CCS lengths 
sqanti_files()              # SQANTI2 filtered classification files (Reference genome alignment)
junc_files()                # SQATNI2 filtered junction files
lncrna_class_files()        # SQANTI2 filtered classification files (Reference lncRNA alignment)
ERCC_sqanti_files()         # SQANTI2 ERCC mouse 
rarefaction_files()        # Rarefaction curves plotted from Liz's scripts
rnaseq_sqanti_files()       # rnaseq transcriptome stringtie


# Homologous genes between human and mouse 
# mouse genome informatics syntenic gene list (http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt)
homologs_genes <- read.csv(paste0(reference_dir, "mousehumangeneconversion.csv"))
homologs_genes$mouse <- as.factor(toupper(homologs_genes$mouse))

# read in mouse and human reference genome (adjusted in table)
# datawrangle to create separate new column for gene names taken taken from "Genesymbol" (i.e. remove white space, capitalise)
mm10 <- read.csv(paste0(reference_dir, "gencode.vM22_gene_annotation_table.csv")) %>% mutate(Gene = toupper(word(.$GeneSymbol,c(3),  sep = fixed (' '))))
hg38 <- read.table(paste0(aaron_reference_dir, "gencode.v31_gene_annotation_table.txt"), sep = "\t", header = T) %>% 
  mutate(GeneSymbol = gsub(" ", "", GeneSymbol)) %>% mutate(Gene = GeneSymbol)


### QC Plots ########################################################################################
### CCS 
#ccs_length_distribution ("NA","NA","Human_Mouse_Merge") # main
ccssuppA <- ccs_length_distribution ("NA","NA","Fetal_Adult_Merge")
ccssuppB <- ccs_length_distribution(ccs_input_dir$Adult,"Adult", "Individual")
ccssuppC <- ccs_length_distribution(ccs_input_dir$Fetal,"Fetal", "Individual")
ccssuppD <- ccs_length_distribution(ccs_input_dir$Mouse,"Mouse", "Individual")

### Rarefaction
rarefaction_plots <- rarefaction_distribution()
rrsuppA <- rarefaction_distribution()[[2]]
rrsuppB <- rarefaction_distribution()[[3]][["Human"]]
rrsuppC <- rarefaction_distribution()[[3]][["Adult"]]
rrsuppD <- rarefaction_distribution()[[3]][["Fetal"]]
rrsuppE <- rarefaction_distribution()[[3]][["Mouse"]]

### GO results from top 500 most abundant genes 
#most_abundant500Genes()

### RNASeq & Junction Support###################################################################
# Correlation
# p1 <- Gene Level Correlation
# p2 <- Gene Level Correlation with cutoff expression (Log Isoseq Gene TPM > 2.5)
# p3 <- Isoform Level Correlation
# Mouse correlation
#rnaseq_isoseq_correlation(class.files$Mouse)
# Fetal correlation
fcrnaseqsuppA <- rnaseq_isoseq_correlation(class.files$Fetal)[[1]]
fcrnaseqsuppB <- rnaseq_isoseq_correlation(class.files$Fetal)[[3]]

# Junction Support - main figure 
# rnaseq_junction_support()

# Junction Support novel transcripts of annotated genes 
junction_support_subset_plot <- junction_support_subset("rnaseq")

# Splice sites of all transcripts in annotated genes
splice_site_frequency_output <- splice_site_frequency()

# RNASeq coverage of novel transcripts in annotated genes
#pdf(paste0(output_plot_dir,"/rnaseq_coverage_noveltranscripts_annotated_genes.pdf"), width = 11, height = 8.5)
#fetal_mouse_rnaseq_distribution()
#dev.off()

rnaseq_isoseq_countsplot <- rnaseq_isoseq_counts(class.files$Mouse)

### Isoforms Descriptive Plots ######################################################################
# Number of Isoforms Group
num_iso <- no_of_isoforms_group()

# Number of Isoforms Group - Boxplot
bpnumisosuppA <- no_of_isoforms_gene()[[1]]
bpnumisosuppB <- no_of_isoforms_gene()[[2]]

# Structural categories 
#structural_categories_distribution("Individual")
#structural_categories_distribution("Merged")


# Length vs Number of Isoform correlation
length_iso_corrsuppA <- exon_length_isoform_correlation(class.files$Human)[[4]]
length_iso_corrsuppB <- exon_length_isoform_correlation(class.files$Mouse)[[4]]
length_iso_corrsuppC <- exon_length_isoform_correlation(class.files$Human)[[5]]
length_iso_corrsuppD <- exon_length_isoform_correlation(class.files$Mouse)[[5]]


# Exon vs Number of Isoform correlation
exon_iso_corrsuppA <- exon_length_isoform_correlation(class.files$Human)[[6]]
exon_iso_corrsuppB <- exon_length_isoform_correlation(class.files$Mouse)[[6]]
exon_iso_corrsuppC <- exon_length_isoform_correlation(class.files$Human)[[7]]
exon_iso_corrsuppD <- exon_length_isoform_correlation(class.files$Mouse)[[7]]

# Distance to cage peak 
all_hist_peaks_merged <- all_hist_peaks("All_Transcripts")

# novel transcripts of annotated genes
# List of Plots: cage peaks, tss, tts, and then: "Fetal","Mouse","Adult","Human"
novel_transcripts_peaks <- all_hist_peaks("Novel_Transcripts_Annotated_Genes")
NovIso_AnnoGene_cagesuppA <- novel_transcripts_peaks[[4]]
NovIso_AnnoGene_cagesuppB <- novel_transcripts_peaks[[5]]
NovIso_AnnoGene_cagesuppC <- novel_transcripts_peaks[[6]]

### Annotated Genes Plots ##################################################################
# Percentage of annotated trancripts (annotated genes) across categories 
AnnoGene_NovelIso_catesuppA <- novel_transcripts_annotated_genes_cate()[[1]]
AnnoGene_NovelIso_catesuppB <- novel_transcripts_annotated_genes_cate()[[2]]

# Novel vs Annotated Transcript Expression and Length in Annotated Gene
AnnvsNovIso_AnnoGenesuppA <- novel_annotated_transcript_expression("Mann_Whitney_Transcript_Expression_Novel_Transcripts.txt")
AnnvsNovIso_AnnoGenesuppB <- novel_annotated_transcript_length("Mann_Whitney_Transcript_Length_Novel_Transcripts.txt")
AnnvsNovIso_AnnoGenesuppC <- novel_annotated_transcript_exon("Mann_Whitney_Transcript_Exons_Novel_Transcripts.txt")


### Numbers of Isoforms Comparisons ##############################################
# Human and mouse comparison numbers, # Adult and Fetal comparison numbers
human_mouse <- human_mouse_comparison_num("homology_yes")
adult_fetal <- adult_fetal_comparison_num()
dataset_comp_suppA <- human_mouse[[1]]
dataset_comp_suppB <- adult_fetal[[1]]
dataset_comp_suppC <- human_mouse[[4]]
dataset_comp_suppD <- adult_fetal[[4]]

### Venn Diagram of Genes ##################################################################
# Venn diagram of all genes
vdgenessuppA <- venn_diagram_plot_twocircles(human_mouse[[5]]$associated_gene, human_mouse[[6]]$associated_gene, "Human","Mouse")
vdgenessuppB <- venn_diagram_plot_twocircles(adult_fetal[[5]]$associated_gene, adult_fetal[[6]]$associated_gene, "Human (Adult)","Human (Fetal)")

# Venn Diagram of novel transcripts of annotated genes in adult vs fetal 
#pdf(paste0(output_plot_dir,"/adultfetal_noveltranscripts_annotated_genes.pdf"), width = 11, height = 8.5)
#overlap_noveltranscripts_annotated_genes()
#dev.off()

# Venn Diagram of annotated genes in adult vs fetal 
#pdf(paste0(output_plot_dir,"/adultfetal_annotated_genes.pdf"), width = 11, height = 8.5)
#overlap_annotated_genes ()
#dev.off()

### lncRNA #################################################################
lncRNA_plots <- lncRNA()

# Length
lncRNA_length_suppA <- lncRNA_plots$Human[[1]]
lncRNA_length_suppB <- lncRNA_plots$Mouse[[1]]
lncRNA_length_suppC <- lncRNA_plots$Human[[2]]
lncRNA_length_suppD <- lncRNA_plots$Mouse[[2]]
plot_grid(lncRNA_length_suppA, lncRNA_length_suppB, lncRNA_length_suppC, lncRNA_length_suppD, 
          labels = "auto", label_size = 30, label_fontfamily = "ArialMT", ncol = 2)

# Expression
lncRNA_exp_suppA <- lncRNA_plots$Human[[5]]
lncRNA_exp_suppB <- lncRNA_plots$Mouse[[5]]
lncRNA_exp_suppC <- lncRNA_plots$Human[[7]]
lncRNA_exp_suppD <- lncRNA_plots$Mouse[[7]]

# ORF
lncRNA_orf_suppA <- lncRNA_plots$Human[[6]]
lncRNA_orf_suppB <- lncRNA_plots$Mouse[[6]]


### IR and NMD #################################################################
# Intron Retention Venn Diagram
# run function to output genes with intron retention transcripts
IR_Genes <- venn_diagram_IR()
# draw venn diagram for human vs mouse (homologs genes only)
human_mouse_IR_venn <- venn_diagram_plot_twocircles(IR_Genes$human_homologs$associated_gene, IR_Genes$mouse_homologs$associated_gene,"Human","Mouse")
# draw venn diagram for human adult vs human cortex
adult_fetal_IR_venn <- venn_diagram_plot_twocircles(IR_Genes$adult$associated_gene, IR_Genes$fetal$associated_gene, "Human (Adult)","Human (Fetal)")

#  NMD_vs_NonNMD Expression
nmd_expression <- NMD_vs_NonNMD(class.files)

# IR_NMD_run Venndiagram of NMD in IR transcripts
nmd_plots <- IR_NMD_run()
IR_nmd_suppA <- nmd_plots$Human
IR_nmd_suppB <- nmd_plots$Adult
IR_nmd_suppC <- nmd_plots$Fetal
IR_nmd_suppD <- nmd_plots$Mouse

# RNASeq expression of intron-retained transcripts 
#IR_rnaseq_support()

# IR rate 
IR_rate_plots <- IR_rate()
#suppA <- IR_rate_plots[[2]]
#suppB <- IR_rate_plots[[3]]
#suppC <- IR_rate_plots[[4]]
#plot_grid(suppA, labels = c("a"), label_size = 30, label_fontfamily = "ArialMT")
#plot_grid(suppB, labels = c("b"), label_size = 30, label_fontfamily = "ArialMT")
#plot_grid(suppC, labels = c("c"), label_size = 30, label_fontfamily = "ArialMT")



### PolyA ##################################################################
# PolyA 
#pdf(paste0(output_plot_dir,"/PolyA_Distribution_Motif.pdf"), width = 11, height = 8.5)
#polyA(class.names.files)
#polyA_freq(class.names.files)
#dev.off()

### Novel Genes ##################################################################
# Length of Novel vs Annotated genes 
NovelGenes_lengthsuppA <- length_novel_annotated_genes()[[1]]
# Expression of Novel vs Annotated genes 
NovelGenes_lengthsuppB <- expression_novel_annotated_genes()[[2]]

### Alternative Splicing ##################################################################
# SUPPA2 
fetal_adult_suppa2_output <- fetal_adult_suppa2()
human_mouse_suppa2_output <- human_mouse_suppa2("homology")
human_mouse_suppa2_output_no_homology <- human_mouse_suppa2("not_homology")

# overlap of genes with SUPPA2
AS_venn_suppA <- human_mouse_suppa2_output[[1]]
AS_venn_suppB <- fetal_adult_suppa2_output[[1]]

# common overlap of splicing events between datasets
commonASEvents_suppA <- human_mouse_suppa2_output[[2]]
commonASEvents_suppB <- fetal_adult_suppa2_output[[2]]

# Number of Splicing Events per AS Gene 
numASEvents <- SUPPA2_NumAS_PerGene()

HumanMouse_AS <- AS_genes_events("Human","Mouse")
AdultFetal_AS <- AS_genes_events("Human (Adult)","Human (Fetal)")
#SUPPA2_events_genes_plot()

### ONT Validation ###############################################################
#ONT_validation()

### ERCC ###############################################################
# run_ERCC(Ercc.class.file)


### Main figures ################################################################### 
empty_plot <- plot.new()
pdf(paste0(output_plot_dir,"/Figure1_R.pdf"), width = 11, height = 8.5)
ccs_length_distribution ("NA","NA","Human_Mouse_Merge")
rarefaction_plots[[1]]
dev.off()

pdf(paste0(output_plot_dir,"/Figure2_R.pdf"), width = 11, height = 8.5)
rnaseq_isoseq_correlation(class.files$Mouse)[[1]] + mymaintheme 
rnaseq_isoseq_correlation(class.files$Mouse)[[3]] + mymaintheme 
rnaseq_junction_support()[[2]] + mymaintheme + theme(legend.position = c(0.90,0.75))
dev.off()

pdf(paste0(output_plot_dir,"/Figure3_R.pdf"), width = 11, height = 8.5)
all_hist_peaks_merged[[4]] + mymaintheme + theme(legend.position = c(0.85,0.85))
all_hist_peaks_merged[[5]] + mymaintheme + theme(legend.position = c(0.85,0.85))
all_hist_peaks_merged[[6]] + mymaintheme + theme(legend.position = c(0.85,0.85))
num_iso[[1]] + mymaintheme  + theme(legend.position = c(0.85,0.85))
structural_categories_distribution("Merged")[[1]] + mymaintheme  + theme(legend.position = c(0.85,0.85))
AnnvsNovIso_AnnoGenesuppA[[1]] + mymaintheme + theme(legend.position = c(0.90,0.90))
dev.off()

pdf(paste0(output_plot_dir,"/Figure5_R.pdf"), width = 11, height = 8.5)
HumanMouse_AS[[2]] + mymaintheme + theme(legend.position = "bottom")
plot_grid(grobTree(AS_venn_suppA), grobTree(AS_venn_suppB), ncol = 1)
HumanMouse_AS[[3]] + mymaintheme + theme(legend.position = c(0.85,0.85))
plot_grid(grobTree(human_mouse_IR_venn),grobTree(adult_fetal_IR_venn), ncol = 1)
IR_rate_plots[[1]] + mymaintheme + theme(legend.position = c(0.85,0.85))
dev.off()


### Supplementary figures ###################################################################
pdf(paste0(output_plot_dir,"/Supplementary_Plots.pdf"), width = 11, height = 15)
# CCS
plot_grid(ccssuppA,ccssuppB,ccssuppC,ccssuppD + theme(legend.position = c(0.8,0.7)),labels = c("a", "b","c","d"),label_size = 30,label_fontfamily = "ArialMT") 
# Rarefaction
plot_grid(rrsuppA + theme(legend.position = c(0.8,0.2)),NULL,rrsuppB,rrsuppC,rrsuppD,rrsuppE, labels = c("a","","b","c","d","e"), label_size = 30, label_fontfamily = "ArialMT", ncol = 2, scale = 0.9)
# Venn diagram of overlap of genes
plot_grid(grobTree(vdgenessuppA),grobTree(vdgenessuppB),labels = c("a","b"),label_size = 30,label_fontfamily = "ArialMT",ncol = 1)
# GO results from top 500 most abundant genes
plot_grid(most_abundant500Genes(),empty_plot,ncol = 1,scale = 0.9)
# Fetal Correlation RNASeq
plot_grid(fcrnaseqsuppA, fcrnaseqsuppB,empty_plot,empty_plot,labels = c("a","b"),label_size = 30,label_fontfamily = "ArialMT", nrow = 3)
# Number of isoforms - boxplot
plot_grid(bpnumisosuppA, bpnumisosuppB, empty_plot,empty_plot,empty_plot,empty_plot, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
# Length vs number of isoforms correlation
plot_grid(length_iso_corrsuppA, length_iso_corrsuppB, length_iso_corrsuppC, length_iso_corrsuppD, empty_plot, empty_plot,labels = c("a","b","c","d"), label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
# Exon vs number of isoforms correlation
plot_grid(exon_iso_corrsuppA, exon_iso_corrsuppB, exon_iso_corrsuppC, exon_iso_corrsuppD,empty_plot, empty_plot, labels = c("a","b","c","d"), label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
# Percentage of annotated trancripts (annotated genes) across categories 
plot_grid(AnnoGene_NovelIso_catesuppA, AnnoGene_NovelIso_catesuppB, empty_plot,empty_plot,empty_plot,empty_plot, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", nrow = 3)
# Novel vs Annotated Transcript Expression and Length in Annotated Gene
plot_grid(AnnvsNovIso_AnnoGenesuppA[[1]], AnnvsNovIso_AnnoGenesuppA[[3]], empty_plot,empty_plot,labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", nrow = 2)
plot_grid(AnnvsNovIso_AnnoGenesuppB[[1]],AnnvsNovIso_AnnoGenesuppB[[3]],AnnvsNovIso_AnnoGenesuppC[[1]],AnnvsNovIso_AnnoGenesuppC[[3]], labels = c("a","b","c","d"), label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
# Cage Peak: Novel transcripts of annotated genes
plot_grid(NovIso_AnnoGene_cagesuppA, NovIso_AnnoGene_cagesuppB, empty_plot,empty_plot,empty_plot,empty_plot, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", nrow = 3)
# Junction Support novel transcripts of annotated genes 
plot_grid(junction_support_subset_plot[[1]][[2]],empty_plot, ncol = 1, scale = 0.9)
# ONT
plot_grid(ONT_validation(),empty_plot, ncol = 1, scale = 0.9)
# Length of Novel vs Annotated genes 
# Expression of Novel vs Annotated genes 
plot_grid(NovelGenes_lengthsuppA, NovelGenes_lengthsuppB,empty_plot,empty_plot,empty_plot,empty_plot,labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT",ncol = 2)
# LncRNA length
plot_grid(lncRNA_length_suppA, lncRNA_length_suppB, lncRNA_length_suppC, lncRNA_length_suppD, empty_plot,empty_plot,
          labels = c("a","b","c","d"), label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
# LncRNA expression
plot_grid(lncRNA_exp_suppA, lncRNA_exp_suppB,lncRNA_exp_suppC, lncRNA_exp_suppD, empty_plot,empty_plot, labels = c('a', 'b','c', 'd'), label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
# LncRNA ORF
plot_grid(lncRNA_orf_suppA, lncRNA_orf_suppB, empty_plot,empty_plot, empty_plot,empty_plot, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT",ncol = 2)
# AS Proportion of events
plot_grid(HumanMouse_AS[[1]], AdultFetal_AS[[1]], empty_plot,empty_plot, empty_plot,empty_plot, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
# common overlap of splicing events between datasets
plot_grid(commonASEvents_suppA,commonASEvents_suppB, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", ncol = 1)
# Number of Splicing Events per AS Gene 
plot_grid(numASEvents[[1]],empty_plot, ncol = 1, scale = 0.9)
#  NMD_vs_NonNMD Expression
plot_grid(nmd_expression$Human,nmd_expression$Mouse,empty_plot,empty_plot,empty_plot,empty_plot,labels = c("a","b"),label_size = 30,label_fontfamily = "ArialMT",ncol = 2)
# IR_NMD_run Venndiagram of NMD in IR transcripts
plot_grid(grobTree(IR_nmd_suppA),grobTree(IR_nmd_suppB),grobTree(IR_nmd_suppC),grobTree(IR_nmd_suppD),labels = c("a","b","c","d"),label_size = 30,label_fontfamily = "ArialMT")
# Human and mouse comparison numbers, # Adult and Fetal comparison numbers
plot_grid(dataset_comp_suppA, dataset_comp_suppC,dataset_comp_suppB, dataset_comp_suppD,empty_plot,empty_plot,labels = c("a","b","c","d"),label_size = 30,label_fontfamily = "ArialMT", ncol = 2)
# Splice sites of all transcripts in annotated genes
plot_grid(splice_site_frequency_output[[1]],splice_site_frequency_output[[2]],empty_plot, labels = c("a","b"),label_size = 30,label_fontfamily = "ArialMT",ncol = 1)
# Number of Isoforms in adult and fetal
plot_grid(num_iso[[2]],empty_plot,ncol = 1,scale = 0.9)
# Intron Retention Rate
plot_grid(IR_rate_plots[[5]],empty_plot,ncol = 1,scale = 0.9)
dev.off()

#   detach("package:plyr")
pdf(paste0(output_plot_dir,"/New_Supplementary_Plots.pdf"), width = 11, height = 15)
plot_grid(NovIso_AnnoGene_cagesuppA,NovIso_AnnoGene_cagesuppB,NovIso_AnnoGene_cagesuppC,labels = c("a","b","c"),label_size = 30,label_fontfamily = "ArialMT",ncol = 1)
rnaseq_isoseq_transcriptome(cuffrefmap_input,cufftmap_input)
whole_vs_targeted_plots()
dev.off()

## Addditional Plots ###################################################################
# iso_length(class.files$Mouse)
### RNASeq vs IsoSeq Defined transcriptome
#ggsave(file=paste0(output_plot_dir,"/RNAseq.pdf"), width = 210, height = 297, units = "mm")
#rnaseq_isoseq_transcriptome(cuffrefmap_input,cufftmap_input)
#dev.off()
## Number of isoforms per sample
#no_of_isoforms_sample(class.files$Mouse)[1]
#no_of_isoforms_sample(class.files$Mouse)[2]
