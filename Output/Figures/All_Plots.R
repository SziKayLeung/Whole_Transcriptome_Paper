# Szi Kay Leung
# Script to call functions for plots for IsoSeq Paper (2020)


source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Output/Figures/All_Plots_Functions.R")

# do not output log files for venn diagrams
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

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
human_mouse_ccs <- ccs_length_distribution ("NA","NA","Human_Mouse_Merge") # main
ccssuppA <- ccs_length_distribution ("NA","NA","Fetal_Adult_Merge")
ccssuppB <- ccs_length_distribution(ccs_input_dir$Adult,"Adult", "Individual")
ccssuppC <- ccs_length_distribution(ccs_input_dir$Fetal,"Fetal", "Individual")
ccssuppD <- ccs_length_distribution(ccs_input_dir$Mouse,"Mouse", "Individual")

### Rarefaction
rarefaction_plots <- rarefaction_distribution()
rrsupp <- rarefaction_distribution()[[1]]
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
mcrnaseqsuppA <- rnaseq_isoseq_correlation(class.files$Mouse)[[1]]
mcrnaseqsuppB <- rnaseq_isoseq_correlation(class.files$Mouse)[[3]]
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

# RNASeq defined transcriptome 
RNAseq_defined <- rnaseq_isoseq_transcriptome(cuffrefmap_input,cufftmap_input)

### Isoforms Descriptive Plots ######################################################################
# Number of Isoforms Group
num_iso <- no_of_isoforms_group()

# Number of Isoforms Group - Boxplot
bpnumisosuppA <- no_of_isoforms_gene()[[1]]
bpnumisosuppB <- no_of_isoforms_gene()[[2]]

# CPAT scores 
Human_CPAT_plots <- cpat_plots(Human_CPAT,"Human")
Human_CPAT_Plots_Cate <- CPAT_plots_cate(class.files$Human,Human_CPAT,"Human")

# Structural categories 
Human_structural <- structural_categories_distribution("Human")
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
noveltranscripts_hist_peaks_merged <- all_hist_peaks("Novel_Transcripts_Annotated_Genes")

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

Human_AS <- AS_genes_events("Human")
Mouse_AS <- AS_genes_events("Mouse")
Human_Adult_AS <- AS_genes_events("Human (Adult)")
Human_Fetal_AS <- AS_genes_events("Human (Fetal)")

#AdultFetal_AS <- AS_genes_events("Human (Adult)","Human (Fetal)")
#SUPPA2_events_genes_plot()

### ONT Validation ###############################################################
#ONT_validation()

### ERCC ###############################################################
ERCC <- run_ERCC(Ercc.class.file)


### Main figures ################################################################### 
# Figure 1
pdf(paste0(output_plot_dir,"/Figure1_R_plots.pdf"), width = 10, height = 10)
top_row <- plot_grid(human_mouse_ccs, all_hist_peaks_merged[[4]], nrow = 1, labels = c("A","B"),label_fontfamily = "ArialMT")
bottom_row <- plot_grid(Human_CPAT_plots[[1]],Human_CPAT_plots[[2]],num_iso[[1]], nrow = 1,  labels = c("C","D","E"), rel_widths = c(1.1,1,1.3),label_fontfamily = "ArialMT")
plot_grid(top_row,bottom_row, ncol = 1, scale = 0.9) 
dev.off()

# Figure 2
pdf(paste0(output_plot_dir,"/Figure2_R_plots.pdf"), width = 8.5, height = 11)
top_figure <- plot_grid(rasterGrob(Figure2a),Human_structural, rel_widths = c(1,0.6), scale = 0.9, labels = c("A","B"))
mid_figure <- plot_grid(rasterGrob(Figure2c),rasterGrob(Figure2d), nrow = 1,scale = 0.95,labels = c("C","D"))
bottom_figure <- plot_grid(rasterGrob(Figure2e),NULL,rel_widths = c(1,0),scale = 0.95, labels = c("E"))
plot_grid(top_figure,mid_figure,bottom_figure, nrow = 3)
dev.off()

# Figure 3; Tracks 
# Figure 4; Tracks

pdf(paste0(output_plot_dir,"/Figure5_R.pdf"), width = 8.5, height = 11)
top_row = plot_grid(NULL,Human_AS[[1]],Human_AS[[2]] + theme(legend.position="none"), nrow = 1, scale = 0.9, rel_widths = c(0.5,0.25,0.25), labels = c("a","b","c"))
mid_row = plot_grid(NULL,numASEvents[[2]], nrow = 1, scale = 0.9, rel_widths = c(0.6,0.4), labels = c("d","e"))
bottom_row = plot_grid(NULL,labels = c("f"))
bottom_row2 = plot_grid(NULL,labels = c("g"))
plot_grid(top_row,mid_row,bottom_row,bottom_row2,nrow = 4)
dev.off()


### Supplementary figures ###################################################################
pdf(paste0(output_plot_dir,"/Supplementary_Plots.pdf"), width = 11, height = 15)
# Rarefaction
plot_grid(rrsupp, rrsuppA + theme(legend.position = c(0.8,0.2)),rrsuppB,rrsuppC,rrsuppD,rrsuppE, labels = c("a","b","c","d","e","f"), label_size = 30, label_fontfamily = "ArialMT", ncol = 2, scale = 0.9)

# CCS
plot_grid(ccssuppA,ccssuppB,ccssuppC,ccssuppD + theme(legend.position = c(0.8,0.7)),NULL,NULL,labels = c("a", "b","c","d"),label_size = 30,label_fontfamily = "ArialMT", ncol = 2) 

# GO results from top 500 most abundant genes
#top_row <-plot_grid(most_abundant500Genes(),labels = c("a"),scale = 0.9,label_size = 30, label_fontfamily = "ArialMT")
#bottom_row <-plot_grid(ERCC[[2]],ERCC[[3]],labels = c("b","c"),scale = 0.9,label_size = 30, label_fontfamily = "ArialMT")
#plot_grid(top_row,bottom_row, nrow = 2)
# Fetal and Mouse Correlation RNASeq
plot_grid(fcrnaseqsuppA, fcrnaseqsuppB,mcrnaseqsuppA, mcrnaseqsuppB,ERCC[[2]],ERCC[[3]],labels = "auto",label_size = 30,label_fontfamily = "ArialMT", nrow = 3)

# cage peaks, TTS, TSS
plot_grid(all_hist_peaks_merged[[4]],noveltranscripts_hist_peaks_merged[[4]],all_hist_peaks_merged[[5]],noveltranscripts_hist_peaks_merged[[5]],all_hist_peaks_merged[[6]],noveltranscripts_hist_peaks_merged[[6]],labels = "auto",label_size = 30,label_fontfamily = "ArialMT", ncol = 2, scale = 0.9)

# Length vs number of isoforms correlation
# Exon vs number of isoforms correlation
plot_grid(length_iso_corrsuppA, length_iso_corrsuppB, length_iso_corrsuppC, length_iso_corrsuppD,exon_iso_corrsuppA, exon_iso_corrsuppB, exon_iso_corrsuppC, exon_iso_corrsuppD, labels = "auto", label_size = 30, label_fontfamily = "ArialMT", ncol = 2, scale = 0.9)

# Percentage of annotated trancripts (annotated genes) across categories 
#plot_grid(AnnoGene_NovelIso_catesuppA, AnnoGene_NovelIso_catesuppB, NULL,NULL,NULL,NULL, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", nrow = 3)
# Novel vs Annotated Transcript Expression and Length in Annotated Gene
plot_grid(AnnoGene_NovelIso_catesuppA, AnnoGene_NovelIso_catesuppB,AnnvsNovIso_AnnoGenesuppA[[1]], AnnvsNovIso_AnnoGenesuppA[[3]], AnnvsNovIso_AnnoGenesuppB[[1]], AnnvsNovIso_AnnoGenesuppB[[3]],AnnvsNovIso_AnnoGenesuppC[[1]],AnnvsNovIso_AnnoGenesuppC[[3]], labels = "auto", label_size = 30, label_fontfamily = "ArialMT", ncol = 2, scale = 0.9)

# Venn diagram of overlap of genes
plot_grid(grobTree(vdgenessuppA),grobTree(vdgenessuppB),dataset_comp_suppA,dataset_comp_suppB,dataset_comp_suppC,dataset_comp_suppD,labels = "auto",label_size = 30,label_fontfamily = "ArialMT",ncol = 2, scale = 0.9)

# RNASeq 
# Cage Peak: Novel transcripts of annotated genes
#plot_grid(NovIso_AnnoGene_cagesuppA, NovIso_AnnoGene_cagesuppB, NULL,NULL,NULL,NULL, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", nrow = 3)
# Junction Support novel transcripts of annotated genes 
#plot_grid(junction_support_subset_plot[[1]][[2]],NULL, ncol = 1, scale = 0.9)
# ONT
#plot_grid(ONT_validation(),NULL, ncol = 1, scale = 0.9)
# Length of Novel vs Annotated genes 
# Expression of Novel vs Annotated genes 
#plot_grid(NovelGenes_lengthsuppA, NovelGenes_lengthsuppB,NULL,NULL,NULL,NULL,labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT",ncol = 2, scale = 0.9)

# RNASeq defined transcriptome 
plot_grid(RNAseq_defined[[1]],RNAseq_defined[[2]],RNAseq_defined[[3]],labels = c("a","b","c"),scale = 0.9,label_size = 30, label_fontfamily = "ArialMT", ncol = 2)

# LncRNA length# LncRNA expression
plot_grid(lncRNA_length_suppA, lncRNA_length_suppB, lncRNA_length_suppC, lncRNA_length_suppD, lncRNA_exp_suppA, lncRNA_exp_suppB,lncRNA_exp_suppC, lncRNA_exp_suppD, labels = "auto", label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
# LncRNA ORF
#plot_grid(lncRNA_orf_suppA, lncRNA_orf_suppB, NULL,NULL, NULL,NULL, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT",ncol = 2)

# AS Proportion of events
# common overlap of splicing events between datasets
#plot_grid(HumanMouse_AS[[1]], AdultFetal_AS[[1]], NULL,NULL, NULL,NULL, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
plot_grid(commonASEvents_suppA,commonASEvents_suppB,NULL,NULL, labels = "auto",label_size = 30, label_fontfamily = "ArialMT", ncol = 2)
#plot_grid(commonASEvents_suppA,commonASEvents_suppB, labels = c("a","b"), label_size = 30, label_fontfamily = "ArialMT", ncol = 1)

# Number of Splicing Events per AS Gene 
#plot_grid(numASEvents[[1]],NULL, ncol = 1, scale = 0.9)

#  NMD_vs_NonNMD Expression
# IR_NMD_run Venndiagram of NMD in IR transcripts
# Intron Retention Rate
# nmd_expression$Human,nmd_expression$Mouse,
plot_grid(human_mouse_IR_venn,adult_fetal_IR_venn,grobTree(IR_nmd_suppA),grobTree(IR_nmd_suppB),grobTree(IR_nmd_suppC),grobTree(IR_nmd_suppD),IR_rate_plots[[1]],labels = "auto",label_size = 30,label_fontfamily = "ArialMT", ncol = 2, scale = 0.9)


# Splice sites of all transcripts in annotated genes
#plot_grid(splice_site_frequency_output[[1]],splice_site_frequency_output[[2]],NULL, labels = c("a","b"),label_size = 30,label_fontfamily = "ArialMT",ncol = 1)
# Number of Isoforms in adult and fetal
#plot_grid(num_iso[[2]],NULL,ncol = 1,scale = 0.9)
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
