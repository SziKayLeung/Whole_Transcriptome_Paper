#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=All_Demultiplex_Final.o
#SBATCH --error=All_Demultiplex_Final.e

# 23/02/2020: Finalised alternative pipeline on whole transcriptome 12 samples
  # Prerequisite: Run ccs,lima,refine on individual samples

#************************************* DEFINE GLOBAL VARIABLES
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/Post_IsoSeq/
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/RNASeq
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Isoseq_Paper/Fetal/raw
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019

#mkdir $Isoseq3_WKD $PostIsoseq3_WKD $RNASeq_WKD
#cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER MERGED_CLUSTER
#cd $PostIsoseq3_WKD; mkdir MAP TOFU SQANTI2 KALLISTO TAMA SQANTI_TAMA_FILTER
#cd $PostIsoseq3_WKD/SQANTI2; mkdir GENOME CHESS LNCRNA
#cd $PostIsoseq3_WKD/TAMA; mkdir GENOME CHESS LNCRNA
#cd $PostIsoseq3_WKD/SQANTI_TAMA_FILTER; mkdir GENOME CHESS LNCRNA
#cd $RNASeq_WKD; mkdir MAPPED

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/IsoSeq_Processing/Human
source $FUNCTIONS/All_Human_Functions.sh

HUMANCTX_SAMPLES_NAMES=(AdultCTX1 AdultCTX2 AdultCTX3 AdultCTX4 AdultCTX5 FetalCTX1 FetalCTX2 FetalCTX3 FetalCTX4 FetalCTX5)
ADULTCTX_SAMPLES_NAMES=(AdultCTX1 AdultCTX2 AdultCTX3 AdultCTX4 AdultCTX5)
FETALCTX_SAMPLES_NAMES=(FetalCTX1 FetalCTX2 FetalCTX3 FetalCTX4 FetalCTX5)
FETALHIP_SAMPLES_NAMES=(FetalHIP1 FetalHIP2)
FETALSTR_SAMPLES_NAMES=(FetalSTR1 FetalSTR2)

GROUP_NAMES=(HumanCTX AdultCTX FetalCTX FetalHIP FetalSTR)
GROUP_SAMPLES=($HUMANCTX_SAMPLES_NAMES $ADULTCTX_SAMPLES_NAMES $FETALCTX_SAMPLES_NAMES $FETALHIP_SAMPLES_NAMES $FETALSTR_SAMPLES_NAMES)

################################################################################################
#************************************* Isoseq3 and Post_Isoseq3 [Function 1,2,3]
# 1) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
for i in {0..4}; do merging_at_refine $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER ${GROUP_NAMES[$i]} ${GROUP_SAMPLES[$i]}; done

# 2) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
for i in {0..4}; do run_map_cupcakecollapse ${GROUP_SAMPLES[$i]} $Isoseq3_WKD/MERGED_CLUSTER $PostIsoseq3_WKD/MAP $PostIsoseq3_WKD/TOFU; done

# 3) demux <input path read.stat file> <input path of samples file> <path of output>
for i in {0..4}; do demux $PostIsoseq3_WKD/TOFU/${GROUP_SAMPLES[$i]}.collapsed.read_stat.txt $SAMPLES_LIST/All_Human_Demultiplex.csv $PostIsoseq3_WKD/TOFU/${GROUP_SAMPLES[$i]}.Demultiplexed_Abundance.txt; done

################################################################################################
#************************************* RNAseq [Function 4, 5]
sample=(3005_11921 3005_12971FL 3005_12972FL)
for i in ${sample[@]}; do run_STAR_human $i $RNASeq_Filtered $RNASeq_WKD/MAPPED $REFERENCE; post_STAR_process $RNASeq_WKD/MAPPED; done

# Fetal_merge_fastq <Fetal_RNASEQ_input_dir> <Fetal_Kallisto_output_dir>
Fetal_merge_fastq $RNASeq_Filtered $PostIsoseq3_WKD/KALLISTO

################################################################################################
#************************************* RNASeq & IsoSeq [Function 6]
## 6) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
run_kallisto FetalCTX $PostIsoseq3_WKD/TOFU/FetalCTX.collapsed.rep.fa $PostIsoseq3_WKD/KALLISTO $PostIsoseq3_WKD/KALLISTO

################################################################################################
#************************************* SQANTI2 [Function 7]
# run_sqanti2 <input_tofu_prefix> <inout_gtf> <input_tofu_dir> <coverage/genome=own_rnaseq/introplis/chess/lncrna> <input_RNASEQ_dir> <input_kallisto_file> <input_abundance> <output_dir>
for i in HumanCTX AdultCTX FetalSTR FetalHIP; do run_sqanti2 $i.collapsed $i.collapsed.gff $PostIsoseq3_WKD/TOFU introplis NA NA $PostIsoseq3_WKD/TOFU/$i".Demultiplexed_Abundance.txt" $PostIsoseq3_WKD/SQANTI2/GENOME; done

# fetal with own rnaseq
run_sqanti2 FetalCTX.collapsed FetalCTX.collapsed.gff $PostIsoseq3_WKD/TOFU own_rnaseq $PostIsoseq3_WKD/KALLISTO $PostIsoseq3_WKD/KALLISTO/FetalCTX.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/HumanCTX.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2

# chess and lncrna on human dataset
run_sqanti2 HumanCTX.collapsed HumanCTX.collapsed.gff $PostIsoseq3_WKD/TOFU chess NA NA $PostIsoseq3_WKD/TOFU/HumanCTX.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2/CHESS
run_sqanti2 HumanCTX.collapsed HumanCTX.collapsed.gff $PostIsoseq3_WKD/TOFU lncrna NA NA $PostIsoseq3_WKD/TOFU/HumanCTX.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2/LNCRNA

################################################################################################
#************************************* TAMA filter [Function 8,9]
# 8) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
for i in HumanCTX FetalCTX FetalCTX FetalSTR FetalHIP; do TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2/GENOME/$i".collapsed_classification.filtered_lite.gtf" $i $PostIsoseq3_WKD/TAMA/GENOME; done

# chess and lncrna on human dataset
for dir in CHESS LNCRNA; do TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2/CHESS/HumanCTX.collapsed_classification.filtered_lite.gtf HumanCTX $PostIsoseq3_WKD"/TAMA/"$dir; done

# 9) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
sqname=.collapsed_classification.filtered_lite
for i in HumanCTX FetalCTX FetalCTX FetalSTR FetalHIP; do TAMA_sqanti_filter $PostIsoseq3_WKD/TAMA/GENOME/$i".bed" $PostIsoseq3_WKD/SQANTI2/GENOME/ $i$sqname"_classification.txt" $i$sqname".gtf" $i$sqname".fasta" $i$sqname"_junctions.txt" $i $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/GENOME; done

for dir in CHESS LNCRNA; do TAMA_sqanti_filter $PostIsoseq3_WKD"/TAMA/"$dir"/HumanCTX.bed" $PostIsoseq3_WKD"/SQANTI2/"$dir HumanCTX$sqname"_classification.txt" HumanCTX$sqname".gtf" HumanCTX$sqname".fasta" HumanCTX$sqname"_junctions.txt" HumanCTX $PostIsoseq3_WKD"/SQANTI_TAMA_FILTER/"$dir; done
