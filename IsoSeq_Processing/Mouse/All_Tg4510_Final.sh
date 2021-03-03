#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=75:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=All_Tg4510_Final.o
#SBATCH --error=All_Tg4510_Final.e

# 02/03/2020: Finalised alternative pipeline on whole transcriptome 12 samples

#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq
ERCC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/ERCC

#cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER MERGED_CLUSTER
#cd $PostIsoseq3_WKD; mkdir MAP TOFU SQANTI2 KALLISTO TAMA SQANTI_TAMA_FILTER RAREFACTION
#cd $PostIsoseq3_WKD/SQANTI2; mkdir GENOME LNCRNA
#cd $PostIsoseq3_WKD/TAMA; mkdir GENOME LNCRNA
#cd $PostIsoseq3_WKD/SQANTI_TAMA_FILTER; mkdir GENOME LNCRNA
#cd $RNASeq_WKD; mkdir MAPPED STRINGTIE KALLISTO SQANTI2

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/IsoSeq_Processing/Mouse
SAMPLES_LIST=$FUNCTIONS/All_Tg4510_Demultiplex.csv
source $FUNCTIONS/All_Tg4510_Functions.sh

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
SAMPLES_NAMES=(Q21 O18 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24)

module load Miniconda2/4.3.21
################################################################################################
#************************************* Isoseq3 and Post_Isoseq3 [Function 1,2,3]
# 1) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
merging_at_refine $PostIsoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER WholeIsoSeq Q21 O18 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24

# 2) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
run_map_cupcakecollapse WholeIsoSeq $Isoseq3_WKD/MERGED_CLUSTER $PostIsoseq3_WKD/MAP $PostIsoseq3_WKD/TOFU

# 3) demux <input path read.stat file> <input path of samples file> <path of output>
demux $PostIsoseq3_WKD/TOFU/WholeIsoSeq.collapsed.read_stat.txt $SAMPLES_LIST $PostIsoseq3_WKD/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt

################################################################################################
#************************************* RNAseq [Function 4, 5]
### Processed as batches in All_Batch_Processing.sh
## 4) run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
#for i in "${SAMPLES_NAMES[@]}"; do run_star $i $RNASeq_Filtered $STAR_dir; done

## 5) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
mouse_merge_fastq $RNASeq_Filtered $PostIsoseq3_WKD/KALLISTO WholeIsoSeq

################################################################################################
#************************************* RNASeq & IsoSeq [Function 6]
## 6) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
run_kallisto WholeIsoSeq $PostIsoseq3_WKD/TOFU/WholeIsoSeq.collapsed.rep.fa $PostIsoseq3_WKD/KALLISTO $PostIsoseq3_WKD/KALLISTO

################################################################################################
#************************************* SQANTI2 [Function 7]
## 7) run_sqanti2 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
run_sqanti2 WholeIsoSeq.collapsed WholeIsoSeq.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/WholeIsoSeq.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2/GENOME genome

run_sqanti2 WholeIsoSeq.collapsed WholeIsoSeq.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/WholeIsoSeq.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/WholeIsoSeq.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2/LNCRNA lncrna

################################################################################################
#************************************* TAMA filter [Function 8,9]
# 8) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2/GENOME/WholeIsoSeq.collapsed_classification.filtered_lite.gtf WholeIsoSeq $PostIsoseq3_WKD/TAMA/GENOME
TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2/LNCRNA/WholeIsoSeq.collapsed_classification.filtered_lite.gtf WholeIsoSeq $PostIsoseq3_WKD/TAMA/LNCRNA

# 9) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
sqname=WholeIsoSeq.collapsed_classification.filtered_lite
TAMA_sqanti_filter $PostIsoseq3_WKD/TAMA/GENOME/WholeIsoSeq.bed $PostIsoseq3_WKD/SQANTI2/GENOME $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" WholeIsoSeq $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/GENOME

TAMA_sqanti_filter $PostIsoseq3_WKD/TAMA/LNCRNA/WholeIsoSeq.bed $PostIsoseq3_WKD/SQANTI2/LNCRNA $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" WholeIsoSeq $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/LNCRNA

###############################################################################################
# run_ERCC_analysis <sample_prefix_input/output_name> <isoseq3_input_directory> <output_dir>
run_ERCC_analysis WholeIsoSeq $Isoseq3_WKD $ERCC

# make_file_for_rarefaction <sample_name_prefix> <input_tofu_directory> <input_sqanti_tama_directory> <working_directory>
make_file_for_rarefaction WholeIsoSeq $PostIsoseq3_WKD/TOFU $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/GENOME $PostIsoseq3_WKD/RAREFACTION

# CCS and LIMA stats
# parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name>
parse_stats_per_sample $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA WholeIsoSeqAll

# ccs_bam2fasta <sample_name> <input_ccs.bam_dir> <output_dir>
mkdir $Isoseq3_WKD/CCS/Lengths
for i in "${SAMPLES_NAMES[@]}"; do ccs_bam2fasta $i $Isoseq3_WKD/CCS $Isoseq3_WKD/CCS/Lengths; done

################################################################
# run_stringtie <input_reference_gtf> <input_mapped_dir> <output_stringtie_dir> <output_prefix>
#run_stringtie $REFERENCE/gencode.vM22.annotation.gtf $STAR_dir $RNASeq_WKD/STRINGTIE rnaseq_stringtie_merged

# run_featurecounts_transcript_specified <input_dir> <input_reference_dir> <output_prefix_name> <output_dir>
#run_featurecounts_transcript_specified $STAR_dir $RNASeq_WKD/STRINGTIE/rnaseq_stringtie_merged.gtf WholeIsoSeq $RNASeq_WKD/FEATURECOUNTS

# run_sqanti2 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir>
#awk '{ if ($7 == 6) { print } }' pos_cut1-5.txt | head
#run_sqanti2 rnaseq_stringtie_merged rnaseq_stringtie_merged_mod.gtf $RNASeq_WKD/STRINGTIE $STAR_dir NA $RNASeq_WKD/FEATURECOUNTS/WholeIsoSeq.mod.transcript_id.tsv $RNASeq_WKD/SQANTI2