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
#SBATCH --array=0-11 # 12 samples
#SBATCH --output=All_Batch_Processing-%A_%a.o
#SBATCH --error=All_Batch_Processing-%A_%a.e

# 01/03/2021: ccs to cluster and rnaseq star map individual Tg4510 samples (n = 12)

#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/RNASeq
ERCC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/ERCC

#mkdir $Isoseq3_WKD $PostIsoseq3_WKD $RNASeq_WKD $ERCC
#cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER MERGED_CLUSTER
#cd $PostIsoseq3_WKD; mkdir MAP TOFU SQANTI2 KALLISTO TAMA SQANTI_TAMA_FILTER
#cd $RNASeq_WKD; mkdir MAPPED

# Input directory
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/IsoSeq_Processing/Mouse
SAMPLE_LIST=$FUNCTIONS/All_Tg4510_Demultiplex.csv
source $FUNCTIONS/All_Tg4510_Functions.sh

# ENSURE ORDER OF SAMPLE NAMES AND BAM_FILES IS THE SAME
SAMPLES_NAMES=(Q21 O18 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24)
cd $FUNCTIONS
cat All_Tg4510_RawData.txt
# remove comments in raw.txt (https://kvz.io/blog/2007/07/11/cat-a-file-without-the-comments/)
BAM_FILES=(`cat "All_Tg4510_RawData.txt" | egrep -v "^\s*(#|$)"`)

# Batch run
SAMPLE=${SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}


################################################################################################
#*************************************  Isoseq3 [Function 1, 2, 3]
# Isoseq3.4.0
    # run_CCS_batch <input_ccs_bam> <prefix_output_name> <Output_directory>
    # run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
    # run_REFINE $Sample $Input_LIMA_directory $Output_directory
    # run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
run_CCS ${BAM_FILE} ${SAMPLE} $Isoseq3_WKD/CCS
run_LIMA ${SAMPLE} $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA "no_multiplex"
run_REFINE ${SAMPLE} $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE
run_CLUSTER ${SAMPLE} $Isoseq3_WKD/REFINE $Isoseq3_WKD/CLUSTER

################################################################################################
#************************************* RNAseq [Function 7]
## 4) run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
run_star ${SAMPLE} $RNASeq_Filtered $RNASeq_WKD/MAPPED
