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
#SBATCH --array=0-7 # 8 samples
#SBATCH --output=All_Batch_Processing_NewChem-%A_%a.o
#SBATCH --error=All_Batch_Processing_NewChem-%A_%a.e

# 02/03/2021: ccs to cluster human samples (n = 14)

module load Miniconda2/4.3.21
#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/RNASeq

#mkdir $Isoseq3_WKD $PostIsoseq3_WKD $RNASeq_WKD
#cd $Isoseq3_WKD; mkdir CCS LIMA REFINE CLUSTER MERGED_CLUSTER
#cd $PostIsoseq3_WKD; mkdir MAP TOFU SQANTI2 KALLISTO TAMA SQANTI_TAMA_FILTER
#cd $RNASeq_WKD; mkdir MAPPED

# Input directory
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Isoseq_Paper/Fetal/raw
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/IsoSeq_Processing/Human
source $FUNCTIONS/All_Human_Functions.sh

# ENSURE ORDER OF SAMPLE NAMES AND BAM_FILES IS THE SAME
SAMPLES_NAMES=(AdultCTX1 AdultCTX2 AdultCTX3 AdultCTX4 AdultCTX5 FetalCTX1 FetalCTX3 FetalCTX5)
cd $FUNCTIONS
cat All_Human_RawData_NewChem.txt
# remove comments in raw.txt (https://kvz.io/blog/2007/07/11/cat-a-file-without-the-comments/)
BAM_FILES=(`cat "All_Human_RawData_NewChem.txt" | egrep -v "^\s*(#|$)"`)

# Batch run
SAMPLE=${SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}

################################################################################################
#*************************************  Isoseq3 [Function 1, 2, 3]
# Isoseq3.4.0
    # run_CCS <input_ccs_bam> <prefix_output_name> <Output_directory> <chem=oldchem/newchem>
    # run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
    # run_REFINE $Sample $Input_LIMA_directory $Output_directory
    # run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
run_CCS ${BAM_FILE} ${SAMPLE} $Isoseq3_WKD/CCS newchem
run_LIMA ${SAMPLE} $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA "no_multiplex"
run_REFINE ${SAMPLE} $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE
run_CLUSTER ${SAMPLE} $Isoseq3_WKD/REFINE $Isoseq3_WKD/CLUSTER
