#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=90:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-2 # 3 samples
#SBATCH --reservation=research_project-mrc148213_5
#SBATCH --output=All_Tg4510_Rarefaction-%A_%a.o
#SBATCH --error=All_Tg4510_Rarefaction-%A_%a.o

# 04/03/2020: Finalised alternative pipeline on whole transcriptome 12 samples

#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/RNASeq
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Isoseq_Paper/Fetal/raw
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019

# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/IsoSeq_Processing/Human
source $FUNCTIONS/All_Human_Functions.sh

SAMPLES=(HumanCTX FetalCTX AdultCTX)
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}

module load Miniconda2/4.3.21

################################################################################################
echo "#************************************* QC [Function 14,15,16]"
## 14) make_file_for_rarefaction <sample_name_prefix> <input_tofu_directory> <input_sqanti_tama_directory> <working_directory>
make_file_for_rarefaction ${SAMPLE} $PostIsoseq3_WKD/TOFU $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/GENOME $PostIsoseq3_WKD/RAREFACTION


