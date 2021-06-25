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
#SBATCH --output=All_Tg4510_Rarefaction.o
#SBATCH --error=All_Tg4510_Rarefaction.e

# 04/03/2020: Finalised alternative pipeline on whole transcriptome 12 samples

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
#cd $RNASeq_WKD; mkdir MAPPED 
#cd $RNASeq_WKD; mkdir STRINGTIE KALLISTO SQANTI2


# sourcing functions script and input directories
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/IsoSeq_Processing/Mouse
SAMPLES_LIST=$FUNCTIONS/All_Tg4510_Demultiplex.csv
source $FUNCTIONS/All_Tg4510_Functions.sh

REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
RNASeq_Filtered=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RNASeq/all_filtered
SAMPLES_NAMES=(Q21 O18 L22 K18 O23 S23 S18 K17 M21 K23 Q20 K24)

module load Miniconda2/4.3.21

################################################################################################
echo "#************************************* QC [Function 14,15,16]"
## 14) make_file_for_rarefaction <sample_name_prefix> <input_tofu_directory> <input_sqanti_tama_directory> <working_directory>
make_file_for_rarefaction WholeIsoSeq $PostIsoseq3_WKD/TOFU $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/GENOME $PostIsoseq3_WKD/RAREFACTION


