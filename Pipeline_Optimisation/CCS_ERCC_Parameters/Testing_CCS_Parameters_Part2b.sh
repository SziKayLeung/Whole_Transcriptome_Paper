#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=4-5 #6 samples
#SBATCH --output=PostTesting_Part2_CCS_Parameters-%A_%a.o
#SBATCH --error=PostTesting_Part2_CCS_Parameters-%A_%a.e

# 23/09/2020: Post Testing_CCS_Parameters.sh with Sample M21 and K23, to run Iso-Seq3 pipeline and align with ERCC 
# 24/09/2020: Post Testing_CCS_Parameters.sh with Sample K18

#************************************* DEFINE GLOBAL VARIABLES
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/IsoSeq
CCS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/CCS
LIMA=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/LIMA
REFINE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/REFINE
CLUSTER=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/CLUSTER
MAPPING=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/MAPPING
TOFU=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/TOFU
SQANTI2_output_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing/SQANTI2

# Assign job to each test (Part2) of CCS Parameters
SAMPLES=(M21_1_passes_0.9 M21_3_passes_0.9 K23_1_passes_0.9 K23_3_passes_0.9 K18_1_passes_0.9 K18_3_passes_0.9)
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}

#************************************* Iso-Seq3 Pipeline 
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
run_LIMA ${SAMPLE} $CCS $LIMA "no_multiplex"
run_REFINE ${SAMPLE} $LIMA $REFINE 
run_CLUSTER ${SAMPLE} $REFINE $CLUSTER


#************************************* Extract Stats of CCS and Lima for downstream analysis 
module load Miniconda2 
source activate sqanti2_py3
python $FUNCTIONS/Run_Stats/CCS.py $CCS $CCS Testing_Part2
python $FUNCTIONS/Run_Stats/LIMA.py $LIMA $LIMA Testing_Part2


#************************************* Post Iso-Seq3 Pipeline
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh
convert_fa2fq ${SAMPLE} $CLUSTER
run_minimap2 ${SAMPLE} $CLUSTER ERCC 
tofu ${SAMPLE} $CLUSTER