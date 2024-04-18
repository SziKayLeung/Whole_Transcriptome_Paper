#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-2 # 3 samples
#SBATCH --output=Testing_Part2_CCS_Parameters-%A_%a.o
#SBATCH --error=Testing_Part2_CCS_Parameters-%A_%a.e


# 22/09/2020: Having further established two useful sets of parameters, extend analyses to two additional samples and all of CCS 

#************************************* DEFINE GLOBAL VARIABLES
# File directories 
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/IsoSeq
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing

module load Miniconda2/4.3.21

if [ -d $Isoseq3_WKD/CCS ]; then 
	echo "CCS directory already present"
else 
	mkdir CCS	
fi

# ENSURE ORDER OF SAMPLE NAMES AND BAM_FILES IS THE SAME
SAMPLES_NAMES=(M21 K23 K18)
SAMPLE=${SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}

RAW_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/rawdata
BAM_FILES=($RAW_DIR/m54082_190430_163756.subreads.bam $RAW_DIR/m54082_190524_145911.subreads.bam $RAW_DIR/m54082_190307_045507.subreads.bam)
BAM_FILE=${BAM_FILES[${SLURM_ARRAY_TASK_ID}]}

#************************************* DEFINE FUNCTION

# test_CCS <--min-passes=num> <--min-rq=num>
test_CCS(){
    source activate isoseq3
    ccs --version #ccs 4.0.0 (commit v4.0.0)

    cd $Isoseq3_WKD/CCS
    echo "Processing in $Isoseq3_WKD/CCS"
    echo Processing CCS for Sample ${SAMPLE} with ${BAM_FILE}
    echo $1
    echo $2

    num1=$(echo $1 |cut -d'=' -f 2)
    num2=$(echo $2 |cut -d'=' -f 2)

    # ccs <input.subreads.bam> <output.ccs.bam>
    time ccs ${BAM_FILE} ${SAMPLE}"_ccs.bam" $1 $2 --num-threads=20 --log-level=INFO --log-file=${SAMPLE}"_"$num1"_passes_"$num2"_ccs.log"
    echo CCS for Sample ${SAMPLE} successful

    # renaming files 
    mv ccs_report.txt ${SAMPLE}"_"$num1"_passes_"$num2"_rq_ccs_report.txt"
    mv ${SAMPLE}"_ccs.bam" ${SAMPLE}"_"$num1"_passes_"$num2"_ccs.bam"
    mv ${SAMPLE}"_ccs.bam.pbi" ${SAMPLE}"_"$num1"_passes_"$num2"_ccs.bam.pbi"
    
    source deactivate
}


#************************************* APPLY FUNCTION

test_CCS --min-passes=1 --min-rq=0.9
test_CCS --min-passes=3 --min-rq=0.9