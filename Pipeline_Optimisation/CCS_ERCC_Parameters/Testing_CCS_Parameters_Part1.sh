#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

# 23/04/2020: Created Testing_CCS_Parameters.sh to test different parameters on Q21 sample
# 04/09/2020: Rerun script on K18 sample with ERCC and keep bam files
# 05/09/2020: Post Testing_CCS_Parameters.sh with Sample K18, to run Iso-Seq3 pipeline and align with ERCC 
  # 29/09/2020: Merge this script (Testing_CCS_Parameters.sh) with Post_Testing_CCS_Parameters.sh

#************************************* DEFINE GLOBAL VARIABLES
# File directories 
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/IsoSeq
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/Testing
CCS=$Isoseq3_WKD/CCS
LIMA=$Isoseq3_WKD/LIMA
REFINE=$Isoseq3_WKD/REFINE
CLUSTER=$Isoseq3_WKD/CLUSTER
MAPPING=$Isoseq3_WKD/MAPPING
TOFU=$Isoseq3_WKD/TOFU
SQANTI2_output_dir=$Isoseq3_WKD/SQANTI2

module load Miniconda2/4.3.21

if [ -d $Isoseq3_WKD/CCS ]; then 
	echo "CCS directory already present"
else 
	mkdir CCS	
fi

# ENSURE ORDER OF SAMPLE NAMES AND BAM_FILES IS THE SAME
SAMPLES_NAMES=(K18)
total_sample_num=0 	# number should be number of samples - 1, therefore 0 = running 1 sample

# remove comments in raw.txt (https://kvz.io/blog/2007/07/11/cat-a-file-without-the-comments/)
BAM_FILES=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/WholeTranscriptome/rawdata/m54082_190307_045507.subreads.bam


#************************************* DEFINE FUNCTION

# test_CCS <--min-passes=num> <--min-rq=num>
test_CCS(){
  source activate isoseq3
  ccs --version #ccs 4.0.0 (commit v4.0.0)

  cd $Isoseq3_WKD/CCS
  echo "Processing in $Isoseq3_WKD/CCS"
  sample_num=0
  until [ $sample_num -gt $total_sample_num ]; do 
  
  	output=${SAMPLES_NAMES[sample_num]}
  	input_subreads_bam=(${BAM_FILES[sample_num]})
  	
  	echo Processing Sample: $sample_num
  	echo Processing CCS for Sample $output with $input_subreads_bam
  	echo $1
  	echo $2
  	
  	num1=$(echo $1 |cut -d'=' -f 2)
  	num2=$(echo $2 |cut -d'=' -f 2)
	
  	# ccs <input.subreads.bam> <output.ccs.bam>
  	time ccs $input_subreads_bam $output"_ccs.bam" $1 $2 --chunk=1/10 --num-threads=20 --log-level=INFO --log-file=$output"_"$num1"_passes_"$num2"_ccs.log"
  	echo CCS for Sample $output successful
	
	# renaming files 
	mv ccs_report.txt $output"_"$num1"_passes_"$num2"_rq_ccs_report.txt"
    mv $output"_ccs.bam" $output"_"$num1"_passes_"$num2"_ccs.bam"
    mv $output"_ccs.bam.pbi" $output"_"$num1"_passes_"$num2"_ccs.bam.pbi"
  	#ls *$output"_"$num1"_passes_"$num2*
    ((sample_num++))
  done

  source deactivate
}


#************************************* CCS across different parameters

rq_parameters=(0.8 0.9 0.95 0.99)

# passes = 1 
for i in ${rq_parameters[@]}; do test_CCS --min-passes=1 --min-rq=$i; done 
# passes = 2 
for i in ${rq_parameters[@]}; do test_CCS --min-passes=2 --min-rq=$i; done 
# passes = 3 
for i in ${rq_parameters[@]}; do test_CCS --min-passes=3 --min-rq=$i; done 

# mv Part1 Testing to new folder 
cd $Isoseq3_WKD/CCS; mkdir Part1_K18 ; mv *K18* Part1_K18
cd $Isoseq3_WKD/CLUSTER; mkdir Part1_K18 ; mv *K18* Part1_K18
cd $Isoseq3_WKD/LIMA; mkdir Part1_K18 ; mv *K18* Part1_K18
cd $Isoseq3_WKD/REFINE; mkdir Part1_K18 ; mv *K18* Part1_K18
cd $Isoseq3_WKD/MAPPING; mkdir Part1_K18 ; mv *K18* Part1_K18
cd $Isoseq3_WKD/TOFU; mkdir Part1_K18 ; mv *K18* Part1_K18


#************************************* Run Iso-Seq3 Pipeline across different parameters  
source $FUNCTIONS/Isoseq3.2.2_Functions.sh
for i in $CCS/*ccs.bam; do 
	echo "Processing $i for Lima" 
	sample=$(basename "$i" | cut -d "_" -f 1,2,3,4 )
	run_LIMA $sample $CCS $LIMA "no_multiplex"
	run_REFINE $sample $LIMA $REFINE 
	run_CLUSTER $sample $REFINE $CLUSTER
done


#************************************* Extract Stats of CCS and Lima for downstream analysis 
source activate sqanti2_py3
python $FUNCTIONS/Run_Stats/CCS.py $CCS $CCS Testing
python $FUNCTIONS/Run_Stats/LIMA.py $LIMA $LIMA Testing


#************************************* Post Iso-Seq3 Pipeline
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/general/Post_IsoSeq/Post_Isoseq3_Functions.sh
for i in $CLUSTER/*hq.fasta; do 
	echo "Processing $i for Mapping" 
	sample=$(basename "$i" | cut -d "_" -f 1,2,3,4 | cut -d "." -f 1,2)
	echo $sample
	convert_fa2fq $sample $CLUSTER
	run_minimap2 $sample $CLUSTER ERCC 
	tofu $sample $CLUSTER
done 