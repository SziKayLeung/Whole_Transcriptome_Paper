#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/Whole_Transcriptome_Paper/IsoSeq_Processing/ERCC
source $SC_ROOT/ERCC.config


# minimap on the refine files to reference genome (ERCC92.fa) and collapse files
for basename in ${SAMPLES[@]}; do 
  echo "Processing ${basename}"
  source activate nanopore
  
  cd ${WKD_ROOT}/1_minimap2
  minimap2 -ax splice -t 30 -uf --secondary=no -C5 ${GENOME_FASTA} ${REFINE_FASTA}/${basename}.fa > ${basename}.sam 
  samtools view -S -b ${basename}.sam > ${basename}.bam
  samtools sort ${basename}.bam -o ${basename}.sorted.bam
  
  # collapse 
  source activate sqanti2
  cd ${WKD_ROOT}/2_tama_collapse
  python ${TAMA_SC} -s ${WKD_ROOT}/1_minimap2/${basename}.sorted.bam -f ${GENOME_FASTA} -p ${basename} -x no_cap -b BAM -rm low_mem &> ${basename}.collapse.log
done


# merge data
cd ${WKD_ROOT}/3_tama_merge
python ${SOFTWAREPATH}/tama/tama_merge.py -f tama_merge.txt -p $NAME -d merge_dup -a 100 -z 100
python ${SOFTWAREPATH}/tama/tama_go/read_support/tama_read_support_levels.py -f tama_read_support.txt -o ${NAME}_counts -m ${NAME}_merge.txt


# tama filter 
python ${SOFTWAREPATH}/tama/tama_go/filter_transcript_models/tama_remove_single_read_models_levels.py -b ${NAME}.bed -r ${NAME}_counts_read_support.txt -o ${NAME}_nRead2
python ${SOFTWAREPATH}/tama/tama_go/filter_transcript_models/tama_remove_fragment_models.py -f ${NAME}_nRead2.bed -o ${NAME}_filtFrag


# create gtf file
python ${SOFTWAREPATH}/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py ${NAME}_filtFrag.bed ${NAME}_filtFrag.gtf
python ${SOFTWAREPATH}/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py ${NAME}.bed ${NAME}.gtf


# create new abundance file 
awk '{if ($6 != "removed_transcript") print $5,$6,$7,$3,$4}' ${NAME}_nRead2_singleton_report.txt > ${NAME}_nRead2_counts.txt


# run sqanti3
source activate sqanti2_py3
cd ${WKD_ROOT}/4_sqanti/tama_filtered
python $SQANTI2_dir/sqanti_qc2.py -t 30 --gtf ${WKD_ROOT}/3_tama_merge/${NAME}_filtFrag.gtf $REFERENCE_ERCC/ERCC92.gtf $REFERENCE_ERCC/ERCC92.fa &>> $NAME.sqanti.qc.log

cd ${WKD_ROOT}/4_sqanti/tama_collapsed
python $SQANTI2_dir/sqanti_qc2.py -t 30 --gtf ${WKD_ROOT}/3_tama_merge/${NAME}.gtf $REFERENCE_ERCC/ERCC92.gtf $REFERENCE_ERCC/ERCC92.fa &>> $NAME.sqanti.qc.log
