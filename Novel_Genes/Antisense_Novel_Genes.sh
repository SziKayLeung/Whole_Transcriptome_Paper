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

# 14/08/2020: Identify novel antisense genes in human and mouse, and the subset of which overlap with exonic regions of other genes using bedtools
# Note Bedtools intersect requires input as bed file, thus subset gtf and then convert gtf to bed 

#************************************* DEFINE GLOBAL VARIABLES
Human=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME
Mouse=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME
ANTISENSE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Revised_Paper/Novel_Genes/Antisense
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019


#************************************* 1. Subset gtf to only antisense novel genes 
#************************************* 2. Extract gtf of other detected genes  
module load Miniconda2 
source activate nanopore 

cd $ANTISENSE
# Manually create Rscript external to this bash script to generate the gtf and classification files of antisense novel and non-antisense genes 
Rscript Create_Novel_Genes_Antisense_Classification.R

#************************************* 3. Convert gtf to bed
# gtf output from Rscript does not contain quotes around gene names etc
modify_gtf(){
  sed 's/transcript_id \([^;]\+\)/transcript_id \"\1\"/g' $1".gtf" | \
  sed 's/gene_id \([^;]\+\)/gene_id \"\1\"/g' | sed 's/gene_name \([^;]\+\)/gene_name \"\1\"/g' | \
  sed 's/ref_gene_id \([^;]\+\)/ref_gene_id \"\1\"/g' > $1"_final.gtf"
}

# conda install -c bioconda bedops
for f in combined.NovelGenes.Antisense WholeIsoSeq.NovelGenes.Antisense combined.No_Antisense WholeIsoSeq.No_Antisense; do 
  echo "Processing: $f"
  modify_gtf $f
  gtf2bed < $f"_final.gtf" >> $f"_final.bed"
done

#************************************* 4. Use bedtools to find overlapping regions   
## Bedtools parameters
#-split = report only exons overlap 
#-S = overlaps be found on opposite strands 
#-wo = amount of overlap 
#-wa = report the original "A" feature, -wb = report the original "B" feature 
#-loj = Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B.
#-header = Print the header from the A file prior to results.
#-C Reporting the number of overlapping features for each database file'

# Report depends on the bed file assigned first (-a), therefore switch between bed file containing novel, antsense genes and all other genes 
# human 
bedtools intersect -a combined.NovelGenes.Antisense_final.bed -b combined.No_Antisense_final.bed -split -S -wo > combined.NovelGenes.Antisense.output 
bedtools intersect -a combined.No_Antisense_final.bed -b combined.NovelGenes.Antisense_final.bed -split -S -wo > combined.NovelGenes.Antisense.output_collapsed

# mouse 
bedtools intersect -a WholeIsoSeq.NovelGenes.Antisense_final.bed -b WholeIsoSeq.No_Antisense_final.bed -split -S -wo > WholeIsoSeq.NovelGenes.Antisense.output 
bedtools intersect -a WholeIsoSeq.No_Antisense_final.bed -b WholeIsoSeq.NovelGenes.Antisense_final.bed -split -S -wo > WholeIsoSeq.NovelGenes.Antisense.output_collapsed

# Antisense, novel genes may have overlap to exonic regions to known genes that are not detected in Iso-Seq dataset
# therefore background with genome.gtf  
# bedtools gtf does not work directly with genome.gtf (https://www.biostars.org/p/56280/)
cat $REFERENCE/gencode.v31.annotation.gtf | grep transcript_id | grep gene_id | convert2bed --do-not-sort --input=gtf - > gencode.v31.annotation.bed
cat $REFERENCE/gencode.vM22.annotation.gtf | grep transcript_id | grep gene_id | convert2bed --do-not-sort --input=gtf - > gencode.M22.annotation.bed

bedtools intersect -a combined.NovelGenes.Antisense_final.bed -b gencode.v31.annotation.bed -split -S -wo > combined.NovelGenes.Antisense.genome.output 
bedtools intersect -a WholeIsoSeq.NovelGenes.Antisense_final.bed -b gencode.M22.annotation.bed -split -S -wo > WholeIsoSeq.NovelGenes.Antisense.genome.output 

mkdir gtf_bed; mv *gtf* gtf_bed; mv *bed* gtf_bed
