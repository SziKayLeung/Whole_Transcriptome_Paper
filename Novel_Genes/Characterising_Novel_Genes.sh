#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrcq # submit to the serial queue
#PBS -l walltime=144:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC148213
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M sl693@exeter.ac.uk # email me at job completion



#************************************************************Define Global Variables

Human=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME
Mouse=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/Post_IsoSeq/SQANTI_TAMA_FILTER/GENOME
Mouse_Refine=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/All_Tg4510/IsoSeq/REFINE/fasta
Human_Refine=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Human/IsoSeq/REFINE/fasta
FUNCTIONS=//gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Whole_Transcriptome_Paper/Novel_Genes
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
HUMAN2MOUSE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Revised_Paper/Novel_Genes/Human2Mouse
BLAST_ACROSS_GENOME=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Whole_Transcriptome/Revised_Paper/Novel_Genes/AcrossGenome

#************************************************************Define functions
# Extract_PBIDs <input_classification_file>.txt <input_dir_path> <output_dir_path>
Extract_PBIDs(){
  cd $2
  
  # print PacBio Id and remove first row 
  awk '{ print $1 }' $1.txt >> $1_PBIds.txt 
}

# Grep_Fasta <prefix_PBIds.txt> <SQANTI2_filtered.fasta> <output_name.fasta> <output_dir_path>
Grep_Fasta(){
  # https://askubuntu.com/questions/976414/how-to-grep-sequence-of-fasta-using-list-of-ids-in-another-file
  #grep -x -F -A 1 -f $1_PBIds.txt $2 >> $4/$3
  #sed '/--/d' $4/$3 >> $4/Final_$3 #remove hyphens in-between lines
  #rm $4/$3
  source activate sqanti2_py3
  GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
  python $GENERALFUNC/TAMA/tama_sqanti_fastasubset.py $2 $1 $4/$3
}

Grep_Gtf(){  
  while read p; do
  grep "$p" $SQANTI2/WT8IsoSeq.collapsed.filtered.rep.renamed_corrected.gtf
  done <$WKD/Mouse_Novel_PBIds.txt >> $WKD/WT8IsoSeq_NOVEL.collapsed.filtered.rep.renamed_corrected.gtf
}

#Blast_seq2seq <working directory> <sample.fasta.for.db> <output_db_name> <sample.fasta.for.blast> <blast.output.name> <build> <threshold>
Blast_seq2seq(){
  
  ### Prepared Reference Genome for BLAST
  module load Miniconda2
  source activate nanopore
  cd $1
  if [[ $6 == "build" ]]; then 
  echo "Create blast database with fasta sequence: $2"
  makeblastdb -in $2 -dbtype nucl -out $3
  echo "Output file successfully generated: $2"
  
  ### Blast Probes to Reference Genome
  # https://angus.readthedocs.io/en/2016/running-command-line-blast.html
  # http://envgen.nox.ac.uk/bioinformatics/documentation/blast+/user_manual.pdf
  echo "Blast Probes to indexed database with fasta sequence: $3"
  if [[ $7 == "high" ]]; then
  blastn -query $4 -db $3 -out $5 -outfmt 6 -evalue 1e-5 
  #  -outfmt "7 qacc sacc evalue qstart qend sstart send" -evalue 1e-5  
  elif [[ $7 == "low" ]]; then
  echo "Using low threshold for evalue"
  blastn -query $4 -db $3 -out $5 -outfmt 6 -evalue 1e-2
  else
    echo "Threshold required as high or low"
  fi
  
  
  # column headers:https://molevol.mbl.edu/index.php/BLAST_UNIX_Tutorial
  #echo -e "Query_ID\Subject_ID\%_Identity\alignment_length\mismatches\gap\q.start\q.end\s.start\s.end\e_value\bit_score" | cat - $5 > Final_$5
  else 
     echo "Blast Probes to indexed database with fasta sequence: $3"
    if [[ $7 == "high" ]]; then
      blastn -query $4 -db $3 -out $5 -outfmt 6 -evalue 1e-5 
    #  -outfmt "7 qacc sacc evalue qstart qend sstart send" -evalue 1e-5  
    elif [[ $7 == "low" ]]; then
      echo "Using low threshold for evalue"
      blastn -query $4 -db $3 -out $5 -outfmt 6 -evalue 10
    else
      echo "Threshold required as high or low"
    fi
  fi  
  
}

#*********************************************************** Run functions
# Extract_PBIDs <input_classification_file>txt <input_dir_path> <output_dir_path>
# Grep_Fasta <prefix_PBIds.txt> <SQANTI2_filtered.fasta> <output_name.fasta> <output_dir_path>

cd $HUMAN2MOUSE
module load Miniconda2 
source activate nanopore 
Rscript $FUNCTIONS/Create_Novel_Genes_Classification.R
source deactivate 

Extract_PBIDs Combined_Human_Novel_classification . .
Grep_Fasta Combined_Human_Novel_classification_PBIds.txt $Human/HumanCTX_sqantifiltered_tamafiltered_classification.fasta Combined_Human_Novel.fasta .

Extract_PBIDs Mouse_Novel_classification . . 
Grep_Fasta Mouse_Novel_classification_PBIds.txt $Mouse/WholeIsoSeq_sqantifiltered_tamafiltered_classification.fasta Mouse_Novel.fasta .

# grep PB.15851.3 -A 1 $SQANTI2/WT8IsoSeq.collapsed.filtered.rep_classification.filtered_lite.fasta
# grep PB.15851.3 -A 1 $SQANTI2/WT8IsoSeq.collapsed.filtered.rep.renamed_corrected.fasta

#Blast_seq2seq <working directory> <sample.fasta.for.db> <output_db_name> <sample.fasta.for.blast> <blast.output.name>
Blast_seq2seq . Mouse_Novel.fasta Mouse_Novel.db Combined_Human_Novel.fasta Mouse2HumanSqanti.blast.txt build high
Blast_seq2seq . Mouse_Novel.fasta Mouse_Novel.db Combined_Human_Novel.fasta Mouse2HumanSqanti_Low.blast.txt build low
Blast_seq2seq . Combined_Human_Novel.fasta Combined_Human_Novel.db $Mouse_Refine/All_flnc.fasta Human2MouseFlnc.blast.txt build high
Blast_seq2seq . Mouse_Novel.fasta Mouse_Novel.db $Human_Refine/HumanCTX_flnc.fasta Mouse2HumanFlnc.blast.txt build high

Blast_seq2seq . $Mouse/WholeIsoSeq_sqantifiltered_tamafiltered_classification.fasta Mouse_All.db Combined_Human_Novel.fasta Human2MouseAllSqanti.blast.txt build high
Blast_seq2seq . $Human/HumanCTX_sqantifiltered_tamafiltered_classification.fasta Human_All.db Mouse_Novel.fasta Mouse2HumanAllSqanti.blast.txt build high

# For figures after further R analysis
#grep -A 50 PB.14504.1 Combined_Human_Novel.fasta
#grep -A 2 m54082_180607_173058/57147923/ccs $Mouse_Refine/All_flnc.fasta
#grep -A 50 PB.21588.1 Combined_Human_Novel.fasta 
#grep -A 2 m54082_190430_163756/72418134/ccs $Mouse_Refine/All_flnc.fasta

# number of flnc reads
# grep ccs $Human_Refine/HumanCTX_flnc.fasta | wc -l
grep ccs $Mouse_Refine/All_flnc.fasta | wc -l



#*********************************************************** BLAST across genome 
# novelGene_E2F3_AS (human) matched with mouse
# note using the corrected.fasta for mouse rather than SQANTI2 filtered fasta as cut fasta sequence (shouldn't make a difference whether corrected.fasta or filtered however) 
cd $BLAST_ACROSS_GENOME
grep PB.21939.1 -A 1 /gpfs/mrc0/projects/Research_Project-MRC148213/ISOSEQ/SQ2_v7_4/combinedFL.sample.collapsed.filtered.rep_classification.filtered_lite.renamed.fasta > E2F3_Human.fasta 
grep PB.3680.1 -A 1 $Mouse/WT8IsoSeq.collapsed.filtered.rep_corrected.fasta > E2F3_Mouse.fasta 


Blast_seq2seq . $REFERENCE/mm10.fa mm10.db E2F3_Mouse.fasta E2F3_Mouse_Genome_blast.txt build
Blast_seq2seq . $REFERENCE/hg38.fa hg38.db E2F3_Human.fasta E2F3_Human_Genome_blast.txt build

Blast_seq2seq . $REFERENCE/mm10.fa mm10.db $HUMAN2MOUSE/Mouse_Novel.fasta Mouse_Novel_Mouse_Genome_blast.txt build 
Blast_seq2seq . $REFERENCE/hg38.fa hg38.db $HUMAN2MOUSE/Combined_Human_Novel.fasta Human_Novel_Human_Genome_blast.txt build 

