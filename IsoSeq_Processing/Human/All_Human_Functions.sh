#!/bin/bash
################################################################################################
## Define functions for All_Human_Functions.sh
# 1) run_CCS <input_ccs_bam> <prefix_output_name> <Output_directory>
# 2) run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
# 3) run_REFINE $Sample $Input_LIMA_directory $Output_directory
# 4) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
# 5) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
# 6) demux <input path read.stat file> <input path of samples file> <path of output>
# 7) run_STAR_human <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
# 8) fetal_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
# 9) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
# 10) run_sqanti2 <input_tofu_prefix> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir>
# 11) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
# 12) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
################################################################################################

#************************************* DEFINE VARIABLES
module load Miniconda2/4.3.21
source activate isoseq3

# Listing versions
ccs --version #ccs 5.0.0 (commit v5.0.0)
lima --version #lima 2.0.0 (commit v2.0.0)
isoseq3 --version #isoseq3 3.4.0 (commit v3.4.0)

echo "FASTA SEQUENCE (CLONTECH PRIMERS) FOR NON-MULTIPLEXING"
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
FASTA=$REFERENCE/primer.fasta
cat $FASTA

head /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/README.md
################################################################################################
#*************************************  Isoseq3 [Function 1, 2, 3, 4]
# run_CCS <input_ccs_bam> <prefix_output_name> <Output_directory> <chem=oldchem/newchem>
run_CCS(){
  cd $3
  echo "Processing Sample $2 from Functions script"
  echo "Processing $1"
  
  if [ $4 == "oldchem" ]; then
    source activate isoseq3_paper
    ccs --version
    # ccs <input.subreads.bam> <output.ccs.bam>
    # --minPredictedAccuracy 0.9 = default
    time ccs $1 $2.ccs.bam --minPasses=1 --reportFile $2_ccs_report.txt
    echo "CCS for Sample $2 successful"
    ls *$2*    
  elif [ $4 == "newchem" ]; then
    source activate isoseq3
    ccs --version
    # ccs <input.subreads.bam> <output.ccs.bam>
    time ccs $1 $2.ccs.bam --minPasses=1 --min-rq 0.9 --reportFile $2_ccs_report.txt
    echo "CCS for Sample $2 successful"
    ls *$2*
  else
    echo "4th argument required"
  fi

  source deactivate
}

# run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
run_LIMA(){
  source activate isoseq3

  cd $3
  echo "Processing $1 file for demultiplexing"
    if [ -f $2/$1.fl.json ]; then
      echo "$1.fl.bam file already exists; LIMA no need to be processed on Sample $1"
    elif [ $4 = "multiplex" ]; then
      echo "Multiplex: use Targeted_FASTA"
      #lima <input.ccs.merged.consensusreadset.xml> <input.primerfasta> <output.fl.bam>
      time lima $2/$1.ccs.bam $TARGETED_FASTA $1.fl.bam --isoseq --dump-clips --dump-removed --peek-guess
      echo "lima $1 successful"
      ls $1.fl*
    else
      echo "No-Multiplex: use FASTA"
      time lima $2/$1.ccs.bam $FASTA $1.fl.bam --isoseq --dump-clips --dump-removed
      echo "lima $1 successful"
      ls $1.fl*
    fi
    source deactivate
}

# run_REFINE $Sample $Input_LIMA_directory $Output_directory
run_REFINE(){
  source activate isoseq3

  cd $3
  echo "Processing $1 file for refine"
  if [ -f $1.flnc.bam ]; then
    echo "$1.flnc bam file already exists; Refine no need to be processed on Sample $1"
  else
    #refine --require-polya <input.lima.consensusreadset.xml> <input.primer.fasta> <output.flnc.bam>
    time isoseq3 refine $2/$1.fl.primer_5p--primer_3p.bam $FASTA $1.flnc.bam --require-polya
    echo "refine $1 successful"
    ls $1.flnc*
  fi

  source deactivate
}

# run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
run_CLUSTER(){
  source activate isoseq3

  cd $3
  echo "Processing $1 file for cluster"
  if [ -f $1.clustered.bam ]; then
    echo "$1.clustered.bam file already exists; Cluster no need to be processed on Sample $1"
  else
    # cluster <input.flnc.bam> <output.unpolished.bam>
    time isoseq3 cluster $2/$1.flnc.bam $1.clustered.bam --verbose --use-qvs 2> $1.cluster.log
    echo "cluster $1 successful"
    ls $1.clustered*
  fi

  source deactivate
}


# merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
# aim: merging bam files from refine onwards (similar to run_isoseq3_2_1_merge, but no need to rerun from ccs)
# prerequisite: run ccs, lima, refine for samples
# <input_flnc_bam_dir> = path of flnc bam files generated from refine
# <output_dir> = path of output files from merging
# <output_name> = output name for merged files
# samples.... = list of the sample names
merging_at_refine(){

    source activate isoseq3
    isoseq3 --version

    ###********************* Merging at REFINE
    # Define variable "Merge_Samples" as a list of all samples, in order to find the specified flnc.bam (for dataset create ConsensusReadSet)
    # Define variable "all_flnc_bams" for merged path-directory of all flnc samples (for dataset create ConsensusReadSet)
    Merge_Samples=$(echo "${@:4}")

    echo "Merging flnc of samples $Merge_Samples"
    all_flnc_bams=$(
        for i in ${Merge_Samples[@]}; do
            flnc_bam_name=$(find $1 -name "*.flnc.bam" -exec basename \{} \; | grep ^$i )
            flnc_bam=$(find $1 -name "*.flnc.bam" | grep "$flnc_bam_name" )
            echo $flnc_bam
        done
    )

    cd $2
    printf '%s\n' "${all_flnc_bams[@]}" > $3.flnc.fofn
    cat $3.flnc.fofn
    
    echo "Merge Output: $3.clustered.bam"

    ###*********************

    isoseq3 cluster $3.flnc.fofn $3.clustered.bam --verbose --use-qvs 
    gunzip *.gz*
    source deactivate
}

################################################################################################
#************************************* Post_Isoseq3 (Minimap2, Cupcake, Demultiplex) [Function 5, 6]
# run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
# Prerequisite: mm10 cage peak
run_map_cupcakecollapse(){

    CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake
    SQANTI2_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    ANNOTATION=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/annotation
    SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence

    source activate sqanti2_py3

    # convert between fasta and fastq for downstream process
    echo "fasta2fastq conversion"
    python $SEQUENCE/fa2fq.py $2/$1.clustered.hq.fasta

    samtools --version # echo version
    echo "Minimap2 version:" $(minimap2 --version) # echo version

    echo "Processing Sample $1 for Minimap2 and sort"
    cd $3 #cd to $MAP directory for output
    minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $REFERENCE/hg38.fa $2/$1.clustered.hq.fastq > $1.sam 2> $1.map.log
    samtools sort -O SAM $1.sam > $1.sorted.sam

    echo "Processing Sample $1 for TOFU, with coverage 85% and identity 95%"
    source activate cupcake
    cd $4 #cd to $TOFU directory for output
    collapse_isoforms_by_sam.py -c 0.85 -i 0.95 --input $2/$1.clustered.hq.fastq --fq -s $3/$1.sorted.sam --dun-merge-5-shorter -o $1 &>> $1.collapse.log
    get_abundance_post_collapse.py $1.collapsed $2/$1.clustered.cluster_report.csv &>> $1.abundance.log

    source activate sqanti2_py3
    # convert rep.fq to rep.fa for SQANTI2 input
    seqtk seq -a $1.collapsed.rep.fq > $1.collapsed.rep.fa
    echo "Processing Sample $1 for TOFU successful"
    source deactivate
}


# demux <input path read.stat file> <input path of samples file> <path of output>
# run Cupcake_Demultiplex.R, read in read.stat file from cupcake collapse output and count abundance of each sample based on CCS_ID
demux(){
    source activate sqanti2_py3
    DEMUXFUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation

    # Rscript.R <input path read.stat file> <input path of samples file> <path of output>
    Rscript $DEMUXFUNCTIONS/Cupcake_Demultiplex.R  $1 $2 $3
}

################################################################################################
#************************************* RNAseq [Function 7,8]
# run_STAR_human <input_rnaseq_sample_name> <input_rnaseq_dir> <output_dir> <STAR_reference_input_directory>
# Aim: Individually align RNASeq filtered fastq files to STAR-indexed genome using STAR, and subsequently input SJ.out.bed file into SQANTI2
# Note: run_STAR_human and run_star (WT8_Isoseq3.1.2.sh) has the exact same command for STAR, but just different input of F_File and R_File
run_STAR_human(){

    source activate sqanti2
    cd $2
    echo "Finding reads for Sample $1 for STAR"
    # find reverse and forward file, trimming "./" therefore only printing file name
    F=$(find . | grep "fq" | grep $1 | grep "r1" | sed 's|^./||' )
    R=$(find . | grep "fq" | grep $1 | grep "r2" | sed 's|^./||' )
    # save path directory of files as variable for later mapping
    F_File=$(realpath $F)
    R_File=$(realpath $R)
    echo "Processing Forward Reads: $F_File"
    echo "Processing Reverse Reads: $R_File"

    cd $3
    mkdir $1 
    cd $1
    
    # Parameters according to https://github.com/jennylsmith/Isoseq3_workflow/blob/master/shell_scripts/8__STAR_Junctions_ShortReads.sh
    # 2-pass mode alignment with chimeric read detection
    # at least 25 bp of one read of a pair maps to a different loci than the rest of the read pair
    # require 20 pb on either side of chimeric junction
    # include chimeric reads in the output BAM
    # don't include chimeras that map to reference genome contig with Ns
    # --outSAMtype BAM SortedByCoordinate, output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command.
    # note STAR indexed with genocde vM22 primary assembly annotation gtf therefore no need to include into command line, otherwise create additional folders
    STAR --runMode alignReads --runThreadN 32 --genomeDir $4/hg38 \
    --readFilesIn $F_File $R_File \
    --outSAMtype BAM SortedByCoordinate \
    --chimSegmentMin 25 \
  	--chimJunctionOverhangMin 20 \
  	--chimOutType WithinBAM \
  	--chimFilter banGenomicN \
  	--chimOutJunctionFormat 1 \
  	--twopassMode Basic \
  	--twopass1readsN -1 &>> $1_star.log #use all reads
      #--readFilesCommand zcat \
      #--sjdbGTFfile $4/gencode.vM22.primary_assembly.annotation.gtf \
  
    # to remove duplicates between samples
    picard MarkDuplicates INPUT=$3/$1".sorted.bam" OUTPUT=$1".sorted.nodup.bam" METRICS_FILE=$1".dup.txt" VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> $1.PicardDup.log
  
    # rename output files
    mv Aligned.sortedByCoord.out.bam $1.sorted.bam
    mv Log.final.out $1.Log.final.out
    mv Log.out $1.Log.out
    mv Log.progress.out $1.Log.progress.out
    mv SJ.out.tab $1.SJ.out.bed
    mv $1.SJ.out.bed ../

  source deactivate
}


# fetal_merge_fastq <Fetal_RNASEQ_input_dir> <Fetal_Kallisto_output_dir>
# Aim: Merge forward and reverse of all RNASeq filtered fastq, for further downstream alignment to IsoSeq Tofu output using Kallisto
# PreRequisite: Ensure specified raw fastq files defined from function are in RNASEq_input_dir
# Output: Fetal_R1.fq, Fetal_R2.fq
Fetal_merge_fastq(){
    cd $1
    echo "Processing Fetal R1 files"
    ls *r1.fq
    cat 3005_11921_trimmed_r1.fq  3005_12971FL_trimmed_r1.fq  3005_12972FL_trimmed_r1.fq  > $2/Fetal_R1.fq

    echo "Processing Fetal R2 files"
    ls *r2.fq
    cat 3005_11921_trimmed_r2.fq  3005_12971FL_trimmed_r2.fq  3005_12972FL_trimmed_r2.fq  > $2/Fetal_R2.fq
}



################################################################################################
#************************************* RNASeq & IsoSeq [Function 9]
# run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
# Aim: Align merged RNASeq fastq files to IsoSeq Tofu output fasta using Kallisto
# Prerequisite:
    # run mouse_merge_fastq/any equivalent merging of all RNASeq fastq files (note sample_prefix_output name is the same)
# Output: <output_prefix>.mod.abundance.tsv for input into Sqanti_qc.py (not modified script to take in space delimited rather than tab file)
run_kallisto(){
    source activate sqanti2

    echo "Processing Kallisto for $1"
    cd $4
    kallisto version
    time kallisto index -i $1_Kallisto.idx $2 2> $1_Kallisto.index.log
    time kallisto quant -i $4/$1_Kallisto.idx --fr-stranded $3/$1_R1.fq --rf-stranded $3/$1_R2.fq -o $4 2> $1_Kallisto.quant.log
    mv abundance.tsv $1.abundance.tsv

    # problem: retained co-ordinates, which does not input well into SQANTI2
    echo "Kallisto original $1.abundance.tsv"
    head $1.abundance.tsv
    # solution: retain the PB.ID
    while read line ; do
	    first=$( echo "$line" |cut -d\| -f1 ) # each read line, and remove the first part i.e. PacBio ID
	    rest=$( echo "$line" | cut -d$'\t' -f2-5 ) #save the remaining columns
	    echo $first $rest # concatenate
    done < $4/$1.abundance.tsv > $4/$1.temp.abundance.tsv

    header=$(head -n 1 $4/$1.abundance.tsv)
	sed -i '1d' $4/$1.temp.abundance.tsv # remove header of temp.file to be replaced
    echo $header > foo
	cat foo $4/$1.temp.abundance.tsv > $4/$1.mod.abundance.tsv

    echo "Kallisto $1.mod.abundance.tsv"
    head $4/$1.mod.abundance.tsv
    rm $1.temp.abundance.tsv
	rm foo

    source deactivate
}


################################################################################################
#************************************* SQANTI2 [Function 10]
# run_sqanti2 <input_tofu_prefix> <inout_gtf> <input_tofu_dir> <coverage/genome=own_rnaseq/introplis/chess> <input_RNASEQ_dir> <input_kallisto_file> <input_abundance> <output_dir>
# PreRequisite:
    # if running "own_rnaseq", define $MAPPED directory containing SJ.out.bed files in working script
    # Define $REFERENCE directory containing gencode.v31.annotation.gtf
# Note: Function specific to human by reference genome input, and also have to state whether using "own_rnaseq" or "introplis" data as input; "chess" for usig intropolis dataset but chess.gtf
run_sqanti2(){
    #cp /gpfs/mrc0/projects/Research_Project-MRC148213/ISOSEQ/reference/gencode.v31.long_noncoding_RNAs.gtf /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    # variables
    input_tofu_prefix=$1
    input_gtf=$2
    input_dir=$3
    coverage_decision=$4
    input_RNASEQ_dir=$5
    input_kallisto_file=$6
    input_abundance=$7
    output_dir=$8

    source activate sqanti2_py3
    echo "Processing Sample $input_tofu_prefix for SQANTI2 QC"
    cd $output_dir

    SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
    introplis=$(echo /gpfs/mrc0/projects/Research_Project-MRC148213/ISOSEQ/SQANTI2addons/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified)

    export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
    python $SQANTI2_DIR/sqanti_qc2.py --version

    if [ $coverage_decision == "own_rnaseq" ]; then
        echo "Processing with own RNASeq dataset and genocode.gtf for genome annotation "
        echo "Collapsing with the following bed files"
        echo $input_RNASEQ_dir/*SJ.out.bed
        
        cp $input_RNASEQ_dir/*SJ.out.bed .

        python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf $input_dir/$input_gtf $REFERENCE/gencode.v31.annotation.gtf $REFERENCE/hg38.fa --cage_peak $REFERENCE/hg38.cage_peak_phase1and2combined_coord.bed --coverage "./*SJ.out.bed" --expression $input_kallisto_file --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --fl_count $input_abundance &>> $input_tofu_prefix.sqanti.qc.log
        
        rm *SJ.out.bed

    elif [ $coverage_decision == "intropolis" ]; then
        echo "Processing with own intropolis and genocode.gtf for genome annotation, no expression data"
        echo "Collapsing with the following coverage expression file"
        echo $introplis

        python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf $input_dir/$input_gtf $REFERENCE/gencode.v31.annotation.gtf $REFERENCE/hg38.fa --cage_peak $REFERENCE/hg38.cage_peak_phase1and2combined_coord.bed --coverage $introplis --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --fl_count $input_abundance &>> $input_tofu_prefix.sqanti.qc.log

    elif [ $coverage_decision == "chess" ]; then
        echo "Processing with intropolis and chess.gtf for genome annotation"
        echo "Collapsing with the following coverage expression file"
        echo $introplis

        python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf $input_dir/$input_gtf $REFERENCE/chess2.2.gtf $REFERENCE/hg38.fa --cage_peak $REFERENCE/hg38.cage_peak_phase1and2combined_coord.bed --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --coverage $introplis --fl_count $input_abundance &>> $input_tofu_prefix.sqanti.qc.log

    elif [ $coverage_decision == "lncrna" ]; then
        echo "Processing with intropolis and gencode.v31.long_noncoding_RNAs.gtf for genome annotation"
        echo "Collapsing with the following coverage expression file"
        echo $introplis

        python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf $input_dir/$input_gtf $REFERENCE/gencode.v31.long_noncoding_RNAs.gtf $REFERENCE/hg38.fa --cage_peak $REFERENCE/hg38.cage_peak_phase1and2combined_coord.bed --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --coverage $introplis --fl_count $input_abundance &>> $input_tofu_prefix.sqanti.qc.log

    else
        echo "ERROR: Require 4th argument to be either: own_rnaseq OR intropolis OR chess"
        return
    fi
    
    mkdir $1"_QC_PNG" $1"_FILTER_PNG"
    mv *png* $1"_QC_PNG/"    

    # sqanti filter
    echo "Processing Sample $input_tofu_prefix for SQANTI2 filter"
    python $SQANTI2_DIR/sqanti_filter2.py $input_tofu_prefix"_classification.txt" $input_tofu_prefix"_corrected.fasta" $input_tofu_prefix"_corrected.gtf" -a 0.6 -c 3 &>> $input_tofu_prefix.sqanti.filter.log

    mv *png* $1"_FILTER_PNG/"
}

################################################################################################
#************************************* TAMA [Function 11, 12]
# TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
# remove short fragments from post tofu
# Prerequisite: Require TAMA_prepare.R to change column 4 of bed12 to gene_name: transcript_name for correct file TAMA format
TAMA_remove_fragments(){

    TAMAFUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/TAMA
    TAMA_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/tama/tama_go/filter_transcript_models
    source activate sqanti2
    cd $3

    # convert gtf to bed12
	  gtfToGenePred $1 $2.genepred
	  genePredToBed $2.genepred > $2.bed12
    awk -F'\t' '{print $1,$2,$3,$4,"40",$6,$7,$8,"255,0,0",$10,$11,$12}' $2.bed12| sed s/" "/"\t"/g|sed s/",\t"/"\t"/g|sed s/",$"/""/g > Tama_$2.bed12

    # Rscript script.R <name>.bed12 <input_dir>
	  Rscript $TAMAFUNCTIONS/TAMA_Merge_Prepare.R Tama_$2 $3
	  python $TAMA_DIR/tama_remove_fragment_models.py -f Tama_$2_mod.bed12 -o $2

	  rm *Tama*

    echo "Number of isoforms filtered by TAMA:"
    wc -l $2"_discarded.txt"

	  source deactivate
}

# TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <sqanti_output_junc> <output_prefix_name> <output_dir>
TAMA_sqanti_filter(){
  source activate sqanti2_py3
  GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
# Rscript .R <path/tama_filtered_output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_junc_txt> <output_prefix_name> <output_dir>
  Rscript $GENERALFUNC/TAMA/tama_sqanti_classgtfsubset.R $1 $2 $3 $4 $6 $7 $8
  # extract fasta sequence based on the pbid (https://www.biostars.org/p/319099/)
  # script.py <path/sqanti_filtered.fasta> <path/retained_pbid_tama.txt> <path/output.fasta>
  cd $8
  awk '{ print $4 }' $1| cut -d ";" -f 2  > $7_tama_retained_pbid.txt
  python $GENERALFUNC/TAMA/tama_sqanti_fastasubset.py $2/$5 $8/$7_tama_retained_pbid.txt $8/$7"_sqantifiltered_tamafiltered_classification.fasta"
  
  # generate file 
  SQ_Report=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2/utilities/SQANTI_report2.R
  Rscript $SQ_Report $8/$7"_sqantitamafiltered.classification.txt" $8/$7"_sqantitamafiltered.junction.txt"
  
  mkdir $7"_PNG"
  mv *png* $7"_PNG"
}

################################################################################################
#************************************* QC [Function 13]
# ccs_bam2fasta: convert ccs_bam to fasta
# ccs_bam2fasta <sample_name> <input_ccs.bam_dir> <output_dir>
# output file = <sample_name>.ccs.fasta
ccs_bam2fasta(){

    source activate sqanti2_py3
    bam2fastq --version

    echo "Converting CCS of sample $1 from bam to fasta"
    ccs_file=$(find $2 -name "*.bam" | grep $1)
    echo $ccs_file
    cd $3
    bam2fasta -u -o $1".ccs" $ccs_file &>> $1.bam2fasta.log

    echo "Extracting length of Sample $1 CCS"
    CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake
    export PYTHONPATH=$PYTHONPATH:$CUPCAKE/sequence
    export PATH=$PATH:$CUPCAKE/sequence
    get_seq_stats.py $1".ccs.fasta" &>> $1.ccs.get_seq.log

    source deactivate
}

# make_file_for_rarefaction <sample_name_prefix> <input_tofu_directory> <input_sqanti_tama_directory> <working_directory>
make_file_for_rarefaction(){
    source activate cupcake
    CUPCAKE_ANNOTATION=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/annotation

  	cd $4
  	echo "Working with $1"
    prefix=$1.collapsed
  	# make_file_for_subsampling_from_collapsed.py <sample_name_prefix>.input.file <sample_name_prefix>.output.file <sample_name_prefix>.classification.txt
  	python $CUPCAKE_ANNOTATION/make_file_for_subsampling_from_collapsed.py -i $2/$prefix -o $1.subsampling -m2 $3/$1_sqantitamafiltered.classification.txt &>> $1.makefile.log
    python $CUPCAKE_ANNOTATION/subsample.py --by refgene --min_fl_count 2 --step 1000 $1.subsampling.all.txt > $1.rarefaction.by_refgene.min_fl_2.txt 
    python $CUPCAKE_ANNOTATION/subsample.py --by refisoform --min_fl_count 2 --step 1000 $1.subsampling.all.txt > $1.rarefaction.by_refisoform.min_fl_2.txt
  	python $CUPCAKE_ANNOTATION/subsample_with_category.py --by refisoform --min_fl_count 2 --step 1000 $1.subsampling.all.txt > $1.rarefaction.by_refisoform.min_fl_2.by_category.txt

    source deactivate
}

# parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name> <oldchem/newchem>
parse_stats_per_sample(){

    # variable
    CCS_dir=$1
    LIMA_dir=$2
    output_name=$3
    chem=$4
    samples=$(echo "${@:5}")
   
   # ccs report of samples
    all_report=$(
        for i in ${samples[@]}; do
            name=$(find $CCS_dir -name "*ccs_report.txt" -exec basename \{} \; | grep ^$i )
            report=$(find $CCS_dir -name "*ccs_report.txt" | grep "$name" )
            echo $report
        done
    )
    
    # use different ccs python script for different ccs versions (due to different seqeuncing chemistries)

    FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
    
    if [ $chem == "oldchem" ]; then
      source activate sqanti2
      echo "Stats of: $samples"
      mkdir $CCS_dir/oldchem_ccs 
      for f in ${all_report[@]}; do cp $f $CCS_dir/oldchem_ccs ; done 
      python $FUNCTIONS/IsoSeq_QC/CCS_OldChem.py $CCS_dir/oldchem_ccs $CCS_dir/oldchem_ccs/ $output_name"_oldchem"
    elif [ $chem == "newchem" ]; then 
      source activate sqanti2_py3
      echo "Stats of: $samples"
      mkdir $CCS_dir/newchem_ccs 
      for f in ${all_report[@]}; do cp $f $CCS_dir/newchem_ccs ; done 
      python $FUNCTIONS/IsoSeq_QC/CCS.py $CCS_dir/newchem_ccs "" $output_name"_newchem"
    else 
      echo "4th argument required as newchem/oldchem"
    fi 
    
    python $FUNCTIONS/IsoSeq_QC/LIMA.py $LIMA_dir "" $output_name
    source deactivate
}


## Merged 
# parse_ccs_merged_samples <sample_ccs_dir/*ccs.fasta> <sample_output_prefix> <output_dir>
# Aim: To merge all the ccs.fasta files to see the merged distribution per sample 
# Prerequisite: run parse_ccs_per_sample as require ccs.fasta as input 
parse_ccs_merged_samples(){

    source activate sqanti2_py3    
    cd $3

    echo "Merging the following fasta files to $2_merged.fasta"
    ls $1/*fasta 

    cat $1/*fasta > $2_merged.fasta
    
    CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake
    export PYTHONPATH=$PYTHONPATH:$CUPCAKE/sequence
    export PATH=$PATH:$CUPCAKE/sequence
    get_seq_stats.py $2_merged.fasta 2> $2_merged.ccs.get_seq.log
    mv $2_merged.fasta.seqlengths.txt $2_merged.fastq.seqlengths.txt  # for consistent file naming as fastq for downstream processing in rscript for plots

    source deactivate
}

################################################################################################
#************************************* Alternative Splicing [Function 14]

# run_suppa2 <input_gtf> <input_class> <output_dir> <output_name>
# input_gtf: sqanti tama filtered gtf
# input_class: sqanti tama filtered classification file
run_suppa2(){
  # variable
  input_gtf=$1
  input_class=$2
  output_dir=$3
  output_name=$4

  echo "#******************* $output_name"
  echo "Processing gtf: $input_gtf"
  echo "Processing classification file: $input_class"

  FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation

  source activate sqanti2_py3
  cd $output_dir
  suppa.py generateEvents -i $input_gtf -o $output_name --pool-genes -f ioe -e SE MX FL SS RI &>> $output_name"_suppa2.log"
  python $FUNCTIONS/Suppa2_output_mod_updated.py $output_name $output_dir $input_class &>> $output_name"_suppa2_mod.log"
}
