#!/bin/bash
################################################################################################
## Define functions for All_Tg4510_Functions.sh
# 1) run_CCS <input_ccs_bam> <prefix_output_name> <Output_directory>
# 2) run_LIMA $Sample $Input_CCS_directory $Output_directory <"no_multiplex"/"multiplex">
# 3) run_REFINE $Sample $Input_LIMA_directory $Output_directory
# 4) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
# 5) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
# 6) demux <input path read.stat file> <input path of samples file> <path of output>
# 7) run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
# 8) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
# 9) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>
# 10) run_sqanti2 <input_tofu_prefix> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir>
# 11) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
# 12) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
### ERCC
# 13) run_ERCC_analysis <sample_prefix_input/output_name> <isoseq3_input_directory> <output_dir>
### QC
# 14) ccs_bam2fasta <sample_name> <input_ccs.bam_dir> <output_dir>
# 15) lenghts <sample> <prefix.fasta> <input_dir> <output_dir>
# 16) make_file_for_rarefaction <sample_name_prefix> <input_tofu_directory> <input_sqanti_tama_directory> <working_directory
# 17) parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name>
### RNA-Seq defined transcriptome
# 18) run_stringtie <input_reference_gtf> <input_mapped_dir> <output_stringtie_dir> <output_prefix>
# 19) run_longshort_gffcompare <longread_gtf> <shortread_gtf> <output_dir> <output_name>
# 20) run tama merge <isoseq_gtf> <rnaseq_gtf> <isoseq_gtf_name> <rnaseq_gtf_name> <output_dir> <tama_output_name>
### Iso-Seq vs RNA-Seq defined transcriptome
# 21) run_counts_isoseqrnaseqtranscriptome <alignmentgtf> <fwd_rnaseq> <rev_rnaseq> <output_name> <output_dir>
### ALternative Splicing
# 22) run_suppa2 <input_gtf> <input_class> <output_dir> <output_name>
### Human MAPT
# 23) find_humanMAPT <cluster_dir> <output_dir> <merged_cluster.fa>
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

################################################################################################
#*************************************  Isoseq3 [Function 1, 2, 3, 4]
# run_CCS <input_ccs_bam> <prefix_output_name> <Output_directory>
run_CCS(){
  source activate isoseq3

  cd $3
  echo "Processing Sample $2 from Functions script"
  echo "Processing $1"
  echo "Checking if $2.ccs.bam exists"

    if [ -f $2.ccs.bam ]; then
      echo "$2.ccs.bam file already exists; CCS no need to be processed on Sample $2"
    else
      echo "$2.ccs.bam file does not exist"
      echo "Processing CCS for sample $2"
      # ccs <input.subreads.bam> <output.ccs.bam>
      time ccs $1 $2.ccs.bam --minPasses=1 --min-rq 0.9 --reportFile $2_ccs_report.txt
      echo "CCS for Sample $2 successful"
      ls *$2*
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

refine2fasta(){
  input_dir=$1
  output_dir=$2

  Samples=$(echo "${@:3}")

  REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
  CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
  cd $output_dir
  source activate sqanti2_py3
  for i in ${Samples[@]}; do echo "Processing: $i"; samtools bam2fq $input_dir/$i.flnc.bam | seqtk seq -A > $i.fa; done
  cat *fa* > All_flnc.fasta
  minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $REFERENCE/mm10.fa All_flnc.fasta > All.sam 2> All.map.log
  samtools sort -O SAM All.sam > All.sorted.sam
  python $CUPCAKE/sam_to_gff3.py All.sorted.sam -i All_flnc.fasta -s mm10
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
  module load Miniconda2/4.3.21
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
    minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $REFERENCE/mm10.fa $2/$1.clustered.hq.fastq > $1.sam 2> $1.map.log
    samtools sort -O SAM $1.sam > $1.sorted.sam

    # Alignment stats
    mkdir PAF; cd PAF
    source activate nanopore
    htsbox samview -pS $3/$1.sorted.sam > $1.paf
    awk -F'\t' '{if ($6="*") {print $0}}' $1.paf > $1.allread.paf # all reads
    awk -F'\t' '{if ($6=="*") {print $0}}' $1.paf > $1.notread.paf
    awk -F'\t' '{if ($6!="*") {print $0}}' $1.paf > $1.filtered.paf
    awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $1.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $1"_reads_with_alignment_statistics.txt"
    echo "Number of mapped transcripts to genome:"
    wc -l $1.filtered.paf
    echo "Number of ummapped transcripts to genome:"
    wc -l $1.notread.paf

    echo "Processing Sample $1 for TOFU, with coverage 85% and identity 95%"
    source activate cupcake
    head $CUPCAKE/README.md
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
# run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>
# Aim: Individually align RNASeq filtered fastq files to STAR-indexed genome using STAR, and subsequently input SJ.out.bed file into SQANTI2
run_star(){
    source activate sqanti2
    echo "STAR Version"
    STAR --version
    # samtools version: Version: 1.9 (using htslib 1.9)

    # extract official sample name from all_filtered directory
    # extract only files with "fastq.filtered", beginning with sample name, and R1/R2

    # cd to $J20/Tg4510_input_directory
    F_name=$(find $2 -name "*fastq.filtered" -exec basename \{} \; | grep ^$1 | grep "R1" )
    R_name=$(find $2 -name "*fastq.filtered" -exec basename \{} \; | grep ^$1 | grep "R2" )
    # save path directory of files as variable for later mapping
    F_File=$(find $2 -name "$F_name")
    R_File=$(find $2 -name "$R_name")

    # cd to $WKD
    cd $3
    if [ -f $1.SJ.out.bed ]; then
        echo "$1.SJ.out.bedalready exists; STAR Mapping not needed on Sample $1"
    else
        mkdir $1
        cd $1
        echo "Processing Sample $1 for STAR"
        echo "Processing Forward Reads: $F_File"
        echo "Processing Reverse Reads: $R_File"

        STAR_reference_input_directory=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/STAR_main


      # Parameters according to https://github.com/jennylsmith/Isoseq3_workflow/blob/master/shell_scripts/8__STAR_Junctions_ShortReads.sh
      # 2-pass mode alignment with chimeric read detection
      # at least 25 bp of one read of a pair maps to a different loci than the rest of the read pair
      # require 20 pb on either side of chimeric junction
      # include chimeric reads in the output BAM
      # don't include chimeras that map to reference genome contig with Ns
      # --outSAMtype BAM SortedByCoordinate, output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command.
      # note STAR indexed with genocde vM22 primary assembly annotation gtf therefore no need to include into command line, otherwise create additional folders
        STAR --runMode alignReads --runThreadN 32 --genomeDir $STAR_reference_input_directory \
        --readFilesIn $F_File $R_File \
        --outSAMtype BAM SortedByCoordinate \
        --chimSegmentMin 25 \
    		--chimJunctionOverhangMin 20 \
    		--chimOutType WithinBAM \
    		--chimFilter banGenomicN \
    		--chimOutJunctionFormat 1 \
    		--twopassMode Basic \
    		--twopass1readsN -1 #use all reads
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
    fi
    source deactivate
}

# mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>
# Aim: Merge forward and reverse of all RNASeq filtered fastq, for further downstream alignment to IsoSeq Tofu output using Kallisto
mouse_merge_fastq(){
    R1_READS=($(
        for i in ${SAMPLES_NAMES[@]}; do
            F_name=$(find $1 -name "*fastq.filtered" -exec basename \{} \; | grep ^$i | grep "R1" )
            F_File=$(find $1 -name "$F_name")
            echo $F_File
        done
    ))

    R2_READS=($(
        for i in ${SAMPLES_NAMES[@]}; do
            R_name=$(find $1 -name "*fastq.filtered" -exec basename \{} \; | grep ^$i | grep "R2" )
            R_File=$(find $1 -name "$R_name")
            echo $R_File
        done
    ))

    R1_READS_MERGE=$(echo "${R1_READS[@]}")
    R2_READS_MERGE=$(echo "${R2_READS[@]}")

    echo "Processing R1 READS"
    echo $R1_READS_MERGE
    echo "Processing R2 READS"
    echo $R2_READS_MERGE

    cd $2
    cat $R1_READS_MERGE > $3_R1.fq
    cat $R2_READS_MERGE > $3_R2.fq
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
# run_sqanti2 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/rnaseq/lncrna>
run_sqanti2(){

    source activate sqanti2_py3

    # variables
    SQANTI2_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    CUPCAKE_SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence
    REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
    CAGE_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/CAGE

    # copy STAR output SJ.out.bed files
    SJ_OUT_BED=($(
        for rnaseq in ${SAMPLES_NAMES[@]}; do
            name=$(find $4 -name "*SJ.out.bed" -exec basename \{} \; | grep ^$rnaseq)
            File=$(find $4 -name "$name")
            echo $File
        done
        ))
    for file in ${SJ_OUT_BED[@]}; do cp $file $7 ;done

    # prepare sqanti
    cd $7
    export PYTHONPATH=$PYTHONPATH:$CUPCAKE_SEQUENCE
    python $SQANTI2_DIR/sqanti_qc2.py -v

    # sqanti qc
    echo "Processing Sample $1 for SQANTI2 QC"

    # no kalliso file
    if [ $8 == "rnaseq" ]; then
      python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed --coverage "./*SJ.out.bed" --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt &> $1.sqanti.qc.log

    elif [ $8 == "lncrna" ]; then
      echo "Processing with lncRNA.gtf for genome annotation "
      python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM25.long_noncoding_RNAs.gtf $REFERENCE/mm10.fa --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed --coverage "./*SJ.out.bed" --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --skipORF --fl_count $6 &> $1.sqanti.qc.log

    elif [ $8 == "genome" ]; then
      echo "Processing with gencode.vM22.annotation.gtf for genome annotation "
      python $SQANTI2_DIR/sqanti_qc2.py -t 30 --gtf $3/$2 $REFERENCE/gencode.vM22.annotation.gtf $REFERENCE/mm10.fa --cage_peak $CAGE_DIR/mm10.cage_peak_phase1and2combined_coord.bed --coverage ./*SJ.out.bed --polyA_motif_list $SQANTI2_DIR/../human.polyA.list.txt --expression $5 --fl_count $6 &> $1.sqanti.qc.log
    else
      echo "8th argument required"
    fi

    echo "Processing Sample $1 for SQANTI2 filter"
    python $SQANTI2_DIR/sqanti_filter2.py $1"_classification.txt" $1"_corrected.fasta" $1"_corrected.gtf" -a 0.6 -c 3 &> $1.sqanti.filter.log

    # remove temp SJ.out bed files
    rm *SJ.out.bed
    source deactivate
}



################################################################################################
#************************************* TAMA [Function 11, 12]
# TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
# remove short fragments from post tofu
# Prerequisite: Require TAMA_prepare.R to change column 4 of bed12 to gene_name: transcript_name for correct file TAMA format
TAMA_remove_fragments(){

    TAMAFUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/TAMA
    TAMA_DIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/tama/tama_go/filter_transcript_models
    source activate sqanti2_py3
    cd $3

    # convert gtf to bed12
	  gtfToGenePred $1 $2.genepred
	  genePredToBed $2.genepred > $2.bed12
    awk -F'\t' '{print $1,$2,$3,$4,"40",$6,$7,$8,"255,0,0",$10,$11,$12}' $2.bed12| sed s/" "/"\t"/g|sed s/",\t"/"\t"/g|sed s/",$"/""/g > Tama_$2.bed12

    # Rscript script.R <name>.bed12 <input_dir>
	  Rscript $TAMAFUNCTIONS/TAMA_Merge_Prepare.R Tama_$2 $3
    source activate sqanti2
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
  awk '{ print $4 }' $1| cut -d ";" -f 2  > tama_retained_pbid.txt
  python $GENERALFUNC/TAMA/tama_sqanti_fastasubset.py $2/$5 $8/tama_retained_pbid.txt $8/$7"_sqantifiltered_tamafiltered_classification.fasta"

  # reinsert the quotation marks around the transcript id, gene id etc as required for recognition by SQANTI2
  # note, these are removed after processing in R
  sed 's/transcript_id \([^;]\+\)/transcript_id \"\1\"/g' $7"_sqantitamafiltered.classification.gtf" | sed 's/gene_id \([^;]\+\)/gene_id \"\1\"/g' | sed 's/gene_name \([^;]\+\)/gene_name \"\1\"/g' | sed 's/ref_gene_id \([^;]\+\)/ref_gene_id \"\1\"/g' > $7"_sqantitamafiltered.final.classification.gtf"
}


################################################################################################
#************************************* ERCC [Function 13]
# run_ERCC_analysis <sample_prefix_input/output_name> <isoseq3_input_directory>
run_ERCC_analysis(){

    CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake
    SQANTI2_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/SQANTI2
    REFERENCE_ERCC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/ERCC
    ANNOTATION=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/annotation
    SEQUENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake/sequence

    cd $3 #cd to output

    # Mapping
    source activate sqanti2_py3
    echo "Processing Sample $1 for Minimap2 and sort"
    minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $REFERENCE_ERCC/ERCC92.fa $2/$1.clustered.hq.fastq > $1.sam 2> $1.map.log
    samtools sort -O SAM $1.sam > $1.sorted.sam

    htsbox samview -pS $1.sorted.sam > $1.paf
    awk -F'\t' '{if ($6="*") {print $0}}' $1.paf > $1.allread.paf # all reads
    awk -F'\t' '{if ($6=="*") {print $0}}' $1.paf > $1.notread.paf
    awk -F'\t' '{if ($6!="*") {print $0}}' $1.paf > $1.filtered.paf
    awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $1.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $1"_reads_with_alignment_statistics.txt"

    # Cupcake
    coverage_threshold=0.85
    identity_threshold=0.95
    echo "Processing Sample $1 for TOFU, with coverage $coverage_threshold and identity $identity_threshold"
    source activate cupcake
    head $CUPCAKE/README.md
    collapse_isoforms_by_sam.py -c $coverage_threshold -i $identity_threshold --input $2/$1.clustered.hq.fastq --fq -s $1.sorted.sam --dun-merge-5-shorter -o $1 &>> $1.collapse.log
    get_abundance_post_collapse.py $1.collapsed $2/$1.clustered.cluster_report.csv &>> $1.abundance.log

    # SQANTI
    source activate sqanti2_py3
    export PYTHONPATH=$PYTHONPATH:$SEQUENCE
    prefix=$1.collapsed
    python $SQANTI2_dir/sqanti_qc2.py -v
    python $SQANTI2_dir/sqanti_qc2.py -t 30 --gtf $prefix.gff $REFERENCE_ERCC/ERCC92.gtf $REFERENCE_ERCC/ERCC92.fa --fl_count $prefix.abundance.txt &>> $1.sqanti.qc.log
    python $SQANTI2_dir/sqanti_filter2.py $prefix"_classification.txt" $prefix"_corrected.fasta" $prefix"_corrected.gtf" -a 0.6 -c 3 &>> $1.sqanti.filter.log
    TAMA_remove_fragments WholeIsoSeq.collapsed_classification.filtered_lite.gtf WholeIsoSeq $3
    TAMA_remove_fragments WholeIsoSeq.collapsed_corrected.gtf WholeIsoSeq.collapsed $3 # from SQANTI classification not filter
    TAMA_sqanti_filter WholeIsoSeq.bed $3 WholeIsoSeq.collapsed_classification.filtered_lite_classification.txt WholeIsoSeq.collapsed_classification.filtered_lite.gtf WholeIsoSeq.collapsed_classification.filtered_lite.fasta WholeIsoSeq.collapsed_classification.filtered_lite_junctions.txt WholeIsoSeq $3

    source deactivate
}

################################################################################################
#************************************* QC [Function 14 - 17]
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

# lenghts <sample> <prefix.fasta> <input_dir> <output_dir>
lengths(){
  # variables
  sample=$1
  suffix=$2
  input_dir=$3
  output_dir=$4

  echo "Extracting length of Sample $sample"
  source activate sqanti2_py3
  CUPCAKE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/Post_Isoseq3/cDNA_Cupcake
  export PYTHONPATH=$PYTHONPATH:$CUPCAKE/sequence
  export PATH=$PATH:$CUPCAKE/sequence
  cd $output_dir
  get_seq_stats.py $input_dir/$sample$suffix &>> $sample.get_seq.log
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

# parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name>
parse_stats_per_sample(){

    # variable
    CCS_dir=$1
    LIMA_dir=$2
    output_name=$3

    FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
    source activate sqanti2_py3
    python $FUNCTIONS/IsoSeq_QC/CCS.py $CCS_dir "" $3
    python $FUNCTIONS/IsoSeq_QC/LIMA.py $LIMA_dir "" $3
    source deactivate
}


################################################################################################
#************************************* RNA-Seq defined transcriptome [Function 18]
# RNASeq alone
# run_stringtie <input_reference_gtf> <input_mapped_dir> <output_stringtie_dir> <output_prefix>
# input: Sample names taken from working script
# prerequisite: run star to generated sorted_bam
run_stringtie(){

    ### variables ###
    input_reference_gtf=$1
    input_mapped_dir=$2
    output_stringtie_dir=$3
    output_prefix=$4

    export PATH=$PATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/stringtie
    echo "Stringtie version:" stringtie --version
    cd $output_stringtie_dir

    # run stringtie individually on rnaseq sorted bam file
    for i in ${SAMPLES_NAMES[@]}; do
       echo "Processing stringtie for sample $i against $input_reference_gtf"
       #stringtie -G <reference_gtf> -o <output_gtf> <input mapped rnaseq bam>
       stringtie -G $input_reference_gtf -o $i.out.gtf $input_mapped_dir/$i/$i.sorted.bam -v &> $i.stringtie.log
    done

    # merge stringtie gtf output to one merged gtf
    all_stringtie_gtf=(*gtf*)
    all_stringtie_gtf_names=$(echo "${all_stringtie_gtf[@]}")
    echo "Merging: $all_stringtie_gtf_names"
    stringtie --merge -G $input_reference_gtf -p 8 -o $output_prefix.gtf $all_stringtie_gtf_names

    # modify stringtie gtf for input into SQANTI2 due to some strands with no strand annotation (blips)
    # https://github.com/gpertea/stringtie/issues/322
    source activate sqanti2_py3
    FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation
    # Rscript .R <path; input gtf> <path; output gtf>
    Rscript $FUNCTIONS/RNASeq/Subset_Stringtiegtf.R $output_stringtie_dir/$output_prefix".gtf" $output_stringtie_dir/$output_prefix"_mod.gtf"

    # reinsert the quotation marks around the transcript id, gene id etc as required for recognition by SQANTI2
    # note, these are removed after processing in R
    sed 's/transcript_id \([^;]\+\)/transcript_id \"\1\"/g' $output_prefix"_mod.gtf" | sed 's/gene_id \([^;]\+\)/gene_id \"\1\"/g' | sed 's/gene_name \([^;]\+\)/gene_name \"\1\"/g' | sed 's/ref_gene_id \([^;]\+\)/ref_gene_id \"\1\"/g' > $output_prefix"_final.gtf"

}

# run_longshort_gffcompare <longread_gtf> <shortread_gtf> <output_dir> <output_name>
run_longshort_gffcompare(){
  # variable
  longread_gtf=$1
  shortread_gtf=$2
  output_dir=$3
  output_name=$4

  cd $output_dir

  # reinsert the quotation marks around the transcript id, gene id etc as required for recognition by SQANTI2
  # note, these are removed after processing in R
  echo "Processing: $longread_gtf"
  echo "Against: $shortread_gtf"
  sed 's/transcript_id \([^;]\+\)/transcript_id \"\1\"/g' $longread_gtf | sed 's/gene_id \([^;]\+\)/gene_id \"\1\"/g' | sed 's/gene_name \([^;]\+\)/gene_name \"\1\"/g' | sed 's/ref_gene_id \([^;]\+\)/ref_gene_id \"\1\"/g' > WholeIsoSeq_longread.gtf

  export PATH=$PATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/gffcompare
  gffcompare -R -r WholeIsoSeq_longread.gtf -o $output_name $shortread_gtf
}


# run tama merge <isoseq_gtf> <rnaseq_gtf> <isoseq_gtf_name> <rnaseq_gtf_name> <output_dir> <tama_output_name>
run_tama_merge(){

  source activate sqanti2

  # variables
  input_gtf1=$1
  input_gtf2=$2
  input_gtf1_name=$3
  input_gtf2_name=$4
  output_dir=$5
  tama_output_name=$6


  GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/TAMA
  TAMA=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/tama

  # Modify_genomegtf_TAMAinput <input_gtf> <output_name>
  Modify_genomegtf_TAMAinput(){

    # variables
    input_gtf=$1
    output_name=$2
    output_dir=$3
    typegtf=$4

    source activate sqanti2

    echo "Preparing $input_gtf for TAMA merge"

    gtfToGenePred $input_gtf $output_dir/$output_name"_annotation.genepred"
    genePredToBed $output_dir/$output_name"_annotation.genepred" $output_dir/$output_name"_annotation.bed12"


    ###################################
    # change format specifically for tama merge (tab format)
    ###################################
    awk -F'\t' '{print $1,$2,$3,$4,"40",$6,$7,$8,"255,0,0",$10,$11,$12}' $output_dir/$output_name"_annotation.bed12"|\
    sed s/" "/"\t"/g|sed s/",\t"/"\t"/g|sed s/",$"/""/g > $output_dir/Tama_$output_name"_annotation.bed12"

    ###################################
    # Modify column 4 of tab format to include "gene_id; transcript_id" as specified for TAMA
    ## WTAC_scripts/assign_gene_id_gene_name_to_isoform_list_for_TAMA_merge.py <list_of_gene_transcript_id_path> <input_bed12_format> <output_bed12_format>
    ###################################
    if [ $typegtf == "RNASeq_mousegenome" ]; then
      TAMA_ref=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/TAMA_REFERENCE/gene_id_transcript_id_gene_name_mm10vM22_ERCC_SIRV
      WTAC_scripts=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/Nanopore/WTAC_scripts
      echo "Converting bed12 with usable bed12 using the genome"
      echo "$TAMA_ref"
      # Rscript script.R <input_gtf> <output_file>
      Rscript $GENERALFUNC/MSTRG_format.R $input_gtf $output_dir/rnaseq_gene_id_MSTRG
      cat $TAMA_ref $output_dir/rnaseq_gene_id_MSTRG > $output_dir/rnaseq_finalgeneid
      python $WTAC_scripts/assign_gene_id_gene_name_to_isoform_list_for_TAMA_merge.py $output_dir/rnaseq_finalgeneid $output_dir/Tama_$output_name"_annotation.bed12" $output_dir/Tama_formatted_$output_name"_annotation.bed12"
    elif [ $typegtf == "Isoseq" ]; then
      echo "Converting bed12 with usable bed12 using the info from isoseq gtf"
      #Rscript script.R <name>.bed12 <input_dir>
      Rscript $GENERALFUNC/TAMA_Merge_Prepare.R Tama_$output_name"_annotation" $output_dir
    else
      echo "4th argument required as RNASeq_mousegenome or Isoseq"
    fi
  }


  # Modify gtf for TAMA input
  cd $output_dir
  #Modify_genomegtf_TAMAinput $input_gtf1 $input_gtf1_name $output_dir Isoseq
  #Modify_genomegtf_TAMAinput $input_gtf2 $input_gtf2_name $output_dir RNASeq_mousegenome
  #mv Tama_$input_gtf1_name"_annotation_mod.bed12" Tama_formatted_$input_gtf1_name"_annotation.bed12"

  echo "Merging $input_gtf1_name and $input_gtf2_name"
  ###################################
  # Create input file.list for Tama Merge input
  # 1st line: refers to genome input (and therefore genome gff)
  # 2nd line: refers to sample input (and therefore gff input)
  # sample_gff = sample name with path directory removed user-defined 2nd argument
  ###################################
  echo "$output_dir/Tama_formatted_$input_gtf1_name"_annotation.bed12":no_cap:1,1,1:$input_gtf1_name"|sed s/":"/"\t"/g> $output_dir/file.list
  echo "$output_dir/Tama_formatted_$input_gtf2_name"_annotation.bed12":no_cap:1,1,1:$input_gtf2_name"|sed s/":"/"\t"/g>>$output_dir/file.list
  cat $output_dir/file.list
  python $TAMA/tama_merge.py -f file.list -a 50 -z 50 -m 20 -p $tama_output_name &> $tama_output_name"_tama_merge.log"

  source deactivate
}

################################################################################################
echo "#************************************* Iso-Seq vs RNA-Seq defined transcriptome [Function 21]"

# run_counts_isoseqrnaseqtranscriptome <alignmentgtf> <fwd_rnaseq> <rev_rnaseq> <output_name> <output_dir>
# aim: to align merged rnaseq reads to iso-seq defined or rna-seq defined transcriptome using STAR, followed by featureCounts with the transcript mode
run_counts_isoseqrnaseqtranscriptome(){
  # variable
  mousegtf=$1
  fwd_rnaseq=$2
  rev_rnaseq=$3
  output_name=$4
  output_dir=$5

  REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019
  STAR_reference_input_directory=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/STAR_main
  source activate sqanti2
  echo "STAR Version"
  STAR --version

  cd $output_dir
  echo "STAR alignment with $fwd_rnaseq and $rev_rnaseq"
  STAR --runMode alignReads --runThreadN 32 --genomeDir $STAR_reference_input_directory --readFilesIn $fwd_rnaseq $rev_rnaseq --outSAMtype BAM SortedByCoordinate --chimSegmentMin 25 --chimJunctionOverhangMin 20 --chimOutType WithinBAM --chimFilter banGenomicN --chimOutJunctionFormat 1 --twopassMode Basic --twopass1readsN -1
  mv Aligned.sortedByCoord.out.bam $output_name.sorted.bam

  echo "Featurecounts with transcript mode"
  featureCounts -T 8 -O -g transcript_id -p -a  $mousegtf -o $output_name.transcript_id.tsv $output_name.sorted.bam 2> $1.featurecounts.log
}

################################################################################################
#************************************* Alternative Splicing [Function 22]

# run_suppa2 <input_gtf> <input_class> <output_dir> <output_name>
# input_gtf: sqanti tama filtered gtf
# input_class: sqanti tama filtered classification file
run_suppa2(){
  # variable
  input_gtf=$1
  input_class=$2
  output_dir=$3
  output_name=$4

  FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation

  source activate sqanti2_py3
  cd $output_dir
  suppa.py generateEvents -i $input_gtf -o $output_name --pool-genes -f ioe -e SE MX FL SS RI &>> $output_name"_suppa2.log"
  python $FUNCTIONS/Suppa2_output_mod_updated.py $output_name $output_dir $input_class &>> $output_name"_suppa2_mod.log"
}

################################################################################################
#************************************* Find human MAPT [Function 23]

# find_humanMAPT <cluster_dir> <output_dir> <merged_cluster.fa>
# cluster_dir = directory containing hq.fasta
# aim: to grep only the clustered hq reads with human MAPT sequence in all the samples in the clustered directory
find_humanMAPT(){
  hMAPT_1=TGGTTAATCACTTAACCTGCTTTTGTCACTCGGCTTTGGCTCGGGACTTCAAAATCAGTGATGGGAGTAAGAGCAAATTTCATCTTTCCAAATTGATGGGTGGGCTAGTAATAAAATATTTAAAAAAAAACATTCAAAAACATGGCCACATCCAACATTTCCTCAGGCAATTCCTTTTGATTCTTTTTTCTTCCCCCTCCATGTA
  hMAPT_2=AAAATCAGTGATGGGAGTAAGAGCAAATTTCATCTTTCCAAATTGATGGGTGGGCTAGTAATAAAATATTTAAAAAAAAACATTCAAAAACATGGCCACATCCAACATTTCCTCAGGCAATTCCTTTTGATTCTTTTTTCTTCCCCCTCCATGTAGAAGAGGGAGAAGGAGAGGCTCTGAAAGCTGCTTCTGGGGGATTT
  mMAPT_1=GGGGGGTGGTATTCTGGGATGTGGGTCCCAGGCCTCCCATCCCTCACACAGCCACTGTATCCCCTCTCTCTGTCCTATCATGCCCACGTCTGCCACGAGAGCTAGTCACTGCCGTCCGTACATCACGTCTCACTGTCCTGAGTGCCATGC

  # variables
  cluster_dir=$1
  output_dir=$2
  merged_cluster=$3

  for i in $cluster_dir/*fasta; do
    echo "Processing with $i"
    sample=$(basename "$i" | cut -d "." -f 1)

    # grep the reads with the human MAPT1 and human MAPT2 sequences
    # include just the headers for a count of the number of reads
    grep -B1 $hMAPT_1 $i > $output_dir/$sample.hMAPT1.clustered.hq.fasta
    grep "^>" $output_dir/$sample.hMAPT1.clustered.hq.fasta > $output_dir/$sample.hMAPT1.header

    grep -B1 $mMAPT_1 $i > $output_dir/$sample.mMAPT1.clustered.hq.fasta
    grep "^>" $output_dir/$sample.mMAPT1.clustered.hq.fasta > $output_dir/$sample.mMAPT1.header

    grep -B1 $hMAPT_2 $i > $output_dir/$sample.hMAPT2.clustered.hq.fasta
    grep "^>" $output_dir/$sample.hMAPT2.clustered.hq.fasta > $output_dir/$sample.hMAPT2.header
  done

  # Merged Cluster
  source activate sqanti2_py3
  GENERALFUNC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/General/2_Transcriptome_Annotation/IsoSeq_QC
  #script.py <path/input.fasta> <outputname> <path/outputdir>
  merged_sample=$(basename "$merged_cluster" | cut -d "." -f 1)
  python $GENERALFUNC/Find_Human_Mapt.py $merged_cluster $merged_sample $output_dir > $output_dir/$merged_sample.out
  mv $output_dir/$merged_sample $output_dir/$merged_sample.clustered.hq.fasta
}
