# Aaron Jeffries: A.R.Jeffries@exeter.ac.uk

~/hub/genePredToBed adult.sample.collapsed.filtered.rep_corrected.genePred adult.sample.collapsed.filtered.rep_corrected.bed
refAnnotation_adult.sample.collapsed.filtered.rep.genePred
LC_COLLATE=C
sort -k1,1 -k2,2n adult.sample.collapsed.filtered.rep_corrected.bed > adult.bed
# Colourify script
~/hub/bedToBigBed -extraIndex=name adult.bed https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes adult.sample.bb

~/hub/genePredToBed refAnnotation_adult.sample.collapsed.filtered.rep.genePred adult.2.bed

~/hub/gtfToGenePred adult.sample.collapsed.filtered.rep_classification.filtered_lite.gtf adult.genePred

#Convert to genePred
~/hub/gtfToGenePred WT8IsoSeq.collapsed.filtered.rep_classification.filtered_lite.gtf mouse_lite.genePred
~/hub/gtfToGenePred WT8IsoSeq.collapsed.filtered.rep_corrected.gtf mouse.genePred
~/hub/gtfToGenePred adult.sample.collapsed.filtered.rep_classification.filtered_lite.gtf adult_lite.genePred
~/hub/gtfToGenePred adult.sample.collapsed.filtered.rep_corrected.gtf adult.genePred
~/hub/gtfToGenePred combinedFL.sample.collapsed.filtered.rep_classification.filtered_lite.gtf human_lite.genePred
~/hub/gtfToGenePred combinedFL.sample.collapsed.filtered.rep_corrected.gtf human.genePred
~/hub/gtfToGenePred fetalFL.sample.collapsed.filtered.rep_classification.filtered_lite.gtf fetalCTX_lite.genePred
~/hub/gtfToGenePred fetalFL.sample.collapsed.filtered.rep_corrected.gtf fetalCTX.genePred
~/hub/gtfToGenePred fetalHIP.sample.collapsed.filtered.rep_classification.filtered_lite.gtf fetalHIP_lite.genePred
~/hub/gtfToGenePred fetalHIP.sample.collapsed.filtered.rep_corrected.gtf fetalHIP.genePred
~/hub/gtfToGenePred fetalSTR.sample.collapsed.filtered.rep_classification.filtered_lite.gtf fetalSTR_lite.genePred
~/hub/gtfToGenePred fetalSTR.sample.collapsed.filtered.rep_corrected.gtf fetalSTR.genePred

#Convert to bed12 format and sort
~/hub/genePredToBed mouse_lite.genePred mouse_lite.tmp
~/hub/genePredToBed mouse.genePred mouse.tmp
~/hub/genePredToBed adult_lite.genePred adult_lite.tmp
~/hub/genePredToBed adult.genePred adult.tmp
~/hub/genePredToBed human_lite.genePred human_lite.tmp
~/hub/genePredToBed human.genePred human.tmp
~/hub/genePredToBed fetalCTX_lite.genePred fetalCTX_lite.tmp
~/hub/genePredToBed fetalCTX.genePred fetalCTX.tmp
~/hub/genePredToBed fetalHIP_lite.genePred fetalHIP_lite.tmp
~/hub/genePredToBed fetalHIP.genePred fetalHIP.tmp
~/hub/genePredToBed fetalSTR_lite.genePred fetalSTR_lite.tmp
~/hub/genePredToBed fetalSTR.genePred fetalSTR.tmp
LC_COLLATE=C
sort -k1,1 -k2,2n adult.tmp > adult.bed12
sort -k1,1 -k2,2n adult_lite.tmp > adult_lite.bed12
sort -k1,1 -k2,2n human.tmp > human.bed12
sort -k1,1 -k2,2n human_lite.tmp > human_lite.bed12
sort -k1,1 -k2,2n fetalCTX.tmp > fetalCTX.bed12
sort -k1,1 -k2,2n fetalCTX_lite.tmp > fetalCTX_lite.bed12
sort -k1,1 -k2,2n fetalHIP.tmp > fetalHIP.bed12
sort -k1,1 -k2,2n fetalHIP_lite.tmp > fetalHIP_lite.bed12
sort -k1,1 -k2,2n fetalSTR.tmp > fetalSTR.bed12
sort -k1,1 -k2,2n fetalSTR_lite.tmp > fetalSTR_lite.bed12
sort -k1,1 -k2,2n mouse.tmp > mouse.bed12
sort -k1,1 -k2,2n mouse_lite.tmp > mouse_lite.bed12

#Run R-script for colouring bed12 file
#Then convert to bigBed (binary) for hub integration

~/hub/bedToBigBed -extraIndex=name adult.bed12 https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes adult.bb
~/hub/bedToBigBed -extraIndex=name adult_lite.bed12 https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes adult_lite.bb
~/hub/bedToBigBed -extraIndex=name human.bed12 https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes human.bb
~/hub/bedToBigBed -extraIndex=name human_lite.bed12 https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes human_lite.bb
~/hub/bedToBigBed -extraIndex=name fetalCTX.bed12 https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes fetalCTX.bb
~/hub/bedToBigBed -extraIndex=name fetalCTX_lite.bed12 https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes fetalCTX_lite.bb
~/hub/bedToBigBed -extraIndex=name fetalSTR.bed12 https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes fetalSTR.bb
~/hub/bedToBigBed -extraIndex=name fetalSTR_lite.bed12 https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes fetalSTR_lite.bb
~/hub/bedToBigBed -extraIndex=name fetalHIP.bed12 https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes fetalHIP.bb
~/hub/bedToBigBed -extraIndex=name fetalHIP_lite.bed12 https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes fetalHIP_lite.bb
~/hub/bedToBigBed -extraIndex=name mouse.bed12 http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes mouse.bb
~/hub/bedToBigBed -extraIndex=name mouse_lite.bed12 http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes mouse_lite.bb