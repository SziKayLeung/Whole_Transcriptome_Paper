# Aaron Jeffries: A.R.Jeffries@exeter.ac.uk
setwd("C:/Users/arj207/Desktop/hubproc")

prefix <- c("adult_lite","adult","human_lite","human","fetalCTX_lite","fetalCTX","fetalHIP_lite","fetalHIP","fetalSTR_lite","fetalSTR","mouse","mouse_lite")
bed12 <- paste(prefix,".bed12",sep="")

sqantifiles <-c("adult.sample.collapsed.filtered.rep_classification.filtered_lite_classification.txt","adult.sample.collapsed.filtered.rep_classification.txt","combinedFL.sample.collapsed.filtered.rep_classification.filtered_lite_classification.txt","combinedFL.sample.collapsed.filtered.rep_classification.txt","fetalFL.sample.collapsed.filtered.rep_classification.filtered_lite_classification.txt","fetalFL.sample.collapsed.filtered.rep_classification.txt","fetalHIP.sample.collapsed.filtered.rep_classification.filtered_lite_classification.txt","fetalHIP.sample.collapsed.filtered.rep_classification.txt","fetalSTR.sample.collapsed.filtered.rep_classification.filtered_lite_classification.txt","fetalSTR.sample.collapsed.filtered.rep_classification.txt","WT8IsoSeq.collapsed.filtered.rep_classification.txt", "WT8IsoSeq.collapsed.filtered.rep_classification.filtered_lite_classification.txt")

#location of SQANTI2 classification files
class_filepath <- "C:/Users/arj207/Dropbox/Isoseq paper/Data/SQANTI and FUSION and GENE SUMMARY/sq2v7_4/"

#location of bed12 files
bed12_filepath <- "C:/Users/arj207/Desktop/hubproc/"

for(sq in 1:length(bed12)) {
bedfile <- paste(bed12_filepath,bed12[sq],sep="")
fileprefix <- prefix[sq]
sqantifile <- paste(class_filepath,sqantifiles[sq],sep="")

#first make a BED12 file externally

bed <- read.table(bedfile, stringsAsFactors = F) 
SQANTI <- read.table(sqantifile, sep="\t", stringsAsFactors = F, header=T)
x <- match(bed$V4, SQANTI$isoform) #return a lookup in SQANTI to match bed
  
  #Loop through bed and annotate
  for (b in 1:nrow(bed)) {
  fl <- SQANTI$FL[x[b]]
  type <- SQANTI$structural_category[x[b]]
  associated_gene <- paste(SQANTI$associated_gene[x[b]],".",unlist(strsplit(SQANTI$isoform[x[b]],".", fixed = T))[3],"/",SQANTI$FL[x[b]],sep="") 
    
  if(fl == 2) {intensity <- 200 }
  if(fl == 3) {intensity <- 300 }
  if(fl == 4) {intensity <- 400 }
  if(fl == 5) {intensity <- 500 }
  if(fl == 6) {intensity <- 600 }
  if(fl == 7) {intensity <- 700 }
  if(fl == 8) {intensity <- 800 }
  if(fl == 9) {intensity <- 900 }
  if(fl > 9) {intensity <- 1000 }
    
  if(type == "incomplete-splice_match") { genecolour <- "0,255,255" } #cyan
  if(type == "full-splice_match") { genecolour <- "0,0,255" } #blue
  if(type == "novel_not_in_catalog") { genecolour <- "255,0,0" } #red
  if(type == "novel_in_catalog") { genecolour <- "255,128,0" } #orange
  if(type == "antisense") { genecolour <- "255,0,255" } #purple
  if(type == "genic") { genecolour <- "128,128,128" } #grey
  if(type == "intergenic") { genecolour <- "128,128,128" } #grey
  

bed$V5[b] <- intensity
bed$V9[b] <- genecolour
bed$V4[b] <- associated_gene #rename geneID - consider adding to classification too?

  }
  #add header (dont use when making hub as can be encoded into TrackDb.txt and referred to same BED file for two annotations)
  #header <- c("track name='Test Track' description='Test Track' useScore=1 itemRgb=On","","","","","","","","","","","")
  #bed <- rbind(header,bed)
  write.table(bed,paste(fileprefix,".bed12",sep=""),row.names=F, quote=F,sep="\t", col.names=F)

}
  #use bedtools bamtobed -bed12 -i <input.sam> 