#!/usr/bin/env Rscript

pacman::p_load(data.table,doParallel,tidyverse,ff,foreach,optparse,Rcpp,zoo)

gc()
rm(list=ls())
# no_cores <- detectCores()-1
# registerDoParallel(no_cores)
args <- commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("-n", "--nnpl"), type="character", default=NULL,
              help="nanopolish", metavar="character"),
  make_option(c("-t", "--tsv"), type="character", default=NULL,
              help="tsv", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="outname.inAE",
              help="output file name", metavar="character")			  
)

opt_parser <- OptionParser(option_list=option_list,add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$nnpl)|is.null(opt$tsv)){
  print_help(opt_parser)
  stop("Input files must be supplied (-n,-t).", call.=FALSE)
}

getdinodir <- function(){
    commandArgs() %>%
       tibble::enframe(name=NULL) %>%
       tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
       dplyr::filter(key == "--file") %>%
       dplyr::pull(value) %>%
       word(., start=1, end=-3, sep="/")
}
dinodir <- getdinodir()

Rcpp::sourceCpp(paste0(dinodir,"/code/all_functions.cpp"))

print('Start')

tmp.bam <- fread(file=paste("./",opt$tsv,sep=""),quote = "", header=TRUE, sep="\t")
print('Done loading')

names(tmp.bam)[2]="strand"
tmp.bam$REF=toupper(tmp.bam$REF)
#Add quality coulumn
tmp.bam$REF_POS=as.numeric(tmp.bam$REF_POS)
tmp.bam$READ_POS=as.numeric(tmp.bam$READ_POS)
tmp.bam$QUAL=as.factor(tmp.bam$QUAL)

tmp.bam[,"qual"] <- 0
myascii(tmp.bam$QUAL,tmp.bam$qual)
#Add mis/match + deletion columns
tmp.bam[,"mat"] <- 0
tmp.bam[,"del"] <- 0
tmp.bam[,"mis"] <- 0
mymatmis(tmp.bam$OP,tmp.bam$REF,tmp.bam$BASE,tmp.bam$mat,tmp.bam$mis,tmp.bam$del)

tmp.bam[,"ins"] <- 0
tmp.bam.negative <- tmp.bam[tmp.bam$strand == "n",]
tmp.bam.positive <- tmp.bam[tmp.bam$strand == "p",]
rm(tmp.bam)
#Add insertion column
insertfunP(tmp.bam.positive$OP,tmp.bam.positive$ins)
insertfunN(tmp.bam.negative$OP,tmp.bam.negative$ins)
print("Done features")

#Filter file
tmp.bam2.positive <- tmp.bam.positive[(tmp.bam.positive$OP!='I') & (tmp.bam.positive$REF_POS>0),]
rm(tmp.bam.positive)
readname.positive <- unique(tmp.bam2.positive$`#READ_NAME`)
merge.bam.positive <- tmp.bam2.positive[,c(1,2,3,5,7,8,10:14)]

rm(tmp.bam2.positive)

tmp.bam2.negative <- tmp.bam.negative[(tmp.bam.negative$OP!='I') & (tmp.bam.negative$REF_POS>0),]
rm(tmp.bam.negative)
readname.negative <- unique(tmp.bam2.negative$`#READ_NAME`)
merge.bam.negative <- tmp.bam2.negative[,c(1,2,3,5,7,8,10:14)]
rm(tmp.bam2.negative)
merge.bam.negative$BASE <- as.character(merge.bam.negative$BASE)
merge.bam.negative$BASE1=""
complement_char(merge.bam.negative$BASE,merge.bam.negative$BASE1)
merge.bam.negative=merge.bam.negative[,-c("BASE")]
names(merge.bam.negative)[11]="BASE"
print("Done filter bam file")


#S2:Nanopolish file
data.nnpl <- fread(file=paste("./",opt$nnpl,sep=""),header=TRUE, sep="\t", colClasses=c(rep("factor", 2),rep("numeric", 4),rep("factor", 1)))
merge.nnpl.postive <- data.nnpl[data.nnpl$read_name %in% readname.positive,]
merge.nnpl.postive$position <- merge.nnpl.postive$position+2#set coordinate to second position in kmer
merge.nnpl.negative <- data.nnpl[data.nnpl$read_name %in% readname.negative,]
rm(data.nnpl)
print("done positive strand for nanopolish")
merge.nnpl.negative$position <- merge.nnpl.negative$position+4#set coordinate to fourth position in kmer
merge.nnpl.negative$reference_kmer <- as.character(merge.nnpl.negative$reference_kmer)
merge.nnpl.negative$reference_kmer1 <- ""
rev_comp(merge.nnpl.negative$reference_kmer,merge.nnpl.negative$reference_kmer1)
merge.nnpl.negative=merge.nnpl.negative[,-c("reference_kmer")]
names(merge.nnpl.negative)[7]="reference_kmer"
print("done negative strand for nanopolish")

#S3:merge file
combine.postive <- merge.nnpl.postive[merge.bam.positive, on = c("read_name"="#READ_NAME","position"="REF_POS","contig"="CHROM")] #left join
rm(merge.nnpl.postive,merge.bam.positive)
combine.negative <- merge.nnpl.negative[merge.bam.negative, on = c("read_name"="#READ_NAME","position"="REF_POS","contig"="CHROM")]
rm(merge.nnpl.negative,merge.bam.negative)
print ("done merging files")

#S4:Add neighboring
data.sorted.positive <- combine.postive[order(combine.postive$read_name,combine.postive$position),]
rm(combine.postive)
data.sorted.negative <- combine.negative[order(combine.negative$read_name,combine.negative$position),]
rm(combine.negative)
print ("done adding neighboring positions")
gc()

#### compile features table for positive-stranded sequences ####
dat.split.positive <- NROW(unique(data.sorted.positive$read_name))/(20)
name.list <- split(unique(data.sorted.positive$read_name), trunc(0:(NROW(unique(data.sorted.positive$read_name))-1 )/dat.split.positive))
dat.com.positive <- foreach(b=1:length(name.list),.combine=rbind, .errorhandling="pass") %do%
  { 
    data.sorted.positive[data.sorted.positive$read_name %in% unlist(name.list[b]),][,c("qual1","qual2","qual4","qual5",
                                                                                       "del1","del2","del4","del5",
                                                                                       "mis1","mis2","mis4","mis5",
                                                                                       "mat1","mat2","mat4","mat5",
                                                                                       "ins1","ins2","ins4","ins5"):=.(lag(qual,n=2L),lag(qual,n=1L),lead(qual,n=1L),lead(qual,n=2L),
                                                                                                                      lag(del,n=2L),lag(del,n=1L),lead(del,n=1L),lead(del,n=2L),
                                                                                                                      lag(mis,n=2L),lag(mis,n=1L),lead(mis,n=1L),lead(mis,n=2L),
                                                                                                                      lag(mat,n=2L),lag(mat,n=1L),lead(mat,n=1L),lead(mat,n=2L),
                                                                                                                      lag(ins,n=2L),lag(ins,n=1L),lead(ins,n=1L),lead(ins,n=2L)),by=.(read_name)]
  }
rm(data.sorted.positive);gc()
dat.com.positive$posi_id=paste(dat.com.positive$contig,dat.com.positive$position,dat.com.positive$strand,sep=":")
fwrite(dat.com.positive,file=paste0(opt$out, ".positive"),,quote= F, sep="\t")
col=colnames(dat.com.positive)
rm(dat.com.positive);gc()
print ("Done positive strand")

#### compile features table for negative-stranded sequences ####
dat.split.negative <- NROW(unique(data.sorted.negative$read_name))/(20)
name.list <- split(unique(data.sorted.negative$read_name), trunc(0:(NROW(unique(data.sorted.negative$read_name))-1 )/dat.split.negative))
dat.com.negative <- foreach(b=1:length(name.list),.combine=rbind, .errorhandling="pass") %do%
  { 
    data.sorted.negative[data.sorted.negative$read_name %in% unlist(name.list[b]),][,c("qual1","qual2","qual4","qual5",
                                                                                       "del1","del2","del4","del5",
                                                                                       "mis1","mis2","mis4","mis5",
                                                                                       "mat1","mat2","mat4","mat5",
                                                                                       "ins1","ins2","ins4","ins5"):=.(lag(qual,n=2L),lag(qual,n=1L),lead(qual,n=1L),lead(qual,n=2L),
                                                                                                                       lag(del,n=2L),lag(del,n=1L),lead(del,n=1L),lead(del,n=2L),
                                                                                                                       lag(mis,n=2L),lag(mis,n=1L),lead(mis,n=1L),lead(mis,n=2L),
                                                                                                                       lag(mat,n=2L),lag(mat,n=1L),lead(mat,n=1L),lead(mat,n=2L),
                                                                                                                       lag(ins,n=2L),lag(ins,n=1L),lead(ins,n=1L),lead(ins,n=2L)),by=.(read_name)]
  }
rm(data.sorted.negative)
print ("done negative strand")
dat.com.negative$posi_id=paste(dat.com.negative$contig,dat.com.negative$position,dat.com.negative$strand,sep=":")
dat.com.negative=as.data.frame(dat.com.negative)
dat.com.negative=dat.com.negative[,col]
fwrite(dat.com.negative,file=paste0(opt$out, ".negative"),quote= F, sep="\t",col.names=F,append=T)
print("Done!")
