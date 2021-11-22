#!/usr/bin/env Rscript

#This script is used to combine signal from output of nanopolish. This script must be run under script s2.Combine_raw_nnpl.sh, don't run alone :)

pacman::p_load(foreach,doParallel,tidyverse,Rcpp,pracma,optparse,data.table)

no_cores <- detectCores()-1

args = commandArgs(trailingOnly=TRUE)
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Event file name from Nanopolish", metavar="character"),
  make_option(c("-t", "--thread"), type="integer", default=no_cores,
              help="Number of cores allocated", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="combined",
              help="output file name", metavar="character"),
  make_option(c("-s", "--size"), type="numeric", default=NULL,
              help="number of line in input", metavar="integer"),
  make_option(c("-n", "--num"), type="numeric", default=NULL,
              help="order number of line in input", metavar="integer"),
  make_option(c("-c", "--chunk"), type="numeric", default=1000000,
              help="chunk size", metavar="integer")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file) | is.null(opt$size)|is.null(opt$num)){
  print_help(opt_parser)
  stop("Input file and size of input must be supplied.", call.=FALSE)
}
no_cores=opt$thread
registerDoParallel(no_cores)

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

head=fread("raw_nanopolish.header")
n=opt$num
ch=opt$chunk
ind=c()
opt$size <- as.numeric(opt$size)

dat=fread(opt$file,header=F,nThread=no_cores)
colnames(dat)=colnames(head)
if((n+ch)< opt$size){	
    read=dat$read_name[nrow(dat)]
    pos=dat$position[nrow(dat)]
    ind=which(dat$position==pos & dat$read_name==read)
    dat=dat[-ind,]
 }
    print(nrow(dat))
  
  dat=dat[dat$strand=="t" & dat$reference_kmer != "NNNNN" & dat$event_stdv<50,]
  dat[,"count"]=round(3012*dat$event_length)
  
  dat.name=unique(dat$read_name)
  dat.split=NROW(dat.name)/(no_cores+1)
  name.list= split(dat.name, trunc(0:(NROW(dat.name)-1 )/dat.split))
  dat=as.data.frame(dat)
  
  dat.com= foreach(c = 1:length(name.list), .combine=rbind, .errorhandling="pass") %dopar%
    {
      dat[dat$read_name %in% unlist(name.list[c]),] %>% group_by(contig,read_name,position) %>% 
        summarise(event_stdv=sd_d(event_stdv,event_level_mean,count),
                  event_level_mean=mean_c(event_level_mean,count), 
                  count=sum(count),
                  reference_kmer=unique(reference_kmer))
    }

  write_tsv(dat.com,opt$out,col_names = FALSE,append = FALSE, quote_escape = "none")
  if (length(ind)>0){
  eli=length(ind)
  write(length(ind), file="tmp.eli")
}

