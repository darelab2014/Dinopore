pacman::p_load(data.table,tidyverse,Matrix,caret,Rcpp,usefun,scales,keras,abind,plyr,optparse)

rm(list=ls())
gc()

options(scipen=999)

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-t", "--thread"), type="integer", default=1,
              help="Number of cores allocated", metavar="integer"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file name", metavar="character"),
  make_option(c("-c", "--classref"), type="character", default=NULL,
              help="ground truth of class for coordinates", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list,add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$input)|is.null(opt$out)|is.null(opt$classref)){
  print_help(opt_parser)
  stop("Input files must be supplied (-i,-o, -c).", call.=FALSE)
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

load(paste0(dinodir, "/code/misc/minmax_Xen_5pos_012_10bin_morefts_combine.RData"))

generate_before_after <- function(x){
  x=as.data.frame(x)
  x=x %>% separate(x,c("chr","position","strand"),sep=":")
  x$po.le1=as.numeric(x$position)+1
  x$po.le2=as.numeric(x$position)+2
  x$po.la1=as.numeric(x$position)-1
  x$po.la2=as.numeric(x$position)-2
  x$po.le.id1=paste(x$chr,x$po.le1,x$strand,sep = ":")
  x$po.la.id1=paste(x$chr,x$po.la1,x$strand,sep = ":")
  x$po.le.id2=paste(x$chr,x$po.le2,x$strand,sep = ":")
  x$po.la.id2=paste(x$chr,x$po.la2,x$strand,sep = ":")
  y=list('ld1'=x$po.le.id1, 'lg1'=x$po.la.id1,'ld2'=x$po.le.id2, 'lg2'=x$po.la.id2)
  return(y)
}


replace0=function(x){
  x$event_stdv.m[x$event_stdv.m==0] <- 0.00001
  x$event_stdv.s[x$event_stdv.s==0] <- 0.00001
  x$event_stdv.A[x$event_stdv.A==0] <- 0.00001
  x$count.m[x$count.m==0] <- 0.00001
  x$count.s[x$count.s==0] <- 0.00001
  x$ins3.m[x$ins3.m==0] <- 0.00001
  x$ins1.m[x$ins1.m==0] <- 0.00001
  x$ins2.m[x$ins2.m==0] <- 0.00001
  x$ins4.m[x$ins4.m==0] <- 0.00001
  x$ins5.m[x$ins5.m==0] <- 0.00001
  x$ins3.s[x$ins3.s==0] <- 0.00001
  x$ins1.s[x$ins1.s==0] <- 0.00001
  x$ins2.s[x$ins2.s==0] <- 0.00001
  x$ins4.s[x$ins4.s==0] <- 0.00001
  x$ins5.s[x$ins5.s==0] <- 0.00001
  return(x)
}

decode=function(x){
  setDT(x)
  x$contig=replace_list(chroms,chromid,x$contig)
  x$contig=as.numeric(x$contig)
  
  x$reference_kmer=gsub("A", "1", x$reference_kmer) 
  x$reference_kmer=gsub("C", "2", x$reference_kmer) 
  x$reference_kmer=gsub("G", "3", x$reference_kmer) 
  x$reference_kmer=gsub("T", "4", x$reference_kmer) 
  x$reference_kmer=as.numeric(x$reference_kmer)
 
  patterns = c("A","C","G","T")
  replacement = as.character(c(1:4))
  x$REF=replace_list(patterns,replacement,x$REF)
  x$REF=as.numeric(x$REF)
  
  x$strand=gsub("p", "1", x$strand) 
  x$strand=gsub("n", "2", x$strand) 
  x$strand=as.numeric(x$strand)
  
  x$chr.str=paste0(x$contig,x$strand)
  x=x %>% select(chr.str,everything())
  x$chr.str=as.numeric(x$chr.str)
  x=x[,-c("contig","strand","posi_id")]
  
  return(x)
}

check_mono <- function(x,n){
  xt=x %>% select(c(1:5,n))
  xt=xt[order(xt$strand,xt$contig,xt$position),]
  p=x[((x$REF=="A" & x$strand=="p")| (x$REF=="T" & x$strand=="n")),]$posi_id
  xt$cld1=0
  xt$clg1=0
  xt$cld2=0
  xt$clg2=0
  countnb5(xt$REF,xt$strand,xt$position,xt$cld1,xt$clg1,xt$cld2,xt$clg2)
  check=xt %>% group_by(contig,position,strand,posi_id,REF) %>% dplyr::summarise(ld1=max(cld1),lg1=max(clg1),ld2=max(cld2),lg2=max(clg2),n=n())
  check=check[((check$REF=="A" & check$strand=="p") | (check$REF=="T" & check$strand=="n")),]
  check=check[check$posi_id %in% p,]
  check$t=check$ld1*check$lg1*check$n*check$ld2*check$lg2
  return(check)
}


transform_dat_big <- function(checkdf,df, ncols = 70){
  set.seed(1999)
  posm=checkdf[checkdf$t>10 & checkdf$n<10,]$posi_id
  dfm=df[df$posi_id %in% posm & ((df$REF=="A" & df$strand=="p")| (df$REF=="T" & df$strand=="n")),]
  dfm2=decode(dfm)[,c(1,2)] %>% distinct()
  dfm2=as.matrix(dfm2)
  df1=decode(df)
  df1=as.matrix(df1)
  
  lt.mt=list()
  for (i in 1:nrow(dfm2)){
    lt.mt[[i]]=transform_each_pos5_regular(dfm2[i,],df1, ncols)
  }

  lt.mt1=do.call(rbind, lt.mt)
  lt.mt2 <-  aperm(`dim<-`(t(lt.mt1), c(ncols,5, nrow(lt.mt1)/5)), c(2, 1,3))
  x_train=lt.mt2[,6:ncols,]
  x_train=aperm(x_train,c(3,1,2))
  x_train=array_reshape(x_train,dim=c(nrow(lt.mt1)/5,5,65,1))
  info_train=lt.mt2[,1:5,]
  
  return(list('x'=x_train, 'info'=info_train))
}

transform_dat <- function(checkdf,df, ncols = 70){
  posm=checkdf[checkdf$t>0 & checkdf$t <= 10,]$posi_id
  dfm=df[df$posi_id %in% posm & ((df$REF=="A" & df$strand=="p")| (df$REF=="T" & df$strand=="n")),]
  dfm=dfm[sample(nrow(dfm)),]
  dfm=decode(dfm)
  dfm=as.matrix(dfm)
  df1=decode(df)
  df1=as.matrix(df1)
  
  
  lt.mt=list()
  for (i in 1:nrow(dfm)){
    lt.mt[[i]]=transform_each_pos5(dfm[i,],df1, ncols)
  }
  
  lt.mt1=do.call(rbind, lt.mt)
  lt.mt2 <-  aperm(`dim<-`(t(lt.mt1), c(ncols,5, nrow(lt.mt1)/5)), c(2, 1,3))
  x_train=lt.mt2[,6:ncols,]
  x_train=aperm(x_train,c(3,1,2))
  x_train=array_reshape(x_train,dim=c(nrow(lt.mt1)/5,5,65,1))
  info_train=lt.mt2[,1:5,]
  
  return(list('x'=x_train, 'info'=info_train))
}


transform_dat_bb <- function(checkdf,df, ncols = 70){
  posm=checkdf[checkdf$t>10 & checkdf$n==10,]$posi_id
  dfm=df[df$posi_id %in% posm & ((df$REF=="A" & df$strand=="p")| (df$REF=="T" & df$strand=="n")),]
  dfm2=decode(dfm)[,c(1,2)] %>% distinct()
  dfm2=as.matrix(dfm2)
  set.seed(1999)
  df1=decode(df)
  df1=df1[sample(nrow(df1)),]
  df1=as.matrix(df1)
  
  lt.mt=list()
  set.seed(1999)
  for (i in 1:nrow(dfm2)){
    lt.mt[[i]]=transform_each_pos5_big(dfm2[i,],df1, ncols)
  }

  lt.mt1=do.call(rbind, lt.mt)
  lt.mt2 <-  aperm(`dim<-`(t(lt.mt1), c(ncols,5, nrow(lt.mt1)/5)), c(2, 1,3))
  x_train=lt.mt2[,6:ncols,]
  x_train=aperm(x_train,c(3,1,2))
  x_train=array_reshape(x_train,dim=c(nrow(lt.mt1)/5,5,65,1))
  info_train=lt.mt2[,1:5,]
  
  return(list('x'=x_train, 'info'=info_train))
}

print(opt$input)
file.df=fread(file=opt$input, nThread=opt$thread) %>% distinct()

								
dat=file.df %>% distinct()%>% filter(reference_kmer !="")
dat=replace0(dat)
dat[,c("event_stdv.m","event_stdv.s","event_stdv.A","count.m","count.s","ins3.m","ins1.m","ins2.m","ins4.m","ins5.m","ins3.s","ins1.s","ins2.s","ins4.s","ins5.s")]=
   lapply(dat[,c("event_stdv.m","event_stdv.s","event_stdv.A","count.m","count.s","ins3.m","ins1.m","ins2.m","ins4.m","ins5.m","ins3.s","ins1.s","ins2.s","ins4.s","ins5.s")],log)
dat=dat[,-c(30:34,51:55)]
dat=dat[complete.cases(dat),]
dat=as.data.frame(dat)

col=colnames(dat)[8:72]
for (i in 1:length(col)){
  dat[,col[i]]=((dat[,col[i]]-mm.df[mm.df$col==col[i],"min"])*100)/(mm.df[mm.df$col==col[i],"diff"])
}

checkdf=check_mono(dat,72)
chroms=sort(unique(checkdf$contig)) 
chromid=as.character(seq(length(chroms)))
chrmapping <- data.frame(chroms, chromid); saveRDS(chrmapping, file = "chrmapping.rds")

test.s=transform_dat(checkdf,dat)
test.r=transform_dat_big(checkdf,dat)
test.b=transform_dat_bb(checkdf,dat)
  
test.tmp=list()
if(exists("test.r")==F){
  if(exists("test.b")==F){
    test.tmp=test.s
  } else {
    test.tmp$x=abind(test.s$x,test.b$x,along = 1)
  test.tmp$y=abind(test.s$y,test.b$y,along = 1)
    test.tmp$info=abind(test.s$info,test.b$info,along = 3)
  }
} else{
  if(exists("test.b")==F){
    test.tmp$x=abind(test.s$x,test.r$x,along = 1)
  test.tmp$y=abind(test.s$y,test.r$y,along = 1)
    test.tmp$info=abind(test.s$info,test.r$info,along = 3)
  } else {
    test.tmp$x=abind(test.s$x,test.r$x,test.b$x,along = 1)
  test.tmp$y=abind(test.s$y,test.r$y,test.b$y,along = 1)
    test.tmp$info=abind(test.s$info,test.r$info,test.b$info,along = 3)
  }
}

inf=as.data.frame(t(test.tmp$info[3,,]))
colnames(inf)=c("chr.str","pos","REFbase","kmer","cov")
inf$ind=(1:nrow(inf))
inf$id=paste0(inf$chr.str,":",inf$pos)

df1 <- fread(opt$classref, nThread=opt$thread)

df2=left_join(inf,df1,by="id")
df2=df2[complete.cases(df2),]
df2=df2[order(df2$ind),]

test=list()
test$y=as.vector(df2$ref)
test$x=test.tmp$x[df2$ind,,,]
test$x=array_reshape(test$x,dim=c(length(test$y),5,65,1))
test$info=as.matrix(df2[,c(1:7)])
test$y2=as.vector(df2$rate)

save(test,file=opt$out)
print("Done")
