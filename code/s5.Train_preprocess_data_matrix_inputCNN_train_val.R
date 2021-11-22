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
df=fread(file=opt$input, nThread=opt$thread) %>% distinct()%>% filter(reference_kmer !="")		

truth= fread(opt$classref, nThread=opt$thread)###truth needs to have 5 columns: contig, position,strand, edit, rate	

chroms=sort(unique(truth$contig)) 
chromid=as.character(seq(length(chroms)))
chrmapping <- data.frame(chroms, chromid); saveRDS(chrmapping, file = "chrmapping.rds")

truth$posi_id=paste(truth$contig,truth$position,truth$strand,sep=":")
truth$REF="A"
truth2=decode(truth)
truth2$id=paste0(truth2$chr.str,":",truth2$position)
truth=bind_cols(truth,truth2[,"id"])
rm(truth2);gc()
pos=unique(truth$posi_id)
posab=generate_before_after(pos)
truth=truth[order(truth$id),]
cut=nrow(truth)*0.9
test.pos=truth[round(cut):nrow(truth),]
train.pos=setdiff(truth,test.pos)

dat=df[df$posi_id %in% unique(c(pos,posab$ld1,posab$lg1,posab$ld2,posab$lg2)),] %>% distinct()
rm(df);gc()
dat=replace0(dat)
dat[,c("event_stdv.m","event_stdv.s","event_stdv.A","count.m","count.s","ins3.m","ins1.m","ins2.m","ins4.m","ins5.m","ins3.s","ins1.s","ins2.s","ins4.s","ins5.s")]=
   lapply(dat[,c("event_stdv.m","event_stdv.s","event_stdv.A","count.m","count.s","ins3.m","ins1.m","ins2.m","ins4.m","ins5.m","ins3.s","ins1.s","ins2.s","ins4.s","ins5.s")],log)
dat=dat[,-c(30:34,51:55)]
dat=dat[complete.cases(dat),]

train.ab=generate_before_after(train.pos$posi_id)
train.dat=dat[dat$posi_id %in% c(train.pos$posi_id,train.ab$ld1,train.ab$lg1,train.ab$ld2,train.ab$lg2),]
test.ab=generate_before_after(test.pos$posi_id)
test.dat=dat[dat$posi_id %in% c(test.pos$posi_id,test.ab$ld1,test.ab$lg1,test.ab$ld2,test.ab$lg2),]

col=colnames(train.dat)[8:72]
max=train.dat %>% dplyr::summarise_each_(funs(max(.,na.rm = TRUE)), col)
max=unlist(unname((max)))
min=train.dat %>% dplyr::summarise_each_(funs(min(.,na.rm = TRUE)), col)
min=unlist(unname((min)))
mm.df=data.frame(col,max,min)
mm.df$diff=mm.df$max-mm.df$min

train.dat=as.data.frame(train.dat)
test.dat=as.data.frame(test.dat)
for (i in 1:length(col)){
  test.dat[,col[i]]=((test.dat[,col[i]]-mm.df[mm.df$col==col[i],"min"])*100)/(mm.df[mm.df$col==col[i],"diff"])
  train.dat[,col[i]]=((train.dat[,col[i]]-mm.df[mm.df$col==col[i],"min"])*100)/(mm.df[mm.df$col==col[i],"diff"])
}

checktrain=check_mono(train.dat,72)
checktest=check_mono(test.dat,72)

train.s=transform_dat(checktrain,train.dat)
train.r=transform_dat_big(checktrain,train.dat)
train.b=transform_dat_bb(checktrain,train.dat)

test.s=transform_dat(checktest,test.dat)
test.r=transform_dat_big(checktest,test.dat)
test.b=transform_dat_bb(checktest,test.dat)

train.tmp=list()
if(exists("train.r")==F){
  if(exists("train.b")==F){
    train.tmp=train.s
  } else {
    train.tmp$x=abind(train.s$x,train.b$x,along = 1)
    train.tmp$info=abind(train.s$info,train.b$info,along = 3)
  }
} else{
  if(exists("train.b")==F){
    train.tmp$x=abind(train.s$x,train.r$x,along = 1)
    train.tmp$info=abind(train.s$info,train.r$info,along = 3)
  } else {
    train.tmp$x=abind(train.s$x,train.r$x,train.b$x,along = 1)
    train.tmp$info=abind(train.s$info,train.r$info,train.b$info,along = 3)
  }
}  

ineli=which(is.na(train.tmp$x))
if(length(ineli)>0){
	train.tmp$x=train.tmp$x[-c(ineli),,,]
	train.tmp$info=train.tmp$info[,,-c(ineli)]
	train.tmp$x=array_reshape(train.tmp$x,dim=c(nrow(train.tmp$x),5,65,1))
 }
 
test.tmp=list()
if(exists("test.r")==F){
  if(exists("test.b")==F){
    test.tmp=test.s
  } else {
    test.tmp$x=abind(test.s$x,test.b$x,along = 1)
    test.tmp$info=abind(test.s$info,test.b$info,along = 3)
  }
} else{
  if(exists("test.b")==F){
    test.tmp$x=abind(test.s$x,test.r$x,along = 1)
    test.tmp$info=abind(test.s$info,test.r$info,along = 3)
  } else {
    test.tmp$x=abind(test.s$x,test.r$x,test.b$x,along = 1)
    test.tmp$info=abind(test.s$info,test.r$info,test.b$info,along = 3)
  }
}

###train matrix
inf=as.data.frame(t(train.tmp$info[3,,]))
colnames(inf)=c("chr.str","pos","REFbase","kmer","cov")
inf$ind=(1:nrow(inf))
inf$id=paste0(inf$chr.str,":",inf$pos)

df1=left_join(inf,train.pos,by="id")
df1=df1[complete.cases(df1),]
df1=df1[order(df1$ind),]

train.matrix=list()
train.matrix$y=as.vector(df1$edit)
train.matrix$x=train.tmp$x[df1$ind,,,]
train.matrix$x=array_reshape(train.matrix$x,dim=c(length(train.matrix$y),5,65,1))
train.matrix$info=as.matrix(df1[,c("chr.str","pos","REFbase","kmer","cov","ind","id")])
train.matrix$y2=as.vector(df1$rate)

###test matrix
inf=as.data.frame(t(test.tmp$info[3,,]))
colnames(inf)=c("chr.str","pos","REFbase","kmer","cov")
inf$ind=(1:nrow(inf))
inf$id=paste0(inf$chr.str,":",inf$pos)

df1=left_join(inf,test.pos,by="id")
df1=df1[complete.cases(df1),]
df1=df1[order(df1$ind),]

test.matrix=list()
test.matrix$y=as.vector(df1$edit)
test.matrix$x=test.tmp$x[df1$ind,,,]
test.matrix$x=array_reshape(test.matrix$x,dim=c(length(test.matrix$y),5,65,1))
test.matrix$info=as.matrix(df1[,c("chr.str","pos","REFbase","kmer","cov","ind","id")])
test.matrix$y2=as.vector(df1$rate)

saveRDS(test.matrix, file = paste0(opt$out,".validation_matrix.rds"))
saveRDS(train.matrix, file = paste0(opt$out,".training_matrix.rds"))
print("Done")
