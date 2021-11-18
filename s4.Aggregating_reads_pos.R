pacman::p_load(data.table,tidyverse,R.utils,plyr,doParallel,Rcpp,optparse)

rm(list=ls())
gc()

Rcpp::sourceCpp("/code/all_functions.cpp")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-t", "--thread"), type="integer", default=1,
              help="Number of cores allocated", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file name", metavar="character"),
  make_option(c("-r", "--regex"), type="character", default=NULL,
              help="output file name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list,add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$out)|is.null(opt$regex)){
  print_help(opt_parser)
  stop("Input files must be supplied (-o,-r).", call.=FALSE)
}

###Functions###
impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
entr=function(x){
  y=x*log2(x)
  return(y)
}

###Start###
files=list.files(pattern = paste0("_grp",opt$regex))
list.df=list()
for (i in 1:length(files)){
	list.df[[i]]=fread(files[i],nThread=opt$thread,header=T)
}
df=rbindlist(list.df) %>% distinct()
rm(list.df);gc()

df$REF=toupper(df$REF)
df$event_level_mean=as.numeric(df$event_level_mean)
df$event_stdv=as.numeric(df$event_stdv)
df$count=as.numeric(df$count)
print("Done loading")

######calculate entropy####
A.df=df[df$BASE=="A",][, .(cov.A = .N), by = posi_id]
G.df=df[df$BASE=="G",][, .(cov.G = .N), by = posi_id]
C.df=df[df$BASE=="C",][, .(cov.C = .N), by = posi_id]
T.df=df[df$BASE=="T",][, .(cov.T = .N), by = posi_id]
df.Total=df[, .(total = .N), by = posi_id]

dat.com=join_all(list(df.Total,A.df,G.df,T.df,C.df), by='posi_id', type='full')
dat.com[is.na(dat.com)==T]=0
dat.com$pct.A=dat.com$cov.A/dat.com$total
dat.com$pct.G=dat.com$cov.G/dat.com$total
dat.com$pct.C=dat.com$cov.C/dat.com$total
dat.com$pct.T=dat.com$cov.T/dat.com$total

dat.com$entr.A=entr(dat.com$pct.A)
dat.com$entr.G=entr(dat.com$pct.G)
dat.com$entr.C=entr(dat.com$pct.C)
dat.com$entr.T=entr(dat.com$pct.T)
dat.com$entropy=-rowSums(dat.com[,c(11:14)],na.rm = T)
dat.com=dat.com[(dat.com$total>=6),]
rm(A.df,C.df,G.df,T.df); gc()
print("Done entropy")
##########################

####aggregating####
cov=df[, .(n = .N), by = posi_id]
cov=cov[cov$n>=6,]
under.posi=cov[cov$n<10,]
over.posi1=(cov[cov$n>=20 & cov$n < 110,])
over.posi2=(cov[cov$n>=110,])

dat.under=df[df$posi_id %in% under.posi$posi_id,]
dat.under2=ddply(dat.under,.(posi_id),function(x) x[sample(nrow(x),(10-nrow(x))),])
rm(dat.under);gc()
print("Done under")

dat.over=df[df$posi_id %in% c(over.posi1$posi_id,over.posi2$posi_id),]
dat.over=dat.over[sample(nrow(dat.over)),]
dat.over2=left_join(dat.over[dat.over$posi_id %in% over.posi1$posi_id,],over.posi1,by="posi_id")
dat.over3=left_join(dat.over[dat.over$posi_id %in% over.posi2$posi_id,],over.posi2,by="posi_id")
rm(dat.over);gc()
dat.over2$index=0
dat.over3$index=0
setDT(dat.over2)[,c("index"):=.(indexgroup10(index,unique(n))),by=.(posi_id)]
setDT(dat.over3)[,c("index"):=.(indexgroup_big(index,unique(n))),by=.(posi_id)]
print("Done over")

df=df[df$posi_id %in% c(over.posi1$posi_id,over.posi2$posi_id) == F,]
df=df[df$posi_id %in% cov$posi_id==T,]
df=bind_rows(df,dat.under2)

dat.split1 <- uniqueN(df$posi_id)/(20)
name.list1 <- split(unique(df$posi_id), trunc(0:(NROW(unique(df$posi_id))-1 )/dat.split1))
df1 <- foreach(b=1:length(name.list1),.combine=rbind, .errorhandling="pass") %do%
{
  df[df$posi_id %in% unlist(name.list1[b]),][,c("event_level_mean","event_stdv","count"):=.(impute.mean(event_level_mean),impute.mean(event_stdv),impute.mean(count)),
											 by=.(posi_id)][,by=.(contig,position,strand,REF,posi_id),
															.(reference_kmer=unique(na.omit(reference_kmer)),
															  cov=.N,
															  event_stdv.m=mean(event_stdv,na.rm=T),
															  event_level_mean.m=mean(event_level_mean,na.rm=T),
															  event_stdv.s=sd(event_stdv,na.rm=T),
															  event_level_mean.s=sd(event_level_mean,na.rm=T),
															  event_stdv.A=sd_d(event_stdv,event_level_mean,count),
															  event_level_mean.A=mean_c(event_level_mean,count),
															  count.m=mean(count,na.rm=T),
															  qual1.m=mean(qual1,na.rm=T),
															  qual2.m=mean(qual2,na.rm=T),
															  qual3.m=mean(qual,na.rm=T),
															  qual4.m=mean(qual4,na.rm=T),
															  qual5.m=mean(qual5,na.rm=T),
															  del1.m=mean(del1,na.rm=T),
															  del2.m=mean(del2,na.rm=T),
															  del3.m=mean(del,na.rm=T),
															  del4.m=mean(del4,na.rm=T),
															  del5.m=mean(del5,na.rm=T),
															  mis1.m=mean(mis1,na.rm=T),
															  mis2.m=mean(mis2,na.rm=T),
															  mis3.m=mean(mis,na.rm=T),
															  mis4.m=mean(mis4,na.rm=T),
															  mis5.m=mean(mis5,na.rm=T),
															  ins1.m=mean(ins1,na.rm=T),
															  ins2.m=mean(ins2,na.rm=T),
															  ins3.m=mean(ins,na.rm=T),
															  ins4.m=mean(ins4,na.rm=T),
															  ins5.m=mean(ins5,na.rm=T),
															  count.s=sd(count,na.rm=T),
															  qual1.s=sd(qual1,na.rm=T),
															  qual2.s=sd(qual2,na.rm=T),
															  qual3.s=sd(qual,na.rm=T),
															  qual4.s=sd(qual4,na.rm=T),
															  qual5.s=sd(qual5,na.rm=T),
															  del1.s=sd(del1,na.rm=T),
															  del2.s=sd(del2,na.rm=T),
															  del3.s=sd(del,na.rm=T),
															  del4.s=sd(del4,na.rm=T),
															  del5.s=sd(del5,na.rm=T),
															  mis1.s=sd(mis1,na.rm=T),
															  mis2.s=sd(mis2,na.rm=T),
															  mis3.s=sd(mis,na.rm=T),
															  mis4.s=sd(mis4,na.rm=T),
															  mis5.s=sd(mis5,na.rm=T),
															  ins1.s=sd(ins1,na.rm=T),
															  ins2.s=sd(ins2,na.rm=T),
															  ins3.s=sd(ins,na.rm=T),
															  ins4.s=sd(ins4,na.rm=T),
															  ins5.s=sd(ins5,na.rm=T),
																							event_stdv.ms=mean(event_stdv**2,na.rm=T),
																							event_level_mean.ms=mean(event_level_mean**2,na.rm=T),
																							count.ms=mean(count**2,na.rm=T),
																							qual1.ms=mean(qual1**2,na.rm=T),
																							qual2.ms=mean(qual2**2,na.rm=T),
																							qual3.ms=mean(qual**2,na.rm=T),
																							qual4.ms=mean(qual4**2,na.rm=T),
																							qual5.ms=mean(qual5**2,na.rm=T),
																							event_stdv.rms=sqrt(mean(event_stdv**2,na.rm=T)),
																							event_level_mean.rms=sqrt(mean(event_level_mean**2,na.rm=T)),
																							count.rms=sqrt(mean(count**2,na.rm=T)),
																							qual1.rms=sqrt(mean(qual1**2,na.rm=T)),
																							qual2.rms=sqrt(mean(qual2**2,na.rm=T)),
																							qual3.rms=sqrt(mean(qual**2,na.rm=T)),
																							qual4.rms=sqrt(mean(qual4**2,na.rm=T)),
																							qual5.rms=sqrt(mean(qual5**2,na.rm=T)))]
  
}
rm(df);gc()
print("done aggregate undersampling")

dat.over=bind_rows(dat.over2,dat.over3)
rm(dat.over2,dat.over3)
dat.split2 <- uniqueN(dat.over$posi_id)/(20)
name.list2 <- split(unique(dat.over$posi_id), trunc(0:(NROW(unique(dat.over$posi_id))-1 )/dat.split2))
df2 <- foreach(b=1:length(name.list2),.combine=rbind, .errorhandling="pass") %do%
{
  dat.over[dat.over$posi_id %in% unlist(name.list2[b]),][,c("event_level_mean","event_stdv","count"):=.(impute.mean(event_level_mean),impute.mean(event_stdv),impute.mean(count)),
														   by=.(posi_id)][,by=.(contig,position,strand,REF,posi_id,index),
																		  .(reference_kmer=unique(na.omit(reference_kmer)),
																			cov=.N,
																			event_stdv.m=mean(event_stdv,na.rm=T),
																			event_level_mean.m=mean(event_level_mean,na.rm=T),
																			event_stdv.s=sd(event_stdv,na.rm=T),
																			event_level_mean.s=sd(event_level_mean,na.rm=T),
																			event_stdv.A=sd_d(event_stdv,event_level_mean,count),
																			event_level_mean.A=mean_c(event_level_mean,count),
																			count.m=mean(count,na.rm=T),
																			qual1.m=mean(qual1,na.rm=T),
																			qual2.m=mean(qual2,na.rm=T),
																			qual3.m=mean(qual,na.rm=T),
																			qual4.m=mean(qual4,na.rm=T),
																			qual5.m=mean(qual5,na.rm=T),
																			del1.m=mean(del1,na.rm=T),
																			del2.m=mean(del2,na.rm=T),
																			del3.m=mean(del,na.rm=T),
																			del4.m=mean(del4,na.rm=T),
																			del5.m=mean(del5,na.rm=T),
																			mis1.m=mean(mis1,na.rm=T),
																			mis2.m=mean(mis2,na.rm=T),
																			mis3.m=mean(mis,na.rm=T),
																			mis4.m=mean(mis4,na.rm=T),
																			mis5.m=mean(mis5,na.rm=T),
																			ins1.m=mean(ins1,na.rm=T),
																			ins2.m=mean(ins2,na.rm=T),
																			ins3.m=mean(ins,na.rm=T),
																			ins4.m=mean(ins4,na.rm=T),
																			ins5.m=mean(ins5,na.rm=T),
																			count.s=sd(count,na.rm=T),
																			qual1.s=sd(qual1,na.rm=T),
																			qual2.s=sd(qual2,na.rm=T),
																			qual3.s=sd(qual,na.rm=T),
																			qual4.s=sd(qual4,na.rm=T),
																			qual5.s=sd(qual5,na.rm=T),
																			del1.s=sd(del1,na.rm=T),
																			del2.s=sd(del2,na.rm=T),
																			del3.s=sd(del,na.rm=T),
																			del4.s=sd(del4,na.rm=T),
																			del5.s=sd(del5,na.rm=T),
																			mis1.s=sd(mis1,na.rm=T),
																			mis2.s=sd(mis2,na.rm=T),
																			mis3.s=sd(mis,na.rm=T),
																			mis4.s=sd(mis4,na.rm=T),
																			mis5.s=sd(mis5,na.rm=T),
																			ins1.s=sd(ins1,na.rm=T),
																			ins2.s=sd(ins2,na.rm=T),
																			ins3.s=sd(ins,na.rm=T),
																			ins4.s=sd(ins4,na.rm=T),
																			ins5.s=sd(ins5,na.rm=T),
																			event_stdv.ms=mean(event_stdv**2,na.rm=T),
																			event_level_mean.ms=mean(event_level_mean**2,na.rm=T),
																			count.ms=mean(count**2,na.rm=T),
																			qual1.ms=mean(qual1**2,na.rm=T),
																			qual2.ms=mean(qual2**2,na.rm=T),
																			qual3.ms=mean(qual**2,na.rm=T),
																			qual4.ms=mean(qual4**2,na.rm=T),
																			qual5.ms=mean(qual5**2,na.rm=T),
																			event_stdv.rms=sqrt(mean(event_stdv**2,na.rm=T)),
																			event_level_mean.rms=sqrt(mean(event_level_mean**2,na.rm=T)),
																			count.rms=sqrt(mean(count**2,na.rm=T)),
																			qual1.rms=sqrt(mean(qual1**2,na.rm=T)),
																			qual2.rms=sqrt(mean(qual2**2,na.rm=T)),
																			qual3.rms=sqrt(mean(qual**2,na.rm=T)),
																			qual4.rms=sqrt(mean(qual4**2,na.rm=T)),
																			qual5.rms=sqrt(mean(qual5**2,na.rm=T)))]
}
rm(dat.over);gc()
df2$index = NULL
print("Done aggregating")
df.com=rbindlist(list(df1,df2))
rm(df1,df2);gc()
df.com$cov=NULL
names(dat.com)[2]="cov"

all=df.com[dat.com[,c(1,2,7:10,15)],on="posi_id",nomatch = 0]
all=all %>% relocate(cov,.after = reference_kmer)
all=all[,c(1:55,72:76,56:71)]
all$ms.pct=rowSums(all[,c(56:59)],na.rm = T)/4
all$rms.pct=sqrt(rowSums(all[,c(56:59)],na.rm = T)/4)
all$ms.mis=rowSums(all[,c(25:29)],na.rm = T)/5
all$rms.mis=sqrt(rowSums(all[,c(25:29)],na.rm = T)/5)
all$ms.del=rowSums(all[,c(20:24)],na.rm = T)/5
all$rms.del=sqrt(rowSums(all[,c(20:24)],na.rm = T)/5)
all=all[all$reference_kmer != "",]
				
print("Start writing")
fwrite(all,file=opt$out,quote=F,sep="\t")
print("Done")
