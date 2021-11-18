pacman::p_load(data.table,tidyverse,caret,stringr,keras,tensorflow, optparse, multiROC)

rm(list=ls())
gc()
Rcpp::sourceCpp("/code/all_functions.cpp")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input file name", metavar="character"),
  make_option(c("-t", "--thread"), type="integer", default=1,
              help="Number of cores allocated", metavar="integer"),
  make_option(c("-c", "--classref"), type="character", default=NULL,
              help="ground truth of class for coordinates", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list,add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input files must be supplied (-i).", call.=FALSE)
}


chrmapping <- readRDS("chrmapping.rds")
# replacement=chrmapping$chroms
model=load_model_hdf5("/data/model/best_pos5_mix_3class_resnet_1992.h5")
model1=load_model_hdf5("/data/model/best_pos5_mix_3c_1vs1_3010_resnet_10.h5")

load(opt$input)

agggrp <- word(opt$input,sep = "\\.", 1)
  
k_set_learning_phase(0)
pred <- predict(model,test$x[,,c(1:43),])
print("3 class")
pred=as.data.frame(pred)
pred_test=colnames(pred[,c(1:3)])[apply(pred[,c(1:3)],1,which.max)]
pred_test=gsub("V1",0,pred_test)
pred_test=gsub("V2",1,pred_test)
pred_test=gsub("V3",2,pred_test)
pred_test=as.numeric(pred_test)
c01=which(pred_test!=2)
pred$n2=1-pred$V3
pred$c1=0
pred$c2=0
  
pred2 <- predict(model1,test$x[c01,,c(1:43),])
print("2class")
pred2=as.data.frame(pred2)
  
pred[c01,c(5,6)]=pred2
pred[c01,"V1"]=pred[c01,"c1"]*pred[c01,"n2"]
pred[c01,"V2"]=pred[c01,"c2"]*pred[c01,"n2"]
  
# pre=data.frame(pred[,c(1:3)],'ref'=test$y,t(test$info[3,c(1,2),]))
pre=data.frame(pred[,c(1:3)],'id'=test$info[,"id"],'cov'=test$info[,"cov"])

# pre2=pre %>% group_by(id,ref) %>% dplyr::summarise(n=n(),'0'=mean(V1),'1'=mean(V2),'2'=mean(V3))
pre2=pre %>% group_by(id) %>% dplyr::summarise(cov=unique(cov),'0'=mean(V1),'1'=mean(V2),'2'=mean(V3))

# pre2$pred=colnames(pre2[,c(4:6)])[apply(pre2[,c(4:6)],1,which.max)]
pre2$pred=colnames(pre2[,c(3:5)])[apply(pre2[,c(3:5)],1,which.max)]

truth=fread(opt$classref, nThread=opt$thread)
pre2=left_join(pre2,truth,by="id")
pre2$cov=as.numeric(pre2$cov)

pre3=pre2 %>% separate(id,into=c("chr_str","position"),sep=":")
pre3$strand=str_sub(pre3$chr_str,-1)
pre3$strand[pre3$strand==1]="+"
pre3$strand[pre3$strand==2]="-"
  
pre3$chr=str_sub(pre3$chr_str,1,-2)
pre3$contig=replace_list(chrmapping$chromid,chrmapping$chroms,pre3$chr)
class0 <- pre3[pre3$pred == 0,c("contig","position","strand","cov","0","1","2","pred","ref")]
class1 <- pre3[pre3$pred == 1,c("contig","position","strand","cov","0","1","2","pred","ref")]
class2 <- pre3[pre3$pred == 2,c("contig","position","strand","cov","0","1","2","pred","ref")]
  
fwrite(class0, paste0(agggrp, ".output_prediction_CNN_class0.txt"),sep = "\t")
fwrite(class1, paste0(agggrp, ".output_prediction_CNN_class1.txt"),sep = "\t")
fwrite(class2, paste0(agggrp, ".output_prediction_CNN_class2.txt"),sep = "\t")


######ROC PR curve#####

y_ref=to_categorical(pre3$ref,3)
perform.df <- data.frame(y_ref,pre3[,c(4:6)])
colnames(perform.df)=c("S0_true","S1_true","S2_true","S0_pred_n","S1_pred_n","S2_pred_n")

roc_df <- plot_roc_data(multi_roc(perform.df))
roc_df$FPR <- 1-roc_df$Specificity
names(roc_df)[2] <- "TPR"
auroc <- roc_df %>% group_by(Group) %>% summarise(Group=unique(Group),AUC=unique(AUC))
print(auroc)

mic.roc=roc_df[roc_df$Group=="Micro",]
mic.roc[nrow(mic.roc) + 1,c("TPR","FPR","AUC")] = c(0,0,unique(mic.roc$AUC))

roc=ggplot(data=mic.roc,aes(x=FPR,y=TPR))+
  geom_line(size=1)+theme_bw()+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1, alpha=0.5, linetype = "dashed")+
  xlab("False Positive Rate")+
  ylab("True Positive Rate")+
  annotate("text", label=paste("AUC =",round(unique(mic.roc$AUC),3)), 
           x=0.36, y=0, hjust=0, size=10, vjust=0,color="blue")+
  theme(axis.text.y =element_text(size=20, color="black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x =element_text(size=20, color="black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title=element_text(size=25,color="black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.ticks.length = unit(1, "line"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        plot.title = element_text(size=25,hjust=0.5))
ggsave(file=paste0(agggrp, ".ROC_curve.png"), plot=roc, width=6, height=6,units = "in")


pr_df=plot_pr_data(multi_pr(perform.df))
names(pr_df)[c(1,2)]=c("recall","precision")
aupr=pr_df %>% group_by(Group) %>% summarise(Group=unique(Group),AUC=unique(AUC))
print(aupr)

mic.pr=pr_df[pr_df$Group=="Micro",]
# mic.pr[nrow(mic.pr) + 1,c("precision","recall","AUC")] = c(0,1,unique(mic.pr$AUC))

pr=ggplot(data=mic.pr,aes(x=recall,y=precision))+
  geom_line(size=1)+theme_bw()+
  ylim(0,1)+
  xlab("Recall")+
  ylab("Precision")+
  annotate("text", label=paste("AUC =",round(unique(mic.pr$AUC),3)), 
           x=0, y=0, hjust=0, size=10, vjust=0,color="blue")+
  theme(axis.text.y =element_text(size=20, color="black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x =element_text(size=20, color="black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title=element_text(size=25,color="black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.key.width = unit(5, "line"),
        axis.ticks.length = unit(1, "line"),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        plot.title = element_text(size=25,hjust=0.5))
ggsave(file=paste0(agggrp, ".PR_curve.png"), plot=pr, width=6, height=6,units = "in")

