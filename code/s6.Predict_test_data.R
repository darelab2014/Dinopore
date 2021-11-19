pacman::p_load(data.table,tidyverse,caret,stringr,keras,tensorflow, optparse, multiROC)

rm(list=ls())
gc()
Rcpp::sourceCpp("/code/all_functions.cpp")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input file name", metavar="character"),
  make_option(c("-t", "--thread"), type="integer", default=1,
              help="Number of cores allocated", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list,add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input file must be supplied (-i).", call.=FALSE)
}

chrmapping <- readRDS("chrmapping.rds")

model1=load_model_hdf5("/data/model/best_pos5_mix_3class_resnet_1992.h5")
model2=load_model_hdf5("/data/model/best_pos5_mix_3c_1vs1_3010_resnet_10.h5")
model3=load_model_hdf5("/data/model/best_regression_mixnHEK_morefts_16384_1024_b1024_init_650k_XHe0_Mus_asin06.h5")

load(opt$input)

agggrp <- word(opt$input,sep = "\\.", 1)
  
k_set_learning_phase(0)

pred <- predict(model1,test$x[,,c(1:43),])
pred <- as.data.frame(pred)

pred_test <- colnames(pred[,c(1:3)])[apply(pred[,c(1:3)],1,which.max)]
pred_test <- gsub("V1",0,pred_test)
pred_test <- gsub("V2",1,pred_test)
pred_test <- gsub("V3",2,pred_test)
pred_test <- as.numeric(pred_test)

c01=which(pred_test!=2)
pred$n2=1-pred$V3
pred$c1=0
pred$c2=0
  
pred2 <- predict(model2,test$x[c01,,c(1:43),])
pred2 <- as.data.frame(pred2)
  
pred[c01,c(5,6)] <- pred2
pred[c01,"V1"] <- pred[c01,"c1"]*pred[c01,"n2"]
pred[c01,"V2"] <- pred[c01,"c2"]*pred[c01,"n2"]

pred.2model <- data.frame(pred[,c(1:3)], 'ref'=test$y, 'id'=test$info[,"id"], 'cov'=test$info[,"cov"]) %>% 
  group_by(id, cov, ref) %>% 
  dplyr::summarise('0'=mean(V1),'1'=mean(V2),'2'=mean(V3))

pred.2model$pred <- colnames(pred.2model[,c(4:6)])[apply(pred.2model[,c(4:6)],1,which.max)]
pred.2model <- pred.2model %>% 
  separate(id,into=c("chr_str","position"), sep=":", remove = F) %>% 
  mutate(strand = str_sub(chr_str,-1),
         strand = ifelse(strand == 1, "+", "-"),
         chr = str_sub(chr_str, 1, -2),
         cov = as.numeric(cov))

pred.2model$contig=replace_list(chrmapping$chromid,chrmapping$chroms,pred.2model$chr)

#####Regression model to predict editing rate####
pred.2modelcov40 <- pred.2model %>% filter(cov >= 40)
id=pred.2modelcov40$id[pred.2modelcov40$pred %in% c("1")]
ind=which(test$info[,"id"] %in% id)
pred3 <- predict(model3,0.01*test$x[ind,,,])

pred.regression <- data.frame('id'=(test$info[ind,"id"]),
                              'cov'=test$info[ind,"cov"],
                              'ref'=test$y[ind],
                              'rate'=test$y2[ind],
                              'pred.rate'=pred3[,1]) %>% 
  group_by(id,cov,ref, rate) %>% 
  dplyr::summarise(pred.rate=mean(pred.rate)) %>% 
  ungroup()

#Combine classification and regression model together
pred.all <- left_join(pred.2model, pred.regression[,c(1,4,5)],by="id") %>% 
  filter(cov >= 20) %>% 
  mutate(pred.rate = sin(pred.rate**(5/3)))

class0 <- pred.all[pred.all$pred == 0,c("contig","position","strand","cov","0","1","2","ref","pred")]
class1 <- pred.all[pred.all$pred == 1,c("contig","position","strand","cov","0","1","2","ref","pred","rate","pred.rate")]
class2 <- pred.all[pred.all$pred == 2,c("contig","position","strand","cov","0","1","2","ref","pred")]

fwrite(class0, paste0(agggrp, ".output_prediction_CNN_class0.txt"),sep = "\t")
fwrite(class1, paste0(agggrp, ".output_prediction_CNN_class1.txt"),sep = "\t")
fwrite(class2, paste0(agggrp, ".output_prediction_CNN_class2.txt"),sep = "\t")
