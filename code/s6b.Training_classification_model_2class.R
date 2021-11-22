pacman::p_load(data.table,tidyverse,caret,stringr,keras,tensorflow, optparse, multiROC)

rm(list=ls())
gc()

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-v", "--vali"), type="character", default=NULL,
              help="validation file name", metavar="character"),
  make_option(c("-t", "--train"), type="character", default=NULL,
              help="training file name", metavar="character"),
  make_option(c("-o", "--name"), type="character", default="out",
              help="common name for model", metavar="character"),
  make_option(c("-e", "--epoch"), type="integer", default=500,
              help="Number of epochs", metavar="integer"),
  make_option(c("-b", "--batch"), type="integer", default=256,
              help="Batch size", metavar="integer"),
  make_option(c("-s", "--seed"), type="integer", default=9999,
              help="Seed number", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$vali)|is.null(opt$train)){
  print_help(opt_parser)
  stop("Validation/testing and training files must be supplied (-v & -t).", call.=FALSE)
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

train=readRDS(opt$train)
train$x=train$x[,,c(1:43),]
train$x=array_reshape(train$x,dim=c(nrow(train$x),5,43,1))
vali=readRDS(opt$vali)
vali$x=vali$x[,,c(1:43),]
vali$x=array_reshape(vali$x,dim=c(nrow(vali$x),5,43,1))

reset_keras=function(){
  sess=tf$compat$v1$keras$backend$get_session()
  backend()$clear_session
  sess$close()
  sess=tf$compat$v1$keras$backend$get_session()
  set.seed(opt$seed)
}

gpus = tf$config$experimental$list_physical_devices('GPU')
tf$config$experimental$set_visible_devices(gpus[1], 'GPU')

reset_keras()
tf$random$set_seed(opt$seed)

residual_module <- function(layer_in,number_class){
  X5=layer_in  %>% layer_conv_2d(filters = 60,
                               kernel_size = c(5,5),
                               activation = "relu",
                               kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.001))
  X4=layer_in  %>% layer_conv_2d(filters = 60,
                                 kernel_size = c(4,4),
                                 activation = "relu",
                                 kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.001))
  X3=layer_in  %>% layer_conv_2d(filters = 60,       
                                 kernel_size = c(3,3),             
                                 activation = "relu",
                                 kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.001)) 
  X=layer_in  %>% layer_conv_2d(filters = 60,       
                                 kernel_size = c(2,2),             
                                 activation = "relu",
                                 kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.001)) %>% 
                  layer_locally_connected_2d(filters =60,
                                             kernel_size = c(2,2),
                                             activation = "relu",
                                             kernel_regularizer = regularizer_l2(l = 0.001))
  Y=layer_add(c(X3,X))
  Y3=Y %>% layer_locally_connected_2d(filters = 60,
                                     kernel_size = c(3,3),
                                     activation = "relu",
                                     kernel_regularizer = regularizer_l2(l = 0.001))
  Y2=Y %>% layer_locally_connected_2d(filters = 60,
                               kernel_size = c(2,2),
                               activation = "relu",
                               kernel_regularizer = regularizer_l2(l = 0.001))
  Z=layer_add(c(Y2,X4))
  Z=Z %>% layer_conv_2d(filters = 60,       
                  kernel_size = c(2,2),             
                  activation = "relu",
                  kernel_regularizer = regularizer_l1_l2(l1 = 0.01, l2 = 0.001))
  final=layer_add(c(Z,Y3,X5))
  layer_out=final %>%   layer_dropout(rate = 0.2) %>%
    layer_flatten()%>%
    layer_dense(units = 256, activation = "relu")%>% 
    layer_dropout(rate= 0.3)%>%
    layer_dense(units = 512, activation = "relu")%>%
    layer_dropout(rate= 0.2)%>%
    layer_dense(units = number_class, activation = "softmax")
  
  return(layer_out)
}

###2-class model
reset_keras()
tf$random$set_seed(opt$seed)

cl2=which(train$y==2)
train2=list()
train2$x=train$x[-cl2,,,]
train2$y=train$y[-cl2]
train2$x=array_reshape(train2$x,dim=c(length(train2$y),5,43,1))
train2$info=train$info[-cl2,]

cl2=which(vali$y==2)
vali2=list()
vali2$x=vali$x[-cl2,,,]
vali2$y=vali$y[-cl2]
vali2$x=array_reshape(vali2$x,dim=c(length(vali2$y),5,43,1))
vali2$info=vali$info[-cl2,]

train_lab <-to_categorical(train2$y,2) #Catagorical vector for training classes
vali_lab <-to_categorical(vali2$y,2)

visible=layer_input(shape=c(5,43,1))
layer=residual_module(visible,2)
backend()$clear_session
model1=keras_model(inputs = visible,outputs = layer)

model1 %>% compile(
  loss = 'categorical_crossentropy',
  optimizer =  'adam',
  metrics = c('accuracy')
)

save.checkpoint=callback_model_checkpoint(paste0(opt$name,'.best_pos5_mix_3c_1vs1_resnet.h5'),monitor = "val_accuracy",save_best_only = T,mode="auto")

history <- model1 %>% fit(
  x=train2$x, y=train_lab, 
  epochs = opt$epoch, batch_size = opt$batch,  
  callbacks = save.checkpoint,
  validation_data=list(vali2$x,vali_lab),
  class_weight=list('0'=0.8,'1'=1.8)
)

