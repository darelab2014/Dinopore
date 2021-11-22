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
  make_option(c("-e", "--epoch"), type="integer", default=900,
              help="Number of epochs", metavar="integer"),
  make_option(c("-b", "--batch"), type="integer", default=1024,
              help="Batch size", metavar="integer"),
  make_option(c("-s", "--seed"), type="integer", default=9999,
              help="Seed number", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list,add_help_option=FALSE)
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
vali=readRDS(opt$vali)

cl1=which(train$y==1)
train.r=list()
train.r$x=train$x[cl1,,,]
train.r$y=train$y2[cl1]
train.r$x=array_reshape(train.r$x,dim=c(length(train.r$y),5,65,1))
train.r$info=train$info[cl1,]

cl1=which(vali$y==1)
vali.r=list()
vali.r$x=vali$x[cl1,,,]
vali.r$y=vali$y2[cl1]
vali.r$x=array_reshape(vali.r$x,dim=c(length(vali.r$y),5,65,1))
vali.r$info=vali$info[cl1,]

train.r$xn=0.01*train.r$x
train.r$y=asin(train.r$y)**0.6
vali.r$xn=0.01*vali.r$x
vali.r$y=asin(vali.r$y)**0.6


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

residual_module <- function(layer_in){
  X5=layer_in  %>% layer_conv_2d(filters = 20,
                                 kernel_size = c(5,5),
                                 kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001)) %>% layer_activation_relu()
  X4=layer_in  %>% layer_conv_2d(filters = 20,
                                 kernel_size = c(4,4),
                                 kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001))%>% layer_activation_relu()
  X3=layer_in  %>% layer_conv_2d(filters = 20,    
                                 kernel_size = c(3,3),
                                 kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001)) %>% layer_activation_relu()
  X=layer_in  %>% layer_conv_2d(filters = 20,       
                                kernel_size = c(2,2),
                                kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001)) %>% layer_activation_relu() %>%
    layer_locally_connected_2d(filters =20,
                               kernel_size = c(2,2),
                               kernel_regularizer = regularizer_l2(l = 0.001))%>% layer_activation_relu()
  Y=layer_add(c(X3,X))
  Y3=Y %>% layer_locally_connected_2d(filters = 20,
                                      kernel_size = c(3,3),
                                      kernel_regularizer = regularizer_l2(l = 0.001))%>% layer_activation_relu()
  Y2=Y %>% layer_locally_connected_2d(filters = 20,
                                      kernel_size = c(2,2),
                                      kernel_regularizer = regularizer_l2(l = 0.001))%>% layer_activation_relu()
  Z=layer_add(c(Y2,X4))
  Z=Z %>% layer_conv_2d(filters = 20,       
                        kernel_size = c(2,2),  
                        kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001))%>% layer_activation_relu()
  final=layer_add(c(Z,Y3,X5))
  layer_out=final %>% 
    layer_dropout(rate = 0.3) %>%
    layer_flatten()%>%
    layer_dense(units = 16384, kernel_initializer="normal")%>% layer_activation_relu() %>%
    layer_dropout(rate= 0.3)%>%
    layer_dense(units = 1024, kernel_initializer="normal")%>% layer_activation_relu() %>%
    layer_dropout(rate= 0.3)%>%
    layer_dense(units = 1, activation = "linear")
  
  return(layer_out)
}
visible=layer_input(shape=c(5,65,1))
layer=residual_module(visible)

backend()$clear_session
model2=keras_model(inputs = visible,outputs = layer)


model2 %>% compile(
  loss = "mean_squared_error",
  optimizer = 'adam',
  metrics = c('mean_absolute_error')
)

save.checkpoint=callback_model_checkpoint(paste0(opt$name,'.best_regression_morefts_16384_1024_asin06.h5'),monitor = "val_mean_absolute_error",save_best_only = T,mode="auto")


history <- model2 %>% fit(
  x=train.r$xn, y=train.r$y, 
  epochs = opt$epoch, batch_size = opt$batch, 
  callbacks = save.checkpoint,
  validation_data=list(vali.r$xn,vali.r$y),
)
