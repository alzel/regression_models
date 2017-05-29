#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()


all_methods = c("lmProfile",
                "rfProfile",
                "glmStepAICModel",
                "fobaModel",
                "plsModel",
                "enetModel",
                "lassoModel",
                "earthModel",
                "svmRModel",
                "nnetModel",
                "rpartModel",
                "ctreeModel",
                "treebagModel",
                "rfModel",
                "gbmModel",
                "rfBorutaModel",
                "rknnBelModel",
                "mtModel",
                "cbModel")

parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")

parser$add_argument("-p", "--preprocess", action="store_true", default=FALSE,
                    help="Apply Box-Cox, centering, scaling [default]")

parser$add_argument("-f", "--force_pca", action="store_true", default=FALSE,
                    help="Force doing PCA on predictors [default]")

parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")

parser$add_argument("-r", "--repeats", type="integer", default=1,
                    help="Number of repeats used for cross-validation [default %(default)s]",
                    metavar="number")

parser$add_argument("-c", "--cores", type="integer", default=1,
                    help="Number of cores to use [default %(default)s]",
                    metavar="number")
# 
# parser$add_argument("-b", "--best", action="store_true", default=FALSE,
#                     help="Use best predictors found in linear regression [default]")

parser$add_argument("-l", "--report", action="store_true", default=TRUE,
                    help="Save report to output_file.pdf [default]")

parser$add_argument("-t", "--threshold", type="double", default=0.85,
                    help="Treshold to remove highly correlated variables [default %(default)s]",
                    metavar="number")

parser$add_argument("-m", "--method", type="character", default="lmProfile",
                    help=c("Choose one of the following methods [default %(default)s]", all_methods))
#, paste0(all_methods, collapse = "\n")))
                    
parser$add_argument("input_file", nargs=1, 
                    help="File with R data.frame to be loaded, first column must be response")

parser$add_argument("output_file", nargs=1, 
                    help="File to save all the models")

args <- parser$parse_args()

input_file <- args$input_file
output_file <- args$output_file
report = args$report
cores = args$cores
repeats = args$repeats
#select.best = args$best
preprocess = args$preprocess
cor_thr = args$threshold
forcePCA = args$force_pca
method = args$method
# # 

# rm(list=ls())
# #
# input_file = "./results/2016-02-24/data.AA/data.AA.alanine.3.0.RData"
# output_file = "test.Rdata"
# method = "lmProfile"
# report = T
# cores = 1
# repeats = 1
# select.best = F
# preprocess = T
# cor_thr = 1
# forcePCA = T
# method = "gbm"


## -- SETTINGS ----

if (cores > 1 ) {
  library("doParallel")
  ncores = detectCores()
  cl <- makeCluster(ncores) # Register cluster
  registerDoParallel(cl)  
}

library(caret)
library(plyr)
library(dplyr)
library(reshape2)


controlObject <- trainControl(method = "repeatedcv", 
                              repeats = repeats,
                              number = 10,
                              verbose = F)

my_models = list()


if( file.access(input_file) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", input_file))
} else {
  input.data <- get(load(input_file))
}

## -- data transformation ---- 
input.data = na.omit(input.data)

input.data.tmp = input.data[,-1]
toRemove = findCorrelation(cor(input.data.tmp), cutoff = 0.95, exact = TRUE)

#always do PCA in the case if number of good (non-correlated) variables is less than n-max(length(k-fold)) of samples
doPCA <- ifelse(ncol(input.data[,-1]) - length(toRemove) > length(na.omit(input.data[,1])) - max(laply(createFolds(na.omit(input.data[,1])), length)), T, F)

#always do PCA if number of highly correlated predictors is very high  
doPCA <- ifelse(ncol(input.data[,-1]) - length(toRemove) <= 2, T, F)

#doPCA anyway if forced
if(forcePCA) {
  doPCA <- TRUE
}


if (!doPCA & length(input.data) >  3 & length(toRemove) > 0) {
  input.data.tmp = as.data.frame(cbind(input.data[,1],input.data.tmp[,-toRemove]))  
  if (length(input.data.tmp) >2) {
    names(input.data.tmp)[1] = names(input.data)[1]
    rownames(input.data.tmp) = rownames(input.data)
  } else if (length(input.data.tmp) == 2) {
    stop(paste("Something wrong in file", input_file) )
    #tmp.idx <- which(!(1:length(input.data) %in% toRemove))
    #names(input.data.tmp)[tmp.idx] = names(input.data)[tmp.idx]
  } else {
    stop(paste("Something wrong in file", input_file) )
  }
  input.data = input.data.tmp
}


if (preprocess) {
  trans = preProcess(x = input.data, method = c("BoxCox", "center", "scale"))  
} else {
  trans = preProcess(x = input.data, method = c("center", "scale"))  
}

input.data.trans = predict(trans, input.data)

y = input.data.trans[,1]
X = input.data.trans[,-1]
trans.x = NULL

if (doPCA) {
  
  if (preprocess) {
    trans.x = preProcess(x = input.data[,-1], method = c("BoxCox", "center", "scale", "pca"), thresh = 0.99)  
  } else {
    trans.x = preProcess(x = input.data[,-1], method = c("center", "scale", "pca"), thresh = 0.99)  
  }
  
  X = predict(trans.x, input.data[,-1])  
}

# seeds for multicores
myseeds = as.list(rep(list(rep(123,ncol(X))),repeats * 10))
myseeds[[length(myseeds) + 1]] = 123


switch(method,
  lmProfile = {
    ## -- linear regression variable selection ----
    
    ctrl <- rfeControl(functions = lmFuncs,
                       method = "repeatedcv",
                       repeats = repeats,
                       number = 10,
                       verbose = F,
                       seeds = myseeds,
                       returnResamp = "all")
    
    SUBS = c(1:(ncol(X)-1)) #sizes of variable subset to search
    result <- rfe(X, y,
                     sizes = SUBS,
                     rfeControl = ctrl)
    my_models[["lmProfile"]] = result  
  },
  rfProfile = {
    ## -- random forest variable selection ----
    ctrl <- rfeControl(functions = rfFuncs,
                       method = "repeatedcv",
                       repeats = repeats,
                       verbose = F,
                       returnResamp = "all")
    
    result <- rfe(X, y,
                     sizes = SUBS,
                     rfeControl = ctrl)
    
    my_models[["rfProfile"]] = result
    
  },
  glmStepAICModel = {
    ## -- glmStepAIC ----
    message("glmStepAIC")
    glmStepAICModel <- train(y = y, x = X,
                             method = "glmStepAIC",
                             trControl = controlObject)
    my_models[["glmStepAICModel"]] = glmStepAICModel
  },
  fobaModel = {
    ## -- foba ----
    message("foba")
    fobaModel <- train(y = y, x = X,
                     method = "foba",
                     tuneGrid = expand.grid(k = 1:(ncol(X)-1),
                                            lambda = c(0, .001, .01, .1)),
                     trControl = controlObject)
    my_models[["fobaModel"]] = fobaModel
  },
  plsModel = {
    ## -- pls ----
    message("pls")
    set.seed(123)
    plsModel <- train(y = y, x = X,
                      method = "pls",
                      tuneLength = 15,
                      trControl = controlObject)
    my_models[["plsModel"]] = plsModel
    
  },
  enetModel = {
    ## -- enet ----
    message("enet")
    
    enetGrid <- expand.grid(.lambda = c(0, .001, .01, .1), 
                            .fraction = seq(0.05, 1, length = 20))
    set.seed(123)
    enetModel <- train(y=y, x = X,
                       method = "enet",
                       tuneGrid = enetGrid,
                       trControl = controlObject)
    
    my_models[["enetModel"]] = enetModel
  },  
  lassoModel = {
    ## -- lasso ----
    message("lasso")
    
    lassoGrid <- expand.grid(.fraction = seq(0.05, 1, length = 20))
    
    set.seed(123)
    lassoModel <- train(y=y, x = X,
                        method = "lasso",
                        tuneGrid = lassoGrid,
                        trControl = controlObject)
    
    my_models[["lassoModel"]] = lassoModel
    
  },
  earthModel = {
    ## -- earth ----
    message("earth")
    
    set.seed(123)
    earthModel <- train(y = y,
                        x = X,
                        method = "earth",
                        tuneGrid = expand.grid(.degree = 1:3,
                                               .nprune = 2:10),    
                        trControl = controlObject)
    my_models[["earthModel"]] = earthModel
    
  },
  svmRModel = {
    ## -- svmRadial ----
    message("svmRadial")
    set.seed(123)
    svmRModel <- train(y = y, x = X,
                       method = "svmRadial",
                       tuneLength = 20,
                       trControl = controlObject)
    my_models[["svmRModel"]] = svmRModel
    
  },
  nnetModel = {
    ## -- avNNet ----
    message("avNNet")
    
    nnetGrid <- expand.grid(.decay = c(0.001, .01, .1), 
                            .size = seq(1, 10, by = 1),
                            .bag = FALSE)
    set.seed(123)
    nnetModel <- train(y = y,
                       x =  X,
                       method = "avNNet",
                       tuneGrid = nnetGrid,
                       linout = TRUE,
                       trace = FALSE,
                       maxit = 1000,
                       trControl = controlObject)
    my_models[["nnetModel"]] = nnetModel
    
  },
  rpartModel = {
    ## -- rpart ----
    message("rpart")
    
    set.seed(123)
    rpartModel <- train(y = y, 
                        x = X,
                        method = "rpart",
                        tuneLength = 30,
                        trControl = controlObject)
    my_models[["rpartModel"]] = rpartModel
  },
  ctreeModel = {
    ## -- ctree ----
    message("ctree")
    
    set.seed(123)
    ctreeModel <- train(y = y,
                        x = X,                          
                        method = "ctree",
                        tuneLength = 10,
                        trControl = controlObject)
    my_models[["ctreeModel"]] = ctreeModel
  },
  treebagModel = {
    ## -- treebag ----
    message("treebag")
    
    set.seed(123)
    treebagModel <- train(y = y,
                          x = X,
                          method = "treebag",
                          trControl = controlObject)
    my_models[["treebagModel"]] = treebagModel
  },
  rfModel = {
    ## -- rf ----
    message("rf")
    set.seed(123)
    rfModel <- train(y = y,
                     x = X,
                     method = "rf",
                     tuneLength = 10,
                     ntrees = 1000,
                     importance = TRUE,
                     trControl = controlObject)
    my_models[["rfModel"]] = rfModel
  },
  gbmModel = {
    ## -- gmb ----
    message("gmb")
    
    gbmGrid <- expand.grid(interaction.depth = seq(1, 7, by = 1),
                           n.minobsinnode = 2:5,
                           n.trees = seq(100, 1000, by = 50),
                           shrinkage = c(0.01, 0.1))
    
    set.seed(123)
    gbmModel <- train(y = y,
                      x = X,
                      method = "gbm",
                      tuneGrid = gbmGrid,
                      verbose = FALSE,
                      trControl = controlObject)
    my_models[["gbmModel"]] = gbmModel
  },
  rfBorutaModel = {
    ## -- boruta ----
    message("boruta")
    rfBorutaModel <<- train(y = y, x = X,
                            method = "Boruta",
                            tuneGrid = expand.grid(.mtry = unique(round(ncol(X)/3:5))),
                            trControl = controlObject)
    my_models[["rfBorutaModel"]] = rfBorutaModel  
  },
  rknnBelModel = {
    ## -- rknnBel ----
    message("rknnBel")
    rknnBelModel <<- train(y = y, x = X,
                           method = "rknnBel",
                           tuneGrid = expand.grid(mtry = unique(2:round(ncol(X)/2)),
                                                  d = 1,
                                                  k = 2:5),
                           trControl = controlObject)
    my_models[["rknnBelModel"]] = rknnBelModel
  },
  mtModel = {
    ## -- M5 ----
    message("M5")
    
    set.seed(123)
    mtModel <- train(y = y,
                     x = X,
                     method = "M5",
                     trControl = controlObject)
    my_models[["mtModel"]] = mtModel
  },
  cbModel = {
    ## -- cubist ----
    message("cubist")
    
    cubistGrid <- expand.grid(.committees = c(1, 5, 10, 50, 75, 100),
                              .neighbors = c(0, 1, 3, 5, 7, 9))
    set.seed(123)
    cbModel <- train(y = y,
                     x = X,
                     method = "cubist",
                     tuneGrid = cubistGrid,
                     trControl = controlObject)
    my_models[["cbModel"]] = cbModel
  },
  {
    stop("No such method")  
  }
)
  
  
  
  
  
  
  ## -- robust regression variable selection ----
  # rlmFuncs = lmFuncs
  # 
  # rlmFuncs$fit = function (x, y, first, last, ...)
  # {
  #   tmp = data.frame()
  #   if (is.data.frame(x)) {
  #     tmp <- x
  #   } else {
  #     tmp <- as.data.frame(x)
  #   }
  #   tmp$y <- y
  #   rlm(y ~ ., data = tmp, maxit = 1000000, method = "M")
  # }
  # 
  # ctrl <- rfeControl(functions = rlmFuncs,
  #                    method = "repeatedcv",
  #                    repeats = repeats,
  #                    verbose = T,
  #                    returnResamp = "all")
  # 
  # rlmProfile <- rfe(X, y,
  #                  rfeControl = ctrl)
  # 
  # my_models[["rlmProfile"]] = rlmProfile

## -- save output ----

if (cores > 1) {
  stopCluster(cl)  
}

results = list(input_data = input.data,
               trans = trans,
               trans.x = trans.x,
               models = my_models)

file_name = paste(output_file)
file_path = file_name
save(results, file=file_path)

if (report) {
  if (length(my_models) == 0) {stop ("Nothing to report, no models found")}
  file_name = paste(output_file, "pdf", sep = ".")
  pdf(file_name, height=247/25.4, width=183/25.4)
    par(pty="s")
    p = ggplot(my_models[[1]], metric = "Rsquared") + 
      ggtitle(output_file) + 
      theme(aspect.ratio = 5/8)
    print(p)
    p = ggplot(my_models[[1]], metric = "RMSE")+ 
      ggtitle(output_file) + 
      theme(aspect.ratio = 5/8)
    print(p)
  dev.off()
}


