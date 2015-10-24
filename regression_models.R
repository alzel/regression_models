#!/usr/bin/env Rscript-3.1.0

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")

parser$add_argument("-p", "--preprocess", action="store_false", default=TRUE,
                    help="Apply Box-Cox, centering, scaling, if false no Box-Cox [default]")

parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")

parser$add_argument("-r", "--repeats", type="integer", default=1,
                    help="Number of repeats used for cross-validation [default %(default)s]",
                    metavar="number")

parser$add_argument("-c", "--cores", type="integer", default=1,
                    help="Number of cores to use [default %(default)s]",
                    metavar="number")

parser$add_argument("-b", "--best", action="store_true", default=FALSE,
                    help="Use best predictors found in linear regression [default]")

parser$add_argument("-l", "--report", action="store_true", default=TRUE,
                    help="Save report to output_file.pdf [default]")

parser$add_argument("-t", "--threshold", type="double", default=0.85,
                    help="Treshold to remove highly correlated variables [default %(default)s]",
                    metavar="number")


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
select.best = args$best
preprocess = args$preprocess
cor_thr = args$threshold
# # 
# input_file = "./results/2015-09-29/data.AA/data.AA.alanine.1.0.RData"
# output_file = "test.Rdata"
# report = T
# cores = 1
# repeats = 1
# select.best = F
# preprocess = F
# cor_thr = 0.99

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
toRemove = findCorrelation(cor(input.data.tmp), cutoff = cor_thr, exact = TRUE)

if (length(toRemove) > 0) {
  input.data.tmp = as.data.frame(cbind(input.data[,1],input.data.tmp[,-toRemove]))  
  if (length(input.data.tmp) >2) {
    names(input.data.tmp)[1] = names(input.data)[1]
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


if (ncol(X) - max(laply(createFolds(na.omit(y)), length)) > length(na.omit(y)) ) { # checking if there is enough data points for CV
  
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



## -- linear regression variable selection ----

ctrl <- rfeControl(functions = lmFuncs,
                   method = "repeatedcv",
                   repeats = repeats,
                   number = 10,
                   verbose = F,
                   seeds = myseeds,
                   returnResamp = "all")

SUBS = c(1:(ncol(X)-1)) #sizes of variable subset to search
lmProfile <- rfe(X, y,
                 sizes = SUBS,
                 rfeControl = ctrl)

my_models[["lmProfile"]] = lmProfile

## -- robust regression variable selection ----
rlmFuncs = lmFuncs

rlmFuncs$fit = function (x, y, first, last, ...)
{
  tmp <- if (is.data.frame(x)) {
    x
  } else {
    as.data.frame(x)
  }
    
  library(MASS)
  tmp$y <- y
  rlm(y ~ ., data = tmp, maxit = 100000, method = "M")
}

ctrl <- rfeControl(functions = rlmFuncs,
                   method = "repeatedcv",
                   repeats = repeats,
                   verbose = F,
                   returnResamp = "all")


rlmProfile <- rfe(X, y,
                 sizes = SUBS,
                 rfeControl = ctrl)

my_models[["rlmProfile"]] = rlmProfile

## -- random forest variable selection ----
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = repeats,
                   verbose = F,
                   returnResamp = "all")

rfProfile <- rfe(X, y,
                sizes = SUBS,
                rfeControl = ctrl)
                
my_models[["rfProfile"]] = rfProfile


# forTraining = createDataPartition(y, p=1)[[1]]
# trainingSet = input.data[forTraining,] 
# 
# trans = preProcess(x = trainingSet, method = c("BoxCox", "center", "scale"))
# input.data.trans = predict(trans, trainingSet)
# 
# y = input.data.trans[,1]
# X = input.data.trans[,-1]


## -- glmStepAIC ----
message("glmStepAIC")
glmStepAICModel <- train(y = y, x = X,
                  method = "glmStepAIC",
                  trControl = controlObject)

my_models[["glmStepAICModel"]] = glmStepAICModel


## -- foba ----
message("foba")
fobaModel <- train(y = y, x = X,
                       method = "foba",
                       tuneGrid = expand.grid(k = 1:(ncol(X)-1),
                                              lambda = c(0, .001, .01, .1)),
                       trControl = controlObject)
my_models[["fobaModel"]] = fobaModel

# ## -- boruta ----
# message("boruta")
# rfBorutaModel <<- train(y = y, x = X,
#                        method = "Boruta",
#                        tuneGrid = expand.grid(.mtry = unique(round(ncol(X)/3:5))),
#                        trControl = controlObject)
# 
# }
#my_models[["rfBorutaModel"]] = rfBorutaModel


## -- rknnBel ----
# message("rknnBel")
# rknnBelModel <<- train(y = y, x = X,
#                        method = "rknnBel",
#                        tuneGrid = expand.grid(mtry = unique(2:round(ncol(X)/2)),
#                                               d = 1,
#                                               k = 2:5),
#                        trControl = controlObject)
# 
#my_models[["rknnBelModel"]] = rknnBelModel


if (select.best) {
  X = X[,names(lmProfile$fit$coefficients)[-1]]
}

## -- pls ----
message("pls")
set.seed(123)
plsModel <- train(y = y, x = X,
                    method = "pls",
                    tuneLength = 15,
                    trControl = controlObject)
my_models[["plsModel"]] = plsModel

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


## -- lasso ----
message("lasso")

lassoGrid <- expand.grid(.fraction = seq(0.05, 1, length = 20))

set.seed(123)
lassoModel <- train(y=y, x = X,
                   method = "lasso",
                   tuneGrid = lassoGrid,
                   trControl = controlObject)

my_models[["lassoModel"]] = lassoModel


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


## -- svmRadial ----
message("svmRadial")

set.seed(123)
svmRModel <- train(y = y, x = X,
                   method = "svmRadial",
                   tuneLength = 20,
                   trControl = controlObject)
my_models[["svmRModel"]] = svmRModel


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


## -- rpart ----
message("rpart")

set.seed(123)
rpartModel <- train(y = y, 
                    x = X,
                    method = "rpart",
                    tuneLength = 30,
                    trControl = controlObject)
my_models[["rpartModel"]] = rpartModel


## -- ctree ----
message("ctree")

set.seed(123)
ctreeModel <- train(y = y,
                    x = X,                          
                    method = "ctree",
                    tuneLength = 10,
                    trControl = controlObject)
my_models[["ctreeModel"]] = ctreeModel


# ## -- M5 ----
# message("M5")
# 
# set.seed(123)
# mtModel <- train(y = y, 
#                   x = X,
#                   method = "M5",
#                   trControl = controlObject)
# my_models[["mtModel"]] = mtModel
# 



## -- treebag ----
message("treebag")

set.seed(123)
treebagModel <- train(y = y,
                      x = X,
                      method = "treebag",
                      trControl = controlObject)
my_models[["treebagModel"]] = treebagModel

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


## -- cubist ----
# message("cubist")
# 
# cubistGrid <- expand.grid(.committees = c(1, 5, 10, 50, 75, 100), 
#                           .neighbors = c(0, 1, 3, 5, 7, 9))
# set.seed(123)
# cbModel <- train(y = y,
#                  x = X,
#                  method = "cubist",
#                  tuneGrid = cubistGrid,
#                  trControl = controlObject)
# my_models[["cbModel"]] = cbModel


## -- save output ----

allResamples = resamples(my_models)

if (cores > 1) {
  stopCluster(cl)  
}


file_name = paste(output_file)
file_path = file_name
save(my_models, file=file_path)

if (report) {
  t = resamples(my_models)
  t.long = melt(t$values[,grep("Rsquared", names(t$values))],)
  t.long.stats = t.long %>% group_by(variable) %>% summarise(meanR2=mean(value, na.rm=T)) %>% arrange(meanR2)
  t.long.stats$variable = factor(t.long.stats$variable, levels=t.long.stats$variable)
  p = ggplot(data=t.long.stats, aes(x=variable, ymax=meanR2, ymin=0) ) + 
    geom_linerange() + 
    geom_point(data=t.long.stats, aes(x=variable,y=meanR2)) +
    coord_flip()
  
  file_name = paste(output_file, "pdf", sep=".")
  ggsave(plot=p, filename=file_name,  width = 210, height = 297, units = "mm")
}



