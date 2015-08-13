#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")

parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")

parser$add_argument("-r", "--repeats", type="integer", default=1,
                    help="Number of repeats used for cross-validation [default %(default)s]",
                    metavar="number")

parser$add_argument("-k", "--random_samples", type="integer", default=100,
                    help="Number of random samples used for AIC approach [default %(default)s]",
                    metavar="number")

parser$add_argument("-c", "--cores", type="integer", default=1,
                    help="Number of cores to use [default %(default)s]",
                    metavar="number")

parser$add_argument("input_file", nargs=1, 
                    help="File with R data.frame to be loaded, first column must be response")

parser$add_argument("output_file", nargs=1, 
                    help="File to save all the models")

args <- parser$parse_args()

input_file <- args$input_file
output_file <- args$output_file

cores = args$cores
repeats = args$repeats
k = args$random_samples

#input_file ="./results/2015-08-03/data.AA/data.AA.alanine.1.0.RData"

if( file.access(input_file) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", input_file))
} else {
  input.data <- get(load(input_file))
}


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
library(leaps)

input.data = na.omit(input.data)
trans = preProcess(x = input.data, method = c("BoxCox", "center", "scale"))
input.data.trans = predict(trans, input.data)

y = input.data.trans[,1]
X = input.data.trans[,-1]

trans.x = NULL
if (ncol(X) - max(laply(createFolds(y), length)) > nrow(X) ) { # checking if there is enough data points for CV
  trans.x = preProcess(x = input.data[,-1], method=c( "BoxCox", "center", "scale","pca"))
  X = predict(trans.x, input.data[,-1])
}

data = cbind(y,X)
tmp.df.scaled = data

sample_size = round(nrow(tmp.df.scaled)*0.9)
total_samples = nrow(tmp.df.scaled)


sample.matrix = matrix(rep(0, k*sample_size), nrow=k)

set.seed(123)
for (ik in 1:k) {
  selected = sample(1:total_samples, size=sample_size)  
  sample.matrix[ik,] = selected
}

formulas = c() #storing models

i = "y"
#choosing best model using exhaustive approach
if (ncol(tmp.df.scaled)-1 <= 30) {
  for(j in 1:nrow(sample.matrix)) {
    
    sample.data = tmp.df.scaled[sample.matrix[j,],]
    
    NVMAX = ncol(tmp.df.scaled) - 1
    if ( NVMAX >= nrow(sample.data[,-1]) - max(laply(createFolds(sample.data[,1]), length)) )  {
      #NVMAX = round(nrow(sample.data)/10)*9 - 1 #to ensure 10-fold cross validation validity
      NVMAX = ncol(sample.data[,-1]) - max(laply(createFolds(sample.data[,1], 10), length))
    }
        
    b <- regsubsets(as.formula(paste0(i," ~ ", ".")), data=sample.data, nbest=1, nvmax=NVMAX)
    rs = summary(b)
    
    n_points = nrow(na.omit(sample.data))
    k_params = apply(rs$which, 1, sum)
    tmp.best = data.frame(n_params = apply(rs$which, 1, sum),
                          cp = rs$cp,
                          adjr2 = rs$adjr2,
                          aic = n_points*log(rs$rss/n_points) + 2*k_params,
                          bic = rs$bic)
    
    best_idx = which(tmp.best$aic == min(tmp.best$aic))
    #best_idx = which(tmp.best$adjr2 == max(tmp.best$adjr2))
    #best_idx = which(tmp.best$cp == min(tmp.best$cp))
    #best_idx = which(tmp.best$bic == min(tmp.best$bic))
    F1 = paste0(i," ~ ", paste(sort(names(which(rs$which[best_idx,-1] == TRUE))), collapse=" + "))
    formulas = c(formulas, F1)
  }
  
} else { #using step approach
  
  for(j in 1:nrow(sample.matrix)) {
    sample.data <- na.omit(tmp.df.scaled[sample.matrix[j,],])
    
    #null.lm <- lm(formula=paste0(i," ~ ", "1"), data=sample.data)
    #full.lm <- lm(formula=paste0(i," ~ ", "."), data=sample.data)
    
    null.lm <- do.call("lm", list(paste0(i," ~ ", "1"),
                                  data = sample.data))
    
    full.lm <- do.call("lm", list(paste0(i," ~ ", "."),
                                  data = sample.data))
    
    result = step(null.lm, scope=list(lower=null.lm, upper=full.lm), direction="both", steps=100000)
    tmp_vars = paste(sort(names(result$coefficients)[-1]), collapse=" + ")
    
    if ( tmp_vars == "" ) {
      tmp_vars = 1
    }
    F1 = paste0(i, " ~ ", tmp_vars)
    formulas = c(formulas, F1)
  }
}


tbl.fomulas = table(formulas)

if (length(tbl.fomulas[tbl.fomulas > 1]) != 0 ) {
  best.models = na.omit(names(sort(-tbl.fomulas[tbl.fomulas > 1])[1:5]))
} else {
  set.seed(123)
  best.models = sample(formulas,size=10)
}
                  

## -- model selection ----
controlObject <- trainControl(method = "repeatedcv", 
                              repeats = 100,
                              number = 10,
                              verbose = F)


my_models = list()
sumarries.df = data.frame()
fits.list = list()

for (m in 1:length(best.models)) {
  current.model = best.models[m]
  lmModel <- train(as.formula(current.model), data=tmp.df.scaled,
                   method = "lm",
                   trControl = controlObject)
  my_models[[m]] = lmModel
}


sumarries.df = data.frame()
fits.list = list()

for(j in 1:nrow(sample.matrix)) {
  sample.data = tmp.df.scaled[sample.matrix[j,],]
  m.list = list()
  for (m in 1:length(best.models)) {
    current.model = best.models[m]
    fit = lm(formula=as.formula(current.model), sample.data)
    fit.s = summary(fit)
    m.list[[m]] = fit
    
    tmp = data.frame(model = m, 
                     data = j,
                     formula = best.models[m],
                     coeficients = rownames(fit.s$coefficients),
                     values = fit.s$coefficients[,1],
                     r.squared = fit.s$r.squared,
                     adj.r.squared = fit.s$adj.r.squared,                 
                     p.value = ifelse(is.null(fit.s$fstatistic[1]), NA, (1 - pf(fit.s$fstatistic[1],fit.s$fstatistic[2],fit.s$fstatistic[3]))),
                     aic = AIC(fit))
    
    sumarries.df = rbind(sumarries.df, tmp)  
  }
  fits.list[[j]] = m.list
}


results = list(caret_models = my_models,
               sample_models = sumarries.df)

if (cores > 1) {
  stopCluster(cl)  
}

file_name = paste(output_file)
file_path = file_name
save(results, file=file_path)



