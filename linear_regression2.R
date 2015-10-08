#!/usr/bin/env Rscript

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

parser$add_argument("-k", "--random_samples", type="integer", default=100,
                    help="Number of random samples used for AIC approach [default %(default)s]",
                    metavar="number")
# 
# parser$add_argument("-c", "--cores", type="integer", default=1,
#                     help="Number of cores to use [default %(default)s]",
#                     metavar="number")

parser$add_argument("input_file", nargs=1, 
                    help="File with R data.frame to be loaded, first column must be response")

parser$add_argument("output_file", nargs=1, 
                    help="File to save all the models")

args <- parser$parse_args()

input_file <- args$input_file
output_file <- args$output_file

# cores = args$cores

# input_file ="./results/2015-08-03/data.AA//data.AA.log.aspartate.1.0.RData"
repeats = args$repeats
k = args$random_samples
preprocess = args$preprocess

# cores = 1
# repeats = 100
# k = 100
# preprocess = T

repeatedCV = function(fit, repeats = 100) {
  
  #cross-valitation of all data
  theta.fit <- function(x,y){lsfit(x,y)}
  theta.predict <- function(fit,x){cbind(1,x)%*%fit$coef} 
  
  X <- fit$model[,-1]
  y <- fit$model[,1]
  
  CVs = c()
  for (tmp.i in 1:repeats) {
    cv.results <- crossval(x=X, y=y, 
                           theta.fit=theta.fit, 
                           theta.predict=theta.predict, ngroup=10)
    cv.r.squared = cor(y,cv.results$cv.fit)**2
    CVs = c(CVs, cv.r.squared)
  }
  
  CVs = data.frame(CVs)  
  return(CVs)
}



if( file.access(input_file) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", input_file))
} else {
  input.data <- get(load(input_file))
}


## -- SETTINGS ----

library(caret)
library(plyr)
library(dplyr)
library(reshape2)
library(leaps)
library(lmtest)
library(bootstrap)
library(car)

input.data = na.omit(input.data)

if (preprocess) {
  trans = preProcess(x = input.data, method = c("BoxCox", "center", "scale"))  
} else {
  trans = preProcess(x = input.data, method = c("center", "scale"))  
}

input.data.trans = predict(trans, input.data)

y = input.data.trans[,1]
X = input.data.trans[,-1]

trans.x = NULL
if (ncol(X) - max(laply(createFolds(y), length)) > nrow(X) ) { # checking if there is enough data points for CV
  if (preprocess) {
    trans.x = preProcess(x = input.data, method = c("BoxCox", "center", "scale", "pca"))  
  } else {
    trans.x = preProcess(x = input.data, method = c("center", "scale", "pca"))  
  }
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

if (ncol(tmp.df.scaled)-1 <= 40) {
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

#if no models which were identified more than once choose randomly any
if (length(tbl.fomulas[tbl.fomulas > 1]) != 0 ) {
  best.models = na.omit(names(sort(-tbl.fomulas[tbl.fomulas > 1])[1:5]))
} else { 
  set.seed(123)
  best.models = sample(formulas,size=10)
}

best.models[which(best.models ==  "y ~ 1")] = NA

best.models = as.vector(na.omit(best.models))




m.list = list()
sumarries.df = data.frame()

#MODEL evaluation
for (m in 1:length(best.models)) {
  
  current.model = best.models[m]
  fit = lm(formula=as.formula(current.model), tmp.df.scaled)
  fit.s = summary(fit)
  
  m.list[[m]] = list()
  m.list[[m]]$before = fit
  
  dw.res = dwtest(formula = as.formula(current.model), data=fit$model)
  bg.res = bgtest(formula = as.formula(current.model), data=fit$model)
  
  all.data = na.omit(tmp.df.scaled)  
  
  set.seed(123)
  CVs = repeatedCV(fit, repeats = repeats)
      
  tmp = data.frame(model = m,
                   type = factor("before", levels=c("before", "after")),
                   formula = best.models[m],
#                    coeficients = rownames(fit.s$coefficients),
#                    coef.values = fit.s$coefficients[,1],
                   r.squared = fit.s$r.squared,
                   adj.r.squared = fit.s$adj.r.squared,                 
                   p.value = ifelse(is.null(fit.s$fstatistic[1]), NA, (1 - pf(fit.s$fstatistic[1],fit.s$fstatistic[2],fit.s$fstatistic[3]))),
                   median.cv.r2 = median(CVs$CVs, na.rm=T),
                   dw.p.value = dw.res$p.value,
                   bg.p.value = bg.res$p.value,
                   datafile = input_file,
                   isPreprocessed = preprocess)
  
  sumarries.df = rbind(sumarries.df, tmp)  
    
  outliers = outlierTest(fit)
  cooks_thr = 4/(nrow(na.omit(fit$model)) - length(fit$coefficients)-2)
  
  clean.data = droplevels(fit$model[!(rownames(fit$model) %in% names(outliers$bonf.p < 0.05)),])
  #cooks_thr = 4/nrow(na.omit(tmp.df.scaled))
  
  if ( (nrow(fit$model) - sum(cooks.distance(fit) > cooks_thr) - ceiling(nrow(fit$model)/10)) < length(fit$coefficients)) {
    cooks_thr  = sort(cooks.distance(fit), decreasing=T)[3]
  }
  
  #all.data = droplevels(all.data[!(rownames(all.data) %in% names(which(cooks.distance(fit) == max(cooks.distance(fit))))),])
  clean.data = droplevels(clean.data[!(rownames(clean.data) %in% names(which((cooks.distance(fit) > cooks_thr) == TRUE))),])
  
  fit.after = lm(formula=current.model, data=clean.data)
  fit.after.s = summary(fit.after)
  
  m.list[[m]]$after = fit.after

  set.seed(123)
  CVs.after = repeatedCV(fit.after, repeats = repeats)  
  dw.res.after = dwtest(formula = as.formula(current.model), data=fit.after$model)
  bg.res.after = bgtest(formula = as.formula(current.model), data=fit.after$model)

  tmp.after = data.frame(model = m,
                         type = factor("after", levels=c("before", "after")),
                         formula = best.models[m],
#                          coeficients = rownames(fit.after.s$coefficients),
#                          coef.values = fit.after.s$coefficients[,1],
                         r.squared = fit.after.s$r.squared,
                         adj.r.squared = fit.after.s$adj.r.squared,                 
                         p.value = ifelse(is.null(fit.after.s$fstatistic[1]), NA, (1 - pf(fit.after.s$fstatistic[1],
                                                                                          fit.after.s$fstatistic[2],
                                                                                          fit.after.s$fstatistic[3]))),
                         median.cv.r2 = median(CVs.after$CVs, na.rm=T),
                         dw.p.value = dw.res.after$p.value,
                         bg.p.value = bg.res.after$p.value,
                         datafile = input_file,
                         isPreprocessed = preprocess)
      
  sumarries.df = rbind(sumarries.df, tmp.after)  
}

results = list(models = m.list,
               summaries = sumarries.df)


file_name = paste(output_file)
file_path = file_name
save(results, file=file_path)



