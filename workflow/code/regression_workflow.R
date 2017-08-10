library(randomForest)
library(ggplot2)
library(ROCR)
library(gplots)
library(reshape)
library(MASS)
library(caTools)

setwd('~/Desktop/Projects/git/2016_summerproject/new_plots/')

load(file='random_state_seed.RData') # Loading the seed ensures the workflow is always reproducible

#---------------------------------------------
# Functions
#---------------------------------------------
seq_rf = function(train, imp.ord, reps) { # Function for building iterative models for feature selection
  dat.t = train
  imp = imp.ord
  genes = row.names(subset(imp, Rank<=1))
  dat.sub = as.data.frame(dat.t[,c(genes)])
  colnames(dat.sub) = genes
  row.names(dat.sub) = row.names(dat.t)
  dat.sub$Activity = dat.t$Activity
  oob = 0
  print(1)
  for (j in c(1:reps)) {
    rf = randomForest(data=dat.sub, Activity~., ntree=1001, type='regression', mtry=1)
    oob = oob + as.numeric(rf$mse[1001])
  }
  oob = oob/reps
  error = data.frame(Genes=1, OOB=oob)
  for (i in c(2:nrow(imp))) {
    print(i)
    genes = row.names(subset(imp, Rank<=i))
    dat.sub = dat.t[,c(genes)]
    dat.sub$Activity = dat.t$Activity
    oob = 0
    for (j in c(1:reps)) {
      rf = randomForest(data=dat.sub, Activity~., ntree=1001, type='regression', mtry=i)
      oob = oob + as.numeric(rf$mse[1001])
    }
    oob = oob/reps
    error = rbind(error, c(i,oob))
  }
  return(error)
}
#---------------------------------------------
# Workflow
#---------------------------------------------
  
### Read in datasets
c.dat = read.delim(file='data/Chari.txt', sep="", header=T, row.names=1)   # Indel Dataset: used for feature selection and training of final model

### Perform feature selection on Indel dataset
chari.full.rf  = randomForest(Activity ~ ., data=c.dat, ntree=10001, importance=T)
re.imp         = as.data.frame(varImpPlot(chari.full.rf))
re.impOrd      = re.imp[order(-re.imp$`%IncMSE`),]
re.impOrd$Rank = c(1:nrow(re.impOrd))

re.error      = seq_rf(c.dat, re.impOrd, 3) # Construct interative models, adding progressively more features in order of improtance

re.min.rf  = as.numeric(row.names(subset(re.error, OOB<=min(re.error$OOB))))[1] # Extract smallest model with lowest OOB error
re.vars.rf = row.names(subset(re.impOrd, Rank<=re.min.rf))

### Build final model
re.rf = randomForest(Activity ~ ., data=subset(c.dat, select=c(re.vars.rf,'Activity')), ntree=10001) 
