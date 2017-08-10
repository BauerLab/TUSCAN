### Workflow for building TUSCAN-Classification model

library(randomForest)
library(ggplot2)
library(ROCR)
library(gplots)
library(reshape)
library(MASS)
library(caTools)

load('random_state_seed.RData') # Loading the seed ensures the workflow is always reproducible

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
  dat.sub$Class = dat.t$Class
  oob = 0
  print(1)
  for (j in c(1:reps)) {
    rf = randomForest(data=dat.sub, Class~., ntree=1001, type='classification', mtry=1)
    oob = oob + as.numeric(rf$err.rate[1001,1])
  }
  oob = oob/reps
  error = data.frame(Genes=1, OOB=oob)
  for (i in c(2:nrow(imp))) {
    print(i)
    genes = row.names(subset(imp, Rank<=i))
    dat.sub = dat.t[,c(genes)]
    dat.sub$Class = dat.t$Class
    oob = 0
    for (j in c(1:reps)) {
      rf = randomForest(data=dat.sub, Class~., ntree=1001, type='classificatin', mtry=i)
      oob = oob + as.numeric(rf$err.rate[1001,1])
    }
    oob = oob/reps
    error = rbind(error, c(i,oob))
  }
  return(error)
}

rf_cvTest = function(df, fold, m) { # Function for performing Cross Validation
  dat = df
  nTest = round(nrow(dat)/fold,0)
  pred = c(0,0)
  class = c(0)
  
  for (k in (c(1:fold))) {
    dat.test = dat[sample(nrow(dat), nTest),]
    class = c(class, as.numeric(as.character(dat.test$Class)))
    dat.train = dat[!(row.names(dat) %in% row.names(dat.test)),]
    cv.rf = randomForest(Class~., data=dat.train, ntree=10001, type='prob')
    predi = predict(cv.rf, dat.test, type='prob')
    pred = rbind(pred, predi)
  }
  pred = pred[-1,]
  class = class[-1]
  pred.to.roc = pred[,2]
  pred.rocr = prediction(pred.to.roc, as.factor(class))
  perf.rocr = performance(pred.rocr, measure='auc', x.measure='cutoff')
  AUC = deparse(round(as.numeric(perf.rocr@y.values),3))
  perf.tpr.rocr = performance(pred.rocr,'tpr','fpr')
  perf.precrec.rocr = performance(pred.rocr,'prec','rec')
  perf.sens.rocr = performance(pred.rocr,'tpr','tnr')
  res = list()
  res[[1]] = AUC
  res[[2]] = perf.tpr.rocr
  res[[3]] = perf.precrec.rocr
  res[[4]] = perf.sens.rocr
  names(res) = c('AUC','ROC','PrecRec','SensSpe')
  return(res)
}


#---------------------------------------------
# Workflow
#---------------------------------------------

### Read in datasets
c.dat = read.delim(file='../data/Chari.txt',         sep="", header=T, row.names=1)   # Indel Dataset: used for initial feature selection
r.dat = read.delim(file='../data/Resistance.txt',    sep="", header=T)                # RES Dataset: used to train final model  
f.dat = read.delim(file='../data/FlowCytometry.txt', sep="", header=T)                # FC Dataset: used to train final model

### Process data
c.dat$Class = 0
for (i in c(1:nrow(c.dat))) { if (c.dat[i,'Activity'] > 1) { c.dat[i,'Class'] = 1 } }
c.dat$Class = as.factor(c.dat$Class)
c.dat$Activity = NULL

r.ord = r.dat[order(r.dat$Activity),]
r.dat = rbind(r.ord[1:(nrow(r.ord)*0.25),],r.ord[(0.75*nrow(r.ord)):nrow(r.ord),])
r.dat$Class = 0
for (i in c(1:nrow(r.dat))) { if (r.dat[i,'Activity'] > -2 ) { r.dat[i,'Class'] = 1 } }
r.dat$Activity = NULL

f.ord = f.dat[order(f.dat$Activity),]
f.dat = rbind(r.ord[1:(nrow(f.ord)*0.25),],f.ord[(0.75*nrow(f.ord)):nrow(f.ord),])
f.dat$Class = 0
for (i in c(1:nrow(f.dat))) { if (f.dat[i,'Activity'] > -2 ) { f.dat[i,'Class'] = 1 } }
f.dat$Activity = NULL

com = rbind(c.dat,f.dat,r.dat) # Combined dataset used to train final model

### Perform feature selection on chari dataset
chari.full.rf  = randomForest(Class ~ ., data=c.dat, ntree=10001, importance=T)
cl.imp         = as.data.frame(varImpPlot(chari.full.rf))
cl.impOrd      = cl.imp[order(-cl.imp$MeanDecreaseGini),]
cl.impOrd$Rank = c(1:nrow(cl.impOrd))

cl.error       = seq_rf(c.dat, cl.impOrd, 3) # Construct iterative models, adding progressively more features in order of importance

cl.min  = as.numeric(row.names(subset(cl.error, OOB<=min(cl.error$OOB))))[1] # Extract smallest model with lowest OOB error
cl.vars = row.names(subset(cl.impOrd, Rank<=cl.min))

### Perform 10-fold Cross Validation
co.cvRes = rf_cvTest(subset(com,   select=c(cl.vars,'Class')),10,cl.min)

plot(co.cvRes[[2]], col='blue', lwd=3, main='10-fold Cross Validation\nROC curve')
plot(co.cvRes[[3]], col='blue', lwd=3, ylim=c(0,1), main='10-fold Cross Validation\nPrecision Recall curve')

### Build final model using full combined dataset
combined.rf = randomForest(Class ~ ., data = com, importance = T, type = 'classification', ntree = 10001) 

