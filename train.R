library(haven)
library(readxl)
library(data.table)
library(magrittr)
library(stringi)
library(ewastools)
library(ppcor)
library(parallel)
library(svd)
library(purrr)
library(robustbase)
library(quadprog)
library(ggplot2)
library(ff)


source('cytometry.R')
source('meth.R')

ffload('intermediate/meth')
  load('intermediate/meth.RData')
  load('intermediate/processed.rda')

# ---------------------------------------------------

pheno = LC[pheno,on='FID']

# ---------------------------------------------------
# Split in test/train

plates = as.character(98:129)
plates = stri_pad(plates,width=3,pad='0')
plates = paste0('Thyagarajan_Sample_',plates)

train = pheno[plate %in% plates]$j
test  = setdiff(pheno$j,train)

cat('# samples in training set:',length(train),'\n')
cat('# samples in test     set:',length(test ),'\n')

rm(plates)

# ---------------------------------------------------
# Model training

cat('Model training\n')

# do not include NE among the covariates
covars = setdiff(celltypes,'NE')

pheno$meth = 0
frml = paste0('meth~1+',paste0(covars,collapse='+')) %>% formula
mm = model.matrix(frml,data=pheno[train])

# Compute partial correlation coefficients
PC = mclapply(1:nrow(meth),function(i){
		y = meth[i,train]
		x = copy(mm)
		x[,1] = y
		x = na.omit(x)
		pcor(x,method='spearman')$estimate[-1,1]
		},mc.cores=10)

PC = do.call('rbind',PC)

# Select for each cell type the 50 markers with the highest absolute correlations
markers = apply(abs(PC),2,order,decreasing=TRUE)
markers = markers[1:50,]
markers = unique(markers[TRUE])
markers = rownames(meth)[markers]

# Learn the regression coefficients
coefs = apply(meth[markers,],1,function(y){
	pheno$meth = y
	lmrob(frml,pheno[train])$coefficients
})

coefs = t(coefs)
coefs[,-1] = coefs[,-1]+coefs[,1]
colnames(coefs)[1] = 'NE'
rownames(coefs) = markers

write.table(coefs,file='intermediate/HRS.txt')
