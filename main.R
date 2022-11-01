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
# source('meth.R')

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

train = pheno[plate %in% plates & rep == 1]$j
test  = setdiff(pheno[rep==1]$j,train)

pheno[train,split:='train']
pheno[test ,split:='test' ]

cat('# samples in training set:',length(train),'\n')
cat('# samples in test     set:',length(test ),'\n')

tmp = melt(pheno,id.vars=c('FID','split'),measure.vars=c('sex','race','hispanic'))
dcast(tmp,variable+value ~ split,fun=length)


rm(plates)

# ---------------------------------------------------
# Model training

source('train.R')
write.table(coefs,file='intermediate/HRS.txt')

# ---------------------------------------------------
# Model evaluation

est_LC = function(beta,coefs){
 
	ib = match(rownames(coefs),names(beta))
	ic = !is.na(ib) & !is.na(beta[ib])
	ib = ib[ic]

	n_celltypes = ncol(coefs)

	props = solve.QP(
		 t(coefs[ic,]) %*% coefs[ic,]
		,t(coefs[ic,]) %*% beta[ib]
		,diag(n_celltypes)
		,rep(0,times=n_celltypes)
		)$sol

	names(props) = colnames(coefs)
	props
}

coefs = read.table('intermediate/HRS.txt',header=TRUE,)
coefs %<>% as.matrix

reinius =
	system.file('data/Reinius.txt',package='ewastools') %>%
	read.table %>%
	as.matrix

table(rownames(reinius) %in% common)
# FALSE  TRUE 
#     3   597 

reinius = reinius[rownames(reinius) %in% common,]

source('eval.R')
source('ewas.R')
