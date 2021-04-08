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

train = pheno[plate %in% plates]$j
test  = setdiff(pheno$j,train)

cat('# samples in training set:',length(train),'\n')
cat('# samples in test     set:',length(test ),'\n')

rm(plates)

# ---------------------------------------------------
# Model training

source('train.R')
write.table(coefs,file='intermediate/HRS.txt')

# ---------------------------------------------------
# Model evaluation

coefs = read.table('intermediate/HRS.txt',header=TRUE,)
coefs %<>% as.matrix

source('evalR')
