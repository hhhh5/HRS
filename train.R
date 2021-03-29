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
