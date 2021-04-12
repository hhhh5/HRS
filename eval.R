LC0 = pheno[,..celltypes]

LC1 = map(pheno$j,~ est_LC(meth[,.x],coefs  ))
LC1 %<>% do.call('rbind',.) %>% data.table

raw_correlations = sapply(celltypes,function(ct){
	cor(LC0[[ct]][test],LC1[[ct]][test])
})

# ---------------------------------------------------
# Estimates based on Reinius

LC2 = map(pheno$j,~ est_LC(meth[,.x],reinius)) 
LC2 %<>% do.call('rbind',.) %>% data.table

# ---------------------------------------------------
# Principal components of control metrics

pcs = t(as.matrix(metrics))
tmp = rowMeans(pcs,na.rm=TRUE)
pcs = pcs - tmp
pcs[is.na(pcs)] = 0
pcs = trlan.svd(t(pcs),neig=5)
pcs = pcs$u

pheno$pc1 = pcs[,1]
pheno$pc2 = pcs[,2]
pheno$pc3 = pcs[,3]
pheno$pc4 = pcs[,4]
pheno$pc5 = pcs[,5]

rm(tmp,pcs)

# ---------------------------------------------------
# How much variance explain the various estimates

exp_var = function(covars){
	frml = paste0("meth~",paste0(names(covars),collapse="+"))
	frml = formula(frml)
	covars[,meth:=runif(.N)]

	m = lm(frml,covars[test])
	mm = model.matrix(m)
	
	f = function(i){
		y = meth[i,test]
		i = !is.na(y)
		m[1:8] = lm.fit(mm[i,],y[i])
		summary(m)$r.squared
	}

	possibly(f,otherwise=NA_real_)
}

features = setdiff(common,rownames(coefs))

BASELINE = pheno[,.(sex,age,plate,pc1,pc2,pc3,pc4,pc5)]


EXPV = list()

f = exp_var(copy(BASELINE))
EXPV$REF = unlist(mclapply(features,f,mc.cores=10))

f = exp_var(cbind(BASELINE,LC0))
EXPV$LC0 = unlist(mclapply(features,f,mc.cores=10))

f = exp_var(cbind(BASELINE,LC1))
EXPV$LC1 = unlist(mclapply(features,f,mc.cores=10))

f = exp_var(cbind(BASELINE,LC2))
EXPV$LC2 = unlist(mclapply(features,f,mc.cores=10))

save(pheno,train,test,LC0,LC1,LC2,EXPV,file="intermediate/results.rda")

# ---------------------------------------------------
# RMSE

tmp0 = cbind(pheno[,.(FID)],LC0)[test]
tmp1 = cbind(pheno[,.(FID)],LC1)[test]

tmp0 = melt(tmp0,id.vars='FID',variable.name='cell_type',value.name='prop0')
tmp1 = melt(tmp1,id.vars='FID',variable.name='cell_type',value.name='prop1')

tmp = tmp0[tmp1,on=c('FID','cell_type')]

tmp[,.(rsme=round(sqrt(mean((prop1-prop0)^2))*100,1)),cell_type]

rm(tmp,tmp0,tmp1)

# ---------------------------------------------------

mean(EXPV$LC0 - EXPV$REF)
mean(EXPV$LC1 - EXPV$REF)
mean(EXPV$LC2 - EXPV$REF)

gains = data.table(
	 probe_id = features
	,measured  = EXPV$LC0-EXPV$REF
	,estimated = EXPV$LC1-EXPV$REF
	)

gains = gains[order(-measured)][1:10000]
gains[,x:=.I]

p =
( ggplot(gains)
+ geom_point(aes(x=x,y=estimated),col='red',alpha=0.4,shape=20,size=0.5,show.legend=FALSE)
+ geom_line (aes(x=x,y=measured),show.legend=FALSE)
+ theme_linedraw()
+ xlab('Index') + ylab('Gain in explained variance')
)

ggsave(p,file='gains.png',width=5,height=5)

