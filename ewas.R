# External dataset (450K)

# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/suppl/GSE87571_RAW.tar
# tar -xf GSE87571_RAW.tar

meta = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/matrix/GSE87571_series_matrix.txt.gz'
meta = fread(meta,sep='\t',skip=29,nrows=35,header=FALSE)

setnames(meta,'V1','variable')
meta = meta[,variable:=stri_sub(variable,9,-1)]
meta = meta[variable %in% c('geo_accession','characteristics_ch1','supplementary_file')]

meta = data.table(t(meta[,-1]))
names(meta) = c('gsm','sex','age','tissue','disease_state','grn','red')

meta[sex == 'gender: Male'  ,sex:='m']
meta[sex == 'gender: Female',sex:='f']
meta[sex == 'gender: NA'    ,sex:=NA]
meta[,age:=as.integer(stri_sub(age,6,-1))]
meta[,idat:=stri_sub(grn,68,-13)]
meta[,idat:=paste0('GSE87571/',idat)]
meta[,c('tissue','disease_state','grn','red'):=NULL]
meta = na.omit(meta)

meth = 
	meta$idat %>%
	read_idats(quiet=FALSE) %>%
	detectionP %>%
	correct_dye_bias %>%
	mask(0.05)

metrics =
	meth %>%
	control_metrics %>%
	as.data.table

meth %<>% dont_normalize

# ---------------------------------------------------
# Principal components of control metrics

pcs = t(as.matrix(metrics))
tmp = rowMeans(pcs,na.rm=TRUE)
pcs = pcs - tmp
pcs[is.na(pcs)] = 0
pcs = trlan.svd(t(pcs),neig=5)
pcs = pcs$u

meta$pc1 = pcs[,1]
meta$pc2 = pcs[,2]
meta$pc3 = pcs[,3]
meta$pc4 = pcs[,4]
meta$pc5 = pcs[,5]

rm(metrics,pcs)

# ---------------------------------------------------
# Estimate cell composition

LC1 = map(1:ncol(meth), ~ est_LC(meth[,.x],coefs) )
LC1 %<>% do.call('rbind',.) %>% as.data.table

LC2 = map(1:ncol(meth), ~ est_LC(meth[,.x],reinius) )
LC2 %<>% do.call('rbind',.) %>% as.data.table

# ---------------------------------------------------
# How much variance explain the various estimates + EWAS of `age`

assoc = function(covars){
	frml = paste0("meth~",paste0(names(covars),collapse="+"))
	frml = formula(frml)
	covars[,meth:=runif(.N)]

	m = lm(frml,covars)
	mm = model.matrix(m)
	
	f = function(i){
		y = meth[i,]
		i = !is.na(y)
		m[1:8] = lm.fit(mm[i,],y[i])
		tmp = summary(m)
		c(tmp$r.squared,tmp$coefficients['age','Pr(>|t|)'])
	}

	possibly(f,otherwise=c(NA_real_,NA_real_))
}

features = setdiff(rownames(meth),rownames(coefs))

BASELINE = meta[,.(sex,age,pc1,pc2,pc3,pc4,pc5)]

EXPV = list()
EWAS = list()

f = assoc(copy(BASELINE))
tmp = mclapply(features,f,mc.cores=10)
EXPV$REF = map_dbl(tmp,1)
EWAS$REF = map_dbl(tmp,2)

f = assoc(cbind(BASELINE,LC1))
tmp = mclapply(features,f,mc.cores=10)
EXPV$LC1 = map_dbl(tmp,1)
EWAS$LC1 = map_dbl(tmp,2)

f = assoc(cbind(BASELINE,LC2))
tmp = mclapply(features,f,mc.cores=10)
EXPV$LC2 = map_dbl(tmp,1)
EWAS$LC2 = map_dbl(tmp,2)

mean(EXPV$LC1 - EXPV$REF)
mean(EXPV$LC2 - EXPV$REF)

length(features)
sum(p.adjust(EWAS$REF,m='b') < 0.05,na.rm=TRUE)
sum(p.adjust(EWAS$LC1,m='b') < 0.05,na.rm=TRUE)
sum(p.adjust(EWAS$LC2,m='b') < 0.05,na.rm=TRUE)

save(meta,LC1,LC2,EXPV,EWAS,file="intermediate/results_geo.rda")
