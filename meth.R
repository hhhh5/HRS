# Import methylation data and perform quality checks

pheno = '/nfs/turbo/bakulski1/Datasets/HRS/jonheiss/sensitive/flow/expanded_pheno.rds'
pheno %<>% readRDS %>% as.data.table

# 4018 subjects
pheno = pheno[,.(
	 FID
	,plate = Sample_Plate # {Thyagarajan_Sample_127}
	,file = paste0('/nfs/turbo/bakulski1/Datasets/HRS/jonheiss/sensitive/idats/',Slide,'_',Array)
	,sex = ifelse(Reported == 'M','m','f')
	,age = PAGE  # patient age  
	,BIRTHMO     # ‘%Y’ Year with century.
	,BIRTHYR     # ‘%m’ Month as decimal number (01-12).
	,COLLECTDATE # {%Y-%m-%d %H:%M:%S}
	)]

## ---------------------------------------------------
## Are all .idat files available?
table(file.exists(pheno$file %s+% '_Grn.idat'))
table(file.exists(pheno$file %s+% '_Red.idat'))

##---------------------------------------------------
## Check the age variable for consistency with collection date/birthday

pheno$COLLECTDATE %<>%
	as.IDate(format='%Y-%m-%d %H:%M:%S') %>%
	round(unit='day')

pheno[,birthday:=paste0(BIRTHYR,'-',BIRTHMO,'-01')]
pheno[,birthday:=as.IDate(birthday,format='%Y-%m-%d')]

pheno[,calculated_age:=COLLECTDATE-birthday]
pheno[,calculated_age:=round(calculated_age/365)]

# No samples where the age deviates by more than 2 years
pheno[abs(age-calculated_age) > 2]

pheno[,c('calculated_age','BIRTHMO','BIRTHYR','COLLECTDATE','birthday'):=NULL]

## ---------------------------------------------------
## Read in .idat files in chunks

common = intersect(
	 ewastools:::manifest_450K[probe_type!="rs"]$probe_id
	,ewastools:::manifest_epic[probe_type!="rs"]$probe_id
	)

snps = ewastools:::manifest_epic[probe_type=="rs"]$probe_id
autosomal = ! ewastools:::manifest_epic$chr %in% c('chrX','chrY')

meth = matrix(NA_real_,nrow=length(common),ncol=pheno[,.N])
rownames(meth) = common

pheno[,j:=.I]
chunks = split(pheno,pheno$j %/% 200)

f = function(pheno){

	tmp = 
		pheno$file %>%
		read_idats(quiet=TRUE) %>%
		detectionP %>%
		correct_dye_bias

	pheno[,c("X","Y"):=check_sex(tmp)]

	pheno$undetected = colSums(tmp$detP[autosomal,] > 0.05,na.rm=TRUE)
	
	metrics =
		tmp %>%
		control_metrics %>%
		as.data.table

	tmp =
		mask(tmp,threshold=0.05) %>%
		dont_normalize

	meth[,pheno$j] <- tmp[common,]

	# Epigenetic age
	pheno$horvath = ewastools:::methylation_score(tmp,model="horvath_clock")

	chunk = list(
		 pheno = pheno
		,metrics = metrics
		,snps = tmp[snps,]
		)

	rm(tmp); gc()

	chunk
}

chunks %<>% map(f)

snps    = chunks %>% map('snps')    %>% do.call('cbind',.)
metrics = chunks %>% map('metrics') %>% rbindlist
pheno   = chunks %>% map('pheno')   %>% rbindlist

## ---------------------------------------------------
## Quality control

pheno[,failed:=FALSE]


## -----------
## Sex check
pheno[,predicted_sex:=predict_sex(X,Y,which(sex=="m"),which(sex=="f"))]
table(pheno$sex == pheno$predicted_sex)
# FALSE  TRUE 
#     1  3994 

png('intermediate/qc1.png')
tmp = pheno[sex==predicted_sex]
plot  (Y ~ X,data=tmp,pch=ifelse(tmp$sex=="f",1,4),asp=1,xlab="Normalized X chromosome intensities",ylab="Normalized Y chromosome intensities")
tmp = pheno[sex!=predicted_sex]
points(Y ~ X,data=tmp,pch=ifelse(tmp$sex=="f",1,4),col=2)
legend("topright",pch=c(1,4),legend=c("female","male"))
dev.off()

pheno[sex!=predicted_sex,failed:=TRUE]

## -----------
## Undetected probes
table(pheno$undetected > 1e5)
# FALSE  TRUE 
#  3534   484
pheno[undetected > 1e5,failed:=TRUE]

## -----------
## SNP outliers
tmp = call_genotypes(snps)
pheno$snp_outlier = snp_outliers(tmp)
table(pheno$snp_outlier > -3)
# FALSE  TRUE 
#  3910   108 

pheno[snp_outlier > -3,failed:=TRUE]

png("intermediate/snps.png")
stripchart(pheno$snp_outlier,m="j")
abline(v=-3,lty=3,col=2)
dev.off()


## -----------
## Control probe metrics
i = sample_failure(metrics)
table(i)
# FALSE  TRUE 
#  3836   182 
pheno[i,failed:=TRUE]


tmp = copy(metrics)
tmp$plate = pheno$plate
tmp = melt(tmp,id.var='plate')

p = ( ggplot(tmp)
+ geom_boxplot(aes(x=plate,y=value))
+ facet_wrap(facets=~variable,ncol=2,scales="free_y")
)

ggsave(p,file="intermediate/metrics.pdf")

## ---------------------------------------------------

keep = which(pheno$failed==FALSE)

pheno   = pheno  [ keep]
meth    = meth   [,keep]
snps    = snps   [,keep]
metrics = metrics[ keep]

save(pheno,snps,metrics,autosomal,common,file='intermediate/processed.rda')

# Store matrix of beta-values on hard drive
meth = as.ff(meth,filename='intermediate/meth.ff')

rm(tmp,f,i,p,chunks,keep)
