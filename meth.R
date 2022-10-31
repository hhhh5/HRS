# Import methylation data and perform quality checks

pheno = '/nfs/turbo/bakulski1/Datasets/HRS/jonheiss/sensitive/flow/expanded_pheno.rds'
pheno %<>% readRDS %>% as.data.table

# 4018 subjects
pheno = pheno[,.(
	 FID
	,plate = Sample_Plate # {Thyagarajan_Sample_127}
	,file  = paste0('/nfs/turbo/bakulski1/Datasets/HRS/jonheiss/sensitive/idats/',Slide,'_',Array)
	,sex   = ifelse(Reported == 'M','m','f')
	,age   = PAGE  # patient age
	,race  = RACE 
	,BIRTHMO     # ‘%Y’ Year with century.
	,BIRTHYR     # ‘%m’ Month as decimal number (01-12).
	,COLLECTDATE # {%Y-%m-%d %H:%M:%S}
	)]

# 39 subjects with technical replicates
dupes = '/nfs/turbo/bakulski1/Datasets/HRS/jonheiss/sensitive/flow/methyl_dupes_pheno.rds'
dupes %<>% readRDS %>% data.table

dupes = dupes[,.(
	 FID   = Original.Sample.ID        # {F0000000}
	,plate = dupe_plate                # {Thyagarajan_Sample_190}
	,file  = paste0('/nfs/turbo/bakulski1/Datasets/HRS/jonheiss/sensitive/idats/',dupe_methid) # {202711530006_R04C01}
	,sex   = ifelse(GENDER==1,'m','f') # {1/2}
	,age   = PAGE                    
	,race  = RACE                      # {1/2/7}
	# ,blindedID = Blind.Dup.ID        # {F0000000}
	# ,orig_methid                     # {202163110034_R05C01}
	# ,HHID                            # Household Identifier {000000}
	# ,PN                              # PERSON NUMBER {000}
	# ,Pair.ID                         # {HRSDNA11}
	# ,methylation.plates              # {onto T13}
	# ,X                               # {12-H9}
	# ,Bharat.s.Lab.plates             # {dup 11}
	# ,X.1                             # {onto T31}
	# ,orig_plate                      # {Thyagarajan_Sample_111}
	# ,HISPANIC                        # {1/2/5}
	# ,ethrace                         # {h/nhb/nhw}
	)]

# Mark original and duplicate samples
pheno$rep = 1L
dupes$rep = 2L

pheno = rbind(pheno,dupes,use.names=TRUE,fill=TRUE); rm(dupes)

pheno[,race:=factor(race,levels=c(1,2,7,0),labels=c('white','black','other','unknown'))]

pheno = pheno[FID %in% LC$FID]

setkey(pheno,'FID')

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

# Reduce to Cpg sites common to EPIC and 450K. Also necessary because full set of sites too large (more than Machine$integer.max)
common = intersect(
	 ewastools:::manifest_450K[probe_type!='rs' & !chr %in% c('X','Y')]$probe_id
	,ewastools:::manifest_epic[probe_type!='rs' & !chr %in% c('X','Y')]$probe_id
	)

snps = ewastools:::manifest_epic[probe_type=='rs']$probe_id
autosomal = ! ewastools:::manifest_epic$chr %in% c('X','Y')

pheno[,j:=.I]
chunks = split(pheno,pheno$j %/% 200)

# Store matrix of beta-values on hard drive
meth = ff(
	 initdata=NA_real_
	,dim=c(length(common),nrow(pheno)))



f = function(pheno){

	tmp = 
		pheno$file %>%
		read_idats(quiet=TRUE) %>%
		detectionP %>%
		correct_dye_bias

	pheno[,c('X','Y'):=check_sex(tmp)]

	pheno$undetected = colSums(tmp$detP[autosomal,] > 0.05,na.rm=TRUE)
	
	metrics =
		tmp %>%
		control_metrics %>%
		as.data.table

	tmp =
		mask(tmp,threshold=0.05) %>%
		dont_normalize

	# Epigenetic age
	pheno$horvath = ewastools:::methylation_score(tmp,model='horvath_clock')

	meth[,pheno$j] <<- tmp[common,]

	chunk = list(
		 pheno = pheno
		,metrics = metrics
		,snps = tmp[snps,]
		)

	chunk
}

x = lapply(chunks,f)

snps    = x %>% map('snps')    %>% do.call('cbind',.)
metrics = x %>% map('metrics') %>% rbindlist
pheno   = x %>% map('pheno')   %>% rbindlist

rm(x,chunks,f)

## ---------------------------------------------------
## Quality control

pheno[,failed:=FALSE]


## -----------
## Sex check
pheno[,predicted_sex:=predict_sex(X,Y,which(sex=='m'),which(sex=='f'))]
table(pheno$sex == pheno$predicted_sex)                           
# FALSE  TRUE                                            
#     1  2908  

png('intermediate/qc1.png')
tmp = pheno[sex==predicted_sex]
plot  (Y ~ X,data=tmp,pch=ifelse(tmp$sex=='f',1,4),asp=1,xlab='Normalized X chromosome intensities',ylab='Normalized Y chromosome intensities')
tmp = pheno[sex!=predicted_sex]
points(Y ~ X,data=tmp,pch=ifelse(tmp$sex=='f',1,4),col=2)
legend('topright',pch=c(1,4),legend=c('female','male'))
dev.off()

pheno[sex!=predicted_sex,failed:=TRUE]

## -----------
## Undetected probes
table(pheno$undetected > 1e5)
# FALSE  TRUE 
#  2571   351 

pheno[undetected > 1e5,failed:=TRUE]

## -----------
## SNP outliers
tmp = call_genotypes(snps)
pheno$snp_outlier = snp_outliers(tmp)
table(pheno$snp_outlier > -3)
# FALSE  TRUE 
#  2846    76 

pheno[snp_outlier > -3,failed:=TRUE]

png('intermediate/snps.png')
stripchart(pheno$snp_outlier,m='j')
abline(v=-3,lty=3,col=2)
dev.off()

## -----------
# Do the genotypes match for dupes

pheno$genotype = enumerate_sample_donors(tmp)

tmp = dcast(pheno,FID ~ rep,value.var='genotype')
table(tmp[,`1` == `2`])

## -----------
## Control probe metrics
i = sample_failure(metrics)
table(i)
# FALSE  TRUE                            
#  2778   144

pheno[i,failed:=TRUE]

tmp = copy(metrics)
tmp$plate = pheno$plate
tmp = melt(tmp,id.var='plate')

p = ( ggplot(tmp)
+ geom_boxplot(aes(x=plate,y=value))
+ facet_wrap(facets=~variable,ncol=2,scales='free_y')
)

ggsave(p,file='intermediate/metrics.pdf')

## ---------------------------------------------------
# Subsetting

keep = which(pheno$failed==FALSE)

pheno   = pheno  [ keep]
meth    = meth   [,keep]
snps    = snps   [,keep]
metrics = metrics[ keep]

# update column indices
pheno[,j:=.I]

# drop probes for which a third of measurements are missing
i = rowSums(is.na(meth))
i = which(i < ncol(meth)/3)
meth   = meth  [i,]
common = common[i]

rownames(meth) = common

# Store matrix of beta-values on hard drive
meth = as.ff(meth)
ffsave(meth,file='intermediate/meth')

save(pheno,snps,metrics,autosomal,common,file='intermediate/processed.rda')

rm(tmp,i,p,keep)
