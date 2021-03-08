# Import methylation data and perform quality checks

pheno = '/nfs/turbo/bakulski1/Datasets/HRS/jonheiss/sensitive/flow/expanded_pheno.rds'
pheno %<>% readRDS %>% as.data.table

# 4018 subjects
pheno = pheno[,.(
	 FID
	,plate = Sample_Plate # {Thyagarajan_Sample_127}
	,file = paste0('/nfs/turbo/bakulski1/Datasets/HRS/jonheiss/sensitive/idats/',Slide,'_',Array)
	,reported_sex = Reported
	,age = PAGE  # patient age  
	,BIRTHMO     # ‘%Y’ Year with century.
	,BIRTHYR     # ‘%m’ Month as decimal number (01-12).
	,COLLECTDATE # {%Y-%m-%d %H:%M:%S}
	)]

## --------
## Are all .idat files available?
table(file.exists(pheno$file %s+% '_Grn.idat'))
table(file.exists(pheno$file %s+% '_Red.idat'))

## --------
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

## --------
# Read in .idat files in chunks

common = intersect(
	 ewastools:::manifest_450K[probe_type!="rs"]$probe_id
	,ewastools:::manifest_epic[probe_type!="rs"]$probe_id
	)

snps = ewastools:::manifest_epic[probe_type=="rs"]$probe_id
autosomal = ! ewastools:::manifest_epic$chr %in% c('chrX','chrY')

# Store matrix of beta-values on hard drive
meth = ff(initdata=NA_real_,dim =c(length(common),pheno[,.N]))
rownames(meth) = common

pheno[,j:=.I]
chunks = split(pheno,pheno$j %/% 200)

f = function(pheno){

	cat('sdf')

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

