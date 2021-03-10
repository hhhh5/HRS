# Flow cytometry measurements

flow = "/nfs/turbo/bakulski1/Datasets/HRS/jonheiss/sensitive/flow/hrsflowdata2016.sas7bdat"
vbs  = "/nfs/turbo/bakulski1/Datasets/HRS/jonheiss/sensitive/flow/hrs2016vbs.sas7bdat"

flow %<>% read_sas %>% as.data.table
vbs  %<>% read_sas %>% as.data.table

# Make sure that the combination of `HHID` and `PN` is unique
# Result should be `0`
anyDuplicated(flow,by=c('HHID','PN'))
anyDuplicated(vbs ,by=c('HHID','PN'))

# Inner join
LC = flow[vbs,on=c('HHID','PN'),nomatch=NULL];

# 9933 subjects
LC[,.N]

# Select the relevant variables
LC = LC[,.(

	 HHID                           # Household Identifier
	,PN                             # PERSON NUMBER

	# ,PANEU                          # NEUTROPHIL COUNT - X10E9/L
	# ,PAEOS                          # EOSINOPHIL COUNT - X10E9/L
	# ,PABAS	                      # BASOPHIL   COUNT - X10E9/L
	# ,PALYM                          # LYMPHOCYTE COUNT - X10E9/L
	# ,PAMON                          # MONOCYTE   COUNT - X10E9/L

	,PWBC                           # WHITE BLOOD CELL COUNT - 10E9/L
		
	# of the following 5 cell types I won't use lymphocytes but the more granular subtypes from `flow`
	,NE = PNEUT / 100              # PERCENT NEUTROPHILS
	,EO = PEOS  / 100              # PERCENT EOSINOPHILS
	,BA = PBASO / 100              # PERCENT BASOPHILSq
	,LY = PLYMP / 100              # PERCENT LYMPHOCYTES
	,MO = PMONO / 100              # PERCENT MONOCYTES

	,T     = PTcell_count / PWBC   # T cells

	,CT    =     PCT_count / PWBC  # cytotoxic T cells
	,EM_CT =  PEM_CT_count / PWBC  # Effector memory cytotoxic T cells
	,CM_CT =  PCM_CT_count / PWBC  # Central memory cytotoxic T cells
	,E_CT  =   PE_CT_count / PWBC  # Effector cytotoxic T cells
	,N_CT  =   PN_CT_count / PWBC  # Naive cytotoxic T cells
	
	,HT    =     PHT_count / PWBC  # Helper T cells
	,EM_HT =  PEM_HT_count / PWBC  # Effector memory helper T cells !!!!!! Should I use the new measurements here?
	,CM_HT =  PCM_HT_count / PWBC  # Central memory helper T cells
	,E_HT  =   PE_HT_count / PWBC  # Effector Memory
	,N_HT  =   PN_HT_count / PWBC  # Naive helper T cells 

	,B              =          PBcell_count / PWBC  # B lymphocytes
	,NaiveB         =         PNaiveB_count / PWBC  # CD27- naive B cells
	,IgD_Plus_MemB  =  PIgD_Plus_MemB_count / PWBC  # IgD+ Memory B cells
	,IgD_Minus_MemB = PIgD_Minus_MemB_count / PWBC  # IgD- Memory B cells

	,MON   =   PMONO_count / PWBC  # Monocytes
	,MONC  =  PMONOc_count / PWBC  # classical monocytes
	,MONNC = PMONOnc_count / PWBC  # non-classical monocytes

	,DC    =     PDC_count / PWBC  # Dendritic cells
	,DCm   =    PDCm_count / PWBC  # CD11c+ myeloid dendritic cells
	,DCp   =    PDCp_count / PWBC  # CD123+ plasmacytoid dendritic cells

	,NK    =     PNK_count / PWBC  # Natural killer cells
	,NKHI  =   PNKHI_count / PWBC  # CD56 high NK cells
	,NKLO  =   PNKLO_count / PWBC  # CD56 low NK cells

   	# Percentages for the same blood cell types

	# ,PDC_pct             
	# ,PDCm_pct
	# ,PDCp_pct
	# ,PNK_pct            
	# ,PNKHI_pct
	# ,PNKLO_pct
	# ,PMONO_pct           
	# ,PMONOc_pct
	# ,PMONOnc_pct    
	# ,PBcell_pct         
	# ,PCM_CT_pct
	# ,PCM_HT_pct
	# ,PCT_pct           
	# ,PEM_CT_pct
	# ,PEM_HT_pct
	# ,PE_CT_pct        
	# ,PE_HT_pct
	# ,PHT_pct
	# ,PIgD_Plus_MemB_pct
	# ,PIgD_Minus_MemB_pct
	# ,PN_CT_pct
	# ,PN_HT_pct          
	# ,PTcell_pct
	# ,PNaiveB_pct      

	## Irrelevant variables

	# ,PALB                         # ALBUMIN - G/DLLB
	# ,PALKP2                       # ALKALINE PHOSPHATASE - U/L
	# ,PALT                         # ALANINE AMINOTRANSFERASE - U/L
	# ,PAST                         # ASPARTATE AMINOTRANSFERASE - U/L
	# ,PBILT                        # BILIRUBIN, TOTAL - MG/DL
	# ,PBUN                         # UREA NITROGEN (BUN) - MG/DL
	# ,PCA                          # CALCIUM - MG/DL
	# ,PCHOL                        # CHOLESTEROL, TOTAL - MG/DL
	# ,PCL                          # CHLORIDE - MMOL/L
	# ,PCMVGE                       # CMV IGG - COI
	# ,PCO2                         # BICARBONATE (CO2) - MMOL/L
	# ,PCR                          # CREATININE - MG/DL
	# ,PCRP                         # C-REACTIVE PROTEIN (HIGH SENSITIVITY) - MG/L
	# ,PCYSC                        # CYSTATIN C -  MG/L
	# ,PDHEASE                      # DEHYDROEPIANDROSTERONE SULFATE (DHEAS) - UMOL/L         
	# ,PFERTN                       # FERRITIN - UG/L
	# ,PGLUFF                       # GLUCOSE, FASTING - MG/DL
	# ,PHCT                         # HEMATOCRIT      
	# ,PHDLD                        # HDL-CHOLESTEROL, DIRECT-MEASURE - MG/DL
	# ,PHGB                         # HEMOGLOBIN - G/DL
	# ,PK                           # POTASSIUM - MMOL/L
	# ,PLDLC                        # LDL-CHOLESTEROL, CALCULATED - MG/DL
	# ,PMCH                         # MEAN CORPUSCULAR HEMOGLOBIN - PG
	# ,PMCHC                        # MEAN CORPUSCULAR HEMOGLOBIN CONCENTRATION - G/DL
	# ,PMCV                         # MEAN CORPUSCULAR VOLUME - FL
 	# ,PMPV                         # MEAN PLATELET VOLUME - FL
	# ,PNA                          # SODIUM - MMOL/L
	# ,PNTBNPE                      # B-TYPE NATRIURETIC PEPTIDE,N-TERMINAL PRO (NT-PROBNP)-PG/ML
	# ,PPDW                         # PLATELET DISTRIBUTION WIDTH - FL
	# ,PPLT                         # PLATELET COUNT - 10E9/L
	# ,PRBC                         # RED BLOOD CELL COUNT - 10E12/L
	# ,PRDW                         # RED CELL DISTRIBUTION WIDTH - %
	# ,PTGF                         # RIGLYCERIDES - MG/DL  
	# ,PTP                          # PROTEIN, TOTAL - G/DL
	# ,PVBSWGTR                     # VBS RESPONDENT WEIGHT      
	# ,PVBSWHY0WGT                  # VBS WHY ZERO WEIGHT
	# ,PVBS_N_DAYS                  # NUMBER OF DAYS BETWEEN HRS IW DATE AND COLLECTION DATE
	# ,PCMVGINT                     # CMV IGG (NON-REACTIVE, REACTIVE, BORDERLINE)    
	# ,PFASTYN                      # FASTING STATUS (Y/N)

	# ,VERSION                      # DATASET VERSION
	)]


# These are the cell types we'd like to estimate
celltypes = c("NE","EO","BA","MO","EM_CT","CM_CT","E_CT","N_CT","EM_HT","CM_HT","E_HT","N_HT","B","DC","NK")

# Drop rows with missing values for the required cell types
i = is.na(LC[,..celltypes])
i = apply(i,1,sum)
i = i == 0
LC = LC[i]
LC[,.N]

# Cell types should sum up to approximately 1 for each sample
# Drop samples for which this does not hold true
LC$total = rowSums(as.matrix(LC[,..celltypes]))
LC = LC[total %between% c(0.9,1.02)]
LC[,total:=NULL]

## Cleanup
rm(flow,vbs,i)