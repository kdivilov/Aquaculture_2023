library(sommer)

A = as.matrix(read.csv("A.csv",row.names = 1,check.names = F))
pheno = read.csv("phenotypes.csv")
pheno$Family = factor(pheno$Family)
pheno$Year = factor(pheno$Year)

#Fit the linear mixed model
modelfit = mmer(Status~Year,
                random=~vsr(Family,Gu=A),
                rcov=~units,
                data=pheno)

#Estimate and SE of heritability on the observed scale
vpredict(modelfit,h_o~V1/(V1+V2))

#Observed to underlying liability scale factor
scalefactor = mean(pheno$Status)*(1-mean(pheno$Status))/dnorm(qnorm(1-mean(pheno$Status)))^2

#Estimate and SE of heritability on the underlying liability scale
vpredict(modelfit,h_u~V1/(V1+V2)*scalefactor)

#Estimated breeding values (EBVs) of families
EBV = modelfit$U$`u:Family`$Status

#Mean EBV of cohorts
mean(EBV[grep("27.",names(EBV))])*100
mean(EBV[grep("28.",names(EBV))])*100
mean(EBV[grep("29.",names(EBV))])*100
mean(EBV[grep("30.",names(EBV))])*100

MAS_families = c(30.001,30.002,30.005,30.006,30.007,30.008,30.010,30.011,30.012,30.013,
                 30.014,30.015,30.016,30.017,30.018,30.019,30.022,30.028,30.029,30.033,
                 30.053,30.054,30.055,30.056,30.058,30.060,30.061,30.062,30.063,30.064,
                 30.065,30.066,30.068,30.071,30.073,30.075,30.076,30.078,30.079)

#EBVs of cohort 30 families
EBV_C30 = EBV[grep("30.",names(EBV))]

#Mean EBV of cohort 30 families selected with marker-assisted selection
mean(EBV_C30[names(EBV_C30) %in% MAS_families])*100

#Mean EBV of cohort 30 families selected without marker-assisted selection
mean(EBV_C30[!names(EBV_C30) %in% MAS_families])*100