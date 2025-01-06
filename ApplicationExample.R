library(LSJM)

data(threeC)
threeC$age.visit65 <- (threeC$age.visit-65)/10
threeC$age.final65 <- (threeC$age.final-65)/10
threeC$age0_65 <- (threeC$age0-65)/10
threeC$age.last65 <- (threeC$age.last-65)/10
threeC$age.first65 <- (threeC$age.first-65)/10
threeC$SBP <- threeC$SBP/10

m1 <- lsmm(formFixed = SBP ~ age.visit65,
           formRandom = ~ age.visit65,
           formGroup = ~ ID,
           timeVar = "age.visit65",
           data.long = threeC,
           formVar = "inter-intra",
           random_inter = T,
           random_intra = T,
           formGroupVisit = ~num.visit,
           correlated_re = FALSE,
           S1 = 500,
           S2 = 1000,
           nproc = 4)

summary(m1)

l1 <- lsjm(m1,
           survival_type = "IDM",
           formSurv_01=~sex,
           formSurv_02=~sex,
           formSurv_12=~sex,
           sharedtype_01 = c("value", "slope", "variability inter", "variability intra"),
           sharedtype_02 = c("value", "slope", "variability inter", "variability intra"),
           sharedtype_12 = c("value", "slope", "variability inter", "variability intra"),
           hazardBase_01 = "Weibull",
           hazardBase_02 = "Weibull",
           hazardBase_12 = "Weibull",
           delta1=~dem,
           delta2=~death,
           Time_T =~age.final65,
           Time_L =~age.last65,
           Time_R =~age.first65,
           Time_T0 =~age0_65,
           formSlopeFixed =~1,
           formSlopeRandom =~1,
           index_beta_slope = c(2),
           index_b_slope = c(2),
           S1 = 1000,
           S2 = 2000,
           nproc = 1,
           print.info = T,
           file = "")

summary(l1)



