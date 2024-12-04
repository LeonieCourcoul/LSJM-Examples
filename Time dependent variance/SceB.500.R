### Script simulations - Scenario B - 500

library(LSJM)
library(dplyr)
library(mvtnorm)
gaussKronrod <-
  function (k = 15) {
    sk <- c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
            0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
            -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
            0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329)
    wk15 <- c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
              0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
              0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
              0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)
    wk7 <- c(0.129484966168869693270611432679082, 0.279705391489276667901467771423780, 0.381830050505118944950369775488975,
             0.417959183673469387755102040816327, 0.381830050505118944950369775488975, 0.279705391489276667901467771423780, 0.129484966168869693270611432679082)
    if (k == 7)
      list(sk = sk[1:7], wk = wk7)
    else
      list(sk = sk, wk = wk15)
  }

# Choices of parameters
n <- 500
## longitudinal parameters
beta0 <- 142
beta1 <- 3
cholesky.globale <- matrix(c(14.4,0,0,0,-1.2,2.8,0,0,0,0,0.01,0,0,0,-0.06,0.11), ncol = 4, byrow = T)
B <- cholesky.globale%*%t(cholesky.globale)
M0 <- 2.4
M1 <- 0.05

## Survival parameters
alpha0 <- -7
alpha_sigma <- 0.07
alpha_y <- 0.02
alpha.slope <-0.01
## Survival parameters Competing Risk
alpha0.CR <- -4
alpha_sigma.CR <- 0.15
alpha_y.CR <- -0.01
alpha.slope.CR <- -0.14

shape_rac <- 1.1
shape.CR_rac <- 1.3

Gene_data <- function(n, beta0, beta1, B, M0,M1,
                      alpha0, shape_rac, alpha_y, alpha_slope, alpha_sigma,
                      alpha0.CR, shape.CR_rac, alpha_y.CR, alpha_slope.CR, alpha_sigma.CR){
  shape <- shape_rac**2
  shape.CR <- shape.CR_rac**2
  data_long <- c()
  sk <- gaussKronrod()$sk
  wk <- gaussKronrod()$wk
  t.max.event <- 5+1
  for(i in 1:n){
    ## longitudinal part
    random.effects <- rmvnorm(n=1, sigma = B)
    b0 <- random.effects[,1]; b1 <- random.effects[,2]
    omega0 <- random.effects[,3]; omega1 <- random.effects[,4]
    visit_i <- c(0)
    visit_i <- c(visit_i, runif(1,2/12,4/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,2/12,4/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,2/12,4/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,2/12,4/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,5/12,7/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,5/12,7/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,5/12,7/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,5/12,7/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,5/12,7/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,5/12,7/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,5/12,7/12))
    visit_i <- c(visit_i, visit_i[length(visit_i)]+runif(1,5/12,7/12))
    data_long_i <- as.data.frame(cbind(rep(i, length(visit_i)), # ID patient
                                       sort(visit_i))) # value of visit
    colnames(data_long_i) <- c("ID", "visit")
    sigma <- exp(M0 + omega0 + visit_i*(M1+omega1))
    error <- rnorm(length(visit_i), sd = sigma)
    data_long_i$y <- beta0 + b0 + visit_i*(beta1+b1) + error

    ## survival part
    u1 <- runif(1); u2 <- runif(1)

    survInvS1 <- function(tstar){
      (tstar/2)*sum(shape*wk*(((tstar/2)*(sk+1))**(shape-1))*exp(alpha0 + alpha.slope*(beta1+b1)+(beta0 +b0)*alpha_y
                                                                 + alpha_y*(beta1+b1)*((tstar/2)*(sk+1))+
                                                                   alpha_sigma*exp((M0+omega0)+(M1+omega1)*((tstar/2)*(sk+1))))) + log(u1)
    }
    Root.1 <- try(expr = uniroot(survInvS1,
                                 interval = c(1e-05, t.max.event))$root,
                  silent = TRUE)
    survInvS2 <- function(tstar){
      (tstar/2)*sum(shape.CR*wk*(((tstar/2)*(sk+1))**(shape.CR-1))*exp(alpha0.CR + alpha.slope.CR*(beta1+b1)+(beta0 +b0)*alpha_y.CR
                                                                       + alpha_y.CR*(beta1+b1)*((tstar/2)*(sk+1))+
                                                                         alpha_sigma.CR*exp((M0+omega0)+(M1+omega1)*((tstar/2)*(sk+1))))) + log(u2)
    }
    Root.2 <- try(expr = uniroot(survInvS2,
                                 interval = c(1e-05, t.max.event))$root,
                  silent = TRUE)

    if(inherits(Root.1, "try-error") & inherits(Root.2, "try-error")){
      data_long_i$time <- t.max.event
      data_long_i$event <- 0
    }
    else{
      if(inherits(Root.1, "try-error")){
        data_long_i$time <- ifelse(Root.2 < t.max.event, Root.2, t.max.event)
        data_long_i$event <- ifelse(Root.2 < t.max.event,2,0)
      }
      else{
        if(inherits(Root.2, "try-error")){
          data_long_i$time <- ifelse(Root.1 < t.max.event, Root.1, t.max.event)
          data_long_i$event <- ifelse(Root.1 < t.max.event,1,0)
        }
        else{
          if(Root.1 <= Root.2){
            data_long_i$time <- ifelse(Root.1 < t.max.event, Root.1, t.max.event)
            data_long_i$event<- ifelse(Root.1 < t.max.event,1,0)
          }
          if(Root.2 < Root.1){
            data_long_i$time <- ifelse(Root.2 < t.max.event, Root.2, t.max.event)
            data_long_i$event <- ifelse(Root.2 < t.max.event,2,0)
          }
        }
      }

    }
    data_long <- rbind(data_long, data_long_i)

  }

  data_long <- data_long[which(data_long$visit <= data_long$time),]
  data_long <- data_long[which(data_long$visit >=0),]
  data_long

}

nb.simu <- 4
nparam <- 26
estimateur.step1 <- matrix(NA, nrow = nb.simu, ncol = nparam)
se.estimateur.step1 <- matrix(NA, nrow = nb.simu, ncol = nparam)
estimateur.step2 <- matrix(NA, nrow = nb.simu, ncol = nparam)
se.estimateur.step2 <- matrix(NA, nrow = nb.simu, ncol = nparam)
non_cvg1 <- 0
list_noncv1 <- c()
reason1 <- c()
non_cvg2 <- 0
list_noncv2 <- c()
reason2 <- c()

for(rep in 1:nb.simu){
  print(rep)
  echantillon <- Gene_data(n, beta0, beta1, B, M0,M1,
                           alpha0, shape_rac, alpha_y, alpha_slope, alpha_sigma,
                           alpha0.CR, shape.CR_rac, alpha_y.CR, alpha_slope.CR, alpha_sigma.CR)
  max.obs.visit <- 6
  echantillon$event <- ifelse(echantillon$time > max.obs.visit, 0, echantillon$event)
  echantillon$time <- ifelse(echantillon$time > max.obs.visit, max.obs.visit, echantillon$time)
  echantillon$event1 <- NA
  echantillon$event2 <- NA
  echantillon$event1 <- ifelse(echantillon$event == 1, 1, 0)
  echantillon$event2 <- ifelse(echantillon$event == 2, 1, 0)

  lsmm.simu <- lsmm(formFixed = y ~ visit , formRandom = ~ visit ,
                    formGroup = ~ ID,  timeVar = "visit", formVar = "cov-dependent",
                    formFixedVar = ~ visit  , formRandomVar = ~ visit, correlated_re = F,
                    data.long = echantillon, nproc = 2, print.info = F,  S1 = 500, S2 = 5000, maxiter = 500)

  lsjm.simu <- lsjm(Objectlsmm = lsmm.simu, survival_type = 'CR',
                    formSurv_01 = ~ 1, formSurv_02 = ~ 1,
                    sharedtype_01 = c("value", "slope", "variability") ,
                    sharedtype_02 = c("value","slope", "variability") ,
                    hazardBase_01 = "Weibull",   hazardBase_02 = "Weibull",
                    formSlopeFixed = ~ 1, formSlopeRandom = ~ 1, index_beta_slope = c(2),
                    index_b_slope = c(2),  delta1 = ~ event1, delta2 = ~ event2,
                    Time_T = ~ time,  nproc = 10, print.info = F, file ="", S1 = 500, S2 = 5000, maxiter = 500,
                    epsa = 0.0001, epsb = 0.0001, epsd = 0.0001)


  #save(lsjm.simu, file = paste("lsjm.simu_SceB.500","_",rep,".RData",sep=""))


  if(lsjm.simu$info_conv_step1$conv !=1){
    non_cvg1 <- non_cvg1+1
    list_noncv1 <- c(list_noncv1,rep)
    reason1 <- c(reason1, lsjm.simu$info_conv_step1$conv)
  }
  else{
    #Results of Step 1
    ## Computation of covariance matrix of random effects
    nb.chol <- 6
    nb.e.a <- 2
    nb.e.a.sigma <- 2
    var_trans <- matrix(rep(0,length(lsjm.simu$result_step1$b)**2),nrow=length(lsjm.simu$result_step1$b),ncol=length(lsjm.simu$result_step1$b))
    var_trans[upper.tri(var_trans, diag=T)] <- lsjm.simu$result_step1$v
    sd.param <- sqrt(diag(var_trans))
    curseur <- length(lsjm.simu$result_step1$b) - nb.chol + 1
    borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
    C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
    C1[lower.tri(C1, diag=T)] <- lsjm.simu$result_step1$b[curseur:borne1]
    C1 <- as.matrix(C1)
    Index.C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
    Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a,2)+nb.e.a)
    Index.C1 <- as.matrix(Index.C1)
    borne3 <- borne1 + choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma
    C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
    C3[lower.tri(C3, diag=T)] <- lsjm.simu$result_step1$b[(borne1+1):borne3]
    C3 <- as.matrix(C3)

    Index.C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
    Index.C3[lower.tri(Index.C3, diag=T)] <- 1:(choose(nb.e.a.sigma,2)+nb.e.a.sigma)
    Index.C3 <- as.matrix(Index.C3)

    MatCovb <- C1%*%t(C1)
    MatCovSig <- C3%*%t(C3)
    param_est <- c(lsjm.simu$result_step1$b,unique(c(t(MatCovb))),unique(c(t(MatCovSig))))

    var_trans <- matrix(rep(0,length(lsjm.simu$result_step1$b)**2),nrow=length(lsjm.simu$result_step1$b),ncol=length(lsjm.simu$result_step1$b))
    var_trans[upper.tri(var_trans, diag=T)] <- lsjm.simu$result_step1$v
    trig.cov <- var_trans[curseur:borne1,curseur:borne1]
    trig.cov <- trig.cov+t(trig.cov)
    diag(trig.cov) <- diag(trig.cov)/2
    cov.cholesky <- trig.cov
    Cov.delta <- matrix(NA, ncol = ncol(cov.cholesky), nrow = ncol(cov.cholesky))

    element.chol <- lsjm.simu$result_step1$b[curseur:length(lsjm.simu$result_step1$b)]
    for(i in 1:ncol(C1)){
      for(j in i:ncol(C1)){
        resultat <- 0
        k <- i
        m <- j
        for(t in 1:min(i,j,ncol(Cov.delta))){
          for(s in 1:min(k,m,ncol(Cov.delta))){
            resultat <- resultat +
              C1[j,t]*C1[m,s]*cov.cholesky[Index.C1[i,t],Index.C1[k,s]] +
              C1[j,t]*C1[k,s]*cov.cholesky[Index.C1[i,t],Index.C1[m,s]] +
              C1[i,t]*C1[m,s]*cov.cholesky[Index.C1[j,t],Index.C1[k,s]] +
              C1[i,t]*C1[k,s]*cov.cholesky[Index.C1[j,t],Index.C1[m,s]]
          }
        }
        sd.param <- c(sd.param,sqrt(resultat))
      }
    }
    trig.cov <- var_trans[(borne1+1):borne3,(borne1+1):borne3]
    trig.cov <- trig.cov+t(trig.cov)
    diag(trig.cov) <- diag(trig.cov)/2
    cov.cholesky <- trig.cov
    Cov.delta <- matrix(NA, ncol = ncol(cov.cholesky), nrow = ncol(cov.cholesky))

    element.chol <- lsjm.simu$result_step1$b[(borne1+1):borne3]
    for(i in 1:ncol(C3)){
      for(j in i:ncol(C3)){
        resultat <- 0
        k <- i
        m <- j
        for(t in 1:min(i,j,ncol(Cov.delta))){
          for(s in 1:min(k,m,ncol(Cov.delta))){
            resultat <- resultat +
              C3[j,t]*C3[m,s]*cov.cholesky[Index.C3[i,t],Index.C3[k,s]] +
              C3[j,t]*C3[k,s]*cov.cholesky[Index.C3[i,t],Index.C3[m,s]] +
              C3[i,t]*C3[m,s]*cov.cholesky[Index.C3[j,t],Index.C3[k,s]] +
              C3[i,t]*C3[k,s]*cov.cholesky[Index.C3[j,t],Index.C3[m,s]]
          }
        }
        sd.param <- c(sd.param,sqrt(resultat))
      }
    }


    table.res1 <- cbind(param_est, sd.param)
    table.res1 <- as.data.frame(table.res1)

    estimateur.step1[rep,] <- table.res1[,1]
    se.estimateur.step1[rep,] <- table.res1[,2]

    if(lsjm.simu$info_conv_step2$conv !=1){
      non_cvg2 <- non_cvg2+1
      list_noncv2 <- c(list_noncv2,rep)
      reason2 <- c(reason2, lsjm.simu$info_conv_step2$conv)
    }
    else{
      # Results of step 2
      estimateur.step2[rep,] <- lsjm.simu$table.res[,1]
      se.estimateur.step2[rep,] <- lsjm.simu$table.res[,2]
    }

  }

}


## Information about convergence :
non_cvg1
list_noncv1
reason1

non_cvg2
list_noncv2
reason2

estimateur.step2 <- na.omit(estimateur.step2)
se.estimateur.step2 <- na.omit(se.estimateur.step2)

estimateur.step1 <- na.omit(estimateur.step1)
se.estimateur.step1 <- na.omit(se.estimateur.step1)

## Results
true_param <- c(1.1,-7,0.02,0.01,0.07,1.3,-4,-0.01,-0.14,0.15,142,3,2.4,0.05,14.4,-1.2,2.8,0.01,-0.06,0.11,
                27.36,-17.28,9.224,0.0001,-0.0006,0.0157)

### Step 1
Simu_analysis <- cbind(true_param,
                       apply(estimateur.step1, 2, mean),
                       apply(estimateur.step1, 2, sd),
                       apply(se.estimateur.step1,2,mean))
Simu_analysis <- as.data.frame(Simu_analysis)
bias <- estimateur.step1- matrix(rep(true_param,nrow(estimateur.step1)), nrow = nrow(estimateur.step1),byrow = T)
bias.step1 <- apply(bias,2,mean)
Simu_analysis <- cbind(Simu_analysis, bias.step1)
CI.inf.estimateur.step1 <- estimateur.step1 -1.96*se.estimateur.step1
CI.sup.estimateur.step1 <- estimateur.step1 +1.96*se.estimateur.step1
nb.simu.cv <- nrow(estimateur.step1)
CR1 <- matrix(ncol = nparam, nrow = nb.simu.cv)
for(i in 1:nb.simu.cv){
  CR1[i,] <- (true_param >= CI.inf.estimateur.step1[i,]) & (true_param <= CI.sup.estimateur.step1[i,])
}

CR1 <- as.data.frame(CR1)
colnames(CR1) <- c("shape","alpha0","alphaCV","alphaSlope","alphaSigma",
                   "shape.CR","alpha0.CR","alphaCV.CR","alphaSlope.CR","alphaSigma.CR",
                   "beta0","beta1","M0","M1",
                   "chol1","chol2","chol3","chol4","chol5","chol6","var_b0","cov_b0b1","var_b1","var_tau0","cov_tau0tau1","var_tau1")
CR_pour1 <- c((table(CR1$shape)/nb.simu.cv)[2],
              (table(CR1$alpha0)/nb.simu.cv)[2],
              (table(CR1$alphaCV)/nb.simu.cv)[2],
              (table(CR1$alphaSlope)/nb.simu.cv)[2],
              (table(CR1$alphaSigma)/nb.simu.cv)[2],
              (table(CR1$shape.CR)/nb.simu.cv)[2],
              (table(CR1$alpha0.CR)/nb.simu.cv)[2],
              (table(CR1$alphaCV.CR)/nb.simu.cv)[2],
              (table(CR1$alphaSlope.CR)/nb.simu.cv)[2],
              (table(CR1$alphaSigma.CR)/nb.simu.cv)[2],
              (table(CR1$beta0)/nb.simu.cv)[2],
              (table(CR1$beta1)/nb.simu.cv)[2],
              (table(CR1$M0)/nb.simu.cv)[2],
              (table(CR1$M1)/nb.simu.cv)[2],
              (table(CR1$chol1)/nb.simu.cv)[2],
              (table(CR1$chol2)/nb.simu.cv)[2],
              (table(CR1$chol3)/nb.simu.cv)[2],
              (table(CR1$chol4)/nb.simu.cv)[2],
              (table(CR1$chol5)/nb.simu.cv)[2],
              (table(CR1$chol6)/nb.simu.cv)[2],
              (table(CR1$var_b0)/nb.simu.cv)[2],
              (table(CR1$cov_b0b1)/nb.simu.cv)[2],
              (table(CR1$var_b1)/nb.simu.cv)[2],
              (table(CR1$var_tau0)/nb.simu.cv)[2],
              (table(CR1$cov_tau0tau1)/nb.simu.cv)[2],
              (table(CR1$var_tau1)/nb.simu.cv)[2]
)

CR_pour1 <- CR_pour1*100
Simu_analysis <- cbind(Simu_analysis,CR_pour1)

### Step 2
Simu_analysis <- cbind(Simu_analysis,
                       apply(estimateur.step2, 2, mean),
                       apply(estimateur.step2, 2, sd),
                       apply(se.estimateur.step2,2,mean))
bias <- estimateur.step2- matrix(rep(true_param,nrow(estimateur.step2)), nrow = nrow(estimateur.step2),byrow = T)
bias.step2 <- apply(bias,2,mean)
Simu_analysis <- cbind(Simu_analysis, bias.step2)
CI.inf.estimateur.step2 <- estimateur.step2 -1.96*se.estimateur.step2
CI.sup.estimateur.step2 <- estimateur.step2 +1.96*se.estimateur.step2
nb.simu.cv <- nrow(estimateur.step2)
CR2 <- matrix(ncol = nparam, nrow = nb.simu.cv)
for(i in 1:nb.simu.cv){
  CR2[i,] <- (true_param >= CI.inf.estimateur.step2[i,]) & (true_param <= CI.sup.estimateur.step2[i,])
}

CR2 <- as.data.frame(CR2)
colnames(CR2) <- colnames(CR1)
CR_pour2 <- c((table(CR2$shape)/nb.simu.cv)[2],
              (table(CR2$alpha0)/nb.simu.cv)[2],
              (table(CR2$alphaCV)/nb.simu.cv)[2],
              (table(CR2$alphaSlope)/nb.simu.cv)[2],
              (table(CR2$alphaSigma)/nb.simu.cv)[2],
              (table(CR2$shape.CR)/nb.simu.cv)[2],
              (table(CR2$alpha0.CR)/nb.simu.cv)[2],
              (table(CR2$alphaCV.CR)/nb.simu.cv)[2],
              (table(CR2$alphaSlope.CR)/nb.simu.cv)[2],
              (table(CR2$alphaSigma.CR)/nb.simu.cv)[2],
              (table(CR2$beta0)/nb.simu.cv)[2],
              (table(CR2$beta1)/nb.simu.cv)[2],
              (table(CR2$M0)/nb.simu.cv)[2],
              (table(CR2$M1)/nb.simu.cv)[2],
              (table(CR2$chol1)/nb.simu.cv)[2],
              (table(CR2$chol2)/nb.simu.cv)[2],
              (table(CR2$chol3)/nb.simu.cv)[2],
              (table(CR2$chol4)/nb.simu.cv)[2],
              (table(CR2$chol5)/nb.simu.cv)[2],
              (table(CR2$chol6)/nb.simu.cv)[2],
              (table(CR2$var_b0)/nb.simu.cv)[2],
              (table(CR2$cov_b0b1)/nb.simu.cv)[2],
              (table(CR2$var_b1)/nb.simu.cv)[2],
              (table(CR2$var_tau0)/nb.simu.cv)[2],
              (table(CR2$cov_tau0tau1)/nb.simu.cv)[2],
              (table(CR2$var_tau1)/nb.simu.cv)[2]
)

CR_pour2 <- CR_pour2*100
Simu_analysis <- cbind(Simu_analysis,CR_pour2)
colnames(Simu_analysis) <- c("True Value","Mean1", "ESE1", "ASE1", "Bias1", "CR1",
                             "Mean2", "ESE2", "ASE2", "Bias2", "CR2")
rownames(Simu_analysis) <- colnames(CR2)
print(Simu_analysis)
View(Simu_analysis[c(1:14,21:26),])




















