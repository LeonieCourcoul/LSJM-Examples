library(LSJM)
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

# Algorithm of data generation:
Gen_data <- function(n, B, mu_sigma, mu_kappa, nb.mesures, beta0, beta1,
                     shape_01,alpha0_01, alpha_y_01, alpha_slope_01, alpha_sigma_01, alpha_kappa_01,
                     shape_02,alpha0_02, alpha_y_02, alpha_slope_02, alpha_sigma_02, alpha_kappa_02,
                     shape_12,alpha0_12, alpha_y_12, alpha_slope_12, alpha_sigma_12, alpha_kappa_12,
                     Censure){
  sk <- gaussKronrod()$sk
  wk <- gaussKronrod()$wk
  shape_01 <- shape_01**2
  shape_02 <- shape_02**2
  shape_12 <- shape_12**2
  data_long <- c()
  i <- 1
  nb.tour <- 0
  while(i<=n){
    add.i <- T
    #On tire les effets aléatoires
    random.effects <- rmvnorm(n=1, sigma = B)
    b0 <- random.effects[1]
    b1 <- random.effects[2]
    tau.sigma <- random.effects[3]
    sigma <- exp(mu_sigma+tau.sigma)
    tau.kappa <- random.effects[4]
    kappa <- exp(mu_kappa+tau.kappa)
    #On génère les temps de visites
    A_0i <- rbeta(1,1.5,3)*65+65
    T1i <- A_0i + runif(1, 0.8, 1.2)
    T2i <- T1i + runif(1, 0.8, 1.2)
    T4i <- T2i + runif(1, 1.8, 2.2)
    T7i <- T4i + runif(1, 2.8, 3.2)
    T10i <- T7i + runif(1, 2.8, 3.2)
    T12i <- T10i + runif(1, 1.8, 2.2)
    T14i <- T12i + runif(1, 1.8, 2.2)
    T17i <- T14i + runif(1, 2.8, 3.2)
    visit <- c(A_0i, T1i,T2i, T4i, T7i, T10i, T12i, T14i, T17i)
    visit <- (visit-65)/10
    visit <- rep(visit, each = nb.mesures)
    #browser()
    data_long_i <- c()
    data_long_i <- as.data.frame(cbind(rep(i, length(visit)),
                                       rep(c("S0","S1","S2","S4","S7","S10","S12","S14","S17"), each = nb.mesures),
                                       sort(visit)))
    colnames(data_long_i) <- c("ID", "num.visit",  "visit")
    data_long_i$num.mesures <- rep(1:nb.mesures,length(visit)/nb.mesures)
    error1 <- rnorm(length(visit)/nb.mesures, mean = 0, sd = sigma)
    error1 <- rep(error1, each = nb.mesures)
    error2 <- rnorm(length(visit), mean = 0, sd = kappa)
    data_long_i$y <- beta0+b0+visit*(beta1+b1) + error1+error2
    data_long_i$Time_T0 <- (A_0i-65)/10
    u_01 <- runif(1); u_02 <- runif(1)

    S_01_inv <- function(tstar){
      (tstar/2)*sum(shape_01*wk*(((tstar/2)*(sk+1))**(shape_01-1))*exp(alpha0_01 +
                                                                         alpha_y_01*(beta0 + b0 + (beta1+b1)*((tstar/2)*(sk+1))) +
                                                                         alpha_slope_01*((beta1+b1)) +
                                                                         alpha_sigma_01*sigma +
                                                                         alpha_kappa_01*kappa
      )) + log(u_01)
    }
    T_01 <- try(expr = uniroot(S_01_inv,
                               interval = c(0, max(visit)))$root,
                silent = TRUE)

    S_02_inv <- function(tstar){
      (tstar/2)*sum(shape_02*wk*(((tstar/2)*(sk+1))**(shape_02-1))*exp(alpha0_02 +
                                                                         alpha_y_02*(beta0 + b0 + (beta1+b1)*((tstar/2)*(sk+1))) +
                                                                         alpha_slope_02*((beta1+b1)) +
                                                                         alpha_sigma_02*sigma +
                                                                         alpha_kappa_02*kappa
      )) + log(u_02)
    }
    T_02 <- try(expr = uniroot(S_02_inv,
                               interval = c(0, max(visit)))$root,
                silent = TRUE)
    if(inherits(T_01, "try-error")){
      T_01 <- 1000000
    }
    if(inherits(T_02, "try-error")){
      T_02 <- 1000000
    }
    if(T_01 < (A_0i -65)/10 || T_02 < (A_0i-65)/10 || A_0i > 85){ #L'individu n'est pas inclu
      add.i <- F
    }
    else{
      if(T_02 < T_01){
        if(T_02 > max(visit) ){ # L'individu est censuré sans démence
          data_long_i$age.evt <- max(visit)
          data_long_i$evt <- 0
          data_long_i$delta1 <- 0
          data_long_i$delta2 <- 0
          data_long_i$delta1_IC <- 0
          data_long_i$delta2_IC <- 0
          data_long_i$Time_T <- max(visit)
          data_long_i$Time_L <- max(visit)
          data_long_i$Time_R <- max(visit)
          data_long_i$Time_T_IC <- max(visit)
          data_long_i$Time_L_IC <- max(visit)
          data_long_i$Time_R_IC <- max(visit)
          i <- i+1
        }
        else{
          data_long_i$age.evt <- T_02
          data_long_i$evt <- 2
          data_long_i$delta1 <- 0
          data_long_i$delta2 <- 1
          data_long_i$Time_T <- T_02
          data_long_i$Time_L <- T_02
          data_long_i$Time_R <- T_02
          data_long_i$delta1_IC <- 0
          data_long_i$delta2_IC <- 1
          data_long_i$Time_T_IC <- T_02
          data_long_i$Time_L_IC <- max(visit[which(visit <= T_02)])
          data_long_i$Time_R_IC <- max(visit[which(visit <= T_02)])
          data_long_i <- data_long_i[which(data_long_i$visit <= T_02),]
          i <- i + 1
        }
      }
      else{
        if(T_01 > max(visit)){# L'individu est censuré sans démence
          data_long_i$age.evt <- max(visit)
          data_long_i$evt <- 0
          data_long_i$delta1 <- 0
          data_long_i$delta2 <- 0
          data_long_i$Time_T <- max(visit)
          data_long_i$Time_L <- max(visit)
          data_long_i$Time_R <- max(visit)
          data_long_i$delta1_IC <- 0
          data_long_i$delta2_IC <- 0
          data_long_i$Time_T_IC <- max(visit)
          data_long_i$Time_L_IC <- max(visit)
          data_long_i$Time_R_IC <- max(visit)
          i <- i+1
        }
        else{ #L'individu est dément, on oublie T_02, on calcule T_12
          u12 <- runif(1)
          S_12 <- function(tps){
            (tps/2)*sum(shape_12*wk*(((tps/2)*(sk+1))**(shape_12-1))*exp(alpha0_12 +
                                                                           alpha_y_12*(beta0 + b0 + (beta1+b1)*((tps/2)*(sk+1))) +
                                                                           alpha_sigma_12*sigma +
                                                                           alpha_slope_12*((beta1+b1)) +
                                                                           alpha_kappa_12*kappa
            ))
          }
          u12_corrige <- u12*exp(-S_12(T_01))
          S_12_inv <- function(tstar){
            (tstar/2)*sum(shape_12*wk*(((tstar/2)*(sk+1))**(shape_12-1))*exp(alpha0_12 +
                                                                               alpha_y_12*(beta0 + b0 + (beta1+b1)*((tstar/2)*(sk+1))) +
                                                                               alpha_slope_12*((beta1+b1)) +
                                                                               alpha_sigma_12*sigma +
                                                                               alpha_kappa_12*kappa
            )) + log(u12_corrige)
          }
          T_12 <- try(expr = uniroot(S_12_inv,
                                     interval = c(0, max(visit)))$root,
                      silent = TRUE)
          if(inherits(T_12, "try-error")){
            T_12 <- 100000000
          }

          if(T_12 > max(visit)){
            if(T_01 > max(visit)){
              print("non")
              data_long_i$age.evt <- max(visit)
              data_long_i$evt <- 0
              data_long_i$delta1 <- 0
              data_long_i$delta2 <- 0
              data_long_i$Time_T <- max(visit)
              data_long_i$Time_L <- max(visit)
              data_long_i$Time_R <- max(visit)
              data_long_i$delta1_IC <- 0
              data_long_i$delta2_IC <- 0
              data_long_i$Time_T_IC <- max(visit)
              data_long_i$Time_L_IC <- max(visit)
              data_long_i$Time_R_IC <- max(visit)
              i <- i+1
            }
            else{
              data_long_i$age.evt <- T_01
              data_long_i$evt <- 1
              data_long_i$delta1 <- 1
              data_long_i$delta2 <- 0
              data_long_i$Time_T <- max(visit)
              data_long_i$Time_L <- T_01
              data_long_i$Time_R <- T_01
              data_long_i$delta1_IC <- 1
              data_long_i$delta2_IC <- 0
              data_long_i$Time_T_IC <- max(visit)
              data_long_i$Time_L_IC <- max(visit[which(visit<=T_01)])
              data_long_i$Time_R_IC <- min(visit[which(visit>=T_01)])
              i <- i+1
            }
          }
          else{
            if(T_01 > max(visit)){
              print("non2")
              data_long_i$age.evt <- T_12
              data_long_i$evt <- 2
              data_long_i$delta1 <- 1
              data_long_i$delta2 <- 1
              data_long_i$Time_T <- T_12
              data_long_i$Time_L <- T_01
              data_long_i$Time_R <- T_01
              data_long_i$delta1_IC <- 0
              data_long_i$delta2_IC <- 1
              data_long_i$Time_T_IC <- T_12
              data_long_i$Time_L_IC <- max(visit)
              data_long_i$Time_R_IC <- max(visit)
              i <- i+1
            }
            else{
              tps_R <- min(visit[which(visit>=T_01)])
              if(tps_R <= T_12){
                data_long_i$age.evt <- T_01
                data_long_i$evt <- 1
                data_long_i$delta1_IC <- 1
                data_long_i$delta2_IC <- 1
                data_long_i$Time_T_IC <- T_12
                data_long_i$Time_L_IC <- max(visit[which(visit<=T_01)])
                data_long_i$Time_R_IC <- min(visit[which(visit>=T_01)])
                data_long_i$delta1 <- 1
                data_long_i$delta2 <- 1
                data_long_i$Time_T <- T_12
                data_long_i$Time_L <- T_01
                data_long_i$Time_R <- T_01
                data_long_i <- data_long_i[which(data_long_i$visit <= T_12),]
                i <- i+1
              }
              else{
                data_long_i$age.evt <- T_12
                data_long_i$evt <- 2
                data_long_i$delta1_IC <- 0
                data_long_i$delta2_IC <- 1
                data_long_i$Time_T_IC <- T_12
                data_long_i$Time_L_IC <- max(visit[which(visit<=T_12)])
                data_long_i$Time_R_IC <- max(visit[which(visit<=T_12)])
                data_long_i$delta1 <- 1
                data_long_i$delta2 <- 1
                data_long_i$Time_T <- T_12
                data_long_i$Time_L <- T_01
                data_long_i$Time_R <- T_01
                data_long_i <- data_long_i[which(data_long_i$visit <= T_12),]
                i <- i+1
              }
            }
          }
        }
        #browser()
      }
    }
    if(add.i){
      data_long <- rbind(data_long,data_long_i)
    }
  }

  data_long
}

# Parameters to generate the data
beta0 <- 14
beta1 <- 0.17
mu_sigma <- 0.30
mu_kappa <- -0.23
alpha0_01 <- -4
alpha_y_01 <- -0.06
alpha_slope_01 <- 1.6e-05
alpha_sigma_01 <- 0.5
alpha_kappa_01 <- 0.01
alpha0_02 <- -2.5
alpha_y_02 <- -0.1
alpha_slope_02 <- -0.4
alpha_sigma_02 <- 0.46
alpha_kappa_02 <- 0.21
alpha0_12 <- -2.2
alpha_y_12 <- 0.04
alpha_slope_12 <- 0.02
alpha_sigma_12 <- -0.12
alpha_kappa_12 <- -0.18
shape_01 <- 2
shape_02 <- 1.7
shape_12 <- 1.7
chol <- matrix(rep(0,16), ncol=4,nrow=4)
chol[1,1] <- 2.2
chol[2,1:2] <- c(-0.86, -0.68)
chol[3,3] <- 0.27
chol[4,3:4] <- c(0.06,0.25)
B <- chol%*%t(chol)

#Initialisation of the result tables
nparam <- 28
nb.simu <- 500
estimateur.step2.IC <- matrix(NA, nrow = nb.simu, ncol = nparam)
se.estimateur.step2.IC <- matrix(NA, nrow = nb.simu, ncol = nparam)
non_cvg1.IC <- 0
non_cvg2.IC <- 0
estimateur.step2.naive <- matrix(NA, nrow = nb.simu, ncol = nparam)
se.estimateur.step2.naive <- matrix(NA, nrow = nb.simu, ncol = nparam)
non_cvg1.naive <- 0
non_cvg2.naive <- 0

# Loop for nb.simu replications:
for(rep in 1:nb.simu){
  donnees <- Gen_data(1000,B = B,mu_sigma=mu_sigma,mu_kappa=mu_kappa,nb.mesures = 2,beta0=beta0,beta1=beta1,
                      shape_01 = shape_01, alpha0_01=alpha0_01, alpha_y_01=alpha_y_01, alpha_slope_01=alpha_slope_01, alpha_sigma_01=alpha_sigma_01, alpha_kappa_01=alpha_kappa_01,
                      shape_02 = shape_02, alpha0_02 =  alpha0_02, alpha_y_02 =  alpha_y_02,  alpha_slope_02=alpha_slope_02, alpha_sigma_02=alpha_sigma_02, alpha_kappa_02 =  alpha_kappa_02,
                      shape_12 = shape_12, alpha0_12 =  alpha0_12, alpha_y_12 = alpha_y_12,  alpha_slope_12=alpha_slope_12, alpha_sigma_12 =  alpha_sigma_12, alpha_kappa_12 =  alpha_kappa_12,
                      Censure = 20)

  # Some statistics about the generated data set
  data.id <- donnees[!duplicated(donnees$ID),]
  table(data.id$delta1_IC) #Nb 0-1
  nrow(data.id[which(data.id$delta1_IC == 1& data.id$delta2_IC == 1),]) #Nb 1-2
  nrow(data.id[which(data.id$delta1_IC == 0& data.id$delta2_IC == 1),]) #Nb 0-2
  donnees$visit <- as.numeric(donnees$visit)
  donnees$Time_Center <- NA
  donnees$Time_Center[which(donnees$delta1_IC == 1)] <- (donnees$Time_L_IC[which(donnees$delta1_IC == 1)] + donnees$Time_R_IC[which(donnees$delta1_IC == 1)])/2
  donnees$Time_Center[which(donnees$delta1_IC == 0)] <- donnees$Time_T_IC[which(donnees$delta1_IC == 0)]

  # Estimation for the Interval censoring model
  lsmm.simu <- lsmm(formFixed = y ~ visit,
                    formRandom = ~ visit,
                    formGroup = ~ ID,
                    timeVar = "visit",
                    data.long = donnees,
                    formVar = "inter-intra",
                    random_inter = T,
                    random_intra = T,
                    formGroupVisit = ~num.visit,
                    correlated_re = FALSE,
                    S1 = 1000,
                    S2 = 5000,
                    nproc = 4)

  lsjm.simu.IC <- lsjm(lsmm.simu,
                    survival_type = "IDM",
                    formSurv_01=~1,
                    formSurv_02=~1,
                    formSurv_12=~1,
                    sharedtype_01 = c("value", "slope","variability inter",  "variability intra"),
                    sharedtype_02 = c("value", "slope","variability inter",  "variability intra"),
                    sharedtype_12 = c("value", "slope","variability inter",  "variability intra"),
                    hazardBase_01 = "Weibull",
                    hazardBase_02 = "Weibull",
                    hazardBase_12 = "Weibull",
                    delta1=~delta1_IC,
                    delta2=~delta2_IC,
                    Time_T =~Time_T_IC,
                    Time_L =~Time_L_IC,
                    Time_R =~Time_R_IC,
                    Time_T0 =~Time_T0,
                    formSlopeFixed =~1,
                    formSlopeRandom =~1,
                    index_beta_slope = c(2),
                    index_b_slope = c(2),
                    S1 = 1000,
                    S2 = 1000,
                    nproc = 15)
  save(lsjm.simu.IC, file = paste("lsjm.simu_SceA.IC","_",rep,".RData",sep=""))

  #Estimation for the naive model
  lsjm.simu.naive <- lsjm(lsmm.simu,
                           survival_type = "IDM",
                           formSurv_01=~1,
                           formSurv_02=~1,
                           formSurv_12=~1,
                           sharedtype_01 = c("value", "slope","variability inter",  "variability intra"),
                           sharedtype_02 = c("value", "slope","variability inter",  "variability intra"),
                           sharedtype_12 = c("value", "slope","variability inter",  "variability intra"),
                           hazardBase_01 = "Weibull",
                           hazardBase_02 = "Weibull",
                           hazardBase_12 = "Weibull",
                           delta1=~delta1_IC,
                           delta2=~delta2_IC,
                           Time_T =~Time_T_IC,
                           Time_L =~Time_Center,
                           Time_R =~Time_Center,
                           Time_T0 =~Time_T0,
                           formSlopeFixed =~1,
                           formSlopeRandom =~1,
                           index_beta_slope = c(2),
                           index_b_slope = c(2),
                           S1 = 1000,
                           S2 = 1000,
                           nproc = 15)
                           
  save(lsjm.simu.naive, file = paste("lsjm.simu_SceA.naive","_",rep,".RData",sep=""))


  # Results
  if(lsjm.simu.IC$result_step1$istop !=1){
    non_cvg1.IC <- non_cvg1.IC+1
    if(lsjm.simu.IC$result_step2$istop !=1){
      non_cvg2.IC <- non_cvg2.IC+1
    }
  }
  else{
    if(lsjm.simu.IC$result_step2$istop !=1){
      non_cvg2.IC <- non_cvg2.IC+1
    }
    else{
      estimateur.step2.IC[rep,] <- lsjm.simu.IC$table.res[c(1:22,29:34),1]
      se.estimateur.step2.IC[rep,] <- lsjm.simu.IC$table.res[c(1:22,29:34),2]
    }
  }

  if(lsjm.simu.naive$result_step1$istop !=1){
    non_cvg1.naive <- non_cvg1.naive+1
    if(lsjm.simu.naive$result_step2$istop !=1){
      non_cvg2.naive <- non_cvg2.naive+1
    }
  }
  else{
    if(lsjm.simu.naive$result_step2$istop !=1){
      non_cvg2.naive <- non_cvg2.naive+1
    }
    else{
      estimateur.step2.naive[rep,] <- lsjm.simu.naive$table.res[c(1:22,29:34),1]
      se.estimateur.step2.naive[rep,] <- lsjm.simu.naive$table.res[c(1:22,29:34),2]
    }
  }


}

## Information about convergence:
### For the interval censoring model:
non_cvg1.IC
non_cvg2.IC
### For the naive model:
non_cvg1.naive
non_cvg2.naive


estimateur.step2.IC <- na.omit(estimateur.step2.IC)
se.estimateur.step2.IC <- na.omit(se.estimateur.step2.IC)
estimateur.step2.naive <- na.omit(estimateur.step2.naive)
se.estimateur.step2.naive <- na.omit(se.estimateur.step2.naive)

true_param <-c(2,-4, -0.06,1.6e-05,0.5,0.01,
               1.7,-2.5, -0.1, -0.4, 0.46,0.21,
               1.7,-2.2, 0.04, 0.02,-0.12,-0.18,
               14, 0.17, 0.30, -0.23, 4.84,-1.892,0.0729,0.0162,0.0661)

Simu_analysis.IC <- cbind(true_param,
                       apply(estimateur.step2.IC, 2, mean),
                       apply(estimateur.step2.IC, 2, sd),
                       apply(se.estimateur.step2.IC,2,mean))

Simu_analysis.IC <- as.data.frame(Simu_analysis.IC)
colnames(Simu_analysis.IC) <- c("True Value", "Mean", "ESE", "ASE")

CI.inf.estimateur.step2.IC <- estimateur.step2.IC -1.96*se.estimateur.step2.IC
CI.sup.estimateur.step2.IC <- estimateur.step2.IC +1.96*se.estimateur.step2.IC
nb.simu.cv <- nrow(estimateur.step2.IC)
CR2.IC <- matrix(ncol = nparam, nrow = nb.simu.cv)

for(i in 1:nb.simu.cv){
  CR2.IC[i,] <- (true_param >= CI.inf.estimateur.step2.IC[i,]) & (true_param <= CI.sup.estimateur.step2.IC[i,])
}

CR2.IC <- as.data.frame(CR2.IC)


colnames(CR2.IC) <- c("shape_01",
                   "alpha0_01", "alpha_y_01", "alpha_slope_01" , "alpha_01_sigma", "alpha_01_kappa",
                   "shape_02",
                   "alpha0_02", "alpha_y_02", "alpha_slope_02" , "alpha_02_sigma", "alpha_02_kappa",
                   "shape_12",
                   "alpha0_12", "alpha_y_12", "alpha_slope_12" , "alpha_12_sigma", "alpha_12_kappa",
                   "beta0", "beta1", "mu.sigma", "mu.kappa",
                   "var_b0", "cov_b0b1", "var_b1", "var_inter", "cov_interintra", "var_intra"
)
CR_pour2.IC <- c((table(CR2.IC$shape_01)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha0_01)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_y_01)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_slope_01)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_01_sigma)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_01_kappa)/nb.simu.cv)[2],
                  (table(CR2.IC$shape_02)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha0_02)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_y_02)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_slope_02)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_02_sigma)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_02_kappa)/nb.simu.cv)[2],
                  (table(CR2.IC$shape_12)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha0_12)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_y_12)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_slope_12)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_12_sigma)/nb.simu.cv)[2],
                  (table(CR2.IC$alpha_12_kappa)/nb.simu.cv)[2],
                  (table(CR2.IC$beta0)/nb.simu.cv)[2],
                  (table(CR2.IC$beta1)/nb.simu.cv)[2],
                  (table(CR2.IC$mu.sigma)/nb.simu.cv)[2],
                  (table(CR2.IC$mu.kappa)/nb.simu.cv)[2],
                  (table(CR2.IC$var_b0)/nb.simu.cv)[2],
                  (table(CR2.IC$cov_b0b1)/nb.simu.cv)[2],
                  (table(CR2.IC$var_b1)/nb.simu.cv)[2],
                  (table(CR2.IC$var_inter)/nb.simu.cv)[2],
                  (table(CR2.IC$cov_interintra)/nb.simu.cv)[2],
                  (table(CR2.IC$var_intra)/nb.simu.cv)[2]
)

CR_pour2.IC <- CR_pour2.IC*100
Simu_analysis.IC$CR <- CR_pour2.IC

# Results for the estimation dealing with interval censoring
print(Simu_analysis.IC)

# Naive
Simu_analysis.naive <- cbind(true_param,
                          apply(estimateur.step2.naive, 2, mean),
                          apply(estimateur.step2.naive, 2, sd),
                          apply(se.estimateur.step2.naive,2,mean))

Simu_analysis.naive <- as.data.frame(Simu_analysis.naive)
colnames(Simu_analysis.naive) <- c("True Value", "Mean", "ESE", "ASE")

CI.inf.estimateur.step2.naive <- estimateur.step2.naive -1.96*se.estimateur.step2.naive
CI.sup.estimateur.step2.naive <- estimateur.step2.naive +1.96*se.estimateur.step2.naive
nb.simu.cv <- nrow(estimateur.step2.naive)
CR2.naive <- matrix(ncol = nparam, nrow = nb.simu.cv)

for(i in 1:nb.simu.cv){
  CR2.naive[i,] <- (true_param >= CI.inf.estimateur.step2.naive[i,]) & (true_param <= CI.sup.estimateur.step2.naive[i,])
}

CR2.naive <- as.data.frame(CR2.naive)


colnames(CR2.naive) <- c("shape_01",
                      "alpha0_01", "alpha_y_01", "alpha_slope_01" , "alpha_01_sigma", "alpha_01_kappa",
                      "shape_02",
                      "alpha0_02", "alpha_y_02", "alpha_slope_02" , "alpha_02_sigma", "alpha_02_kappa",
                      "shape_12",
                      "alpha0_12", "alpha_y_12", "alpha_slope_12" , "alpha_12_sigma", "alpha_12_kappa",
                      "beta0", "beta1", "mu.sigma", "mu.kappa",
                      "var_b0", "cov_b0b1", "var_b1", "var_inter", "cov_interintra", "var_intra"
)
CR_pour2.naive <- c((table(CR2.naive$shape_01)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha0_01)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_y_01)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_slope_01)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_01_sigma)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_01_kappa)/nb.simu.cv)[2],
                 (table(CR2.naive$shape_02)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha0_02)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_y_02)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_slope_02)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_02_sigma)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_02_kappa)/nb.simu.cv)[2],
                 (table(CR2.naive$shape_12)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha0_12)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_y_12)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_slope_12)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_12_sigma)/nb.simu.cv)[2],
                 (table(CR2.naive$alpha_12_kappa)/nb.simu.cv)[2],
                 (table(CR2.naive$beta0)/nb.simu.cv)[2],
                 (table(CR2.naive$beta1)/nb.simu.cv)[2],
                 (table(CR2.naive$mu.sigma)/nb.simu.cv)[2],
                 (table(CR2.naive$mu.kappa)/nb.simu.cv)[2],
                 (table(CR2.naive$var_b0)/nb.simu.cv)[2],
                 (table(CR2.naive$cov_b0b1)/nb.simu.cv)[2],
                 (table(CR2.naive$var_b1)/nb.simu.cv)[2],
                 (table(CR2.naive$var_inter)/nb.simu.cv)[2],
                 (table(CR2.naive$cov_interintra)/nb.simu.cv)[2],
                 (table(CR2.naive$var_intra)/nb.simu.cv)[2]
)

CR_pour2.naive <- CR_pour2.naive*100
Simu_analysis.naive$CR <- CR_pour2.naive

# Results for the naive estimation
print(Simu_analysis.naive)





