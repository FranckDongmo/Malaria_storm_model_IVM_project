#setwd("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion")
setwd("C:/Users/franc/Desktop/PowerPoint PhD/R studio/New folder (2)")
rm(list = ls())#remet à jour tout l'espace de travail

reset <- function() {
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

library(xlsx)
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
# DRC model to estimate dose response link ---- 
library(drc)# load after GAM, package incompatibility 
library(patchwork)#pour montage figure

trapz<-function(x=NULL,y){
  if(is.null(x)){
    x<-1:length(y)
  }
  if(length(y)==1){z<-0}
  else{z<-caTools::trapz(x,y)}
  return(z)
}

P_eta.fun = function(eta, eta50, alpha) {
  p_eta = 1 / (1 + (eta / eta50)^alpha)
  return(p_eta)
}

ModelEeta50Alpha<-function(strategy){
  if(strategy==0){#Long lasting 0.6
    #fitting parameter for LAIF 0.6
    # load data ---- 
    # P4D_PK = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/IDR_tab4J_PK.rds")
    P4D_PK = readRDS("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\dataIVM\\IDR_tab4J_PK.rds")
    # mortality data, proportion of dead at 4 days after exposition 
    # gam_all = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/pred_gam_all.rds")
    gam_all = readRDS("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\dataIVM\\pred_gam_all.rds")
    CritereIVMformulation = "F31-0.6"
    sub_PK_F31 = P4D_PK %>%
      filter(strain == "vk5" & IVM_formulation == CritereIVMformulation)
    sub_PK_F31["DAI"] = as.numeric(sub_PK_F31$DAI)
    #fitting the formula
    SimulDrmLongLasting06= drm(prop_dead ~ DAI, data = sub_PK_F31,
                               fct = LL.4(fixed = c(NA, 0, 1, NA),
                                          names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                               type = "binomial")
    summary(SimulDrmLongLasting06)
    Eta50=SimulDrmLongLasting06$coefficients[2]
    Alpha=SimulDrmLongLasting06$coefficients[1]
    
    #pas très clair ici: juste pour faire de plot visibelement?
    gam_F31 = gam_all %>%
      filter(IVM_formulation == CritereIVMformulation) %>%
      mutate(DAI = Hour / 24, prop_dead = NA)
    Day_seq=length(unique(gam_F31$DAI))
    pred_F31 = predict(SimulDrmLongLasting06, newdata = gam_F31, interval = "confidence")
    predprop_F31 = data.frame(pred_F31) %>%
      cbind(gam_F31) %>%
      mutate(Formulation = "F31-0.6")
    predprop_F31 = predprop_F31[, !names(predprop_F31) %in% 
                                  c("IVM_concentration", "IVM_formulation", "Hour", "prop_dead")]
    predprop=predprop_F31
    
    sub_PK_F31 = sub_PK_F31 %>%mutate(DAI = as.numeric(DAI))
    sub_PK_F31 = sub_PK_F31 %>% mutate(Day = as.numeric(Day))
    sub_PK=sub_PK_F31
    DAI = seq(0,166, 1)
  }
  
  if(strategy==1){#Long lasting 1
    #fitting parameter for LAIF 1
    # Only for Lamidi study 1 mg/kg
    # load data ---- 
    # 1. Data from study of Lamidi (Bapla area?) with cattle injected at 1.0 mg/kg  
    #lamidi_data = read_xlsx("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/Survie_anopheles_sauvage110923.xlsx")
    lamidi_data = read_xlsx("C:\\Users\\franc\\Downloads\\mortality_rate\\mortality_rate\\data\\Survie_anopheles_sauvage110923.xlsx")
    surv_gamb = lamidi_data %>% dplyr::filter(Espece=="gambiae" & Traitement=="IVM") %>% 
      mutate(DAI2=as.numeric(str_extract(DAI, "\\d+")), Status4J=if_else(Temps<=4, 1,0))
    
    # vizu 
    # estimation of proportion of dead mosquitoes at 4 days by cattle and DAI
    # join with GAM_pk estimation of IVM plasma concentration 
    sg_sumdead = surv_gamb %>% group_by(Bovin, DAI2) %>% dplyr::filter(Status4J==1) %>% 
      summarise(ndead=n())
    sg_sum = surv_gamb %>% group_by(Bovin, DAI2) %>%  
      summarise(ntot=n())
    sg_tot = inner_join(sg_sumdead, sg_sum) 
    sg_tot = sg_tot %>% mutate(prop_dead=ndead/ntot, Day=DAI2) %>% droplevels()
    
    #fitting the formula
    SimulDrmLongLasting1=drm(prop_dead~DAI2, data=sg_tot, fct=LL.4(fixed=c(NA, 0, 1, NA),
                                                                   names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
                             type="binomial") 
    summary(SimulDrmLongLasting1)
    Eta50=SimulDrmLongLasting1$coefficients[2]
    Alpha=SimulDrmLongLasting1$coefficients[1]
    
    #pour faire des figures après?
    DAI = seq(0,166, 1)
    ndat = data.frame(DAI2=DAI, prop_dead=NA)
    pred_Pdead_4J = predict(SimulDrmLongLasting1, newdata = ndat, interval="confidence")
    
    pred_Pdead_4J = cbind(pred_Pdead_4J,DAI)
    
    predprop_F1.0 = data.frame(pred_Pdead_4J) %>%
      mutate(Formulation = "Lam-1.0")  # predprop_F1.0, 1.0 pour 1mg/kg
    predprop=predprop_F1.0
    
    sg_tot = sg_tot %>%mutate(DAI2 = as.numeric(DAI2))
    sg_tot = sg_tot %>% mutate(Day = as.numeric(Day))
    sub_PK= sg_tot
  }
  
  #fitting parameter for bohemia and rimdamal
  # Oral_Formulation = read.csv("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/DataIVM/Oral_formulation.csv")
  Oral_Formulation = read.csv("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\dataIVM\\Oral formulation\\Oral_formulation.csv")
  # View(Oral_Formulation)
  #gam_all = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/DataIVM/pred_gam_all.rds")
  gam_all = readRDS("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\dataIVM\\pred_gam_all.rds")
  
  if(strategy==2){#Pour BOHEMIA
    BOHEMIA= Oral_Formulation %>%filter(Article == "BOHEMIA")
    BOHEMIA["DAI"] = as.numeric(BOHEMIA$DAI)
    SimulDrmBohemia= drm(prop_dead ~ DAI, data = BOHEMIA,
                         fct = LL.4(fixed = c(NA, 0, NA, NA),
                                    names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    summary(SimulDrmBohemia)
    Eta50=SimulDrmBohemia$coefficients[2]
    Alpha=SimulDrmBohemia$coefficients[1]
    
    #pour faire des figures après?
    #plot(EDM_BOHEMIA)
    DAI = seq(0,30, 1)
    ndat = data.frame(DAI=DAI, prop_dead=NA)
    
    pred_BOHEMIA = predict(SimulDrmBohemia, newdata = ndat, interval = "confidence", level = 0.45)
    
    pred_BOHEMIA = data.frame(pred_BOHEMIA) %>%
      cbind(ndat) %>%
      mutate(Formulation = "BOHEMIA")
    predprop=pred_BOHEMIA
    
    sub_PK=BOHEMIA %>% mutate(Formulation = "BOHEMIA")
  }
  
  if(strategy==3){#Pour RIMDAMAL
    RIMDAMAL= Oral_Formulation %>%filter(Article == "RIMDAMAL")
    RIMDAMAL["DAI"] = as.numeric(RIMDAMAL$DAI)
    SimulDrmRimdamal= drm(prop_dead ~ DAI, data = RIMDAMAL,
                          fct = LL.4(fixed = c(NA, 0, 0.8632, NA),
                                     names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    summary(SimulDrmRimdamal)
    Eta50=SimulDrmRimdamal$coefficients[2]
    Alpha=SimulDrmRimdamal$coefficients[1]
    
    #plot(EDM_RIMDAMAL)
    DAI = seq(0,30, 1)
    ndat = data.frame(DAI=DAI, prop_dead=NA)
    
    pred_RIMDAMAL=predict(SimulDrmRimdamal,newdata=ndat,interval="confidence",level = 0.7)
    
    pred_RIMDAMAL = data.frame(pred_RIMDAMAL) %>%
      cbind(ndat) %>%
      mutate(Formulation = "RIMDAMAL")
    predprop=pred_RIMDAMAL
    
    sub_PK=RIMDAMAL %>% mutate(Formulation = "RIMDAMAL")
  }
  
  FunList = list("eta50" = Eta50, "alpha"=Alpha, "predprop"=predprop,
                 "sub_PK"=sub_PK,"DAI"=DAI)
  return(FunList)
}

RhoFunction<-function(Nah,Ntau,Vtau,Vah,IVM_age_dependancy,
                      strategy,ProbaDeathThreshold){
  ParamStrategy=ModelEeta50Alpha(strategy)
  eta50=ParamStrategy$eta50 
  alpha=ParamStrategy$alpha
  
  
  rho=matrix(0,nrow = Nah, ncol = Ntau)
  
  # if(strategy==2|strategy==3){#Bohemia or Rimdamal
  #   #On cherche eta0, tel que mum_ivm(eta)=4*mum
  #   y0=1-exp(-3*mu_m)
  #   eta0=eta50*(1/y0-1)^(1/alpha)
  # }
  # 
  # if(strategy==0|strategy==1){#LongLasting 0.6 or 1
  #   #On cherche eta0, tel que PropDeathMosqui(eta)=ProbaDeathThreshold
  #   y0=ProbaDeathThreshold
  #   eta0=eta50*(1/y0-1)^(1/alpha)
  # }
  
  y0=ProbaDeathThreshold
  eta0=eta50*(1/y0-1)^(1/alpha)
  eta_eff=eta0
  
  if(IVM_age_dependancy){
    for (a in 1:Nah) {
      for (j in 1:Ntau) {
        if(a>18){rho[a,j] <- 10*(Vtau[j]>eta_eff)}
        if(a<=18 && a>10){rho[a,j] <- 10*(Vtau[j]>(0.5*eta_eff))}
        if(a<=10 && a>5){rho[a,j] <- 10*(Vtau[j]>(0.25*eta_eff))}
        if(a<=5){rho[a,j] <- 0}
      }
    }
  } else {
    for (a in 1:Nah) {
      for (j in 1:Ntau) {
        rho[a,j] <- 10*(Vtau[j]>eta_eff)
      }
    }
  }
  
  #On a 4 groupes d'ages <5(groupe 1); 5-10(2), 10-18(3) et >18(groupe 4)
  #On calcule les proba d'être encore sous IVM pour chaque group
  PropHuman_in_IVM_Group1=Vtau; IdAgeTemoin=3
  PropHuman_in_IVM_Group1[1]=1
  for (j in 2:Ntau) {
    Id=1:j; xId=Vtau[Id]; yId=rho[IdAgeTemoin,Id]
    PropHuman_in_IVM_Group1[j]=exp(-trapz(xId,yId))
  }
  
  PropHuman_in_IVM_Group2=Vtau; IdAgeTemoin=7
  PropHuman_in_IVM_Group2[1]=1
  for (j in 2:Ntau) {
    Id=1:j; xId=Vtau[Id]; yId=rho[IdAgeTemoin,Id]
    PropHuman_in_IVM_Group2[j]=exp(-trapz(xId,yId))
  }
  
  PropHuman_in_IVM_Group3=Vtau; IdAgeTemoin=12
  PropHuman_in_IVM_Group3[1]=1
  for (j in 2:Ntau) {
    Id=1:j; xId=Vtau[Id]; yId=rho[IdAgeTemoin,Id]
    PropHuman_in_IVM_Group3[j]=exp(-trapz(xId,yId))
  }
  
  PropHuman_in_IVM_Group4=Vtau; IdAgeTemoin=20
  PropHuman_in_IVM_Group4[1]=1
  for (j in 2:Ntau) {
    Id=1:j; xId=Vtau[Id]; yId=rho[IdAgeTemoin,Id]
    PropHuman_in_IVM_Group4[j]=exp(-trapz(xId,yId))
  }
  
  FunList = list("eta_eff" = eta_eff, "rho"=rho,
                 "PropHuman_in_IVM_Group1"=PropHuman_in_IVM_Group1,
                 "PropHuman_in_IVM_Group2"=PropHuman_in_IVM_Group2,
                 "PropHuman_in_IVM_Group3"=PropHuman_in_IVM_Group3,
                 "PropHuman_in_IVM_Group4"=PropHuman_in_IVM_Group4)
  return(FunList)
}

TheIvmStrategy<-function(t_begin_Camp,Number_of_cycle,VectTime_between_cycles,
                         Dur_cycle,Gap,time,Ntime){
  
  phi=rep(0,Ntime)
  if (Number_of_cycle >= 1) {
    for (t in 1:Ntime) {
      for (k in 0:(Number_of_cycle - 1)) { 
        if (k == 0) { ta = t_begin_Camp}
        else {ta = t_begin_Camp + sum(Dur_cycle + VectTime_between_cycles[1:k]) }
        tb = ta + Dur_cycle - 1
        if (ta - Gap <= time[t] & time[t] < ta) { phi[t] = (time[t]-ta+Gap)/Gap } 
        else if (ta <= time[t] & time[t] <= tb) { phi[t] = 1 } 
        else if (tb < time[t] & time[t] <= tb + Gap){phi[t] = (tb + Gap - time[t])/Gap}
      }
    }
  }
  
  FunList = list("PhiIvm" = phi)
  return(FunList)
}

ModelIVM<-function(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                   mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                   ah_max,dah,Nah,tmax,
                   Gap,dt,time,Ntime,
                   d_eta,eta_max,Veta,Neta,
                   dtau,tau_max,Vtau,Ntau,
                   strategy,ProbaDeathThreshold,IVM_age_dependancy){
  Wedge_h=30750 #30750, recrutment rate Bobo dioulasso in 2012
  q=0.48 #proportion of human male
  k_h=5.5*10^(-4)#1.469e-2;
  
  ParamIvmFormulation=ModelEeta50Alpha(strategy)
  p_vec = P_eta.fun(Veta, eta50=ParamIvmFormulation$eta50, alpha=ParamIvmFormulation$alpha)
  mu_m.ivm = mu_m - log(1 - p_vec)
  
  ModelRho=RhoFunction(Nah,Ntau,Vtau,Vah,IVM_age_dependancy,
                       strategy,ProbaDeathThreshold)
  rho=ModelRho$rho
  
  #phi(t,a)
  Modelphi=TheIvmStrategy(t_begin_Camp,Number_of_cycle,VectTime_between_cycles,
                          Dur_cycle,Gap,time,Ntime)
  phi=Modelphi$PhiIvm
  {
    mphi=matrix(0,nrow = Ntime,ncol = Nah);#IVM exposition rate for human male
    fphi=matrix(0,nrow = Ntime,ncol = Nah);#IVM exposition rate for human female
    
    phi=(PropIVM/Dur_cycle)*phi
    for (t in 1:Ntime) {
      for (a in 1:Nah) {
        if(Vah[a]<=5){mphi[t,a]=0; fphi[t,a]=0}
        else if (Vah[a] >= 15 & Vah[a] <= 45){mphi[t,a]=q*phi[t]; fphi[t,a]=0}
        else{mphi[t,a]=q*phi[t]; fphi[t,a]=(1-q)*phi[t]}
      }
    }
  }
  
  # mSh is Susceptible human male and fSh Susceptible human female
  mSh=matrix(0,nrow = Ntime,ncol = Nah);fSh=matrix(0,nrow = Ntime,ncol = Nah)
  mAh=matrix(0,nrow = Ntime,ncol = Nah);fAh=matrix(0,nrow = Ntime,ncol = Nah)
  mIh=matrix(0,nrow = Ntime,ncol = Nah);fIh=matrix(0,nrow = Ntime,ncol = Nah)
  mRh=matrix(0,nrow = Ntime,ncol = Nah);fRh=matrix(0,nrow = Ntime,ncol = Nah)
  Sm=rep(0,Ntime); Im=rep(0,Ntime);
  Sm_ivm=matrix(0,nrow = Ntime,ncol = Neta);Im_ivm=matrix(0,nrow = Ntime,ncol = Neta)
  mSh_ivm=array(0, dim = c(Nah,Ntau,Ntime));fSh_ivm=array(0, dim = c(Nah,Ntau,Ntime))
  mAh_ivm=array(0, dim = c(Nah,Ntau,Ntime));fAh_ivm=array(0, dim = c(Nah,Ntau,Ntime))
  mIh_ivm=array(0, dim = c(Nah,Ntau,Ntime));fIh_ivm=array(0, dim = c(Nah,Ntau,Ntime))
  mRh_ivm=array(0, dim = c(Nah,Ntau,Ntime));fRh_ivm=array(0, dim = c(Nah,Ntau,Ntime))
  Nh=rep(0,Ntime);Nm=rep(0,Ntime)
  
  
  # Initial values (t=0), Ref Quentin Richard et al Bobo Dioulasso 2012
  {
    Int_values = 8136.10 * (c(12.9, 12.5, 11.5, 13.1, 11.9, 9.3, 7.3, 5.5, 4.4,
                              3.2, 2.6, 1.8, 1.4, 0.9, 0.7, 0.3, 0.2, 0.2) / 0.997)
    # apparently the data only count for 99.7 % of the population
    # Age labels corresponding to each range in mSh_values
    age_label_3 = c("0 to 5 years", "5 to 10 years", "10 to 15 years", 
                    "15 to 20 years", "20 to 25 years", "25 to 30 years", 
                    "30 to 35 years", "35 to 40 years", "40 to 45 years", 
                    "45 to 50 years", "50 to 55 years", "55 to 60 years", 
                    "60 to 65 years", "65 to 70 years", "70 to 75 years", 
                    "75 to 80 years", "80 to 85 years", "85 to 90 years")
    
    for (i in seq_along(age_label_3)) {
      if (i == 1) {
        # For the first age range (0 to 5 years)
        mSh[1, Vah <= 5] = 0.45*q*Int_values[i]
        fSh[1, Vah <= 5] = 0.45*(1-q)*Int_values[i]
        mAh[1, Vah <= 5] = 0.3*q*Int_values[i]
        fAh[1, Vah <= 5] = 0.3*(1-q)*Int_values[i]
        mIh[1, Vah <= 5] = 0.1*q*Int_values[i]
        fIh[1, Vah <= 5] = 0.1*(1-q)*Int_values[i]
        mRh[1, Vah <= 5] = 0.15*q*Int_values[i]
        fRh[1, Vah <= 5] = 0.15*(1-q)*Int_values[i]
      } else if (i == length(age_label_3)) {
        # For the last age range (> 85 years)
        mSh[1, Vah > 85] = 0.45*q*Int_values[i]
        fSh[1, Vah > 85] = 0.45*(1-q)*Int_values[i]
        mAh[1, Vah > 85] = 0.3*q*Int_values[i]
        fAh[1, Vah > 85] = 0.3*(1-q)*Int_values[i]
        mIh[1, Vah > 85] = 0.1*q*Int_values[i]
        fIh[1, Vah > 85] = 0.1*(1-q)*Int_values[i]
        mRh[1, Vah > 85] = 0.15*q*Int_values[i]
        fRh[1, Vah > 85] = 0.15*(1-q)*Int_values[i]
      } else {
        # For intermediate age ranges
        age_range <- which(Vah > (i - 1) * 5 & Vah <= i * 5)
        mSh[1, age_range] = 0.45*q*Int_values[i]
        fSh[1, age_range] = 0.45*(1-q)*Int_values[i]
        mAh[1, age_range] = 0.3*q*Int_values[i]
        fAh[1, age_range] = 0.3*(1-q)*Int_values[i]
        mIh[1, age_range] = 0.1*q*Int_values[i]
        fIh[1, age_range] = 0.1*(1-q)*Int_values[i]
        mRh[1, age_range] = 0.15*q*Int_values[i]
        fRh[1, age_range] = 0.15*(1-q)*Int_values[i]
      }
    }
  }
  
  Sm[1]=1.3*10^6#293#378#293 Wedge_m/mu_m;
  Im[1]=0.15*10^6#336#50;
  
  Sm_ivm[1,]=rep(0,Neta); Im_ivm[1,]=rep(0,Neta)
  mSh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fSh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  mAh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fAh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  mIh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fIh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  mRh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fRh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  
  t=1
  Nh[t]= (sum(mSh[t,]+fSh[t,]+mAh[t,]+fAh[t,]+mIh[t,]+fAh[t,]+mRh[t,]+mRh[t,])
          +sum(mSh_ivm[,,t]+fSh_ivm[,,t]+mAh_ivm[,,t]+fAh_ivm[,,t]+mIh_ivm[,,t]
               +fIh_ivm[,,t]+mRh_ivm[,,t]+fRh_ivm[,,t]))
  Nm[t]= Sm[t] + Im[t]+ sum(Sm_ivm[t,] +Im_ivm[t,])
  #Nbrun=0
  eta_eff_age= 0 #initialisation
  
  for (t in 1:(Ntime-1)) {
    lam_m = beta_m*theta[t] * (Im[t] + sum(Im_ivm[t,]))/Nh[t];
    lam_h = theta[t] * (sum( beta_h[]*(mAh[t,]+fAh[t,]) ) 
                     + sum( bar_beta_h[]*(mIh[t,]+fIh[t,])))/Nh[t];
    lam_hI =theta[t]*sum(mAh_ivm[,,t]+fAh_ivm[,,t]+mIh_ivm[,,t]+fIh_ivm[,,t])/Nh[t];
    lam_hS =theta[t]*sum(mSh_ivm[,,t]+fSh_ivm[,,t]+mRh_ivm[,,t]+fRh_ivm[,,t])/Nh[t];
    
    Sm[t+1] = (Wedge_m + Sm[t]/dt )/(1/dt  + mu_m + lam_hI + lam_h + lam_hS);
    Im[t+1] = (Im[t]/dt + lam_h*Sm[t] )/(1/dt  + mu_m + lam_hI +lam_hS);
    
    ### nouvelle naissance
    a=1 
    mSh[t+1,a]=(mSh[t,a]/dt + q*Wedge_h +k_h*mRh[t,a])/
      (1/dt + 1/dah +  mu_h[a] + mphi[t,a] + lam_m); 
    fSh[t+1,a]=(mSh[t,a]/dt + (1-q)*Wedge_h +k_h*mRh[t,a])/
      (1/dt + 1/dah +  mu_h[a] + fphi[t,a] + lam_m); 
    mAh[t+1,a]=0
    fAh[t+1,a]=0
    mIh[t+1,a]=0
    fAh[t+1,a]=0
    mRh[t+1,a]=0
    fRh[t+1,a]=0
    mSh_ivm[a,,t+1]=rep(0,Ntau);fSh_ivm[a,,t+1]=rep(0,Ntau)
    mAh_ivm[a,,t+1]=rep(0,Ntau);fAh_ivm[a,,t+1]=rep(0,Ntau)
    mIh_ivm[a,,t+1]=rep(0,Ntau);fIh_ivm[a,,t+1]=rep(0,Ntau)
    mRh_ivm[a,,t+1]=rep(0,Ntau);fRh_ivm[a,,t+1]=rep(0,Ntau)
    
    #on integre en age
    for (a in 2:Nah) {
      
      mSh[t+1,a]=(mSh[t,a]/dt +mSh[t+1,a-1]/dah +k_h*mRh[t,a]+sum(rho[a,]*mSh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + mphi[t,a] + lam_m); 
      fSh[t+1,a]=(fSh[t,a]/dt +fSh[t+1,a-1]/dah +k_h*fRh[t,a]+sum(rho[a,]*fSh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + fphi[t,a] + lam_m);
      mAh[t+1,a]=(mAh[t,a]/dt +mAh[t+1,a-1]/dah +lam_m*mSh[t,a]+sum(rho[a,]*mAh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + mphi[t,a] + nu_hm[a]);
      fAh[t+1,a]=(fAh[t,a]/dt +fAh[t+1,a-1]/dah +lam_m*fSh[t,a]+sum(rho[a,]*fAh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + fphi[t,a] + nu_hf[a]);
      mIh[t+1,a]=(mIh[t,a]/dt +mIh[t+1,a-1]/dah +nu_hm[a]*mAh[t,a]+sum(rho[a,]*mIh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + delta_h[a]  + gamma_h[a]);
      fIh[t+1,a]=(fIh[t,a]/dt +fIh[t+1,a-1]/dah +nu_hf[a]*fAh[t,a]+sum(rho[a,]*fIh_ivm[a,,t]))/
        (1/dt + 1/dah +  mu_h[a] + delta_h[a]  + gamma_h[a]);
      mRh[t+1,a]=(mRh[t,a]/dt +mRh[t+1,a-1]/dah +gamma_h[a]*mIh[t,a]+sum(rho[a,]*mRh_ivm[a,,t]))/
        (1/dt + 1/dah + mu_h[a] +mphi[t,a] + k_h);
      fRh[t+1,a] = (fRh[t,a]/dt + fRh[t+1,a-1]/dah + gamma_h[a]*fIh[t,a] + sum(rho[a,]*fRh_ivm[a,,t])) /
        (1/dt + 1/dah + mu_h[a] + fphi[t,a] + k_h)
      
      #nouvelle IVM d'age a
      tau=1
      mSh_ivm[a,tau,t+1]=(mSh_ivm[a,tau,t]/dt + mphi[t,a]*mSh[t,a]
                          +mSh_ivm[a-1,tau,t+1]/dah +k_h*mRh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + lam_m + rho[a,tau]);
      fSh_ivm[a,tau,t+1]=(fSh_ivm[a,tau,t]/dt + fphi[t,a]*fSh[t,a] 
                          +fSh_ivm[a-1,tau,t+1]/dah +k_h*fRh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + lam_m + rho[a,tau]);
      mAh_ivm[a,tau,t+1]=(mAh_ivm[a,tau,t]/dt + mphi[t,a]*mAh[t,a] 
                          +mAh_ivm[a-1,tau,t+1]/dah +lam_m*mSh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + nu_hm[a] + rho[a,tau]);
      fAh_ivm[a,tau,t+1]=(fAh_ivm[a,tau,t]/dt + fphi[t,a]*fAh[t,a] 
                          +fAh_ivm[a-1,tau,t+1]/dah +lam_m*fSh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + nu_hf[a] + rho[a,tau]);
      mIh_ivm[a,tau,t+1]=(mIh_ivm[a,tau,t]/dt + 0 +mIh_ivm[a-1,tau,t+1]/dah 
                          +nu_hm[a]*mAh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + gamma_h[a]+delta_h[a] + rho[a,tau]);
      fIh_ivm[a,tau,t+1]=(fIh_ivm[a,tau,t]/dt + 0 +fIh_ivm[a-1,tau,t+1]/dah 
                          +nu_hf[a]*fAh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + gamma_h[a]+delta_h[a] + rho[a,tau]);
      mRh_ivm[a,tau,t+1]=(mRh_ivm[a,tau,t]/dt + mphi[t,a]*mRh[t,a] 
                          +mRh_ivm[a-1,tau,t+1]/dah +gamma_h[a]*mIh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + k_h + rho[a,tau]);
      fRh_ivm[a,tau,t+1]=(fRh_ivm[a,tau,t]/dt + fphi[t,a]*fRh[t,a] 
                          +fRh_ivm[a-1,tau,t+1]/dah +gamma_h[a]*fIh_ivm[a,tau,t])/
        (1/dt + 1/dtau +1/dah + mu_h[a] + k_h + rho[a,tau]);
      
      #on integre les humains sous IVM
      for (tau in 2:Ntau) {
        mSh_ivm[a,tau,t+1]=(mSh_ivm[a,tau,t]/dt + mSh_ivm[a,tau-1,t+1]/dtau 
                            +mSh_ivm[a-1,tau,t+1]/dah +k_h*mRh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + lam_m + rho[a,tau]);
        fSh_ivm[a,tau,t+1]=(fSh_ivm[a,tau,t]/dt + fSh_ivm[a,tau-1,t+1]/dtau 
                            +fSh_ivm[a-1,tau,t+1]/dah +k_h*fRh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + lam_m + rho[a,tau]);
        mAh_ivm[a,tau,t+1]=(mAh_ivm[a,tau,t]/dt + mAh_ivm[a,tau-1,t+1]/dtau 
                            +mAh_ivm[a-1,tau,t+1]/dah +lam_m*mSh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + nu_hm[a] + rho[a,tau]);
        fAh_ivm[a,tau,t+1]=(fAh_ivm[a,tau,t]/dt + fAh_ivm[a,tau-1,t+1]/dtau 
                            +fAh_ivm[a-1,tau,t+1]/dah +lam_m*fSh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + nu_hf[a] + rho[a,tau]);
        mIh_ivm[a,tau,t+1]=(mIh_ivm[a,tau,t]/dt + mIh_ivm[a,tau-1,t+1]/dtau 
                            +mIh_ivm[a-1,tau,t+1]/dah +nu_hm[a]*mAh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + gamma_h[a]+delta_h[a] + rho[a,tau]);
        fIh_ivm[a,tau,t+1]=(fIh_ivm[a,tau,t]/dt + fIh_ivm[a,tau-1,t+1]/dtau 
                            +fIh_ivm[a-1,tau,t+1]/dah +nu_hf[a]*fAh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + gamma_h[a]+delta_h[a] + rho[a,tau]);
        mRh_ivm[a,tau,t+1]=(mRh_ivm[a,tau,t]/dt + mRh_ivm[a,tau-1,t+1]/dtau 
                            +mRh_ivm[a-1,tau,t+1]/dah +gamma_h[a]*mIh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + k_h + rho[a,tau]);
        fRh_ivm[a,tau,t+1]=(fRh_ivm[a,tau,t]/dt + fRh_ivm[a,tau-1,t+1]/dtau 
                            +fRh_ivm[a-1,tau,t+1]/dah +gamma_h[a]*fIh_ivm[a,tau,t])/
          (1/dt + 1/dtau +1/dah + mu_h[a] + k_h + rho[a,tau]);
        
      }
    }
    
    eta=1
    Sm_ivm[t+1,eta]=(Sm_ivm[t,eta]/dt + lam_hS*Sm[t]/d_eta )/
      (1/dt+1/d_eta + mu_m.ivm[eta] + lam_hI +lam_h);
    Im_ivm[t+1,eta]=(Im_ivm[t,eta]/dt + (lam_hI*Sm[t] +lam_hS*Im[t])/d_eta +
                       (lam_hI + lam_h)*Sm_ivm[t,eta])/(1/dt + 1/d_eta + mu_m.ivm[eta] );
    
    for(eta in 2:Neta){
      Sm_ivm[t+1,eta]=(Sm_ivm[t,eta]/dt + Sm_ivm[t+1,eta-1]/d_eta )/
        (1/dt+1/d_eta + mu_m.ivm[eta] + lam_hI +lam_h);
      Im_ivm[t+1,eta]=(Im_ivm[t,eta]/dt + Im_ivm[t+1,eta-1]/d_eta +
                         (lam_hI + lam_h)*Sm_ivm[t,eta])/(1/dt + 1/d_eta + mu_m.ivm[eta] );
    }
    
    
    Nh[t+1]= (sum(mSh[t+1,]+fSh[t+1,]+mAh[t+1,]+fAh[t+1,]+mIh[t+1,]+fAh[t+1,]
                  +mRh[t+1,]+mRh[t+1,]) +sum(mSh_ivm[,,t+1]+fSh_ivm[,,t+1] +mAh_ivm[,,t+1]
                                             +fAh_ivm[,,t+1]+mIh_ivm[,,t+1]+fIh_ivm[,,t+1]+mRh_ivm[,,t+1]+fRh_ivm[,,t+1]))
    
    Nm[t+1]= Sm[t+1] + Im[t+1]+ sum(Sm_ivm[t+1,] +Im_ivm[t+1,])
    
    #Nbrun=Nbrun+1
    #print(Nbrun) #le max etant Ntime
  }  
  
  Sh_ivmTot=time;Ah_ivmTot=time;Ih_ivmTot=time;Rh_ivmTot=time
  Sh_Tot=time;Ah_Tot=time;Ih_Tot=time;Rh_Tot=time;
  Sm_ivmTot=time;Im_ivmTot=time
  
  mIh_a=matrix(0,nrow = Ntime,ncol = Nah); fIh_a=matrix(0,nrow = Ntime,ncol = Nah)
  mAh_a=matrix(0,nrow = Ntime,ncol = Nah); fAh_a=matrix(0,nrow = Ntime,ncol = Nah)
  for (t in 1:Ntime){
    Sh_Tot[t]=sum(mSh[t,] + fSh[t,])
    Ah_Tot[t]=sum(mAh[t,] + fAh[t,])
    Ih_Tot[t]=sum(mIh[t,] + fIh[t,])
    Rh_Tot[t]=sum(mRh[t,] + fRh[t,])
    Sh_ivmTot[t]=sum(mSh_ivm[,,t] +fSh_ivm[,,t])
    Ah_ivmTot[t]=sum(mAh_ivm[,,t] +fAh_ivm[,,t])
    Ih_ivmTot[t]=sum(mIh_ivm[,,t] +fIh_ivm[,,t])
    Rh_ivmTot[t]=sum(mRh_ivm[,,t] +fRh_ivm[,,t])
    
    Sm_ivmTot[t]=sum(Sm_ivm[t,])
    Im_ivmTot[t]=sum(Im_ivm[t,])
    
    for (a in 1:Nah){
      mIh_a[t,a]=mIh[t,a] + sum(mIh_ivm[a,,t]); fIh_a[t,a]=fIh[t,a] + sum(fIh_ivm[a,,t])
      mAh_a[t,a]=mAh[t,a] + sum(mAh_ivm[a,,t]); fAh_a[t,a]=fAh[t,a] + sum(fAh_ivm[a,,t])
    }
  }
  
  IdBegin_Camp= 1+floor(t_begin_Camp/dt)
  IndexTime=IdBegin_Camp:Ntime
  # Delta=trapz(time,Ih+Ih_ivmTot)/(tmax*Ih[1]+ Ih_ivmTot[1])
  # Delta=1-Delta
  Delta=trapz(time[IndexTime],Ih_Tot[IndexTime]+Ih_ivmTot[IndexTime])
  
  PropIhTot=(Ih_Tot +Ih_ivmTot)/Nh
  IdTop=which(PropIhTot[IndexTime]==min(PropIhTot[IndexTime]))
  DurIvmEffect=time[IdTop]
  
  MinPropIhTot=min(PropIhTot[IndexTime])
  PropIhTot0=PropIhTot[t_begin_Camp/dt]
  
  #pour les barplots
  
   cut_age = seq(0, ah_max, by = 5)
  # 
  # #Regrouper mIh_a et fIh_a en classes de 5 ans et le prend a l'instant min(Ih)
  # mIh_a = sapply(1:(length(cut_age) - 1), function(i) {
  #   idx_debut = (cut_age[i] / dah) + 1
  #   idx_fin = (cut_age[i + 1] / dah)
  #   sum(mIh_a[IdBegin_Camp+IdTop, idx_debut:idx_fin])
  # })
  # 
  # fIh_a = sapply(1:(length(cut_age) - 1), function(i) {
  #   idx_debut = (cut_age[i] / dah) + 1
  #   idx_fin = (cut_age[i + 1] / dah)
  #   sum(fIh_a[IdBegin_Camp+IdTop, idx_debut:idx_fin])
  # })
  # 
  
  FunList = list("time" = time, "rho"=rho,"mum.ivm"=mu_m.ivm,
                 "mSh"=mSh,"Sh_Tot"=Sh_Tot,"Ah_Tot"=Ah_Tot,"Ih_Tot"=Ih_Tot,"Rh_Tot"=Rh_Tot,
                 "Nh"=Nh,"Nm"=Nm, "Sh_ivmTot"=Sh_ivmTot,"Ah_ivmTot"=Ah_ivmTot,
                 "Ih_ivmTot"=Ih_ivmTot, "Rh_ivmTot"=Rh_ivmTot,"Sm"=Sm,"Im"=Im,
                 "Sm_ivmTot"=Sm_ivmTot,"Im_ivmTot"=Im_ivmTot,"Delta"=Delta,"lam_h"=lam_h,
                 "lam_m"=lam_m, "lam_hI"=lam_hI, "lam_hS"=lam_hS, "phi"=phi,
                 "mIh_a"=mIh_a, "fIh_a"=fIh_a,"mAh_a"=mAh_a, "fAh_a"=fAh_a, 
                 "PropIhTot"=PropIhTot, "DurIvmEffect"=DurIvmEffect,
                 "MinPropIhTot"=MinPropIhTot, "PropIhTot0"=PropIhTot0, "cut_age"=cut_age)
  return(FunList)
}

#parametres de base
{
  Gap=0.025#0.25;
  dt=Gap;tmax = 540;time=seq(0,tmax,by=dt); Ntime=length(time)
  d_eta=1;eta_max= 200;Veta=seq(0,eta_max,by=d_eta);Neta=length(Veta);
  dah=1;ah_max=90;Vah=seq(0,ah_max,by=dah); Nah=length(Vah)
  dtau=d_eta;tau_max= eta_max;Vtau=seq(0,tau_max,by=dtau);Ntau=length(Vtau);
  mu_m = 0.13  
  
  Wedge_m=2.55*10^5
  beta_m=0.5#0.10;
  Dur_cycle=7 
  #t_begin_Camp=200
  
  #mu_h(a) mortalité naturelle 
  {
    # Values corresponding to age ranges, Ref Quentin Richard et al Bobo Dioulasso
    mu_h_values = c(66.8, 7.6, 1.7, 0.9, 1.3, 1.9, 2.4, 2.8, 3.6, 4.7, 6.3, 8.9,
                    13.2, 19.8, 31.1, 47.7, 71.3, 110.5, 186.7) / 1000 
    # Age labels corresponding to age ranges
    age_labels0=c(0,1,seq(5,90,by=5))
    age_labels = c("0 to 1 year", "1 to 5 years", "5 to 10 years", "10 to 15 years", 
                   "15 to 20 years", "20 to 25 years", "25 to 30 years", 
                   "30 to 35 years", "35 to 40 years", "40 to 45 years", 
                   "45 to 50 years", "50 to 55 years", "55 to 60 years", 
                   "60 to 65 years", "65 to 70 years", "70 to 75 years", 
                   "75 to 80 years", "80 to 85 years", "85 to 90 years")
    
    mu_h = rep(NA, Nah) 
    # Assign mu_h according to age_labels
    for (i in 1:Nah) {
      for (IdAgeGroup in 1:(length(age_labels0)-1)) {
        if(age_labels0[IdAgeGroup]<=Vah[i]&
           Vah[i]<age_labels0[IdAgeGroup+1])
        {mu_h[i]=mu_h_values[IdAgeGroup]}
      }
      if (i>=Nah) {mu_h[i]=mu_h_values[length(mu_h_values)]}
    }
    plot(Vah,mu_h)
  }
  
  #delta_h(a), Ref Quentin Richard et al Bobo Dioulasso. Motalité induite maladie
  {
    # Values of delta_h corresponding to the age ranges
    delta_h_values = c(1.07 * 10^(-3), 7.02 * 10^(-4), 4.55 * 10^(-4), 5.73 * 10^(-5))
    # Age labels for the age ranges
    Age_label_2 = c("0-1 years", "1-5 years", "5-15 years", ">15 years")
    # Create a vector to store delta_h according to age
    delta_h = rep(NA, Nah)
    # Assign nu_h according to Age_label_2
    for (i in 1:Nah) {
      if(0<=Vah[i]&Vah[i]<=1)
      {delta_h[i]=delta_h_values[1]}
      
      if(1<Vah[i]&Vah[i]<=5)
      {delta_h[i]=delta_h_values[2]}
      
      if(5<Vah[i]&Vah[i]<=15)
      {delta_h[i]=delta_h_values[3]}
      
      if(15<Vah[i])
      {delta_h[i]=delta_h_values[4]}
    }
    plot(Vah,delta_h)
  }
  
  #beta_h(a) : infectiosité humaine
  {
    alpha_1 =  0.095 #  0.071[0.023; 0.175]  
    alpha_2 = 0.17 # 0.302[0.16; 0.475] 
    beta_h = rep(NA, Nah)
    
    # Define the function G(a)
    G = function(a) {
      return(22.7 * a * exp(-0.0934 * a))
    }
    # Compute beta_h for each age in Vah
    for (i in 1:Nah) {
      beta_h[i] = alpha_1 * (G(Vah[i])  ^ alpha_2)  
    }
    bar_beta_h=0.8*beta_h
    
    plot(Vah,G(Vah))
    plot(Vah,bar_beta_h)
  }
  
  #nu_hm(a), nu_hf(a): progression Asympto vers cas clinique
  {
    p_f = 0.15 # Fertility probability for women aged 15 to 40 in Bobo Dioulasso,
    nu_0=1/10; nu_1=1/200;nu_2=1/3
    nu_hm=rep(NA,Nah) # disease mortality for men
    nu_hf = rep(NA, Nah) # disease mortality for female
    for (i in 1:Nah) {
      if (Vah[i] <= 5) { nu_hm[i] = nu_2; nu_hf[i] = nu_2 } 
      else if (Vah[i] > 5 && Vah[i] <= 15) {nu_hm[i] = nu_0; nu_hf[i] = nu_0 }
      else {nu_hm[i] = nu_1 }
      
      if (Vah[i] > 15 && Vah[i] <= 45) { nu_hf[i] = p_f * nu_0 + (1 - p_f) * nu_1 } 
      else if (Vah[i] > 40) { nu_hf[i] = nu_1}
    }
    
    
    plot(Vah,nu_hm)
    plot(Vah,nu_hf)
  }
  
  #gamma_h(a): guerison
  {
    # Recovery time ranges for each age group (days)
    recovery_ranges = list(c(163, 345), c(555, 714), c(344, 400), c(181, 204), 
                           c(82, 92), c(56, 61), c(48, 55))
    # Age ranges  corresponding to the recovery time ranges
    age_ranges =list(c(0, 1),c(1, 5),c(5, 8),c(8, 18),c(18, 28),c(28, 43),c(43, ah_max))
    gamma_h = rep(NA, Nah)
    
    #  value of gamma_h
    for (i in 1:Nah) {
      for (j in 1: length(age_ranges)) {
        age_min = age_ranges[[j]][1]
        age_max = age_ranges[[j]][2]
        if(Vah[i] >= age_min & Vah[i] <= age_max)
        {gamma_h[i] = 1 / mean(recovery_ranges[[j]])}
      }
    }
    plot(Vah,gamma_h)
  }
  
  p_theta.fun = function(time, theta_0, q_theta, T_theta, t_begin_Higly_Seas, Gap) {
    ta = t_begin_Higly_Seas
    tb = ta + T_theta
    
    theta = ifelse(time < ta - Gap, (1 - q_theta) * theta_0,
                   ifelse(time < ta, (1 - q_theta) * theta_0 + ((time - (ta - Gap)) / Gap) * (theta_0 - (1 - q_theta) * theta_0),
                          ifelse(time <= tb, theta_0,
                                 ifelse(time <= tb + Gap, theta_0 - ((time - tb) / Gap) * (theta_0 - (1 - q_theta) * theta_0),
                                        (1 - q_theta) * theta_0))))
    
    return(theta)
  }
  
  q_theta=1/3
  theta_0=0.5
  T_theta=120 #dure' de la haute saison de transmission
  t_begin_Higly_Seas=200
  
  Model_theta = function(Saisonality){
    if (Saisonality == 1) {
      theta = p_theta.fun(time, theta_0, q_theta, T_theta, t_begin_Higly_Seas, Gap)
    } else {
      theta = rep(theta_0, Ntime) }
    return(theta)
  }
  
}

#Figure:10 & 11 Scenarios_Formulations LAIF with Seasonality
###########################

Vect_t_begin_Camp=c(170, 200, 230)
t_begin_plot=160
CoefDelay=0.7
IVM_age_dependancy=0
ProbaDeathThreshold=0.75

ColVect = c("#ff3355", "#D98E04", "#5f05f3", "#037153", "#205072")
LineVect <- c(1, 3, 2, 4)
LC <- 2
sink(paste("DurEntreCycle_SeaS_age_dependancy",IVM_age_dependancy,".txt",sep=""))

#without IVM
PropIVM=0;ProbaDeathThreshold=0.75;
strategy=0;Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
Saisonality=1
theta=Model_theta(Saisonality)
ModelZero=list()
for (p in 1:length(Vect_t_begin_Camp)) {
  t_begin_Camp=Vect_t_begin_Camp[p]
ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                       mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                       ah_max,dah,Nah,tmax,
                       Gap,dt,time,Ntime,
                       d_eta,eta_max,Veta,Neta,
                       dtau,tau_max,Vtau,Ntau,
                       strategy,ProbaDeathThreshold,IVM_age_dependancy)
ModelZero[[p]]=ModelOutput
}

VectPropIVM=c(0.5,0.7,0.9)

#LAIF31-0.6
{
FigName = paste0("Scenario_SeaS_Ih_LAIF06_AgeDep",IVM_age_dependancy,"Eff_",
                 gsub("\\.","_",ProbaDeathThreshold) , ".pdf")
pdf(FigName, width = 12, height = 8.2)  
par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.5, 2, .1))
par(mfrow = c(length(VectPropIVM),length(Vect_t_begin_Camp)))

print(paste0("ProbaDeathThreshold=",ProbaDeathThreshold))
print(paste0("CoefDelay=",CoefDelay))

GrNumber=0

for (PropIVM in VectPropIVM) {
  for (p in 1:length(Vect_t_begin_Camp)) {
    
    t_begin_Camp=Vect_t_begin_Camp[p]
    
    print(paste0("Time begin cycle=",t_begin_Camp))
    print(paste0("PropIVM=",PropIVM))

    #On tourne le modèle pour LongLasting 06
    ModelLongLasting06=list()
    
    #On determine les temps entre les cycles (sans seasonality)
    {
      Saisonality=0
    theta=Model_theta(Saisonality)
    strategy=0;NbCycle=3
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,
                           Gap,dt,time,Ntime,
                           d_eta,eta_max,Veta,Neta,
                           dtau,tau_max,Vtau,Ntau,
                           strategy,ProbaDeathThreshold,IVM_age_dependancy)
    
    VectTime_between_cycles=rep(0,NbCycle-1)
    TempsNextCycle0=0
    for (Number_of_cycle in 2:NbCycle) {
      #On récupère le temps pour le cycle suivant
      TempsNextCycle=ModelOutput$DurIvmEffect #temps durée de décoirssance depuis begin camp
      d0=TempsNextCycle-TempsNextCycle0- Dur_cycle
      TempsNextCycle0=TempsNextCycle
      VectTime_between_cycles[Number_of_cycle-1]=CoefDelay*d0
      
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,
                             Gap,dt,time,Ntime,
                             d_eta,eta_max,Veta,Neta,
                             dtau,tau_max,Vtau,Ntau,
                             strategy,ProbaDeathThreshold,IVM_age_dependancy)
    }
    #on récupère les valeurs ici
    print("LongLasting0.6")
    print(VectTime_between_cycles)
    #VectTime_between_cyclesLongLasting06=0.7*c(65, 64)
    }
    
    #on dertmine maintenant ModelLongLasting06
    Saisonality=1
    theta=Model_theta(Saisonality)
    for (Number_of_cycle in 1:NbCycle) {
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,
                             Gap,dt,time,Ntime,
                             d_eta,eta_max,Veta,Neta,
                             dtau,tau_max,Vtau,Ntau,
                             strategy,ProbaDeathThreshold,IVM_age_dependancy)
      ModelLongLasting06[[Number_of_cycle]]=ModelOutput
    }
    
    
    MaxY=max(ModelZero[[p]]$PropIhTot)
    
    GrNumber=GrNumber+1
    plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_plot, tmax), ylab = "",
         ylim = c(0.05, MaxY), cex.lab = 1.2, yaxt = "n", xaxt = "n")
    axis(1, at = seq(t_begin_plot, tmax, by = 50), labels = seq(0, tmax-t_begin_plot, by = 50), 
         las = 1)
    y_vals = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
    axis(2, at = y_vals, labels = y_vals * 100, las = 1)
    par(xpd = NA)
    text(t_begin_plot-0.05, MaxY * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.6)
    par(xpd = FALSE)
    
    lines(time, ModelZero[[p]]$PropIhTot, lwd = LC, col = ColVect[1], lty = LineVect[1])
    
    rect(t_begin_Higly_Seas, par("usr")[3], t_begin_Higly_Seas+T_theta, par("usr")[4], 
         col = adjustcolor("#D3D3D3" , alpha.f = 0.4), border = NA)
    
    abline(v = t_begin_Camp, col=adjustcolor("black", alpha.f = 0.8),lty = 1, lwd=2)
    
    Nbb=1
    for (nc in 1:3) {
      Nbb=Nbb+1
      lines(time, ModelLongLasting06[[nc]]$PropIhTot, lwd = LC, cex.lab = 1.2, 
            lty = LineVect[Nbb], col = ColVect[Nbb])
      # points(time[IdSeq] , ModelLongLasting06[[nc]]$PropIhTot[IdSeq], lwd = 0.7, 
      #        cex.lab = 0.4,   col = ColVect[Nbb], pch = pchVect[1])
    }
    
    Peak_season = ifelse(p == 1, "IVM starts 30 days before peak season", 
                         ifelse(p == 2, "IVM starts at peak season", 
                                "IVM starts 30 days after peak season")) 
    
    if(is.element(GrNumber,c(1,2,3))){
      mtext(Peak_season,side = 3,adj = 0.5,cex = 0.8,line = 1,font = 1)
    }
    
    if(is.element(GrNumber,c(1,4,7))){
      mtext("Clinical cases (%)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8) 
      mtext(bquote(Prop[IVM] == .(PropIVM)),side = 2,adj = 0.5,cex = 0.9,line = 3,font = 1)
    }
    
    if(is.element(GrNumber,c(7,8,9))){
      mtext("Time (Days)", side = 1, adj = 0.5, cex = 0.9, line = 2, font = 0.8)
    }
  }
  
}
#plot legend
#LEGEND=c("Without IVM","BOHEMIA", "RIMDAMAL II")
LegendTex=c("Without IVM",expression(n[c] == 1),
            expression(n[c] == 2),expression(n[c] == 3))

reset()
legend("bottom", legend=c(bquote("LAIF31-0.6" ~ "      " * "Eff" == .(ProbaDeathThreshold*100)*"%")), 
       xpd = NA, horiz = TRUE, inset = c(0,0.045), bty = "n", cex = 1)

legend("bottom", legend=LegendTex,
       xpd = NA, horiz = TRUE, inset = c(0,0.02), 
       bty = "n", lty = LineVect, lw=3,
       col = ColVect, cex = .9)

dev.off()
getwd()
}
#LAIF31-1.0
{
  FigName = paste0("Scenario_SeaS_Ih_LAIF1_0_AgeDep",IVM_age_dependancy,"Eff_",
                   gsub("\\.","_",ProbaDeathThreshold) , ".pdf")
  pdf(FigName, width = 12, height = 8.2)  
  par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.5, 2, .1))
  par(mfrow = c(length(VectPropIVM),length(Vect_t_begin_Camp)))
  
  GrNumber=0
  for (PropIVM in VectPropIVM) {
    for (p in 1:length(Vect_t_begin_Camp)) {
      
      t_begin_Camp=Vect_t_begin_Camp[p]
      
      print(paste0("Time begin cycle=",t_begin_Camp))
      print(paste0("PropIVM=",PropIVM))
      
      #On tourne le modèle pour LongLasting 1
      ModelLongLasting1=list()
      #On determine les temps entre les cycles (sans seasonality)
      {
        Saisonality=0
        theta=Model_theta(Saisonality)
        strategy=1;NbCycle=3
        Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
        ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                               mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                               ah_max,dah,Nah,tmax,
                               Gap,dt,time,Ntime,
                               d_eta,eta_max,Veta,Neta,
                               dtau,tau_max,Vtau,Ntau,
                               strategy,ProbaDeathThreshold,IVM_age_dependancy)
        
        VectTime_between_cycles=rep(0,NbCycle-1)
        TempsNextCycle0=0
        for (Number_of_cycle in 2:NbCycle) {
          #On récupère le temps pour le cycle suivant
          TempsNextCycle=ModelOutput$DurIvmEffect #temps durée de décoirssance depuis begin camp
          d0=TempsNextCycle-TempsNextCycle0- Dur_cycle
          TempsNextCycle0=TempsNextCycle
          VectTime_between_cycles[Number_of_cycle-1]=CoefDelay*d0
          
          ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                                 mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                                 ah_max,dah,Nah,tmax,
                                 Gap,dt,time,Ntime,
                                 d_eta,eta_max,Veta,Neta,
                                 dtau,tau_max,Vtau,Ntau,
                                 strategy,ProbaDeathThreshold,IVM_age_dependancy)
        }
        #on récupère les valeurs ici
        print("LoongLasting1")
        print(VectTime_between_cycles)
      }
      
      #on dertmine maintenant ModelLongLasting1
      Saisonality=1
      theta=Model_theta(Saisonality)
      for (Number_of_cycle in 1:NbCycle) {
        ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                               mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                               ah_max,dah,Nah,tmax,
                               Gap,dt,time,Ntime,
                               d_eta,eta_max,Veta,Neta,
                               dtau,tau_max,Vtau,Ntau,
                               strategy,ProbaDeathThreshold,IVM_age_dependancy)
        ModelLongLasting1[[Number_of_cycle]]=ModelOutput
      }
      
      MaxY=max(ModelZero[[p]]$PropIhTot)
      
      GrNumber=GrNumber+1
      plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_plot, tmax), ylab = "",
           ylim = c(0.05, MaxY), cex.lab = 1.2, yaxt = "n", xaxt = "n")
      axis(1, at = seq(t_begin_plot, tmax, by = 50), labels = seq(0, tmax-t_begin_plot, by = 50), 
           las = 1)
      y_vals = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
      axis(2, at = y_vals, labels = y_vals * 100, las = 1)
      par(xpd = NA)
      text(t_begin_plot-0.05, MaxY * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.6)
      par(xpd = FALSE)
      
      lines(time, ModelZero[[p]]$PropIhTot, lwd = LC, col = ColVect[1], lty = LineVect[1])
      
      rect(t_begin_Higly_Seas, par("usr")[3], t_begin_Higly_Seas+T_theta, par("usr")[4], 
           col = adjustcolor("#D3D3D3" , alpha.f = 0.4), border = NA)
      abline(v = t_begin_Camp, col=adjustcolor("black", alpha.f = 0.8),lty = 1, lwd=2)
      Nbb=1
      for (nc in 1:3) {
        Nbb=Nbb+1
        lines(time, ModelLongLasting1[[nc]]$PropIhTot, lwd = LC, cex.lab = 1.2, 
              lty = LineVect[Nbb], col = ColVect[Nbb])
        # points(time[IdSeq] , ModelLongLasting1[[nc]]$PropIhTot[IdSeq], lwd = 0.7, 
        #        cex.lab = 0.4,   col = ColVect[Nbb], pch = pchVect[2])
      }
      Peak_season = ifelse(p == 1, "IVM starts 30 days before peak season", 
                           ifelse(p == 2, "IVM starts at peak season", 
                                  "IVM starts 30 days after peak season")) 
      if(is.element(GrNumber,c(1,2,3))){
        mtext(Peak_season,side = 3,adj = 0.5,cex = 0.8,line = 1,font = 1)
      }
      
      if(is.element(GrNumber,c(1,4,7))){
        mtext("Clinical cases (%)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8) 
        mtext(bquote(Prop[IVM] == .(PropIVM)),side = 2,adj = 0.5,cex = 0.9,line = 3,font = 1)
      }
      
      if(is.element(GrNumber,c(7,8,9))){
        mtext("Time (Days)", side = 1, adj = 0.5, cex = 0.9, line = 2, font = 0.8)
      }
    }
    
  }
  #plot legend
  #LEGEND=c("Without IVM","BOHEMIA", "RIMDAMAL II")
  LegendTex=c("Without IVM",expression(n[c] == 1),
              expression(n[c] == 2),expression(n[c] == 3))
  
  reset()
  legend("bottom", legend=c(bquote("LAIF31-1.0" ~ "      " * "Eff" == .(ProbaDeathThreshold*100)*"%")), 
         xpd = NA, horiz = TRUE, inset = c(0,0.045), bty = "n", cex = 1)
  
  legend("bottom", legend=LegendTex,
         xpd = NA, horiz = TRUE, inset = c(0,0.02), 
         bty = "n", lty = LineVect, lw=3,
         col = ColVect, cex = .9)
  
  dev.off()
  getwd()
}
#Oral Formulations
{
  FigName = paste0("Scenario_Eff_SeaS_Ih_oralF_AgeDep",IVM_age_dependancy,"Eff_",
                   gsub("\\.","_",ProbaDeathThreshold) , ".pdf")
  pdf(FigName, width = 12, height = 8.2)  
  par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.5, 2, .1))
  par(mfrow = c(length(VectPropIVM),length(Vect_t_begin_Camp)))
  
  GrNumber=0
  for (PropIVM in VectPropIVM) {
    for (p in 1:length(Vect_t_begin_Camp)) {
      t_begin_Camp=Vect_t_begin_Camp[p]
      
      print(paste0("Time begin cycle=",t_begin_Camp))
      print(paste0("PropIVM=",PropIVM))
      
      #On tourne le modèle pour Bohemia
      #On determine les temps entre les cycles (sans seasonality)
      {
        Saisonality=0
        theta=Model_theta(Saisonality)
        strategy=2;NbCycle=3
        Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
        ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                               mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                               ah_max,dah,Nah,tmax,
                               Gap,dt,time,Ntime,
                               d_eta,eta_max,Veta,Neta,
                               dtau,tau_max,Vtau,Ntau,
                               strategy,ProbaDeathThreshold,IVM_age_dependancy)
        VectTime_between_cycles=rep(0,NbCycle-1)
        TempsNextCycle0=0
        for (Number_of_cycle in 2:NbCycle) {
          #On récupère le temps pour le cycle suivant
          TempsNextCycle=ModelOutput$DurIvmEffect #temps durée de décoirssance depuis begin camp
          d0=TempsNextCycle-TempsNextCycle0- Dur_cycle
          TempsNextCycle0=TempsNextCycle
          VectTime_between_cycles[Number_of_cycle-1]=CoefDelay*d0
          
          ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                                 mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                                 ah_max,dah,Nah,tmax,
                                 Gap,dt,time,Ntime,
                                 d_eta,eta_max,Veta,Neta,
                                 dtau,tau_max,Vtau,Ntau,
                                 strategy,ProbaDeathThreshold,IVM_age_dependancy)
        }
        #on récupère les valeurs ici
        print("Bohemia")
        print(VectTime_between_cycles)
      }
      #on dertmine maintenant ModelBohemia
      Saisonality=1
      theta=Model_theta(Saisonality)
      Number_of_cycle=3
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,
                             Gap,dt,time,Ntime,
                             d_eta,eta_max,Veta,Neta,
                             dtau,tau_max,Vtau,Ntau,
                             strategy,ProbaDeathThreshold,IVM_age_dependancy)
      ModelBohemia=ModelOutput
      
      
      #On tourne le modèle pour Rimdamal
      #On determine les temps entre les cycles (sans seasonality)
      {
        Saisonality=0
        theta=Model_theta(Saisonality)
        strategy=3;NbCycle=4
        Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
        ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                               mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                               ah_max,dah,Nah,tmax,
                               Gap,dt,time,Ntime,
                               d_eta,eta_max,Veta,Neta,
                               dtau,tau_max,Vtau,Ntau,
                               strategy,ProbaDeathThreshold,IVM_age_dependancy)
        VectTime_between_cycles=rep(0,NbCycle-1)
        TempsNextCycle0=0
        for (Number_of_cycle in 2:NbCycle) {
          #On récupère le temps pour le cycle suivant
          TempsNextCycle=ModelOutput$DurIvmEffect #temps durée de décoirssance depuis begin camp
          d0=TempsNextCycle-TempsNextCycle0- Dur_cycle
          TempsNextCycle0=TempsNextCycle
          VectTime_between_cycles[Number_of_cycle-1]=CoefDelay*d0
          
          ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                                 mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                                 ah_max,dah,Nah,tmax,
                                 Gap,dt,time,Ntime,
                                 d_eta,eta_max,Veta,Neta,
                                 dtau,tau_max,Vtau,Ntau,
                                 strategy,ProbaDeathThreshold,IVM_age_dependancy)
        }
        #on récupère les valeurs ici
        print("Rimdamal")
        print(VectTime_between_cycles)
      }
      
      #on dertmine maintenant ModelRimdamal
      Saisonality=1
      theta=Model_theta(Saisonality)
      Number_of_cycle=4
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,
                             Gap,dt,time,Ntime,
                             d_eta,eta_max,Veta,Neta,
                             dtau,tau_max,Vtau,Ntau,
                             strategy,ProbaDeathThreshold,IVM_age_dependancy)
      ModelRimdamal=ModelOutput
      
      MaxY=max(ModelZero[[p]]$PropIhTot)
      
      GrNumber=GrNumber+1
      plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_plot, tmax), ylab = "",
           ylim = c(0.05, MaxY), cex.lab = 1.2, yaxt = "n", xaxt = "n")
      axis(1, at = seq(t_begin_plot, tmax, by = 50), labels = seq(0, tmax-t_begin_plot, by = 50), 
           las = 1)
      y_vals = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
      axis(2, at = y_vals, labels = y_vals * 100, las = 1)
      par(xpd = NA)
      text(t_begin_plot-0.05, MaxY * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.6)
      par(xpd = FALSE)
      
      lines(time, ModelZero[[p]]$PropIhTot, lwd = LC, col = ColVect[1], lty = LineVect[1])
      
      rect(t_begin_Higly_Seas, par("usr")[3], t_begin_Higly_Seas+T_theta, par("usr")[4], 
           col = adjustcolor("#D3D3D3" , alpha.f = 0.4), border = NA)
      
      abline(v = t_begin_Camp, col=adjustcolor("black", alpha.f = 0.8),lty = 1, lwd=2)
      
      lines(time, ModelBohemia$PropIhTot, lwd = LC, cex.lab = 1.2, 
            lty = LineVect[3], col = ColVect[2])
      lines(time, ModelRimdamal$PropIhTot, lwd = LC, cex.lab = 1.2, 
            lty = LineVect[4], col = ColVect[5])
      
      Peak_season = ifelse(p == 1, "IVM starts 30 days before peak season", 
                           ifelse(p == 2, "IVM starts at peak season", 
                                  "IVM starts 30 days after peak season")) 
      
      if(is.element(GrNumber,c(1,2,3))){
        mtext(Peak_season,side = 3,adj = 0.5,cex = 0.8,line = 1,font = 1)
      }
      
      if(is.element(GrNumber,c(1,4,7))){
        mtext("Clinical cases (%)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8) 
        mtext(bquote(Prop[IVM] == .(PropIVM)),side = 2,adj = 0.5,cex = 0.9,line = 3,font = 1)
      }
      
      if(is.element(GrNumber,c(7,8,9))){
        mtext("Time (Days)", side = 1, adj = 0.5, cex = 0.9, line = 2, font = 0.8)
      }
    }
  }
  LegendTex=c("Without IVM","BOHEMIA", "RIMDAMAL II")
  
  reset()
  legend("bottom", legend=LegendTex,
         xpd = NA, horiz = TRUE, inset = c(0,0.04), 
         bty = "n", lty =c(1, 2, 4, 3), lw=3,
         col = c("#ff3355", "#D98E04", "#205072", "#5f05f3", "#037153"), cex = .9)
  dev.off()
}
sink()
# sink(file = NULL, type = "output",append = FALSE)
