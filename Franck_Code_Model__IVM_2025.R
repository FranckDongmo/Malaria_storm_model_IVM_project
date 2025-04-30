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
library(gridExtra)
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
    lam_m = beta_m*theta * (Im[t] + sum(Im_ivm[t,]))/Nh[t];
    lam_h = theta * (sum( beta_h[]*(mAh[t,]+fAh[t,]) ) 
                     + sum( bar_beta_h[]*(mIh[t,]+fIh[t,])))/Nh[t];
    lam_hI =theta*sum(mAh_ivm[,,t]+fAh_ivm[,,t]+mIh_ivm[,,t]+fIh_ivm[,,t])/Nh[t];
    lam_hS =theta*sum(mSh_ivm[,,t]+fSh_ivm[,,t]+mRh_ivm[,,t]+fRh_ivm[,,t])/Nh[t];
    
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
  
  cut_age = seq(0, ah_max, by = 5)
  
  #Regrouper mIh_a et fIh_a en classes de 5 ans et le prend a l'instant min(Ih)
  mIh_a = sapply(1:(length(cut_age) - 1), function(i) {
    idx_debut = (cut_age[i] / dah) + 1
    idx_fin = (cut_age[i + 1] / dah)
    sum(mIh_a[IdBegin_Camp+IdTop, idx_debut:idx_fin])
  })
  
  fIh_a = sapply(1:(length(cut_age) - 1), function(i) {
    idx_debut = (cut_age[i] / dah) + 1
    idx_fin = (cut_age[i + 1] / dah)
    sum(fIh_a[IdBegin_Camp+IdTop, idx_debut:idx_fin])
  })
  
  
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
  dt=Gap;tmax = 580;time=seq(0,tmax,by=dt); Ntime=length(time)
  d_eta=1;eta_max= 200;Veta=seq(0,eta_max,by=d_eta);Neta=length(Veta);
  dah=1;ah_max=90;Vah=seq(0,ah_max,by=dah); Nah=length(Vah)
  dtau=d_eta;tau_max= eta_max;Vtau=seq(0,tau_max,by=dtau);Ntau=length(Vtau);
  mu_m = 0.13  
  
  Wedge_m=2.55*10^5
  beta_m=0.5#0.10;
  theta=0.5
  Dur_cycle=7 
  t_begin_Camp=200
  
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
}


#Figure:8  Scenarios_ Impact CoefDelay
###########################
VectCoefDelay=c(0.5, 0.75, 1)
IVM_age_dependancy=0
sink(paste("CoefDelay2_DurEntreCycleIVM_age_dep",IVM_age_dependancy,".txt",sep=""))

VectProbaDeathThreshold=c(0.5,0.6,0.75)
ColVect = c( "#ff3355", "#037153", "#D98E04","#205072" , "#5f05f3")
LineVect <- c(1, 2, 3, 4)
LC <- 2
# IdSeq <- seq(1, length(time), by = 80)
# pchVect<-c(1,8,16)

VectPropIVM=c(0.7)

for (PropIVM in VectPropIVM) {
  print(paste0("PropIVM=",PropIVM))
  
  for (ProbaDeathThreshold in VectProbaDeathThreshold) {
    print(paste0("ProbaDeathThreshold=",ProbaDeathThreshold))
    
    FigName = paste0("Scenario_Ih_Coef_Delay_AgeDep",IVM_age_dependancy,"Eff_",
                     gsub("\\.","_",ProbaDeathThreshold) , ".pdf")
    pdf(FigName, width = 12, height = 3.4)  
    par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.5, 2, .1))
    par(mfrow = c(length(VectPropIVM),3))
    
    GrNumber=0
    Ncp=0
    
    ModelBohemia=list()
    ModelRimdamal=list()
    ModelLongLasting06=list()
    ModelLongLasting1=list()
    for (CoefDelay in VectCoefDelay) {
      Ncp=Ncp+1
      print(paste0("CoefDelay=",CoefDelay))
      
      #On tourne le modèle pour Bohemia
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
      ModelBohemia[[Ncp]] =ModelOutput
      
      
      #On tourne le modèle pour Rimdamal
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
      #VectTime_between_cycles=0.7*c(39.50, 79.25-39.50, 79.25-(79.25-39.50))
      ModelRimdamal[[Ncp]]=ModelOutput
      
      
      #On tourne le modèle pour LongLasting 06
      strategy=0;NbCycle=3
      Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,
                             Gap,dt,time,Ntime,
                             d_eta,eta_max,Veta,Neta,
                             dtau,tau_max,Vtau,Ntau,
                             strategy,ProbaDeathThreshold,IVM_age_dependancy)
      
      while(length(ModelLongLasting06) < Number_of_cycle) {
        ModelLongLasting06[[length(ModelLongLasting06) + 1]] = list()
      }
      
      ModelLongLasting06[[Number_of_cycle]][[Ncp]]=ModelOutput
      
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
        
        while(length(ModelLongLasting06) < Number_of_cycle) {
          ModelLongLasting06[[length(ModelLongLasting06) + 1]] = list()
        }
        
        ModelLongLasting06[[Number_of_cycle]][[Ncp]]=ModelOutput
      }
      #on récupère les valeurs ici
      print("LongLasting0.6")
      print(VectTime_between_cycles)
      #VectTime_between_cyclesLongLasting06=0.7*c(65, 64)
      
      
      #On tourne le modèle pour LongLasting 1
      strategy=1;NbCycle=3
      Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,
                             Gap,dt,time,Ntime,
                             d_eta,eta_max,Veta,Neta,
                             dtau,tau_max,Vtau,Ntau,
                             strategy,ProbaDeathThreshold,IVM_age_dependancy)
      
      while(length(ModelLongLasting1) < Number_of_cycle) {
        ModelLongLasting1[[length(ModelLongLasting1) + 1]] = list()
      }
      
      ModelLongLasting1[[Number_of_cycle]][[Ncp]]=ModelOutput
      
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
        
        while(length(ModelLongLasting1) < Number_of_cycle) {
          ModelLongLasting1[[length(ModelLongLasting1) + 1]] = list()
        }
        
        ModelLongLasting1[[Number_of_cycle]][[Ncp]]=ModelOutput
      }
      #on récupère les valeurs ici
      print("LoongLasting1")
      print(VectTime_between_cycles)
    }
    
    PrevIh =ModelOutput$PropIhTot[t_begin_Camp/dt]
    GainPrev = rep(0.8 * PrevIh, Ntime)
    
    for (nc in 2:3) {
      GrNumber=GrNumber+1
      plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_Camp, tmax), ylab = "",
           ylim = c(0.05, 0.32), cex.lab = 1.2, yaxt = "n", xaxt = "n")
      axis(1, at = seq(t_begin_Camp, tmax, by = 50), labels = seq(0, tmax-t_begin_Camp, by = 50), 
           las = 1)
      y_vals = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
      axis(2, at = y_vals, labels = y_vals * 100, las = 1)
      par(xpd = NA)
      text(t_begin_Camp-0.05, 0.32 * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.6)
      par(xpd = FALSE)
      
      lines(time, GainPrev, lwd = .5, cex.lab = 1.2, lty = LineVect[3])
      
      for (Ncp in 1:length(VectCoefDelay)) { 
        
        lines(time, ModelLongLasting06[[nc]][[Ncp]]$PropIhTot, lwd = LC, cex.lab = 1.2, 
              lty = LineVect[Ncp], col = ColVect[Ncp])
        # points(time[IdSeq] , ModelLongLasting06[[nc]][[Ncp]]$PropIhTot[IdSeq], lwd = 0.7, 
        #           cex.lab = 0.4,   col = ColVect[Ncp], pch = pchVect[Ncp])
        
        
      }
      
      if(is.element(GrNumber,c(1,2))){
        mtext(bquote("LAIF31-0.6 when " ~ n[c] == .(nc)),
              side = 3, adj = 0.5, cex = 0.8, line = 1, font = 1)  
      }
      if(is.element(GrNumber,c(1))){
        mtext("Clinical cases (%)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8) 
        mtext(bquote(Eff == .(ProbaDeathThreshold * 100) * "%"), 
              side = 2, adj = 0.5, cex = 0.9, line = 3, font = 1)
      }
      mtext("Time (Days)", side = 1, adj = 0.5, cex = .9, line = 2, font = .8) 
      
    }
    
    GrNumber=GrNumber+1
    plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_Camp, tmax), ylab = "",
         ylim = c(0.05, 0.32), cex.lab = 1.2, yaxt = "n", xaxt = "n")
    axis(1, at = seq(t_begin_Camp, tmax, by = 50), labels = seq(0, tmax-t_begin_Camp, by = 50),
         las = 1)
    y_vals = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
    axis(2, at = y_vals, labels = y_vals * 100, las = 1)
    par(xpd = NA)
    text(t_begin_Camp-0.05, 0.32 * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.6)
    par(xpd = FALSE)
    
    lines(time, GainPrev, lwd = .5, cex.lab = 1.2, lty = LineVect[3])
    
    for (Ncp in 1:length(VectCoefDelay)) {
      
      # points(time[IdSeq] , ModelRimdamal[[Ncp]]$PropIhTot[IdSeq], lwd = 0.7, 
      #        cex.lab = 0.4,   col = ColVect[Ncp], pch = pchVect[Ncp])
      
      lines(time, ModelRimdamal[[Ncp]]$PropIhTot, lwd = LC, cex.lab = 1.2, 
            lty = LineVect[Ncp], col = ColVect[Ncp])
    }
    
    mtext("RIMDAMAL II", side = 3, adj = 0.5, cex = 0.8, line = 1, font = 1)
    
    # mtext("Clinical cases (%)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8)
    # mtext(bquote(Eff == .(ProbaDeathThreshold * 100) * "%"),
    #       side = 2, adj = 0.5, cex = 0.9, line = 3, font = 1)
    
    mtext("Time (Days)", side = 1, adj = 0.5, cex = .9, line = 2, font = .8)
    
    
    #plot legend
    LEGEND=c(expression(Coef[Delay] == 0.5),
             expression(Coef[Delay] == 0.75),
             expression(Coef[Delay] == 1))
    reset()
    
    # legend("bottom",  
    #        legend = c("LAIF31-0.6","RIMDAMAL II" ), inset = c(0,0.055),
    #        xpd = NA, horiz = TRUE, bty = "n", col = c("black", "black"), 
    #        pch = pchVect,  cex = .8)
    
    legend("bottom", legend=LEGEND, 
           xpd = NA, horiz = TRUE, inset = c(0,0.025), 
           bty = "n", lty = LineVect, lw=LC,
           col = ColVect,   cex = .9)
    
    dev.off() 
  }
}

sink()


#Figure: 4 & 5 Scenarios_Formulations Efficacy # 75%, 60%, 50%
###########################
CoefDelay=0.75
VectIVM_age_dependancy=c(0,1)
for (IVM_age_dependancy in VectIVM_age_dependancy) {
sink(paste("DurEntreCycleIVM_age_dependancy",IVM_age_dependancy,".txt",sep=""))

VectProbaDeathThreshold=c(0.75,0.6,0.5)
ColVect = c("#ff3355", "#037153", "#D98E04", "#205072")
LineVect <- c(1, 3, 2, 4)
LC <- 2

VectPropIVM=c(0.5,0.7,0.9)
for (PropIVM in VectPropIVM) {
  
  FigName = paste0("Scenario_Eff_Ih_AgeDependece",IVM_age_dependancy,"PropIVM_",
                   gsub("\\.","_",PropIVM) , ".pdf")
  pdf(FigName, width = 12, height = 8.2)  
  par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.5, 2, .1))
  par(mfrow = c(length(VectProbaDeathThreshold),3))
  
  print(paste0("CoefDelay=",CoefDelay))
  print(paste0("PropIVM=",PropIVM))
  
  GrNumber=0
  
  for (ProbaDeathThreshold in VectProbaDeathThreshold) {
    
    print(paste0("ProbaDeathThreshold=",ProbaDeathThreshold))
    
    
    
    #On tourne le modèle pour Bohemia
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
    ModelBohemia=ModelOutput
    
    
    #On tourne le modèle pour Rimdamal
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
    #VectTime_between_cycles=0.7*c(39.50, 79.25-39.50, 79.25-(79.25-39.50))
    ModelRimdamal=ModelOutput
    
    
    #On tourne le modèle pour LongLasting 06
    strategy=0;NbCycle=3
    ModelLongLasting06=list()
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,
                           Gap,dt,time,Ntime,
                           d_eta,eta_max,Veta,Neta,
                           dtau,tau_max,Vtau,Ntau,
                           strategy,ProbaDeathThreshold,IVM_age_dependancy)
    ModelLongLasting06[[Number_of_cycle]]=ModelOutput
    
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
      ModelLongLasting06[[Number_of_cycle]]=ModelOutput
    }
    #on récupère les valeurs ici
    print("LongLasting0.6")
    print(VectTime_between_cycles)
    #VectTime_between_cyclesLongLasting06=0.7*c(65, 64)
    
    
    #On tourne le modèle pour LongLasting 1
    strategy=1;NbCycle=3
    ModelLongLasting1=list()
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,
                           Gap,dt,time,Ntime,
                           d_eta,eta_max,Veta,Neta,
                           dtau,tau_max,Vtau,Ntau,
                           strategy,ProbaDeathThreshold,IVM_age_dependancy)
    ModelLongLasting1[[Number_of_cycle]]=ModelOutput
    
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
      ModelLongLasting1[[Number_of_cycle]]=ModelOutput
    }
    #on récupère les valeurs ici
    print("LoongLasting1")
    print(VectTime_between_cycles)
    
    PrevIh =ModelOutput$PropIhTot[t_begin_Camp/dt]
    GainPrev = rep(0.8 * PrevIh, Ntime)
    
    for (nc in 1:3) {
      GrNumber=GrNumber+1
      plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_Camp, tmax), ylab = "",
           ylim = c(0.05, 0.32), cex.lab = 1.2, yaxt = "n", xaxt = "n")
      axis(1, at = seq(t_begin_Camp, tmax, by = 50), labels = seq(0, tmax-t_begin_Camp, by = 50), 
           las = 1)
      y_vals = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
      axis(2, at = y_vals, labels = y_vals * 100, las = 1)
      par(xpd = NA)
      text(t_begin_Camp-0.05, 0.32 * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.6)
      par(xpd = FALSE)
      lines(time, ModelLongLasting06[[nc]]$PropIhTot, lwd = LC, cex.lab = 1.2, 
            lty = LineVect[1], col = ColVect[1])
      lines(time, ModelLongLasting1[[nc]]$PropIhTot, lwd = LC, cex.lab = 1.2, 
            lty = LineVect[2], col = ColVect[2])
      lines(time, ModelBohemia$PropIhTot, lwd = LC, cex.lab = 1.2, 
            lty = LineVect[3], col = ColVect[3])
      lines(time, ModelRimdamal$PropIhTot, lwd = LC, cex.lab = 1.2, 
            lty = LineVect[4], col = ColVect[4])
      lines(time, GainPrev, lwd = .5, cex.lab = 1.2, lty = LineVect[2])
      
      if(is.element(GrNumber,c(1,2,3))){
        mtext(bquote("#cycle for LAIF31: " ~ n[c] == .(nc)),
              side = 3, adj = 0.5, cex = 0.8, line = 1, font = 1)  
      }
      
      if(is.element(GrNumber,c(1,4,7))){
        mtext("Clinical cases (%)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8) 
        mtext(bquote(Eff == .(ProbaDeathThreshold * 100) * "%"), 
              side = 2, adj = 0.5, cex = 0.9, line = 3, font = 1)
      }
      
      if(is.element(GrNumber,c(7,8,9))){
        mtext("Time (Days)", side = 1, adj = 0.5, cex = .9, line = 2, font = .8) 
      }
    }
    
  }
  
  
  #plot legend
  LEGEND=c("LAIF31−0.6", "LAIF31−1.0", "BOHEMIA", "RIMDAMAL II")
  reset()
  legend("bottom", legend=LEGEND, 
         xpd = NA, horiz = TRUE, #inset = c(-3,-0.6), 
         bty = "n", lty = LineVect, lw=LC,
         col = ColVect, cex = .8)
  
  dev.off()
}

sink()
}

#Figure: 6 & 7  Scenarios_Formulations Barplots formulations 
###########################
CoefDelay=0.7
IVM_age_dependancy=0
sink(paste("DurEntreCycleBarplotsIVM_age_dependancy",IVM_age_dependancy,".txt",sep=""))

data_Ih_a_Bohemia = data.frame()
data_Ih_a_Rimdamal = data.frame()
data_Ih_a_LAIF0_6 = data.frame()
data_Ih_a_LAIF1_0 = data.frame()

VectProbaDeathThreshold=c(0.5)
VectDataset=c("No IVM", "PropIVM 0.7")
LegendBarplot=c("Without IVM", "With IVM")
ColBarplot=c("No IVM" = "#148F77", "PropIVM 0.7" = "#F39C12")
VectPropIVM=c(0, 0.7)

for (ProbaDeathThreshold in VectProbaDeathThreshold) {
  
  Nbb=0
  for (PropIVM in VectPropIVM) {
    Nbb= Nbb+1
    
    print(paste0("CoefDelay=",CoefDelay))
    print(paste0("PropIVM=",PropIVM))
    print(paste0("ProbaDeathThreshold=",ProbaDeathThreshold))
    
    #On tourne le modèle pour Bohemia
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
    ModelBohemia=ModelOutput
    
    data_Ih_a_Bohemia = rbind(data_Ih_a_Bohemia, 
                              data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                         mIh_a = ModelBohemia$mIh_a,  fIh_a = ModelBohemia$fIh_a, 
                                         Dataset = VectDataset[Nbb])) 
    
    
    
    #On tourne le modèle pour Rimdamal
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
    #VectTime_between_cycles=0.7*c(39.50, 79.25-39.50, 79.25-(79.25-39.50))
    ModelRimdamal=ModelOutput
    
    data_Ih_a_Rimdamal = rbind(data_Ih_a_Rimdamal, 
                               data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                          mIh_a = ModelRimdamal$mIh_a,  fIh_a = ModelRimdamal$fIh_a, 
                                          Dataset = VectDataset[Nbb]))
    
    
    #On tourne le modèle pour LongLasting 06
    strategy=0;NbCycle=3
    ModelLongLasting06=list()
    Plot_mIh_LAIF06=list()
    Plot_fIh_LAIF06=list()
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,
                           Gap,dt,time,Ntime,
                           d_eta,eta_max,Veta,Neta,
                           dtau,tau_max,Vtau,Ntau,
                           strategy,ProbaDeathThreshold,IVM_age_dependancy)
    ModelLongLasting06[[Number_of_cycle]]=ModelOutput
    
    data_Ih_a_LAIF0_6 = rbind(data_Ih_a_LAIF0_6, 
                              data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                         mIh_a = ModelLongLasting06[[Number_of_cycle]]$mIh_a,
                                         fIh_a = ModelLongLasting06[[Number_of_cycle]]$fIh_a, 
                                         Dataset = VectDataset[Nbb], Ncycle=Number_of_cycle))
    
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
      ModelLongLasting06[[Number_of_cycle]]=ModelOutput
      
      data_Ih_a_LAIF0_6 = rbind(data_Ih_a_LAIF0_6, 
                                data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                           mIh_a = ModelLongLasting06[[Number_of_cycle]]$mIh_a,
                                           fIh_a = ModelLongLasting06[[Number_of_cycle]]$fIh_a, 
                                           Dataset = VectDataset[Nbb], Ncycle=Number_of_cycle))
      
    }
    #on récupère les valeurs ici
    print("LongLasting0.6")
    print(VectTime_between_cycles)
    
    #On tourne le modèle pour LongLasting 1.0
    strategy=1;NbCycle=3
    ModelLongLasting1=list()
    Plot_mIh_LAIF1_0=list()
    Plot_fIh_LAIF1_0=list()
    Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,
                           Gap,dt,time,Ntime,
                           d_eta,eta_max,Veta,Neta,
                           dtau,tau_max,Vtau,Ntau,
                           strategy,ProbaDeathThreshold,IVM_age_dependancy)
    ModelLongLasting1[[Number_of_cycle]]=ModelOutput
    
    data_Ih_a_LAIF1_0 = rbind(data_Ih_a_LAIF1_0, 
                              data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                         mIh_a = ModelLongLasting1[[Number_of_cycle]]$mIh_a,
                                         fIh_a = ModelLongLasting1[[Number_of_cycle]]$fIh_a, 
                                         Dataset = VectDataset[Nbb], Ncycle=Number_of_cycle))
    
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
      ModelLongLasting1[[Number_of_cycle]]=ModelOutput
      
      data_Ih_a_LAIF1_0 = rbind(data_Ih_a_LAIF1_0, 
                                data.frame(Age = ModelOutput$cut_age[-length(ModelOutput$cut_age)], # Bornes inférieures
                                           mIh_a = ModelLongLasting1[[Number_of_cycle]]$mIh_a,
                                           fIh_a = ModelLongLasting1[[Number_of_cycle]]$fIh_a, 
                                           Dataset = VectDataset[Nbb], Ncycle=Number_of_cycle))
      
    }
    #on récupère les valeurs ici
    print("LongLasting1")
    print(VectTime_between_cycles)
    
  }
  
  GrNumber=0
  #plot mIh_a
  GrNumber=GrNumber+1
  #Barplots Bohemia
  {
    maxPop=28000
    PlotmIh_a_Bohemia= ggplot(data = data_Ih_a_Bohemia, aes(x = Age + 2.5, y = mIh_a, fill = Dataset)) +
      geom_bar(
        stat = "identity", 
        color = "black", 
        position = position_dodge(width = 4), # Décalage ajusté
        width = 4                             # Largeur des barres ajustée
      ) +
      annotate(
        "rect", xmin = 5, xmax = 90, ymin = 0, ymax = Inf, fill =   "#76448A",  alpha = 0.15) +
      theme_light() +
      scale_x_continuous(breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
                         limits = c(0, 90))+
      scale_y_continuous(breaks = seq(0, maxPop, by = 5000), 
                         labels = c("0", expression(5%*%10^3), expression(10^4),
                                    expression(1.5%*%10^4),expression(2%*%10^4),expression(2.5%*%10^4)),
                         limits = c(0,maxPop)) +
      labs(
        x = "Age (years)",
        y = "Male clinical cases",
        title = "",
        fill = ""
      ) + 
      ggtitle("BOHEMIA") + 
      labs(subtitle = paste0("(", LETTERS[GrNumber], ")")) +
      scale_fill_manual(
        values = ColBarplot,
        labels = LegendBarplot 
      ) +
      theme(
        legend.position = "none", 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(0.3, "cm"), #dimension des carrés en legende
        plot.title = element_text(hjust = 0.5, size = 12), 
        plot.subtitle =element_text(hjust = 0, size = 14) , 
        panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
      ) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
  }
  GrNumber=GrNumber+1
  #Barplots Rimdamal
  {
    maxPop=28000
    PlotmIh_a_Rimdamal=ggplot(data = data_Ih_a_Rimdamal, aes(x = Age + 2.5, y = mIh_a, fill = Dataset)) +
      geom_bar(
        stat = "identity", 
        color = "black", 
        position = position_dodge(width = 4), # Décalage ajusté
        width = 4                             # Largeur des barres ajustée
      ) +
      annotate(
        "rect", xmin = 5, xmax = 90, ymin = 0, ymax = Inf, fill =   "#76448A",  alpha = 0.15) +
      theme_light() +
      scale_x_continuous(breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
                         limits = c(0, 90))+
      scale_y_continuous(breaks = seq(0, maxPop, by = 5000), 
                         labels = c("0", expression(5%*%10^3), expression(10^4),
                                    expression(1.5%*%10^4),expression(2%*%10^4),expression(2.5%*%10^4)),
                         limits = c(0,maxPop)) +
      labs(
        x = "Age (years)",
        y = "Male clinical cases",
        title = "",
        fill = ""
      ) +
      ggtitle("RIMDAMAL II") + 
      labs(subtitle = paste0("(", LETTERS[GrNumber], ")")) +
      scale_fill_manual(
        values = ColBarplot,
        labels = LegendBarplot 
      ) +
      theme(
        legend.position = "none", 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(0.3, "cm"), #dimension des carrés en legende
        plot.title = element_text(hjust = 0.5, size = 12), 
        plot.subtitle =element_text(hjust = 0, size = 14) , 
        panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
      ) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
    
  }
  
  GrNumber=GrNumber+1
  #plot fIh_a
  #Barplots Bohemia
  {
    maxPop=28000
    PlotfIh_a_Bohemia= ggplot(data = data_Ih_a_Bohemia, aes(x = Age + 2.5, y = fIh_a, fill = Dataset)) +
      geom_bar(
        stat = "identity", 
        color = "black", 
        position = position_dodge(width = 4), # Décalage ajusté
        width = 4                             # Largeur des barres ajustée
      ) +
      annotate(
        "rect", xmin = 5, xmax = 15, ymin = 0, ymax = Inf, fill =  "#76448A",  alpha = 0.15) +
      annotate(
        "rect", xmin = 45, xmax = 90, ymin = 0, ymax = Inf, fill = "#76448A",  alpha = 0.15) +
      theme_light() +
      scale_x_continuous(breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
                         limits = c(0, 90))+
      scale_y_continuous(breaks = seq(0, maxPop, by = 5000), 
                         labels = c("0", expression(5%*%10^3), expression(10^4),
                                    expression(1.5%*%10^4),expression(2%*%10^4),expression(2.5%*%10^4)),
                         limits = c(0,maxPop)) +
      labs(
        x = "Age (years)",
        y = "Female clinical cases",
        title = "",
        fill = ""
      ) +
      labs(subtitle =paste0("(", LETTERS[GrNumber], ")") ) +
      scale_fill_manual(
        values = ColBarplot,
        labels = LegendBarplot 
      ) +
      theme(
        legend.position = "none", 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(0.3, "cm"), #dimension des carrés en legende
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0, size = 14), 
        panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
      ) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
  }
  GrNumber=GrNumber+1
  #Barplots Rimdamal
  {
    maxPop=28000
    PlotfIh_a_Rimdamal= ggplot(data = data_Ih_a_Rimdamal, aes(x = Age + 2.5, y = fIh_a, fill = Dataset)) +
      geom_bar(
        stat = "identity", 
        color = "black", 
        position = position_dodge(width = 4), # Décalage ajusté
        width = 4                             # Largeur des barres ajustée
      ) +
      annotate(
        "rect", xmin = 5, xmax = 15, ymin = 0, ymax = Inf, fill =  "#76448A",  alpha = 0.15) +
      annotate(
        "rect", xmin = 45, xmax = 90, ymin = 0, ymax = Inf, fill = "#76448A",  alpha = 0.15) +
      theme_light() +
      scale_x_continuous(breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
                         limits = c(0, 90))+
      scale_y_continuous(breaks = seq(0, maxPop, by = 5000), 
                         labels = c("0", expression(5%*%10^3), expression(10^4),
                                    expression(1.5%*%10^4),expression(2%*%10^4),expression(2.5%*%10^4)),
                         limits = c(0,maxPop)) +
      labs(
        x = "Age (years)",
        y = "Female clinical cases",
        title = "",
        fill = ""
      ) +
      labs(subtitle =paste0("(", LETTERS[GrNumber], ")") ) +
      scale_fill_manual(
        values = ColBarplot,
        labels = LegendBarplot 
      ) +
      theme(
        legend.position = "none", 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(0.3, "cm"), #dimension des carrés en legende
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0, size = 14), 
        panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
      ) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
  }
  
  
  plot_list_oral=list(PlotmIh_a_Bohemia,PlotmIh_a_Rimdamal,
                      PlotfIh_a_Bohemia,PlotfIh_a_Rimdamal)
  
  
  combined_plot = wrap_plots(plot_list_oral, ncol = 2, nrow = 2) + 
    plot_layout(guides = "collect") + 
    theme(
      legend.position = "bottom",  
      legend.box = "horizontal", 
      legend.key.size = unit(0.3, "cm"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) + 
    plot_annotation(
      theme = theme(legend.position = "bottom")  
    )
  
  ggsave(
    filename = paste0("Barplot_Ih_Oral_Eff_", gsub("\\.", "_", ProbaDeathThreshold), ".pdf"),
    plot = combined_plot, width = 8, height = 6.5, units = "in", dpi = 300 )

  
  # Barplot LAIF0.6
  GrNumber = 0
  for (Number_of_cycle in 1:3) {
    GrNumber = GrNumber + 1
    data_Ih_a_LAIF0_6_filt = data_Ih_a_LAIF0_6 %>% filter(Ncycle == Number_of_cycle)
    
    # plot mIh_a
    {
      maxPop = 28000
      Plot_mIh_LAIF06[[Number_of_cycle]] = ggplot(data = data_Ih_a_LAIF0_6_filt, aes(x = Age + 2.5, y = mIh_a, fill = Dataset)) +
        geom_bar(
          stat = "identity",
          color = "black",
          position = position_dodge(width = 4),
          width = 4
        ) +
        annotate("rect", xmin = 5, xmax = 90, ymin = 0, ymax = Inf, fill = "#76448A", alpha = 0.15) +
        theme_light() +
        scale_x_continuous(
          breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
          limits = c(0, 90)
        ) +
        scale_y_continuous(
          breaks = seq(0, maxPop, by = 5000),
          labels = c("0", expression(5 %*% 10^3), expression(10^4),
                     expression(1.5 %*% 10^4), expression(2 %*% 10^4), expression(2.5 %*% 10^4)),
          limits = c(0, maxPop)
        ) +
        labs(
          x = "Age (years)",
          y = "Male clinical cases",
          title = "",
          fill = ""
        ) +
        ggtitle(bquote("#cycle for LAIF31-0.6: " ~ n[c] == .(Number_of_cycle))) +
        labs(subtitle = paste0("(", LETTERS[GrNumber], ")")) +
        scale_fill_manual(
          values = ColBarplot,
          labels = LegendBarplot
        ) +
        theme(
          legend.position = "none",  
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.3, "cm"),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0, size = 14),
          panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
    }
  }
  
  for (Number_of_cycle in 1:3) {
    GrNumber = GrNumber + 1
    data_Ih_a_LAIF0_6_filt = data_Ih_a_LAIF0_6 %>% filter(Ncycle == Number_of_cycle)
    
    # plot fIh_a
    {maxPop = 28000
      Plot_fIh_LAIF06[[Number_of_cycle]] = ggplot(data = data_Ih_a_LAIF0_6_filt, aes(x = Age + 2.5, y = fIh_a, fill = Dataset)) +
        geom_bar(
          stat = "identity",
          color = "black",
          position = position_dodge(width = 4),
          width = 4
        ) +
        annotate("rect", xmin = 5, xmax = 15, ymin = 0, ymax = Inf, fill = "#76448A", alpha = 0.15) +
        annotate("rect", xmin = 45, xmax = 90, ymin = 0, ymax = Inf, fill = "#76448A", alpha = 0.15) +
        theme_light() +
        scale_x_continuous(
          breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
          limits = c(0, 90)
        ) +
        scale_y_continuous(
          breaks = seq(0, maxPop, by = 5000),
          labels = c("0", expression(5 %*% 10^3), expression(10^4),
                     expression(1.5 %*% 10^4), expression(2 %*% 10^4), expression(2.5 %*% 10^4)),
          limits = c(0, maxPop)
        ) +
        labs(
          x = "Age (years)",
          y = "Female clinical cases",
          title = "",
          fill = ""
        ) +
        labs(subtitle = paste0("(", LETTERS[GrNumber], ")")) +
        scale_fill_manual(
          values = ColBarplot,
          labels = LegendBarplot
        ) +
        theme(
          legend.position = "none",  
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.3, "cm"),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0, size = 14),
          panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
    }
  }
  
  plot_list_LAIF06 = list(Plot_mIh_LAIF06[[1]], Plot_mIh_LAIF06[[2]], Plot_mIh_LAIF06[[3]],
                          Plot_fIh_LAIF06[[1]], Plot_fIh_LAIF06[[2]], Plot_fIh_LAIF06[[3]])
  
  combined_plot_LAIF06 = wrap_plots(plot_list_LAIF06, ncol = 3, nrow = 2) +
    plot_layout(guides = "collect") +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.size = unit(0.3, "cm"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    plot_annotation(
      theme = theme(legend.position = "bottom"))
  
  ggsave(
    filename = paste0("Barplot_Ih_LAIF06_Eff_", gsub("\\.", "_", ProbaDeathThreshold), ".pdf"),
    plot = combined_plot_LAIF06, width = 12, height = 6.5, units = "in", dpi = 300)
  
  
  # Barplot LAIF1_0
  GrNumber = 0
  for (Number_of_cycle in 1:3) {
    GrNumber = GrNumber + 1
    data_Ih_a_LAIF1_0_filt = data_Ih_a_LAIF1_0 %>% filter(Ncycle == Number_of_cycle)
    
    # plot mIh_a
    {
      maxPop = 28000
      Plot_mIh_LAIF1_0[[Number_of_cycle]] = ggplot(data = data_Ih_a_LAIF1_0_filt, aes(x = Age + 2.5, y = mIh_a, fill = Dataset)) +
        geom_bar(
          stat = "identity",
          color = "black",
          position = position_dodge(width = 4),
          width = 4
        ) +
        annotate("rect", xmin = 5, xmax = 90, ymin = 0, ymax = Inf, fill = "#76448A", alpha = 0.15) +
        theme_light() +
        scale_x_continuous(
          breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
          limits = c(0, 90)
        ) +
        scale_y_continuous(
          breaks = seq(0, maxPop, by = 5000),
          labels = c("0", expression(5 %*% 10^3), expression(10^4),
                     expression(1.5 %*% 10^4), expression(2 %*% 10^4), expression(2.5 %*% 10^4)),
          limits = c(0, maxPop)
        ) +
        labs(
          x = "Age (years)",
          y = "Male clinical cases",
          title = "",
          fill = ""
        ) +
        ggtitle(bquote("#cycle for LAIF31-1.0: " ~ n[c] == .(Number_of_cycle))) +
        labs(subtitle = paste0("(", LETTERS[GrNumber], ")")) +
        scale_fill_manual(
          values = ColBarplot,
          labels = LegendBarplot
        ) +
        theme(
          legend.position = "none",  
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.3, "cm"),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0, size = 14),
          panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
    }
  }
  
  for (Number_of_cycle in 1:3) {
    GrNumber = GrNumber + 1
    data_Ih_a_LAIF1_0_filt = data_Ih_a_LAIF1_0 %>% filter(Ncycle == Number_of_cycle)
    
    # plot fIh_a
    {maxPop = 28000
      Plot_fIh_LAIF1_0[[Number_of_cycle]] = ggplot(data = data_Ih_a_LAIF1_0_filt, aes(x = Age + 2.5, y = fIh_a, fill = Dataset)) +
        geom_bar(
          stat = "identity",
          color = "black",
          position = position_dodge(width = 4),
          width = 4
        ) +
        annotate("rect", xmin = 5, xmax = 15, ymin = 0, ymax = Inf, fill = "#76448A", alpha = 0.15) +
        annotate("rect", xmin = 45, xmax = 90, ymin = 0, ymax = Inf, fill = "#76448A", alpha = 0.15) +
        theme_light() +
        scale_x_continuous(
          breaks = seq(0, ah_max, by = 5), minor_breaks = NULL, expand = c(0.02, 0),
          limits = c(0, 90)
        ) +
        scale_y_continuous(
          breaks = seq(0, maxPop, by = 5000),
          labels = c("0", expression(5 %*% 10^3), expression(10^4),
                     expression(1.5 %*% 10^4), expression(2 %*% 10^4), expression(2.5 %*% 10^4)),
          limits = c(0, maxPop)
        ) +
        labs(
          x = "Age (years)",
          y = "Female clinical cases",
          title = "",
          fill = ""
        ) +
        labs(subtitle = paste0("(", LETTERS[GrNumber], ")")) +
        scale_fill_manual(
          values = ColBarplot,
          labels = LegendBarplot
        ) +
        theme(
          legend.position = "none",  
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.3, "cm"),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0, size = 14),
          panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
        ) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
    }
  }
  
  plot_list_LAIF1_0 = list(Plot_mIh_LAIF1_0[[1]], Plot_mIh_LAIF1_0[[2]], Plot_mIh_LAIF1_0[[3]],
                           Plot_fIh_LAIF1_0[[1]], Plot_fIh_LAIF1_0[[2]], Plot_fIh_LAIF1_0[[3]])
  
  combined_plot_LAIF1_0 = wrap_plots(plot_list_LAIF1_0, ncol = 3, nrow = 2) +
    plot_layout(guides = "collect") +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.size = unit(0.3, "cm"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    plot_annotation(
      theme = theme(legend.position = "bottom"))
  
  ggsave(
    filename = paste0("Barplot_Ih_LAIF1_0_Eff_", gsub("\\.", "_", ProbaDeathThreshold), ".pdf"),
    plot = combined_plot_LAIF1_0, width = 12, height = 6.5, units = "in", dpi = 300)
  
}

dev.off()
#dev.new()
sink()


#Figure 9 Scenarios_GainPrev Formulations
#######################
CoefDelay=0.7
Nbrun=0
#sink(paste("Values_GainPrev",".txt",sep=""))

VectPropIVM = c(seq(0.3, 0.95, by = 0.1), 0.95)
VectIVM_age_dependancy=c(0,1)
VectProbaDeathThreshold=c(0.75,0.6,0.5)

GainPrev_array_oral=array(0, dim = c(2,length(VectPropIVM), length(VectProbaDeathThreshold),
                                     length(VectIVM_age_dependancy)))
#GainPrev_matrix_oral[1,,,] pour BOHEMIA, GainPrev_matrix_oral[2,,,] pour RIMDAMAL
Time_max_array_oral=array(0, dim = c(2,length(VectPropIVM),length(VectProbaDeathThreshold),
                                     length(VectIVM_age_dependancy)))

VectNumber_of_cycle=c(1,2,3)
GainPrev_array_LAIF=array(0, dim = c(2,length(VectPropIVM),length(VectNumber_of_cycle),
                                     length(VectProbaDeathThreshold),length(VectIVM_age_dependancy)))
# GainPrev_array_LAIF[1,,,,] pour "LAIF31_0.6", GainPrev_array_LAIF[2,,,,] pour "LAIF31_1.0"
Time_max_array_LAIF=array(0, dim = c(2,length(VectPropIVM),length(VectNumber_of_cycle),
                                     length(VectProbaDeathThreshold),length(VectIVM_age_dependancy)))

for (p in 1:length(VectIVM_age_dependancy)) {
  IVM_age_dependancy=VectIVM_age_dependancy[p]
  print(paste0("IVM_age_dependancy=",IVM_age_dependancy))
  for (j in 1:length(VectPropIVM)) {
    PropIVM = VectPropIVM[j]
    print(paste0("PropIVM=",PropIVM))
    for (k in 1:length(VectProbaDeathThreshold)) {
      ProbaDeathThreshold=VectProbaDeathThreshold[k]
      print(paste0("ProbaDeathThreshold=",ProbaDeathThreshold))
      
      #On tourne le modèle pour Bohemia
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
      ModelBohemia=ModelOutput
      
      GainPrev_array_oral[1,j,k,p]=1-(ModelBohemia$MinPropIhTot/ModelBohemia$PropIhTot0)
      Time_max_array_oral[1,j,k,p]=ModelBohemia$DurIvmEffect
      
      
      #On tourne le modèle pour Rimdamal
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
      #VectTime_between_cycles=0.7*c(39.50, 79.25-39.50, 79.25-(79.25-39.50))
      ModelRimdamal=ModelOutput
      
      GainPrev_array_oral[2,j,k,p]=1-(ModelRimdamal$MinPropIhTot/ModelRimdamal$PropIhTot0)
      Time_max_array_oral[2,j,k,p]=ModelRimdamal$DurIvmEffect
      
      
      #On tourne le modèle pour LongLasting 06
      strategy=0;NbCycle=3
      ModelLongLasting06=list()
      Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,
                             Gap,dt,time,Ntime,
                             d_eta,eta_max,Veta,Neta,
                             dtau,tau_max,Vtau,Ntau,
                             strategy,ProbaDeathThreshold,IVM_age_dependancy)
      ModelLongLasting06[[Number_of_cycle]]=ModelOutput
      
      GainPrev_array_LAIF[1,j,Number_of_cycle,k,p]=1-(ModelLongLasting06[[Number_of_cycle]]$MinPropIhTot/
                                                        ModelLongLasting06[[Number_of_cycle]]$PropIhTot0)
      Time_max_array_LAIF[1,j,Number_of_cycle,k,p]=ModelLongLasting06[[Number_of_cycle]]$DurIvmEffect
      
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
        ModelLongLasting06[[Number_of_cycle]]=ModelOutput
        
        GainPrev_array_LAIF[1,j,Number_of_cycle,k,p]=1-(ModelLongLasting06[[Number_of_cycle]]$MinPropIhTot/
                                                          ModelLongLasting06[[Number_of_cycle]]$PropIhTot0)
        Time_max_array_LAIF[1,j,Number_of_cycle,k,p]=ModelLongLasting06[[Number_of_cycle]]$DurIvmEffect
      }
      #on récupère les valeurs ici
      print("LongLasting0.6")
      print(VectTime_between_cycles)
      
      #On tourne le modèle pour LongLasting 1
      strategy=1;NbCycle=3
      ModelLongLasting1=list()
      Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,
                             Gap,dt,time,Ntime,
                             d_eta,eta_max,Veta,Neta,
                             dtau,tau_max,Vtau,Ntau,
                             strategy,ProbaDeathThreshold,IVM_age_dependancy)
      ModelLongLasting1[[Number_of_cycle]]=ModelOutput
      
      GainPrev_array_LAIF[2,j,Number_of_cycle,k,p]=1-(ModelLongLasting1[[Number_of_cycle]]$MinPropIhTot/
                                                        ModelLongLasting1[[Number_of_cycle]]$PropIhTot0)
      Time_max_array_LAIF[2,j,Number_of_cycle,k,p]=ModelLongLasting1[[Number_of_cycle]]$DurIvmEffect
      
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
        ModelLongLasting1[[Number_of_cycle]]=ModelOutput
        
        GainPrev_array_LAIF[2,j,Number_of_cycle,k,p]=1-(ModelLongLasting1[[Number_of_cycle]]$MinPropIhTot/
                                                          ModelLongLasting1[[Number_of_cycle]]$PropIhTot0)
        Time_max_array_LAIF[2,j,Number_of_cycle,k,p]=ModelLongLasting1[[Number_of_cycle]]$DurIvmEffect
      }
      #on récupère les valeurs ici
      print("LoongLasting1")
      print(VectTime_between_cycles)
      
      Nbrun= Nbrun + 1
      print(Nbrun)
    }
  }
}

#plot Gain Prev 
{
  VectNumber_of_cycle = c(1,2,3)
  
  ColVect = c("#ff3355", "#037153", "#D98E04", "#205072")
  LineVect <- c(1, 3, 2, 4)
  LC <- 2
  GrNumber = 0
  
  FigName=paste0("Scenarios_Gain_Prev",".pdf")
  pdf(FigName, width = 12, height = 8.2)
  
  par(oma = c(2, 1, 1, 1), mar = c(3, 3.5, 2, 1))
  par(mfrow = c(length(VectProbaDeathThreshold),length(VectNumber_of_cycle)))
  
  for (k in 1:length(VectProbaDeathThreshold)) {
    for (i in 1:length(VectNumber_of_cycle)) {
      GrNumber = GrNumber + 1
      plot(-1, 1, type = "l", xlab = "", xlim = c(min(VectPropIVM), 1), ylab = "",
           ylim = c(0, 0.72), cex.lab = 1, yaxt = "n", xaxt = "n")
      y_vals = c(0, 0.1, 0.20, 0.3, 0.40, 0.5, 0.6, 0.7)
      axis(2, at = y_vals, labels = y_vals * 100, las = 1)
      x_vals = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8,  0.9, 1)
      axis(1, at = x_vals, labels = x_vals, las = 1, cex.axis = 1)
      
      par(xpd = NA)
      text(x = min(VectPropIVM)-0.05 , y = 0.72 * 1.15,
           labels = paste0("(", LETTERS[GrNumber], ")"),
           cex = 1.8, font = 1, adj = 0)
      par(xpd = FALSE)
      
      GainPrev = rep(0.2, length(VectPropIVM))
      lines(VectPropIVM, GainPrev, lwd = .7, cex.lab = 1, lty = LineVect[3])
      
      if (is.element(GrNumber,c(1,2,3))) {
        mtext(bquote("#cycle for LAIF31: " ~ n[c] == .(VectNumber_of_cycle[i])),
              side = 3, adj = 0.5, cex = 0.8, line = 1, font = 1)  
      }
      
      if (is.element(GrNumber,c(7,8,9))) {
        mtext("Rate of IVM exposure", side = 1, adj = 0.5, cex = .8, line = 2, font = .8) 
      }
      
      if (i==1) {
        mtext("Relative gain prevalence (%)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8) 
        mtext(bquote(Eff == .(VectProbaDeathThreshold[k] * 100) * "%"), 
              side = 2, adj = 0.5, cex = 0.9, line = 3, font = 1)
      }
      
      for (p in 1:length(VectIVM_age_dependancy)) {
        lines(VectPropIVM, GainPrev_array_LAIF[1,,i,k,p], col = ColVect[1], lwd = LC, lty = LineVect[p])
        lines(VectPropIVM, GainPrev_array_LAIF[2,,i,k,p], col = ColVect[2], lwd = LC, lty = LineVect[p])
        
        lines(VectPropIVM, GainPrev_array_oral[1,,k,p], col = ColVect[3], lwd = LC, lty = LineVect[p])
        lines(VectPropIVM, GainPrev_array_oral[2,,k,p], col = ColVect[4], lwd = LC, lty = LineVect[p])
      }
      
      if (is.element(GrNumber,c(1,2,3))) {
        legend("topleft",   
               legend = c("BOHEMIA","BOHEMIA with AD","RIMDAMAL II","RIMDAMAL II with AD",
                          bquote("LAIF31-0.6, " ~ n[c] == .(VectNumber_of_cycle[i])),
                          bquote("LAIF31-0.6 with AD, " ~ n[c] == .(VectNumber_of_cycle[i])),
                          bquote("LAIF31-1.0, " ~ n[c] == .(VectNumber_of_cycle[i])),
                          bquote("LAIF31-1.0 with AD, " ~ n[c] == .(VectNumber_of_cycle[i]))),
               horiz = FALSE, lwd = 2, bty = "n",    col = c("#D98E04","#D98E04", "#205072", "#205072",
                                                             "#ff3355","#ff3355", "#037153", "#037153"), 
               lty =c(1,3,1,3,1,3,1,3), cex = 0.65)        
      }                
      
    }
  }
  
  dev.off()  
  getwd()
}

# sink(file = NULL, type = "output",append = FALSE)



############################### FIGURE 2 & 3, Dynamics of IVM efficacy loss post-exposure
FigName = "PropHuman_in_IVM_withoutAge.pdf"
pdf(FigName, width = 13, height = 3.2)
{
  
  IVM_age_dependancy=0; 
  VectProbaDeathThreshold=c(0.75,0.6,0.5)
  VectStrategy=c(0,1,2,3)
  Title=c("LAIF31-0.6","LAIF31-1.0","BOHEMIA","RIMDAMAL II")
  XlimVal=c(tau_max,tau_max,30,30)
  XstepLabel=c(30,30,5,5)
  
  par(oma = c(1.2, 2, 1, 1), mar = c(3, 3, 2, 1))
  par(mfrow = c(1, 4)) 
  
  ColVect = c("#037153", "#ff3355", "#5f05f3", "#D98E04", "#205072")
  ltyVect = c(1, 2, 3, 4)
  LC = 2  
  
  cp=0
  for (strategy in VectStrategy) {
    cp=cp+1
    
    
    plot(-1, 1, type = "l", xlab = "", xlim = c(0, XlimVal[cp]), ylab = "", ylim = c(0, 1.05), 
         cex.lab = 0.8, xaxt = "n", yaxt = "n")
    axis(2, at = seq(0, 1.05, by = 0.2), labels = seq(0, 1.1, by = 0.2), cex.axis = 0.65, las = 1)
    axis(1, at = seq(0, XlimVal[cp], by = XstepLabel[cp]), labels = seq(0, XlimVal[cp], by = XstepLabel[cp]), cex.axis = 0.65, las = 1)
    
    ccp=0
    for (ProbaDeathThreshold in VectProbaDeathThreshold) {
      ccp=ccp+1
      ParamRho=RhoFunction(Nah,Ntau,Vtau,Vah,IVM_age_dependancy,
                           strategy,ProbaDeathThreshold)
      lines(Vtau, ParamRho$PropHuman_in_IVM_Group1, col = ColVect[ccp], lty = ltyVect[ccp], lwd = LC)
    }
    
    # if(cp==3){
    #   ParamRho=RhoFunction(Nah,Ntau,Vtau,Vah,IVM_age_dependancy,
    #                        strategy,ProbaDeathThreshold=0.5)
    #   lines(Vtau, ParamRho$PropHuman_in_IVM_Group1, col = ColVect[4], lty = ltyVect[1], lwd = LC)
    # } else if(cp==4){
    #   ParamRho=RhoFunction(Nah,Ntau,Vtau,Vah,IVM_age_dependancy,
    #                        strategy,ProbaDeathThreshold=0.5)
    #   lines(Vtau, ParamRho$PropHuman_in_IVM_Group1, col = ColVect[5], lty = ltyVect[1], lwd = LC)
    # } else {
    #   ccp=0
    #   for (ProbaDeathThreshold in VectProbaDeathThreshold) {
    #     ccp=ccp+1
    #     ParamRho=RhoFunction(Nah,Ntau,Vtau,Vah,IVM_age_dependancy,
    #                          strategy,ProbaDeathThreshold)
    #     lines(Vtau, ParamRho$PropHuman_in_IVM_Group1, col = ColVect[ccp], lty = ltyVect[ccp], lwd = LC)
    #   } 
    # }
    
    
    if(cp==1){
      mtext("Probability of remaining exposed to IVM", side = 2, adj = 0.5, cex = 0.8, line = 2, font = 1) 
    }
    mtext("Time post exposure to IVM (days)", side = 1, adj = 0.5, cex = 0.8, line = 2, font = 1)
    mtext(Title[cp], side = 3, adj = 0.5, cex = 1, line = 1,font = 1.2)
    
    legend_position=ifelse(cp==1, "topright", "bottomleft")
    
    if(cp==1){
      legend(legend_position, legend = c(expression(Eff ==  "75" * "%"),
                                         expression(Eff == "60" * "%"),
                                         expression(Eff == "50" * "%")),  #cex = 1.3 et 1
             col = ColVect, lty = ltyVect, cex = 1.3, horiz = FALSE, bty = "n", lwd = LC) 
    }
    
    par(xpd = NA)
    text(-0.05, 1.05 * 1.15, paste("(", LETTERS[cp], ")", sep = ""), cex = 1.6)
    par(xpd = FALSE)
  }
}
dev.off()

FigName = "PropHuman_in_IVM_withtAge.pdf"
pdf(FigName, width = 13, height = 8.2)
{
  IVM_age_dependancy = 1
  VectProbaDeathThreshold = c(0.75, 0.6, 0.5)
  VectStrategy = c(0, 1, 2, 3)
  Title = c("LAIF31-0.6", "LAIF31-1.0", "BOHEMIA", "RIMDAMAL II")
  XlimVal = c(tau_max, tau_max, 30, 30)
  XstepLabel = c(30, 30, 5, 5)
  
  GrpAgeNames = c("5-10 years", "10-18 years", ">18 years")
  
  par(oma = c(1.2, 2, 1, 1), mar = c(3, 3, 2, 1))
  par(mfrow = c(length(GrpAgeNames),length(VectStrategy))) 
  
  ColVect = c("#037153", "#ff3355", "#5f05f3", "#D98E04", "#205072")
  ltyVect = c(1, 2, 3, 4)
  LC = 2
  
  GrNb = 0
  
  for (age_group in 2:4) {
    for (cp in 1:length(VectStrategy)) {
      strategy = VectStrategy[cp]
      
      GrNb = GrNb + 1
      
      plot(-1, 1, type = "l", xlab = "", xlim = c(0, XlimVal[cp]), ylab = "", ylim = c(0, 1.05), 
           cex.lab = 0.8, xaxt = "n", yaxt = "n")
      axis(2, at = seq(0, 1.05, by = 0.2), labels = seq(0, 1.1, by = 0.2), cex.axis = 0.65, las = 1)
      axis(1, at = seq(0, XlimVal[cp], by = XstepLabel[cp]), labels = seq(0, XlimVal[cp], by = XstepLabel[cp]), 
           cex.axis = 0.65, las = 1)
      
      ccp = 0
      for (ProbaDeathThreshold in VectProbaDeathThreshold) {
        ccp = ccp + 1
        ParamRho = RhoFunction(Nah, Ntau, Vtau, Vah, IVM_age_dependancy,
                               strategy, ProbaDeathThreshold)
        
        if (age_group == 2) {
          prop_data = ParamRho$PropHuman_in_IVM_Group2
        } else if (age_group == 3) {
          prop_data = ParamRho$PropHuman_in_IVM_Group3
        } else {
          prop_data = ParamRho$PropHuman_in_IVM_Group4
        }
        
        lines(Vtau, prop_data, col = ColVect[ccp], lty = ltyVect[ccp], lwd = LC)
      }
      
      
      # if (strategy == 2) {  # BOHEMIA
      #   ParamRho = RhoFunction(Nah, Ntau, Vtau, Vah, IVM_age_dependancy,
      #                          strategy, ProbaDeathThreshold = 0.5)
      #   if (age_group == 2) {
      #     prop_data = ParamRho$PropHuman_in_IVM_Group2
      #   } else if (age_group == 3) {
      #     prop_data = ParamRho$PropHuman_in_IVM_Group3
      #   } else {
      #     prop_data = ParamRho$PropHuman_in_IVM_Group4
      #   }
      #   
      #   lines(Vtau, prop_data, col = ColVect[4], lty = ltyVect[1], lwd = LC)
      # } 
      # else if (strategy == 3) {  # RIMDAMAL II
      #   ParamRho = RhoFunction(Nah, Ntau, Vtau, Vah, IVM_age_dependancy,
      #                          strategy, ProbaDeathThreshold = 0.5)
      #   
      #   if (age_group == 2) {
      #     prop_data = ParamRho$PropHuman_in_IVM_Group2
      #   } else if (age_group == 3) {
      #     prop_data = ParamRho$PropHuman_in_IVM_Group3
      #   } else {
      #     prop_data = ParamRho$PropHuman_in_IVM_Group4
      #   }
      #   
      #   lines(Vtau, prop_data, col = ColVect[5], lty = ltyVect[1], lwd = LC)
      # } 
      # else {  # LAIF31
      #   ccp = 0
      #   for (ProbaDeathThreshold in VectProbaDeathThreshold) {
      #     ccp = ccp + 1
      #     ParamRho = RhoFunction(Nah, Ntau, Vtau, Vah, IVM_age_dependancy,
      #                            strategy, ProbaDeathThreshold)
      #     
      #     if (age_group == 2) {
      #       prop_data = ParamRho$PropHuman_in_IVM_Group2
      #     } else if (age_group == 3) {
      #       prop_data = ParamRho$PropHuman_in_IVM_Group3
      #     } else {
      #       prop_data = ParamRho$PropHuman_in_IVM_Group4
      #     }
      #     
      #     lines(Vtau, prop_data, col = ColVect[ccp], lty = ltyVect[ccp], lwd = LC)
      #   }
      # }
      
      if (age_group == 2) {  # Seulement pour la première ligne
        mtext(Title[cp], side = 3, adj = 0.5, cex = 1, line = 1, font = 1.2)
      }
      
      par(xpd = NA)
      text(-0.05, 1.05 * 1.15, paste("(", LETTERS[GrNb], ")", sep = ""), cex = 1.6)
      par(xpd = FALSE)
      
      if (cp == 1) {
        mtext(GrpAgeNames[age_group - 1], side = 2, adj = 0.5, cex = 1, line = 3, font = 1)
      }
      if (cp == 1) {
        mtext("Probability of remaining exposed to IVM", side = 2, adj = 0.5, cex = 0.8, line = 2, font = 1)
      }
      
      if (age_group == 4) {  # Dernière ligne
        mtext("Time post exposure to IVM (days)", side = 1, adj = 0.5, cex = 0.8, line = 2, font = 1)
      }
      
      if (cp == 1 && age_group == 2) { 
        legend("topright", 
               legend = c(expression(Eff == "75" * "%"),
                          expression(Eff == "60" * "%"),
                          expression(Eff == "50" * "%")),
               col = ColVect, lty = ltyVect, cex = 1.3, horiz = FALSE, bty = "n", lwd = LC)
      }
    }
  }
}
dev.off()

######################FIGURE 1 prop death mosquitoes, mum_IVM, RH pour LOngLasting, Bohemia, Rimdal
{
  strategy=0#Long lasting 0.6
  ParamLongLasting06=ModelEeta50Alpha(strategy)
  strategy=1#Long lasting 1
  ParamLongLasting1=ModelEeta50Alpha(strategy)
  strategy=2#Bohemia
  ParamBohemia=ModelEeta50Alpha(strategy)
  strategy=3#Rimdamal
  ParamRimdamal=ModelEeta50Alpha(strategy)
  
  #Figure parameters LOng Lasting (propd death mosquitoes, mum_ivm, HR)
  # on combine 2 formulations LongLasting
  combined_dataLongLasting = rbind(ParamLongLasting06$predprop, ParamLongLasting1$predprop)
  combined_PkLongLasting = rbind(
    ParamLongLasting06$sub_PK %>% mutate(Formulation = "F31-0.6"),
    ParamLongLasting1$sub_PK %>% mutate(Formulation = "Lam-1.0")
  )
  combined_PkLongLasting$DAI[is.na(combined_PkLongLasting$DAI)] = 
    combined_PkLongLasting$DAI2[is.na(combined_PkLongLasting$DAI)]
  
  # on combine 2 formulations orales
  combined_dataOral = rbind(ParamBohemia$predprop, ParamRimdamal$predprop)
  combined_PkOral = rbind(
    ParamBohemia$sub_PK %>% mutate(Formulation = "BOHEMIA"),
    ParamRimdamal$sub_PK %>% mutate(Formulation = "RIMDAMAL"))
  
  #plot Proportion of dead mosquitoes: LongLasting et Orale
  {
    Plot_Prop_dead_mosquitoesLongLasting=ggplot() +
      geom_point(data = combined_PkLongLasting, 
                 aes(x = DAI, y = prop_dead, colour = Formulation, shape = Formulation), size = 1.8) +
      geom_line(data = combined_dataLongLasting, 
                aes(x = DAI, y = Prediction, colour = Formulation), linewidth = 1) +
      geom_ribbon(data = combined_dataLongLasting, 
                  aes(x = DAI, ymin = Lower, ymax = Upper, fill = Formulation), alpha = 0.3) +
      scale_y_continuous(expand = c(0.02, 0)) +
      scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150),labels = c("0", "30", "60", "90","120","150"),
                         expand = c(0.02, 0),minor_breaks = NULL ) +
      scale_color_manual(values = c("F31-0.6" = "#ff3355", "Lam-1.0" = "#037153"),
                         labels = c("LAIF31-0.6", "LAIF31-1.0")) +
      scale_fill_manual(values = c("F31-0.6" = "#ff3355", "Lam-1.0" = "#037153"),
                        labels = c("LAIF31-0.6", "LAIF31-1.0")) +
      scale_shape_manual(values = c("F31-0.6" = 16, "Lam-1.0" = 18),
                         labels = c("LAIF31-0.6", "LAIF31-1.0")) +
      theme_light() +
      labs(
        x = "Time post exposure to IVM (days)",
        y = "Proportion of dead mosquitoes",
        color = "IVM formulation",
        shape = "IVM formulation",
        fill = "IVM formulation"
      )+
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0, size = 16) #, face = "bold"
      ) +
      ggtitle("(A)")
    
    Plot_Prop_dead_mosquitoesOral=ggplot() +
      geom_point(data = combined_PkOral, 
                 aes(x = DAI, y = prop_dead, colour = Formulation, shape = Formulation), size = 1.8) +
      geom_line(data = combined_dataOral, 
                aes(x = DAI, y = Prediction, colour = Formulation), linewidth = 1) +
      geom_ribbon(data = combined_dataOral, 
                  aes(x = DAI, ymin = Lower, ymax = Upper, fill = Formulation), alpha = 0.3) +
      scale_y_continuous(expand = c(0.02, 0),limits = c(0,1)) +
      scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),labels = c("0", "5", "10", "15","20","25", "30"),
                         expand = c(0.02, 0),minor_breaks = NULL ) +
      scale_color_manual(values = c("BOHEMIA" = "#D98E04", "RIMDAMAL" = "#205072"),
                         labels = c("BOHEMIA", "RIMDAMAL II")) +
      scale_fill_manual(values = c("BOHEMIA" = "#D98E04", "RIMDAMAL" = "#205072"),
                        labels = c("BOHEMIA", "RIMDAMAL II")) +
      scale_shape_manual(values = c("BOHEMIA" = 16, "RIMDAMAL" = 18),
                         labels = c("BOHEMIA", "RIMDAMAL II")) +
      theme_light() +
      labs(
        x = "Time post exposure to IVM (days)",
        y = "Proportion of dead mosquitoes",
        color = "IVM formulation",
        shape = "IVM formulation",
        fill = "IVM formulation" )+
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0, size = 16) #, face = "bold"
      ) +
      ggtitle("(D)")
  }
  
  
  #plot mu.m.ivm and HR: LongLasting and Oral
  Day_plot = 185  
  Veta_plot = seq(0, Day_plot, by = .005)
  LegendTex = c("F31-0.6", "Lam-1.0") 
  p_vecLongLasting06 = P_eta.fun(Veta_plot, eta50=ParamLongLasting06$eta50, alpha=ParamLongLasting06$alpha)
  mu_m_ivmLongLasting06 = mu_m - log(1 - p_vecLongLasting06)
  HRLongLasting06 = mu_m_ivmLongLasting06/mu_m
  p_vecLongLasting1 = P_eta.fun(Veta_plot, eta50=ParamLongLasting1$eta50, alpha=ParamLongLasting1$alpha)
  mu_m_ivmLongLasting1 = mu_m - log(1 - p_vecLongLasting1)
  HRLongLasting1 = mu_m_ivmLongLasting1/mu_m
  mu_m_ivm_dataLongLasting = data.frame()
  mu_m_ivm_dataLongLasting = rbind(mu_m_ivm_dataLongLasting, data.frame(
    eta = Veta_plot,
    mu_m_ivm=mu_m_ivmLongLasting06,
    Formulation = LegendTex[1]
  ))
  mu_m_ivm_dataLongLasting = rbind(mu_m_ivm_dataLongLasting, data.frame(
    eta = Veta_plot,
    mu_m_ivm=mu_m_ivmLongLasting1,
    Formulation = LegendTex[2]
  ))
  
  HR_dataLongLasting = data.frame()
  HR_dataLongLasting = rbind(HR_dataLongLasting, data.frame(
    eta = Veta_plot,
    HR=HRLongLasting06,
    Formulation = LegendTex[1]
  ))
  HR_dataLongLasting = rbind(HR_dataLongLasting, data.frame(
    eta = Veta_plot,
    HR=HRLongLasting1,
    Formulation = LegendTex[2]
  ))
  
  
  
  Day_plot = 31  
  Veta_plot = seq(0, Day_plot, by = .005)
  LegendTex = c("BOHEMIA", "RIMDAMAL") 
  p_vecBohemia = P_eta.fun(Veta_plot, eta50=ParamBohemia$eta50, alpha=ParamBohemia$alpha)
  mu_m_ivmBohemia = mu_m - log(1 - p_vecBohemia)
  HRBohemia = mu_m_ivmBohemia/mu_m
  p_vecRimdamal = P_eta.fun(Veta_plot, eta50=ParamRimdamal$eta50, alpha=ParamRimdamal$alpha)
  mu_m_ivmRimdamal = mu_m - log(1 - p_vecRimdamal)
  HRRimdamal = mu_m_ivmRimdamal/mu_m
  mu_m_ivm_dataOral = data.frame()
  mu_m_ivm_dataOral = rbind(mu_m_ivm_dataOral, data.frame(
    eta = Veta_plot,
    mu_m_ivm=mu_m_ivmBohemia,
    Formulation = LegendTex[1]
  ))
  mu_m_ivm_dataOral = rbind(mu_m_ivm_dataOral, data.frame(
    eta = Veta_plot,
    mu_m_ivm=mu_m_ivmRimdamal,
    Formulation = LegendTex[2]
  ))
  HR_Oral = data.frame()
  HR_Oral = rbind(HR_Oral, data.frame(
    eta = Veta_plot,
    HR=HRBohemia,
    Formulation = LegendTex[1]
  ))
  HR_Oral = rbind(HR_Oral, data.frame(
    eta = Veta_plot,
    HR=HRRimdamal,
    Formulation = LegendTex[2]
  ))
  
  Plot_mu_m_ivmLongLasting=ggplot() +
    geom_line(
      data = filter(mu_m_ivm_dataLongLasting, Formulation == "F31-0.6"), 
      aes(x = eta, y = mu_m_ivm, color = "F31-0.6"), linewidth = 1) +
    geom_line(
      data = filter(mu_m_ivm_dataLongLasting, Formulation == "Lam-1.0"), 
      aes(x = eta, y = mu_m_ivm, color = "Lam-1.0"), linewidth = 1) +
    scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7,8),
                       labels = c("0", "1", "2", "3","4","5","6","7","8"), expand = c(0.01, 0),
                       minor_breaks = NULL, limits = c(0,8.5) )+
    scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150, 180), 
                       labels = c("0", "30","60", "90", "120","150",  "180"),
                       limits = c(0,185),expand = c(0.01, 0),minor_breaks = NULL ) +
    theme_light() +
    scale_color_manual(values = c("F31-0.6" = "#ff3355", "Lam-1.0" = "#037153"),
                       labels = c("LAIF31-0.6", "LAIF31-1.0")) +
    labs(
      x = "Time post exposure to IVM (days)", y = "Mosquitoes mortality rate",
      color = "IVM formulation", shape = "IVM formulation",  fill = "IVM formulation" ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(), 
      plot.title = element_text(hjust = 0, size = 16)
    ) +
    ggtitle("(B)")
  
  Plot_mu_m_ivm_oral=ggplot() +
    geom_line(
      data = filter(mu_m_ivm_dataOral, Formulation == "BOHEMIA"), 
      aes(x = eta, y = mu_m_ivm, color = "BOHEMIA"), linewidth = 1) +
    geom_line(
      data = filter(mu_m_ivm_dataOral, Formulation == "RIMDAMAL"), 
      aes(x = eta, y = mu_m_ivm, color = "RIMDAMAL"), linewidth = 1) +
    scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7,8),
                       labels = c("0", "1", "2", "3","4","5","6","7","8"), expand = c(0.01, 0),
                       minor_breaks = NULL , limits = c(0,7.5))+ 
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),labels = c("0", "5", "10", "15","20","25", "30"),
                       limits = c(0,31),expand = c(0.01, 0),minor_breaks = NULL ) +
    
    theme_light() +
    scale_color_manual(values = c("BOHEMIA" = "#D98E04", "RIMDAMAL" = "#205072"),
                       labels = c("BOHEMIA", "RIMDAMAL II")) +
    labs(
      x = "Time post exposure to IVM (days)", y = "Mosquitoes mortality rate",
      color = "IVM formulation", shape = "IVM formulation",  fill = "IVM formulation" ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(), 
      plot.title = element_text(hjust = 0, size = 16)
    ) +
    ggtitle("(E)")
  View(mu_m_ivm_dataLongLasting)
  #View(HR_dataLongLasting)
  Plot_HRLongLasting=ggplot() +
    geom_line(
      data = filter(HR_dataLongLasting, Formulation == "F31-0.6"), 
      aes(x = eta, y = HR, color = "F31-0.6"), linewidth = 1) +
    geom_line(
      data = filter(HR_dataLongLasting, Formulation == "Lam-1.0"), 
      aes(x = eta, y = HR, color = "Lam-1.0"), linewidth = 1) +
    scale_y_continuous(breaks = c(0, 4, 10, 20, 30, 40, 50, 60),
                       labels = c("0", "4", "10", "20", "30","40","50","60"), 
                       expand = c(0.01, 0),minor_breaks = NULL , limits = c(0,65.25))+ 
    scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150, 180), 
                       labels = c("0", "30","60", "90", "120","150",  "180"),
                       limits = c(0,185),expand = c(0.01, 0),minor_breaks = NULL ) +
    #Pour Ajouter des lignes en pointillé
    geom_segment(aes(x = 0, y = 4, xend =180 , yend = 4), 
                 linetype = "dashed", color = "black", linewidth = 0.6) +
    theme_light() +
    scale_color_manual(values = c("F31-0.6" = "#ff3355", "Lam-1.0" = "#037153"),
                       labels = c("LAIF31-0.6", "LAIF31-1.0")) +
    labs(
      x = "Time post exposure to IVM (days)", y = "Crude hazard ratio",
      color = "IVM formulation", shape = "IVM formulation",  fill = "IVM formulation" ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(), 
      plot.title = element_text(hjust = 0, size = 16)
    ) +
    ggtitle("(C)")
  
  
  
  Plot_HR_oral=ggplot() +
    geom_line(
      data = filter(HR_Oral, Formulation == "BOHEMIA"), 
      aes(x = eta, y = HR, color = "BOHEMIA"), linewidth = 1) +
    geom_line(
      data = filter(HR_Oral, Formulation == "RIMDAMAL"), 
      aes(x = eta, y = HR, color = "RIMDAMAL"), linewidth = 1) +
    scale_y_continuous(breaks = c(0, 4, 10, 20, 30, 40, 50, 60),
                       labels = c("0", "4", "10", "20", "30","40","50","60"), 
                       expand = c(0.02, 0),minor_breaks = NULL , limits = c(0,54.58))+
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30),labels = c("0", "5", "10", "15","20","25", "30"),
                       limits = c(0,31),expand = c(0.01, 0),minor_breaks = NULL ) +
    #Pour Ajouter des lignes en pointillé
    geom_segment(aes(x = 0, y = 4, xend =Day_plot , yend = 4), 
                 linetype = "dashed", color = "black", linewidth = 0.6) +
    theme_light() +
    scale_color_manual(values = c("BOHEMIA" = "#D98E04", "RIMDAMAL" = "#205072"),
                       labels = c("BOHEMIA", "RIMDAMAL II")) +
    labs(
      x = "Time post exposure to IVM (days)", y = "Crude hazard ratio",
      color = "IVM formulation", shape = "IVM formulation",  fill = "IVM formulation" ) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(), 
      plot.title = element_text(hjust = 0, size = 16)
    ) +
    ggtitle("(F)")
  
  
  
  plot_list_data=list(Plot_Prop_dead_mosquitoesLongLasting,Plot_mu_m_ivmLongLasting,Plot_HRLongLasting,
                      Plot_Prop_dead_mosquitoesOral,Plot_mu_m_ivm_oral,Plot_HR_oral)
  
  combined_plot=plot_list_data[[1]]+plot_list_data[[2]]+plot_list_data[[3]]+
    plot_list_data[[4]]+plot_list_data[[5]]+plot_list_data[[6]]+plot_layout(ncol = 3, nrow = 2)
}
combined_plot
ggsave("plot_data_Estimaton_With_HR.pdf", combined_plot, width = 12, height = 7.5, units = "in",dpi = 300)

################# FIGURE S1, the initial population Sh0, natural mortality and nu_h 
{
  plot_list_parms <- list()
  
  
  
  #plot Initial population Sh0
  {
    # Initial human population Sh0
    Int_Humanpop <- 8136.10 * (c(12.9, 12.5, 11.5, 13.1, 11.9, 9.3, 7.3, 5.5, 4.4,
                                 3.2, 2.6, 1.8, 1.4, 0.9, 0.7, 0.3, 0.2, 0.2) / 0.997)
    
    # Create a dataframe for Sh0
    Sh0 <- data.frame(
      Age = seq(2.5, 87.5, by = 5),  # Adjust sequence to match length of Int_Humanpop
      Population = Int_Humanpop)
    
    plot_Sh0=ggplot(Sh0, aes(x = Age, y = Population, fill = "Population S[h0]")) +
      geom_bar(stat = "identity") +
      theme_light() +
      scale_x_continuous(breaks = seq(0, 90, by = 5), minor_breaks = NULL,
                         expand = c(0.02, 0), limits = c(0, 90) ) +
      scale_y_continuous( breaks = c(0, 5000, 15000, 30000, 50000, 75000, 100000),
                          labels = c("0", expression(5 %*% 10^3), expression(1.5 %*% 10^4),
                                     expression(3 %*% 10^4), expression(5 %*% 10^4), 
                                     expression(7.5 %*% 10^4), expression(10^5)),
                          minor_breaks = NULL, limits = c(0, 107000) ) +
      labs( x = "Age (years)",  y = "Population", title = "(A)", fill = "") +
      scale_fill_manual( values = c("Population S[h0]" = "#148F77"),
                         labels = expression("Initial human susceptible population "*S[h0]))+
      theme( legend.position = "bottom", # Position de la légende
             legend.title = element_text(size = 12, face = "bold"), 
             legend.text = element_text(size = 12), 
             legend.key.size = unit(0.3, "cm"), # dimension des carrés
             plot.title = element_text(hjust = 0, size = 16), 
             panel.grid.major.x = element_line(color = "gray", linewidth = 0.1)
      ) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
  }
  
  #plot  mu_h
  {
    # Natural mortality rate of human
    mu_h = c(66.8, 7.6, 1.7, 0.9, 1.3, 1.9, 2.4, 2.8, 3.6, 4.7, 6.3, 8.9, 13.2, 19.8,
             31.1, 47.7, 71.3, 110.5, 186.7) / 1000
    
    # Corresponding ages
    ages <- seq(2.5, 92.5, by = 5)
    
    # Create a dataframe for the mortality rate plot
    mu_h <- data.frame(
      Age = ages,
      MortalityRate = mu_h
    )
    
    
    Age_labels <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 
                    60, 65, 70, 75, 80, 85,90)
    
    plot_mu_h=ggplot(mu_h, aes(x = Age, y = MortalityRate, fill = "mu_h")) +
      geom_bar(stat = "identity") +
      theme_light() +
      labs( x = "Age (years)",  y = "Mortality rate", title = "(B)", fill = "") +
      scale_x_continuous(breaks =  seq(0, 95, by = 5), labels = Age_labels,
                         minor_breaks = NULL, expand = c(0.02, 0),limits = c(0, 95)) +
      scale_y_continuous(breaks = seq(0, 0.2, by = 0.05),
                         minor_breaks = NULL, limits = c(0, 0.21)) +
      scale_fill_manual( values = c("mu_h" = "#148F77"),
                         labels = expression("Natural mortality rate for humans "*mu[h])) +
      theme( legend.position = "bottom", # Position de la légende
             legend.title = element_text(size = 12, face = "bold"), 
             legend.text = element_text(size = 12), 
             legend.key.size = unit(0.3, "cm"), # dimension des carrés
             plot.title = element_text(hjust = 0, size = 16), 
             panel.grid.major.x = element_line(color = "gray", linewidth = 0.1)
      ) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.8)))
  }
  
  #plot  nu_hm nu_hf et gamma_h
  {
    #nu_hm(a), nu_hf(a): progression Asympto vers cas clinique
    {
      p_f = 0.15 # Fertility probability for women aged 15 to 40 in Bobo Dioulasso,
      nu_0=1/10; nu_1=1/50;nu_2=1/5
      eps=0.0005
      nu_hm=rep(NA,Nah) # disease mortality for men
      nu_hf = rep(NA, Nah) # disease mortality for female
      for (i in 1:Nah) {
        if (Vah[i] <= 5) { nu_hm[i] = nu_2; nu_hf[i] = nu_2+eps } 
        else if (Vah[i] > 5 && Vah[i] <= 15) {nu_hm[i] = nu_0; nu_hf[i] = nu_0+eps }
        else {nu_hm[i] = nu_1 }
        
        if (Vah[i] > 15 && Vah[i] <= 40) { nu_hf[i] = p_f * nu_0 + (1 - p_f) * nu_1 } 
        else if (Vah[i] > 40) { nu_hf[i] = nu_1+eps}
      }
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
    }
    
    df_segments <- data.frame(xstart = Vah[-length(Vah)], xend = Vah[-1],
                              nu_hm = nu_hm[-length(nu_hm)], nu_hf= nu_hf[-length(nu_hf)],
                              gamma_h=gamma_h[-length(gamma_h)])
    
    cut_age = seq(0, ah_max, by = 5)
    IdSeq <- seq(1, length(Vah), by = 2)
    
    df_points <- data.frame(Age = Vah[IdSeq],nu_hm = nu_hm[IdSeq], 
                            nu_hf= nu_hf[IdSeq],gamma_h=gamma_h[IdSeq])
    
    coeff=10
    plot_nu_gamma=ggplot() +
      geom_segment(data = df_segments, aes(x = xstart, xend = xend, y = nu_hf, yend = nu_hf, color = "nu_hf"), size = 1) +
      geom_point(data = df_points, aes(x = Age, y = nu_hf, color = "nu_hf", shape = "nu_hf"), size = 1.5) +
      geom_segment(data = df_segments, aes(x = xstart, xend = xend, y = nu_hm, yend = nu_hm, color = "nu_hm"), size = 1) +
      geom_point(data = df_points, aes(x = Age, y = nu_hm, color = "nu_hm", shape = "nu_hm"), size = 1.5) +
      geom_segment(data =df_segments,aes(x=xstart,xend=xend,y=gamma_h*coeff,yend=gamma_h*coeff,color="gamma_h"),size = 1) +
      geom_point(data = df_points, aes(x = Age, y = gamma_h * coeff, color = "gamma_h", shape = "gamma_h"), size = 1.5) +
      theme_light() +
      scale_x_continuous(breaks = cut_age, minor_breaks = NULL, expand = c(0.02, 0),
                         limits = c(0, 90))+
      scale_y_continuous( name = "Asymptomatic progression rate",
                          sec.axis = sec_axis(~ . / coeff, name = "Recovery rate")  ) +
      scale_color_manual(
        values = c("nu_hf" = "#76448A", "nu_hm" = "#148F77", "gamma_h" = "#F39C12"),
        labels = c(expression(gamma[h]), expression(nu[h]^bold(italic("♀"))), expression(nu[h]^bold(italic("♂"))))
      ) +
      scale_shape_manual(values = c("nu_hf" = 8, "nu_hm" = 16, "gamma_h" = 18),
                         labels = "" ) +
      labs( x = "Age (years)",y = "", title = "(C)", fill = "" )+
      theme( 
        legend.position = "bottom",
        legend.title = element_blank(), 
        legend.text = element_text(size = 10), 
        legend.spacing.x = unit(2, 'cm'),
        legend.key.width = unit(1.1, "cm"),
        legend.key.size = unit(0.5, "cm"), #dimension des carrés en legende
        plot.title = element_text(hjust = 0, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        panel.grid.major.x = element_line(color = "gray", linewidth = 0.2)
      ) +
      guides(color = guide_legend(override.aes = list(shape = c(18, 8, 16), size = 3)),
             shape = "none" )
  }
  
  plot_list_parms=list(plot_Sh0, plot_mu_h, plot_nu_gamma)
  
  combined_plot=plot_list_parms[[1]]+ plot_list_parms[[2]]+plot_list_parms[[3]]+plot_layout(ncol = 3)
}
combined_plot
ggsave("plot_parameters_&_Initial_pop.pdf", combined_plot, width = 12, height = 4, units = "in",dpi = 300)