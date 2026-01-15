#setwd("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion")
setwd("C:/Users/franck/Desktop/PhD/R studio/Project IVM5")
#setwd("C:/Users/franc/Desktop/PowerPoint PhD/R studio/New folder (2)")
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
library(future)
library(progressr)
library(future.apply)

handlers(global = TRUE)
plan(multisession, workers = parallel::detectCores() - 2)  #all cores minus 2

trapz<-function(x=NULL,y){
  if(is.null(x)){
    x<-1:length(y)
  }
  if(length(y)==1){z<-0}
  else{z<-caTools::trapz(x,y)}
  return(z)
}

P_tau.fun = function(tau, tau50, alpha, Dmax) {
  p_tau = Dmax / (1 + (tau / tau50)^alpha)
  return(p_tau)
}

ModelEtau50Alpha=function(strategy){
  if(strategy==0){#Long lasting 0.6
    #fitting parameter for LAIF 0.6
    # load data ---- 
    # P4D_PK = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/IDR_tab4J_PK.rds")
    P4D_PK = readRDS("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\dataIVM\\IDR_tab4J_PK.rds")
    # mortality data, proportion of dead at 4 days after exposition 
    # gam_all = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/pred_gam_all.rds")
    gam_all = readRDS("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\dataIVM\\pred_gam_all.rds")
    CritereIVMformulation = "F31-0.6"
    sub_PK_F31 = P4D_PK %>% filter(strain == "vk5" & IVM_formulation == CritereIVMformulation)
    sub_PK_F31["DAI"] = as.numeric(sub_PK_F31$DAI)
    #fitting the formula
    SimulDrmLongLasting06= drm(prop_dead ~ DAI, data = sub_PK_F31,
                               fct = LL.4(fixed = c(NA, 0, 1, NA),
                                          names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                               type = "binomial")
    modeldrm=SimulDrmLongLasting06
    summary(SimulDrmLongLasting06)
    
    coef_summary = summary(SimulDrmLongLasting06)$coefficients
    Dmax = 1; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    
    gam_F31 = gam_all %>%
      filter(IVM_formulation == CritereIVMformulation) %>%
      mutate(DAI = Hour / 24, prop_dead = NA)
    Day_seq=length(unique(gam_F31$DAI))
    pred_F31 = predict(SimulDrmLongLasting06, newdata = gam_F31, interval = "confidence")
    predprop_F31 = data.frame(pred_F31) %>% cbind(gam_F31) %>%  mutate(Formulation = "F31-0.6")
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
    # 1. Data from study of Lamidi with cattle injected at 1.0 mg/kg  
    #lamidi_data = read_xlsx("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/Survie_anopheles_sauvage110923.xlsx")
    lamidi_data = read_xlsx("C:\\Users\\franc\\Downloads\\mortality_rate\\mortality_rate\\data\\Survie_anopheles_sauvage110923.xlsx")
    surv_gamb = lamidi_data %>% dplyr::filter(Espece=="gambiae" & Traitement=="IVM") %>% 
      mutate(DAI2=as.numeric(str_extract(DAI, "\\d+")), Status4J=if_else(Temps<=4, 1,0))
    
    # vizu 
    # estimation of proportion of dead mosquitoes at 4 days by cattle and DAI
    # join with GAM_pk estimation of IVM plasma concentration 
    sg_sumdead = surv_gamb %>% group_by(Bovin, DAI2) %>% dplyr::filter(Status4J==1) %>% 
      summarise(ndead=n())
    sg_sum = surv_gamb %>% group_by(Bovin, DAI2) %>%  summarise(ntot=n())
    sg_tot = inner_join(sg_sumdead, sg_sum) 
    sg_tot = sg_tot %>% mutate(prop_dead=ndead/ntot, Day=DAI2) %>% droplevels()
    
    #fitting the formula
    SimulDrmLongLasting1=drm(prop_dead~DAI2, data=sg_tot, fct=LL.4(fixed=c(NA, 0, 1, NA),
                                                                   names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), type="binomial") 
    modeldrm=SimulDrmLongLasting1
    summary(SimulDrmLongLasting1)
    
    coef_summary = summary(SimulDrmLongLasting1)$coefficients
    Dmax = 1; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    DAI = seq(0,166, 1)
    ndat = data.frame(DAI2=DAI, prop_dead=NA)
    pred_Pdead_4J = predict(SimulDrmLongLasting1, newdata = ndat, interval="confidence")
    
    pred_Pdead_4J = cbind(pred_Pdead_4J,DAI)
    
    predprop_F1.0 = data.frame(pred_Pdead_4J) %>% mutate(Formulation = "Lam-1.0")  # predprop_F1.0, 1.0 pour 1mg/kg
    predprop=predprop_F1.0
    
    sg_tot = sg_tot %>%mutate(DAI2 = as.numeric(DAI2))
    sg_tot = sg_tot %>% mutate(Day = as.numeric(Day))
    sub_PK= sg_tot
  }
  
  if(strategy==5){#Long lasting 0.6 train Kis
    #fitting parameter for LAIF 0.6
    # load data ---- View(gam_all)
    # P4D_PK_kis = read.csv("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/data_mortality_rate_4j_LAIF06_KIS.csv")
    P4D_PK_kis = read.csv("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\dataIVM\\data_mortality_rate_4j_LAIF06_KIS.csv")
    # mortality data, proportion of dead at 4 days after exposition 
    # gam_all = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/dataIVM/pred_gam_all.rds")
    gam_all = readRDS("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\dataIVM\\pred_gam_all.rds")
    CritereIVMformulation = "mdc-STM-001"
    sub_PK_F31 = P4D_PK_kis %>% filter(strain == "KIS" & IVM_formulation == CritereIVMformulation)
    sub_PK_F31["DAI"] = as.numeric(sub_PK_F31$DAI)
    #fitting the formula
    SimulDrmLongLasting06_Kis= drm(prop_dead ~ DAI, data = sub_PK_F31,
                                   fct = LL.4(fixed = c(NA, 0, 1, NA),
                                              names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                   type = "binomial")
    modeldrm=SimulDrmLongLasting06_Kis
    summary(SimulDrmLongLasting06_Kis)
    
    coef_summary = summary(SimulDrmLongLasting06_Kis)$coefficients
    Dmax = 1; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    gam_F31 = gam_all %>% filter(IVM_formulation == "F31-0.6") %>% mutate(DAI = Hour / 24, prop_dead = NA)
    Day_seq=length(unique(gam_F31$DAI))
    pred_F31 = predict(SimulDrmLongLasting06_Kis, newdata = gam_F31, interval = "confidence")
    predprop_F31 = data.frame(pred_F31) %>% cbind(gam_F31) %>% mutate(Formulation = "F31-0.6")
    predprop_F31 = predprop_F31[, !names(predprop_F31) %in% 
                                  c("IVM_concentration", "IVM_formulation", "Hour", "prop_dead")]
    predprop=predprop_F31
    
    sub_PK_F31 = sub_PK_F31 %>%mutate(DAI = as.numeric(DAI))
    # sub_PK_F31 = sub_PK_F31 %>% mutate(Day = as.numeric(Day))
    sub_PK=sub_PK_F31
    DAI = seq(0,166, 1)
  }
  
  #fitting parameter for Bohemia and Rimdamal II
  # Oral_Formulation = read.csv("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/DataIVM/Oral_formulationF.csv")
  Oral_Formulation = read.csv("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\dataIVM\\Oral formulation\\Oral_formulationF.csv")
  # View(Oral_Formulation)
  #gam_all = readRDS("/Users/ramsesdjidjou-demasse/Documents/CodeFranck/CodeFranckIvmMdaGender/LastVersion/DataIVM/pred_gam_all.rds")
  gam_all = readRDS("C:\\Users\\franc\\Desktop\\PowerPoint PhD\\R studio\\dataIVM\\pred_gam_all.rds")
  
  if(strategy==2){#Pour BOHEMIA
    BOHEMIA= Oral_Formulation %>%filter(Article == "BOHEMIA")
    BOHEMIA["DAI"] = as.numeric(BOHEMIA$DAI)
    SimulDrmBohemia= drm(prop_dead ~ DAI, data = BOHEMIA,
                         fct = LL.4(fixed = c(NA, 0, 0.8626, NA),
                                    names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    modeldrm=SimulDrmBohemia
    summary(SimulDrmBohemia)
    coef_summary = summary(SimulDrmBohemia)$coefficients
    Dmax = 0.8626; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    #pour faire des figures après?
    #plot(EDM_BOHEMIA)
    DAI = seq(0,166, 1)
    ndat = data.frame(DAI=DAI, prop_dead=NA)
    
    pred_BOHEMIA = predict(SimulDrmBohemia, newdata = ndat, interval = "confidence")
    
    pred_BOHEMIA = data.frame(pred_BOHEMIA) %>% cbind(ndat) %>% mutate(Formulation = "BOHEMIA")
    predprop=pred_BOHEMIA
    
    sub_PK=BOHEMIA %>% mutate(Formulation = "BOHEMIA")
  }
  
  if(strategy==3){#Pour RIMDAMAL smitETAl 2018
    RIMDAMAL= Oral_Formulation %>%filter(Article == "RIMDAMAL")
    RIMDAMAL["DAI"] = as.numeric(RIMDAMAL$DAI)
    SimulDrmRimdamal= drm(prop_dead ~ DAI, data = RIMDAMAL,
                          fct = LL.4(fixed = c(NA, 0, 0.8632, NA),
                                     names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    modeldrm=SimulDrmRimdamal
    summary(SimulDrmRimdamal)
    
    coef_summary = summary(SimulDrmRimdamal)$coefficients
    Dmax = 0.8632; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    #plot(EDM_RIMDAMAL)
    DAI = seq(0,166, 1)
    ndat = data.frame(DAI=DAI, prop_dead=NA)
    
    pred_RIMDAMAL=predict(SimulDrmRimdamal,newdata=ndat,interval="confidence")
    
    pred_RIMDAMAL = data.frame(pred_RIMDAMAL) %>% cbind(ndat) %>%  mutate(Formulation = "RIMDAMAL")
    predprop=pred_RIMDAMAL
    
    sub_PK=RIMDAMAL %>% mutate(Formulation = "RIMDAMAL")
  }
  
  if(strategy==4){#Pour KamauRIMDAMAL  Kamau et Al 2024
    KamauRIMDAMAL= Oral_Formulation %>%filter(Article == "KamauRIMDAMAL")
    KamauRIMDAMAL["DAI"] = as.numeric(KamauRIMDAMAL$DAI)
    SimulDrmKamauRIMDAMAL= drm(prop_dead ~ DAI, data = KamauRIMDAMAL,
                               fct = LL.4(fixed = c(NA, 0, 0.6744, NA),
                                          names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    modeldrm=SimulDrmKamauRIMDAMAL
    summary(SimulDrmKamauRIMDAMAL)
    
    coef_summary = summary(SimulDrmKamauRIMDAMAL)$coefficients
    Dmax = 0.6744; tau50 = coef_summary[2, "Estimate"]; alpha = coef_summary[1, "Estimate"] 
    tau50_se = coef_summary[2, "Std. Error"]; alpha_se = coef_summary[1, "Std. Error"]
    
    #plot(EDM_KamauRIMDAMAL)
    DAI = seq(0,166, 1)
    ndat = data.frame(DAI=DAI, prop_dead=NA)
    
    pred_KamauRIMDAMAL=predict(SimulDrmKamauRIMDAMAL,newdata=ndat,interval="confidence")
    
    pred_KamauRIMDAMAL = data.frame(pred_KamauRIMDAMAL) %>% cbind(ndat) %>% mutate(Formulation = "KamauRIMDAMAL")
    predprop=pred_KamauRIMDAMAL
    
    sub_PK=KamauRIMDAMAL %>% mutate(Formulation = "KamauRIMDAMAL")
  }
  
  
  FunList = list("tau50" = tau50, "Dmax" = Dmax, 
                 "predprop" = predprop, "sub_PK" = sub_PK, "alpha" = alpha, "DAI" = DAI,
                 "model"=modeldrm)
  return(FunList)
}

calculate_CHR = function(tau_times, Vtau, ParamIvmFormulation, mu_m) {
  p_interp = approx(Vtau, P_tau.fun(Vtau, tau50=ParamIvmFormulation$tau50, 
                                    alpha=ParamIvmFormulation$alpha,
                                    Dmax=ParamIvmFormulation$Dmax), tau_times)$y
  
  mu_m_ivm = mu_m - log(1 - p_interp)
  CHR = mu_m_ivm / mu_m
  
  return(CHR)  
}

RhoFunction<-function(Nah,Ntau,Vtau,Vah,IVM_field_dependancy, strategy){
  ParamStrategy=ModelEtau50Alpha(strategy)
  tau50=ParamStrategy$tau50 
  alpha=ParamStrategy$alpha
  
  
  rho=matrix(0,nrow = Nah, ncol = Ntau)
  
  if(IVM_field_dependancy){
    tau_eff=0.5*tau50
    for (a in 1:Nah) {
      for (j in 1:Ntau) {
        rho[a,j] <- alpha * (Vtau[j] / tau_eff)^(alpha - 1) /
          (tau_eff * (1 + (Vtau[j] / tau_eff)^alpha)) #10*(Vtau[j]>tau_eff)
      }
    }
  } else {
    tau_eff=tau50
    for (a in 1:Nah) {
      for (j in 1:Ntau) {
        rho[a,j] <- alpha * (Vtau[j] / tau_eff)^(alpha - 1) /
          (tau_eff * (1 + (Vtau[j] / tau_eff)^alpha)) #10*(Vtau[j]>tau_eff)
      }
    }
  }
  
  #On calcule les proba d'être encore sous IVM pour chaque group
  PropHuman_in_IVM_Group1=Vtau; IdAgeTemoin=3
  PropHuman_in_IVM_Group1[1]=1
  for (j in 2:Ntau) {
    Id=1:j; xId=Vtau[Id]; yId=rho[IdAgeTemoin,Id]
    PropHuman_in_IVM_Group1[j]=exp(-trapz(xId,yId))
  }
  
  FunList = list( "rho"=rho, "PropHuman_in_IVM_Group1"=PropHuman_in_IVM_Group1)
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
                   ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                   dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum){
  q=0.48 #proportion of human male
  k_h=5.5*10^(-4)#1.469e-2;
  Wedge_h=13250  #recrutment rate Bobo dioulasso in 2012
  
  ParamIvmFormulation=ModelEtau50Alpha(strategy)
  p_vec = P_tau.fun(Vtau, tau50=ParamIvmFormulation$tau50, alpha=ParamIvmFormulation$alpha,
                    Dmax=ParamIvmFormulation$Dmax)
  mu_m.ivm = mu_m - log(1 - p_vec)
  
  ModelRho=RhoFunction(Nah,Ntau,Vtau,Vah,IVM_field_dependancy,strategy)
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
        else if (Vah[a] >= 15 & Vah[a] <= 45){
          if (IVM_Pregnancy) {mphi[t,a] = q * phi[t]; fphi[t,a] = (1 - q)*(1 - p_f)*phi[t]
          } else {
            mphi[t,a] = q * phi[t]; fphi[t,a] = 0
          }
        }
        else{mphi[t,a]=q*phi[t]; fphi[t,a]=(1-q)*phi[t]}
      }
    }
  }
  
  # mSh is Susceptible human male and fSh Susceptible human female
  mSh=fSh=matrix(0,nrow = Ntime,ncol = Nah); mAh=fAh=matrix(0,nrow = Ntime,ncol = Nah)
  mIh=fIh=matrix(0,nrow = Ntime,ncol = Nah); mRh=fRh=matrix(0,nrow = Ntime,ncol = Nah)
  Sm=rep(0,Ntime); Im=matrix(0,nrow = Ntime,ncol = Nsigma);
  Sm_ivm=matrix(0,nrow = Ntime,ncol = Ntau);Im_ivm=array(0, dim = c(Nsigma,Ntau,Ntime))
  mSh_ivm=fSh_ivm=array(0, dim = c(Nah,Ntau,Ntime)); mAh_ivm=fAh_ivm=array(0, dim = c(Nah,Ntau,Ntime))
  mIh_ivm=fIh_ivm=array(0, dim = c(Nah,Ntau,Ntime)); mRh_ivm=fRh_ivm=array(0, dim = c(Nah,Ntau,Ntime))
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
        age_indices <- which(Vah <= 5)
      } else if (i == length(age_label_3)) {
        # For the last age range (> 85 years)
        age_indices <- which(Vah > 85)
      } else {
        # For intermediate age ranges
        age_indices <- which(Vah > (i - 1) * 5 & Vah <= i * 5)
      }
      
      # Calculate n_groups (number of individuals in this age group)
      if (length(age_indices) > 0) {
        n_groups <- length(age_indices)
        
        mSh[1, age_indices] = 0.45 * q * Int_values[i] / n_groups
        fSh[1, age_indices] = 0.45 * (1-q) * Int_values[i] / n_groups
        mAh[1, age_indices] = 0.3 * q * Int_values[i] / n_groups
        fAh[1, age_indices] = 0.3 * (1-q) * Int_values[i] / n_groups
        mIh[1, age_indices] = 0.1 * q * Int_values[i] / n_groups
        fIh[1, age_indices] = 0.1 * (1-q) * Int_values[i] / n_groups
        mRh[1, age_indices] = 0.15 * q * Int_values[i] / n_groups
        fRh[1, age_indices] = 0.15 * (1-q) * Int_values[i] / n_groups
      }
    }
    
  }
  
  Sm[1]=6.5*10^6# Wedge_m/mu_m;
  Im[1,1]=0.75*10^6#rep(0.75*10^6, Nsigma)
  
  Sm_ivm[1,]=rep(0,Ntau); Im_ivm[,,1]=matrix(0,nrow = Nsigma,ncol = Ntau)
  mSh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fSh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  mAh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fAh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  mIh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fIh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  mRh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau);fRh_ivm[,,1]=matrix(0,nrow=Nah,ncol=Ntau)
  
  t=1
  Nh[t]= (sum(mSh[t,]+fSh[t,]+mAh[t,]+fAh[t,]+mIh[t,]+fAh[t,]+mRh[t,]+fRh[t,])
          +sum(mSh_ivm[,,t]+fSh_ivm[,,t]+mAh_ivm[,,t]+fAh_ivm[,,t]+mIh_ivm[,,t]
               +fIh_ivm[,,t]+mRh_ivm[,,t]+fRh_ivm[,,t]))
  Nm[t]= Sm[t] + sum(Im[t,])+ sum(Sm_ivm[t,]) +sum(Im_ivm[,,t])
  #Nbrun=0
  tau_eff_age= 0 #initialisation
  lam_hI=rep(0,Ntau);lam_hS=rep(0,Ntau)
  
  for (t in 1:(Ntime-1)) {
    lam_m = theta[t] * (sum(beta_m[]*Im[t,]) + sum(beta_m[]*Im_ivm[,,t]))/Nh[t];
    lam_h = theta[t] * (sum( beta_h[]*(mAh[t,]+fAh[t,]) ) 
                        + sum( bar_beta_h[]*(mIh[t,]+fIh[t,])))/Nh[t];
    
    for (tau in 1:Ntau) {
      lam_hI[tau] =theta[t]*sum(mAh_ivm[,tau,t]+fAh_ivm[,tau,t]+mIh_ivm[,tau,t]+fIh_ivm[,tau,t])/Nh[t];
      lam_hS[tau] =theta[t]*sum(mSh_ivm[,tau,t]+fSh_ivm[,tau,t]+mRh_ivm[,tau,t]+fRh_ivm[,tau,t])/Nh[t]; 
    }
    
    Sm[t+1] = (Wedge_m[t] + Sm[t]/dt )/(1/dt  + mu_m  + lam_h + sum(lam_hI[] + lam_hS[]));
    
    sigma=1
    Im[t+1,sigma]=(Im[t,sigma]/dt + lam_h*Sm[t]/dsigma)/
      (1/dt + 1/dsigma +  mu_m + sum(lam_hI[] + lam_hS[]))
    
    for (sigma in 2:Nsigma ) {
      Im[t+1,sigma]=(Im[t,sigma]/dt + Im[t+1,sigma-1]/dsigma)/
        (1/dt + 1/dsigma +  mu_m + sum(lam_hI[] + lam_hS[]))
    }
    
    ### nouvelle naissance
    a=1 
    mSh[t+1,a]=(mSh[t,a]/dt + q*Wedge_h +k_h*mRh[t,a])/
      (1/dt + 1/dah +  mu_h[a] + mphi[t,a] + lam_m); 
    fSh[t+1,a]=(fSh[t,a]/dt + (1-q)*Wedge_h +k_h*fRh[t,a])/
      (1/dt + 1/dah +  mu_h[a] + fphi[t,a] + lam_m); 
    
    mAh[t+1,a]=0
    fAh[t+1,a]=0
    mIh[t+1,a]=0
    fIh[t+1,a]=0
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
    
    tau = 1
    Sm_ivm[t+1, tau] = 0
    
    # sigma > 1 à tau=1
    for(sigma in 2:Nsigma) {
      Im_ivm[sigma, tau, t+1] = 0
    }
    # tau > 1  
    for(tau in 2:Ntau) {
      Sm_ivm[t+1, tau] = (lam_hS[tau]*Sm[t]+Sm_ivm[t, tau]/dt + Sm_ivm[t+1, tau-1]/dtau) / 
        (1/dt + 1/dtau + mu_m.ivm[tau] + lam_hI[tau] + lam_h)
      
      sigma = 1
      Im_ivm[sigma, tau, t+1] = (Im_ivm[sigma, tau, t]/dt +  Im_ivm[sigma, tau-1, t+1]/dtau + 
                                   (lam_hI[tau] + lam_h)*Sm_ivm[t, tau]/dsigma + 
                                   lam_hI[tau]*Sm[t]/dsigma ) / 
        (1/dt + 1/dtau + 1/dsigma + mu_m.ivm[tau])
      
      #  sigma > 1
      for(sigma in 2:Nsigma) {
        Im_ivm[sigma, tau, t+1] = (Im_ivm[sigma, tau, t]/dt + (lam_hI[tau] + lam_hS[tau])*Im[t, sigma]+
                                     Im_ivm[sigma, tau-1, t+1]/dtau + 
                                     Im_ivm[sigma-1, tau, t+1]/dsigma) / 
          (1/dt + 1/dtau + 1/dsigma + mu_m.ivm[tau])
      }
    }
    
    
    Nh[t+1]= (sum(mSh[t+1,]+fSh[t+1,]+mAh[t+1,]+fAh[t+1,]+mIh[t+1,]+fAh[t+1,]
                  +mRh[t+1,]+fRh[t+1,]) +sum(mSh_ivm[,,t+1]+fSh_ivm[,,t+1] +mAh_ivm[,,t+1]
                                             +fAh_ivm[,,t+1]+mIh_ivm[,,t+1]+fIh_ivm[,,t+1]+mRh_ivm[,,t+1]+fRh_ivm[,,t+1]))
    
    Nm[t+1]= Sm[t+1] + sum(Im[t+1,])+ sum(Sm_ivm[t+1,]) +sum(Im_ivm[,,t+1])
    
    #Nbrun=Nbrun+1
    #print(Nbrun) #le max etant Ntime
  }  
  
  Sh_Tot=time;Ah_Tot=time;Ih_Tot=time;Rh_Tot=time;Ih_15=time;Nh_15=time
  Sm_Tot=time;Im_Tot=time; Im_TotEIP=time
  
  mIh_a=matrix(0,nrow = Ntime,ncol = Nah); fIh_a=matrix(0,nrow = Ntime,ncol = Nah)
  mAh_a=matrix(0,nrow = Ntime,ncol = Nah); fAh_a=matrix(0,nrow = Ntime,ncol = Nah)
  for (t in 1:Ntime){
    Sh_Tot[t]=sum(mSh[t,] + fSh[t,])+sum(mSh_ivm[,,t] +fSh_ivm[,,t])
    Ah_Tot[t]=sum(mAh[t,] + fAh[t,])+sum(mAh_ivm[,,t] +fAh_ivm[,,t])
    Ih_Tot[t]=sum(mIh[t,] + fIh[t,])+sum(mIh_ivm[,,t] +fIh_ivm[,,t])
    Rh_Tot[t]=sum(mRh[t,] + fRh[t,])+sum(mRh_ivm[,,t] +fRh_ivm[,,t])
    Ih_15[t] = sum(mIh[t,6:15] + fIh[t,6:15]) + sum(mIh_ivm[6:15,,t] + fIh_ivm[6:15,,t])
    Nh_15[t]= (sum(mSh[t,6:15]+fSh[t,6:15]+mAh[t,6:15]+fAh[t,6:15]+mIh[t,6:15]+fAh[t,6:15]+mRh[t,6:15]
                   +fRh[t,6:15]) +sum(mSh_ivm[6:15,,t]+fSh_ivm[6:15,,t]+mAh_ivm[6:15,,t]+fAh_ivm[6:15,,t]
                                     +mIh_ivm[6:15,,t]+fIh_ivm[6:15,,t]+mRh_ivm[6:15,,t]+fRh_ivm[6:15,,t]))
    Sm_Tot[t]=sum(Sm_ivm[t,]) + Sm[t] 
    Im_Tot[t]=sum(Im_ivm[,,t]) + sum(Im[t,])
    Im_TotEIP[t]=sum(Im_ivm[8:Nsigma,,t]) + sum(Im[t,8:Nsigma])
    
    for (a in 1:Nah){
      mIh_a[t,a]=mIh[t,a] + sum(mIh_ivm[a,,t]); fIh_a[t,a]=fIh[t,a] + sum(fIh_ivm[a,,t])
      mAh_a[t,a]=mAh[t,a] + sum(mAh_ivm[a,,t]); fAh_a[t,a]=fAh[t,a] + sum(fAh_ivm[a,,t])
    }
  }
  
  IdBegin_Camp= 1+floor(t_begin_Camp/dt)
  IndexTime=IdBegin_Camp:Ntime
  # Delta=trapz(time,Ih+Ih_ivmTot)/(tmax*Ih[1]+ Ih_ivmTot[1])
  # Delta=1-Delta
  Delta=trapz(time[IndexTime],Ih_Tot[IndexTime])
  
  PropIhTot=Ih_Tot/Nh;PropShTot=Sh_Tot/Nh;PropIh_15=Ih_15/Nh_15
  PropAhTot=Ah_Tot/Nh;PropRhTot=Rh_Tot/Nh
  PropImTot=Im_Tot/Nm;PropSmTot=Sm_Tot/Nm; PropImEIP=Im_TotEIP/Nm;
  
  IdTop=which(PropIhTot[IndexTime]==min(PropIhTot[IndexTime]))
  DurIvmEffect=time[IdTop]
  
  MinPropIhTot=min(PropIhTot[IndexTime])
  PropIhTot0=PropIhTot[t_begin_Camp/dt]
  
  IndexTimeHR = (IdBegin_Camp):Ntime  
  IdTopHR = which(PropIhTot[IndexTimeHR] >= PropIhTot0)[1] 
  
  # IndexTimeHR = (IdBegin_Camp+ 1 ):Ntime 
  # IdTopHR = tail(which(PropIhTot[IndexTimeHR] != PropIhTot0), 1)
  
  if (!is.na(IdTopHR)) {
    IdTopHR = IdTopHR 
    DurIvmEffectHR = time[IdTopHR]
  } else {
    IdTopHR = 0
    DurIvmEffectHR = 0
  }
  
  #Optimisation Timing, et cumulative gain prev, avec le deuxieme gain prev
  # DeltaOptim est le nombre cummule' de Ih (with IVM)
  Nyr=360
  IdBegin_sum= 1+floor(t_begin_sum/dt)
  IdEnd_sum= 1+floor((t_begin_sum+Nyr)/dt)
  IndexOptim=IdBegin_sum:IdEnd_sum
  DeltaOptim=sum(Ih_Tot[IndexOptim])
  DeltaOptim_15=sum(Ih_15[IndexOptim])
  Delta_Tot=sum(Ih_Tot[IndexOptim]+Ah_Tot[IndexOptim])
  
  cut_age = seq(0, ah_max, by = 5)
  
  # time_idx = IdBegin_Camp + IdTop[1]  # Prendre le premier si plusieurs
  # time_idx = min(time_idx, nrow(mIh_a))  # S'assurer qu'on ne dépasse pas
  # 
  # # Regrouper en classes de 5 ans (dah = 1)
  # mIh_age_groups = sapply(1:(length(cut_age) - 1), function(i) {
  #   idx_debut = cut_age[i] + 1
  #   idx_fin = min(cut_age[i + 1], ncol(mIh_a))  
  #   
  #   if (idx_debut <= idx_fin) {
  #     sum(mIh_a[time_idx, idx_debut:idx_fin])
  #   } else {
  #     0
  #   }
  # })
  # 
  # fIh_age_groups = sapply(1:(length(cut_age) - 1), function(i) {
  #   idx_debut = cut_age[i] + 1
  #   idx_fin = min(cut_age[i + 1], ncol(mIh_a))  
  #   
  #   if (idx_debut <= idx_fin) {
  #     sum(fIh_a[time_idx, idx_debut:idx_fin])
  #   } else {
  #     0
  #   }
  # })
  
  FunList = list("time" = time, "rho"=rho,"mum.ivm"=mu_m.ivm, "mSh"=mSh,"Sh_Tot"=Sh_Tot,
                 "Ah_Tot"=Ah_Tot,"Ih_Tot"=Ih_Tot,"Rh_Tot"=Rh_Tot,"Nh"=Nh,"Nm"=Nm, 
                 "Sm_Tot"=Sm_Tot,"Im_Tot"=Im_Tot,"Delta"=Delta,"lam_h"=lam_h,
                 "lam_m"=lam_m, "lam_hI"=lam_hI, "lam_hS"=lam_hS, "phi"=phi,
                 "mIh_a"=mIh_a, "fIh_a"=fIh_a,"mAh_a"=mAh_a, "fAh_a"=fAh_a, "PropRhTot"=PropRhTot,
                 "PropIhTot"=PropIhTot, "PropImTot"=PropImTot, "DurIvmEffect"=DurIvmEffect,
                 "PropShTot"=PropShTot, "PropSmTot"=PropSmTot, "PropAhTot"=PropAhTot,
                 "MinPropIhTot"=MinPropIhTot, "PropIhTot0"=PropIhTot0, "cut_age"=cut_age,
                 "DurIvmEffectHR"=DurIvmEffectHR, "DeltaOptim"=DeltaOptim,"PropIh_15"=PropIh_15,
                 "IndexOptim"=IndexOptim,"DeltaOptim_15"=DeltaOptim_15, "PropImEIP"=PropImEIP,
                 "Delta_Tot"=Delta_Tot)
  return(FunList)
}

# baseline Parameters
{
  Gap=0.25;
  dt=Gap;tmax = 641;time=seq(0,tmax,by=dt); Ntime=length(time)
  dtau=1;tau_max= 300;Vtau=seq(0,tau_max,by=dtau);Ntau=length(Vtau);
  dah=1;ah_max=90;Vah=seq(0,ah_max,by=dah); Nah=length(Vah)
  dtau=dtau;tau_max= tau_max;Vtau=seq(0,tau_max,by=dtau);Ntau=length(Vtau);
  dsigma=dtau;sigma_max=60;Vsigma=seq(0,sigma_max,by=dsigma); Nsigma=length(Vsigma)
  mu_m = 0.13  
  t_begin_sum=250 # par defaut,
  IVM_Pregnancy=0
  EIP=8
  # theta=0.5
  Dur_cycle=7 
  # t_begin_Camp=250
  
  #mu_h(a) natural mortality
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
  
  #delta_h(a), Ref Quentin Richard and al Bobo Dioulasso. Malaria induced mortality
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

  #nu_hm(a), nu_hf(a): progression Asympto to clinical cases
  {
    p_f = 0.15  # Fertility probability for women aged 15 to 40 in Bobo Dioulasso,
    nu_0=1/11; nu_1=1/170;nu_2=1/3
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
  
  #gamma_h(a): Recovery
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
  
  #Name for save
  Name_Prev_mos <- function(Prev_mos10) {
    ifelse(Prev_mos10 == 1, "10", "5")
  }
  # Theta, beta_h(a) et Wedge_M: HBR, human infectiosity et mosquitoes recruitment
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
  
  q_theta = 1/4
  theta_0 = 0.5
  T_theta = 120 # High seasonal transmission duration
  t_begin_Higly_Seas = 250
  ta = t_begin_Higly_Seas
  tb = ta + T_theta
  q_Wedge=0.1
  
  Model_theta = function(Saisonality, Prev_mos10) {
    
    if (Saisonality == 1) {
      theta = p_theta.fun(time, theta_0, q_theta, T_theta, t_begin_Higly_Seas, Gap)
      
      # Calcul of Wedge_m with seasonality
      Wedge_m0 <<- 1.825 * 10^6
      Wedge_m <<- ifelse(time < ta - Gap, Wedge_m0,
                         ifelse(time < ta, Wedge_m0 + ((time - (ta - Gap)) / Gap) * (q_Wedge * Wedge_m0),
                                ifelse(time <= tb, (1 + q_Wedge) * Wedge_m0,
                                       ifelse(time <= tb + Gap, (1 + q_Wedge) * Wedge_m0 - ((time - tb) / Gap) * (q_Wedge * Wedge_m0),
                                              Wedge_m0))))
      
      if (Prev_mos10) {
        beta_m <<- 0.0416 * (Vsigma > EIP)
        alpha_1 <<- 0.122       # 0.071[0.023; 0.175]
        alpha_2 <<- 0.17 
      } else { # Prev_mos5_percent
        beta_m <<- 0.0416 * (Vsigma > EIP)
        alpha_1 <<- 0.122       # 0.071[0.023; 0.175]
        alpha_2 <<- 0.17 
      }
      
    } else { #without seasonality
      theta = rep(theta_0, Ntime)
      
      # Wedge_m constant
      Wedge_m0 <<- 1.33 * 10^6
      Wedge_m <<- rep(Wedge_m0, Ntime)
      
      if (Prev_mos10) {
        beta_m <<- 0.0416 * (Vsigma > EIP)
        alpha_1 <<- 0.122       # 0.071[0.023; 0.175]
        alpha_2 <<- 0.17 
      } else { # Prev_mos5_percent
        beta_m <<- 0.063 * (Vsigma > EIP)
        alpha_1 <<- 0.071        # 0.071[0.023; 0.175]
        alpha_2 <<- 0.15 
      }
    }
    
    # Calcul de beta_h (commun aux deux cas)
    beta_h <<- rep(NA, Nah)
    
    # Define the function G(a)
    G <<- function(a) {
      return(22.7 * a * exp(-0.0934 * a))
    }
    
    # Compute beta_h for each age in Vah
    beta_h <<- alpha_1 * (G(Vah) ^ alpha_2)
    bar_beta_h <<- 0.8 * beta_h
    
    return(theta)
  }
}

#Figure 4: Comparaison seasonality
IVM_Pregnancy=0
Prev_mos10=1 #Prevalence of infectious mosquitoes 10%
VectPropIVM =c(0.5,0.7,0.9)
Vect_t_begin_Camp= c(220, 250, 280)
CoefDelay=0.85
IVM_field_dependancy=0

GainPrev_array_oral_Seas=array(0, dim = c(3,length(VectPropIVM),length(Vect_t_begin_Camp)))
GainPrev_array_oral_Per=array(0, dim = c(3,length(VectPropIVM)))
#GainPrev_matrix_oral_Seas[1,,] for BOHEMIA, GainPrev_matrix_oral_Seas[2,,] for RIMDAMAL
#GainPrev_matrix_oral_Seas[3,,] for Kamau RIMDAMAL
#same for GainPrev_array_oral_Per

VectNumber_of_cycle=c(1,2,3,4)
GainPrev_array_LAIF_Seas=array(0, dim = c(3,length(VectPropIVM),length(VectNumber_of_cycle),
                                          length(Vect_t_begin_Camp)))
GainPrev_array_LAIF_Per=array(0, dim = c(3,length(VectPropIVM),length(VectNumber_of_cycle)))
# GainPrev_array_LAIF_Seas[1,,,] for "mdc-STM-001_0.6-vk5", GainPrev_array_LAIF_Seas[2,,,] for "mdc-STM-001_1.0"
# GainPrev_array_LAIF_Seas[3,,,] for "mdc-STM-001_0.6-kis"
#same GainPrev_array_LAIF_Per

#Prevalence total d'infections
VectNumber_of_cycle=c(1,2,3,4)
GainPrev_Tot_LAIF_Seas=array(0, dim = c(3,length(VectPropIVM),length(VectNumber_of_cycle),
                                          length(Vect_t_begin_Camp)))
GainPrev_Tot_oral_Seas=array(0, dim = c(3,length(VectPropIVM),length(Vect_t_begin_Camp)))

Nbrun=0
#without seasonality
{
  Saisonality=0
  theta=Model_theta(Saisonality, Prev_mos10)
  t_begin_Camp=250
  t_begin_sum=250
  #Without IVM
  {
    strategy=0
    PropIVM=0
    VectTime_between_cycles=c(0,0,0)
    Number_of_cycle=1
    
    OutPut = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                      mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                      ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                      dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
    
    ModelPrevCum= OutPut
    DeltaOptimInit=sum(ModelPrevCum$Ih_Tot[ModelPrevCum$IndexOptim]) 
  }
  
  for (j in 1:length(VectPropIVM)) {
    PropIVM = VectPropIVM[j]
    print(paste0("PropIVM=",PropIVM))
    
    #On tourne le modèle pour Bohemia
    strategy=2;NbCycle=3
    Number_of_cycle=NbCycle;VectTime_between_cycles=c(30,30,30)
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                           dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
    
    #on récupère les valeurs ici
    print("Bohemia")
    print(VectTime_between_cycles)
    VectTime_between_cycles_Bohemia=VectTime_between_cycles
    ModelBohemia=ModelOutput
    
    GainPrev_array_oral_Per[1,j]=1-(ModelBohemia$DeltaOptim/DeltaOptimInit)
    
    #On tourne le modèle pour Rimdamal
    strategy=3;NbCycle=4
    Number_of_cycle=NbCycle;VectTime_between_cycles=c(30,30,30)
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                           dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
    
    #on récupère les valeurs ici
    print("Rimdamal")
    print(VectTime_between_cycles)
    VectTime_between_cycles_Rimdamal=VectTime_between_cycles
    #VectTime_between_cycles=0.7*c(39.50, 79.25-39.50, 79.25-(79.25-39.50))
    ModelRimdamal=ModelOutput
    
    GainPrev_array_oral_Per[2,j]=1-(ModelRimdamal$DeltaOptim/DeltaOptimInit)
    
    #On tourne le modèle pour KamauRimdamal
    strategy=4;NbCycle=4
    Number_of_cycle=NbCycle;VectTime_between_cycles=c(30,30,30)
    ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                           mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                           ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                           dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
    
    #on récupère les valeurs ici
    print("KamauRimdamal")
    print(VectTime_between_cycles)
    VectTime_between_cycles_KamauRimdamal=VectTime_between_cycles
    #VectTime_between_cycles=0.7*c(39.50, 79.25-39.50, 79.25-(79.25-39.50))
    ModelKamauRimdamal=ModelOutput
    
    GainPrev_array_oral_Per[3,j]=1-(ModelKamauRimdamal$DeltaOptim/DeltaOptimInit)
    
    
    #On tourne le modèle pour LongLasting 06
    strategy=0;NbCycle=4
    ModelLongLasting06=list()
    Number_of_cycle=NbCycle;VectTime_between_cycles=c(60,60,60)
   
    for (Number_of_cycle in 1:NbCycle) {
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
      ModelLongLasting06[[Number_of_cycle]]=ModelOutput
      
      GainPrev_array_LAIF_Per[1,j,Number_of_cycle]=1-(ModelLongLasting06[[Number_of_cycle]]$DeltaOptim/DeltaOptimInit)
    }
    #on récupère les valeurs ici
    print("LongLasting0.6")
    print(VectTime_between_cycles)
    VectTime_between_cycles_LongLasting06=VectTime_between_cycles
    
    #On tourne le modèle pour LongLasting 06 kis
    strategy=5;NbCycle=4
    ModelLongLasting06_kis=list()
    Number_of_cycle=NbCycle;VectTime_between_cycles=c(60,60,60)
    for (Number_of_cycle in 1:NbCycle) {
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
      ModelLongLasting06_kis[[Number_of_cycle]]=ModelOutput
      
      GainPrev_array_LAIF_Per[3,j,Number_of_cycle]=1-(ModelLongLasting06_kis[[Number_of_cycle]]$DeltaOptim/DeltaOptimInit)
    }
    #on récupère les valeurs ici
    print("LongLasting0.6 kis")
    print(VectTime_between_cycles)
    VectTime_between_cycles_LongLasting06_kis=VectTime_between_cycles
    
    #On tourne le modèle pour LongLasting 1
    strategy=1;NbCycle=4
    ModelLongLasting1=list()
    Number_of_cycle=NbCycle;VectTime_between_cycles=c(60,60,60)
    
    for (Number_of_cycle in 1:NbCycle) {
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
      ModelLongLasting1[[Number_of_cycle]]=ModelOutput
      
      GainPrev_array_LAIF_Per[2,j,Number_of_cycle]=1-(ModelLongLasting1[[Number_of_cycle]]$DeltaOptim/DeltaOptimInit)
    }
    #on récupère les valeurs ici
    print("LoongLasting1")
    print(VectTime_between_cycles)
    VectTime_between_cycles_LongLasting1=VectTime_between_cycles
    
    Nbrun= Nbrun + 1
    print(Nbrun)
  }
  
  save(GainPrev_array_oral_Per, GainPrev_array_LAIF_Per, VectTime_between_cycles_LongLasting1,
       VectTime_between_cycles_LongLasting06_kis, VectTime_between_cycles_LongLasting06,
       VectTime_between_cycles_KamauRimdamal, VectTime_between_cycles_Rimdamal,
       VectTime_between_cycles_Bohemia,VectPropIVM, VectNumber_of_cycle,Nbrun,
       file = paste0("ResultatsPrev_Perennial_field_dep_", IVM_field_dependancy,"_Pregnant",
                     IVM_Pregnancy,"Prev_mos",Name_Prev_mos(Prev_mos10), ".RData"))
  
}

load(paste0("ResultatsPrev_Perennial_field_dep_", IVM_field_dependancy,"_Pregnant",
            IVM_Pregnancy,"Prev_mos",Name_Prev_mos(Prev_mos10), ".RData"))
ls()

#with Seasonality
{
  Saisonality=1
  theta=Model_theta(Saisonality, Prev_mos10)
  
  Vect_t_begin_Camp= c(220, 250, 280)
  
  
  for (p in 1:length(Vect_t_begin_Camp)) {
      t_begin_Camp= Vect_t_begin_Camp[p]
      
      t_begin_sum= Vect_t_begin_Camp[p]
      
      #Without IVM
      strategy=0; PropIVM=0
      VectTime_between_cycles=c(0,0,0)
      Number_of_cycle=1
      
      OutPut = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                        mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                        ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                        dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
      
      ModelPrevCum= OutPut
      DeltaOptimInit_Seas=sum(ModelPrevCum$Ih_Tot[ModelPrevCum$IndexOptim]) 
      DeltaOptim_Tot_Init_Seas=ModelPrevCum$Delta_Tot
      
      for (j in 1:length(VectPropIVM)) {
        PropIVM = VectPropIVM[j]
      
      print(paste0("Time begin cycle=",t_begin_Camp))
      print(paste0("PropIVM=",PropIVM))
      
      #On tourne le modèle pour LongLasting 06
      ModelLongLasting06=list()
      VectTime_between_cycles=VectTime_between_cycles_LongLasting06
      
      strategy=0;NbCycle=4
      for (Number_of_cycle in 1:NbCycle) {
        ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                               mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                               ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                               dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
        ModelLongLasting06[[Number_of_cycle]]=ModelOutput
        GainPrev_array_LAIF_Seas[1,j,Number_of_cycle,p]=1-(ModelLongLasting06[[Number_of_cycle]]$DeltaOptim/DeltaOptimInit_Seas)
        
        GainPrev_Tot_LAIF_Seas[1,j,Number_of_cycle,p]=1-(ModelLongLasting06[[Number_of_cycle]]$Delta_Tot/DeltaOptim_Tot_Init_Seas)
      }
      
      #On tourne le modèle pour LongLasting 06-kis
      ModelLongLasting06_kis=list()
      
      VectTime_between_cycles=VectTime_between_cycles_LongLasting06_kis
      strategy=5;NbCycle=4
      for (Number_of_cycle in 1:NbCycle) {
        ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                               mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                               ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                               dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
        ModelLongLasting06_kis[[Number_of_cycle]]=ModelOutput
        GainPrev_array_LAIF_Seas[3,j,Number_of_cycle,p]=1-(ModelLongLasting06_kis[[Number_of_cycle]]$DeltaOptim/DeltaOptimInit_Seas)
        GainPrev_Tot_LAIF_Seas[3,j,Number_of_cycle,p]=1-(ModelLongLasting06_kis[[Number_of_cycle]]$Delta_Tot/DeltaOptim_Tot_Init_Seas)
      }
      
      #On tourne le modèle pour LongLasting 1
      ModelLongLasting1=list()
      
      VectTime_between_cycles= VectTime_between_cycles_LongLasting1
      strategy=1;NbCycle=4
      for (Number_of_cycle in 1:NbCycle) {
        ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                               mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                               ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                               dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
        ModelLongLasting1[[Number_of_cycle]]=ModelOutput
        GainPrev_array_LAIF_Seas[2,j,Number_of_cycle,p]=1-(ModelLongLasting1[[Number_of_cycle]]$DeltaOptim/DeltaOptimInit_Seas)
        GainPrev_Tot_LAIF_Seas[2,j,Number_of_cycle,p]= 1-(ModelLongLasting1[[Number_of_cycle]]$Delta_Tot/DeltaOptim_Tot_Init_Seas)
        }
      
      #On tourne le modèle pour Bohemia
      strategy=2; Number_of_cycle=3
      VectTime_between_cycles=VectTime_between_cycles_Bohemia
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
      ModelBohemia=ModelOutput
      GainPrev_array_oral_Seas[1,j,p]=1-(ModelBohemia$DeltaOptim/DeltaOptimInit_Seas)
      GainPrev_Tot_oral_Seas[1,j,p]=1-(ModelBohemia$Delta_Tot/DeltaOptim_Tot_Init_Seas)
      
      #On tourne le modèle pour Rimdamal
      strategy=3; Number_of_cycle=4
      VectTime_between_cycles=VectTime_between_cycles_Rimdamal
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
      ModelRimdamal=ModelOutput
      GainPrev_array_oral_Seas[2,j,p]=1-(ModelRimdamal$DeltaOptim/DeltaOptimInit_Seas)
      GainPrev_Tot_oral_Seas[2,j,p]=1-(ModelRimdamal$Delta_Tot/DeltaOptim_Tot_Init_Seas)
      
      #On tourne le modèle pour KamauRimdamal
      strategy=4; Number_of_cycle=4
      VectTime_between_cycles=VectTime_between_cycles_KamauRimdamal
      ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                             mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                             ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                             dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
      ModelKamauRimdamal=ModelOutput
      GainPrev_array_oral_Seas[3,j,p]=1-(ModelKamauRimdamal$DeltaOptim/DeltaOptimInit_Seas)
      GainPrev_Tot_oral_Seas[3,j,p]=1-(ModelKamauRimdamal$Delta_Tot/DeltaOptim_Tot_Init_Seas)
    }
    Nbrun= Nbrun + 1
    print(Nbrun)
  }
  
  save(GainPrev_array_oral_Seas, GainPrev_array_LAIF_Seas, GainPrev_Tot_oral_Seas, GainPrev_Tot_LAIF_Seas,
       VectPropIVM, VectNumber_of_cycle,Nbrun, 
       file = paste0("ResultatsPrev_Seasonal_field_dep_", IVM_field_dependancy,"_Pregnant",
                     IVM_Pregnancy,"Prev_mos",Name_Prev_mos(Prev_mos10), ".RData"))
  
}

load(paste0("ResultatsPrev_Seasonal_field_dep_", IVM_field_dependancy,"_Pregnant",
            IVM_Pregnancy,"Prev_mos",Name_Prev_mos(Prev_mos10), ".RData"))
ls()

# pdf("Comparisons_Seas_vs_Peren.pdf", width=14, height=9.2)
# par(mfrow=c(3,4), oma=c(8.1,1,2,0.2))
# mar = c(4.1, 2.6, 2, .1)

#PLOT
{ 
  FigName=paste0("Comparisons2_Seas_vs_Peren_field_dep_", IVM_field_dependancy,"_Pregnant",
                 IVM_Pregnancy,"Prev_mos",Name_Prev_mos(Prev_mos10), ".pdf")
  pdf(FigName, width=7.5, height=6.4)
  par(mfrow=c(2,2), oma=c(8.3,1,2,0.2))
  mar = c(2, 3.35, 1.5, .1)
  
  VectPropIVM_target = c(0.5, 0.7, 0.9)
  strategies_oral_names = c("BOH", expression("RII"[S]), expression("RII"[K]))
  strategies_LAIF_names = c("0.6-vk5", "1.0", "0.6-kis")
  
  strategy_colors = c("#D98E04", "#205072", "#8B008B", "#ff3355", "#037153", "#8B4513")
  pch_t = c(1, 15, 16)
  
  GrNumber = 0
  
  for (nc in VectNumber_of_cycle) {
    GrNumber = GrNumber + 1
    par(mar = mar)
    
    # Créer la table des comparaisons pour tous les PropIVM_target
    All_comparisons = data.frame()
    
    for (PropIVM_target in VectPropIVM_target) {
      IndexPropIVM = which.min(abs(VectPropIVM - PropIVM_target))
      
      # Comparaisons pour les stratégies orales
      for (oral_strat in 1:3) {
        for (t_idx in 1:length(Vect_t_begin_Camp)) {
          val = (GainPrev_array_oral_Seas[oral_strat, IndexPropIVM, t_idx] - 
                   GainPrev_array_oral_Per[oral_strat, IndexPropIVM]) * 100
          
          strat_name = c("BOH", "RII[S]", "RII[K]")[oral_strat]
          comp_name = paste0("Delta(", strat_name, ",", strat_name, ")_t", t_idx)
          
          All_comparisons = rbind(All_comparisons,
                                  data.frame(PropIVM=PropIVM_target, 
                                             Strategy=oral_strat, 
                                             t_begin=t_idx,
                                             Comparison=comp_name, 
                                             Value=val, 
                                             Type="Oral"))
        }
      }
      
      # Comparaisons pour les stratégies LAIF
      for (laif_strat in 1:3) {
        for (t_idx in 1:length(Vect_t_begin_Camp)) {
          val = (GainPrev_array_LAIF_Seas[laif_strat, IndexPropIVM, nc, t_idx] - 
                   GainPrev_array_LAIF_Per[laif_strat, IndexPropIVM, nc]) * 100
          
          laif_name = strategies_LAIF_names[laif_strat]
          if (laif_name == "1.0") {
            laif_name = "'1.0'"
          }
          comp_name = paste0("Delta(", laif_name, ",", laif_name, ")_t", t_idx)
          
          All_comparisons = rbind(All_comparisons,
                                  data.frame(PropIVM=PropIVM_target,
                                             Strategy=laif_strat + 3, 
                                             t_begin=t_idx,
                                             Comparison=comp_name, 
                                             Value=val, 
                                             Type="LAIF"))
        }
      }
    }
    
    # Calculer min, max, median pour chaque comparaison unique
    Comparisons_list = unique(All_comparisons$Comparison)
    Summary_stats = data.frame()
    
    for (comp in Comparisons_list) {
      df_comp = subset(All_comparisons, Comparison == comp)
      min_val = min(df_comp$Value)
      max_val = max(df_comp$Value)
      med_val = median(df_comp$Value)
      strategy_num = df_comp$Strategy[1]
      t_begin_num = df_comp$t_begin[1]
      type_comp = df_comp$Type[1]
      
      Summary_stats = rbind(Summary_stats, 
                            data.frame(Comparison=comp, Min=min_val, Max=max_val, Median=med_val, 
                                       Strategy=strategy_num, t_begin=t_begin_num, Type=type_comp))
    }
    
    # Positions sur l'axe 
    n_strategies = 6
    n_t_values = length(Vect_t_begin_Camp)
    spacing = 1
    group_width = n_t_values + 0.5
    
    x_positions = c()
    x_labels = c()
    x_label_positions = c()
    
    for (strat in 1:n_strategies) {
      base_x = (strat - 1) * group_width + 1
      strat_x_positions = base_x + (0:(n_t_values-1)) * spacing
      x_positions = c(x_positions, strat_x_positions)
      
      x_label_positions = c(x_label_positions, base_x + (n_t_values-1) * spacing / 2)
      
      if (strat <= 3) {
        if (strat == 1) {
          x_labels = c(x_labels, "Delta(BOH,BOH)")
        } else if (strat == 2) {
          x_labels = c(x_labels, "Delta(RII[S],RII[S])")
        } else {
          x_labels = c(x_labels, "Delta(RII[K],RII[K])")
        }
      } else {
        laif_idx = strat - 3
        if (laif_idx == 1) {
          x_labels = c(x_labels, "Delta(0.6-vk5,0.6-vk5)")
        } else if (laif_idx == 2) {
          x_labels = c(x_labels, "Delta('1.0','1.0')")
        } else {
          x_labels = c(x_labels, "Delta(0.6-kis,0.6-kis)")
        }
      }
    }
    
    ylim = c(-8.5, 8.5)
    
    plot(-1, 1, type="n", xlim=c(0.5, max(x_positions) + 0.5), 
         ylim=ylim, xaxt="n", yaxt="n", xlab="", ylab="")
    abline(h=0, lty=2, col="gray")
    box()
    axis(2, at = c(-8, -4, 0, 4, 8), cex.axis=0.75)
    
    if(is.element(GrNumber,c(1,3))){
      mtext(expression("Performance" ~ "(" * "seasonal vs perennial" * ")" ~ "(%)"), 
            side = 2, adj = 0.5, cex = .75, line = 2, font = .8) 
    }
    
    axis(1, at = x_label_positions, labels = FALSE)
    if(is.element(GrNumber,c(3,4))){
      for (i in 1:length(x_labels)) {
        label_text = x_labels[i]
        mtext(parse(text=label_text), side=1, line=1, at=x_label_positions[i], cex=0.7, las=2)
      }
    }
    
    # Plot arrows et points
    for (i in 1:nrow(Summary_stats)) {
      row = Summary_stats[i,]
      strat = row$Strategy
      t_val = row$t_begin
      x_pos = x_positions[(strat-1) * n_t_values + t_val]
      
      arrows(x_pos, row$Min, x_pos, row$Max, 
             angle = 90, code = 3, length = 0.015, 
             col = strategy_colors[strat], lwd = 2)
      points(x_pos, row$Median, pch = pch_t[t_val], 
             bg = strategy_colors[strat], col = strategy_colors[strat], cex = 1.6)
    }
    
    par(xpd=NA)
    text(0.2, ylim[2]*1.3, paste("(", LETTERS[GrNumber], ")", sep=""), cex=1.2, adj=0)
    mtext(bquote("#cycle for mdc-STM-001: "~ n[c]==.(nc)), side=3, cex = 0.95, line=1)
    par(xpd=FALSE)
  }
  
  reset <- function() {
    par(mfrow=c(1,1), oma=rep(0,4), mar=rep(0,4), new=TRUE)
    plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
  }
  reset()
  
  LEGEND1 = c( expression("RIMDAMAL II" [S]), expression("RIMDAMAL II" [K]), "BOHEMIA")
  strategy_colors1 = c("#205072", "#8B008B", "#D98E04")
  
  LEGEND2 = c("mdc-STM-001−0.6-vk5", "mdc-STM-001−0.6-kis", "mdc-STM-001−1.0")
  strategy_colors2 = c("#ff3355", "#8B4513", "#037153")
  
  legend("bottom", legend=LEGEND1, col=strategy_colors1, pch=16, xpd=NA, horiz=TRUE,
         inset=c(0,0.027), cex=0.75, bty="n")
  
  legend("bottom", legend=LEGEND2, col=strategy_colors2, pch=16, xpd=NA, horiz=TRUE,
         inset=c(0,0.007), cex=0.75, bty="n")
  
  t_legend = c("IVM starts 30 days before peak season", "IVM starts exactly at peak season", 
               "IVM starts 30 days after peak season") 
  
  legend("bottom", legend=t_legend, pch=pch_t, col="black",
         xpd=NA, horiz=TRUE, inset=c(0,-0.013), cex=0.75, bty="n")
  
  dev.off()
}

dev.off()

#Figure:S7 & S8 Scenarios_Formulations LAIF with Seasonality
###########################
IVM_Pregnancy=0
Prev_mos10=1 
Vect_t_begin_Camp=c(220, 250, 280)
t_begin_plot=210
tmaxmonths=(tmax-t_begin_plot)/30

CoefDelay=0.85
IVM_field_dependancy=0

ColVect = c("#ff3355","#8B4513", "#5f05f3", "#037153", "#205072","#8B008B", "#D98E04")
LineVect <- c(1, 3, 2, 4, 5)
LC=2.5
sink(paste("DurEntreCycle_field_dep",IVM_field_dependancy,"_Pregnant",
           IVM_Pregnancy,"Prev_mos",Name_Prev_mos(Prev_mos10),".txt",sep=""))


#without IVM
PropIVM=0;
strategy=0;Number_of_cycle=1;VectTime_between_cycles=1;#pour le premier cycle, cette variable est inutile
Saisonality=1
theta=Model_theta(Saisonality, Prev_mos10)
ModelZero=list()
for (p in 1:length(Vect_t_begin_Camp)) {
  t_begin_Camp=Vect_t_begin_Camp[p]
  ModelOutput = ModelIVM(PropIVM,VectTime_between_cycles,Number_of_cycle,t_begin_Camp,Dur_cycle,
                         mu_h,delta_h,bar_beta_h,beta_h,nu_hm,nu_hf,gamma_h,theta,beta_m,Wedge_m,mu_m,
                         ah_max,dah,Nah,tmax,Gap,dt,time,Ntime,p_f, IVM_Pregnancy,dtau,tau_max,Vtau,Ntau,
                         dsigma,sigma_max,Vsigma,Nsigma,strategy,IVM_field_dependancy,t_begin_sum)
  ModelZero[[p]]=ModelOutput
}

VectPropIVM=c(0.5,0.7,0.9)

# =============
# FONCTION DE SIMULATION 

run_simulation = function(type) {
  
  simulation_results = list()
  counter = 1
  
  for (PropIVM in VectPropIVM) {
    for (p in 1:length(Vect_t_begin_Camp)) {
      
      t_begin_Camp = Vect_t_begin_Camp[p]
      
      print(paste0("Time begin cycle=", t_begin_Camp))
      print(paste0("PropIVM=", PropIVM))
      
      if (type == "laif") {
        print("Simulation mdc-STM-001-0.6-vk5")
        ModelLongLasting06 = list()
        
          strategy = 0; NbCycle = 4
          VectTime_between_cycles=c(60,60,60)
          print("LongLasting0.6")
          print(VectTime_between_cycles)
          VectTime_between_cycles_LAIF = VectTime_between_cycles
        # Simulation avec seasonalité
        Saisonality = 1
        theta = Model_theta(Saisonality, Prev_mos10)
        
        for (Number_of_cycle in 1:NbCycle) {
          ModelOutput = ModelIVM(PropIVM, VectTime_between_cycles_LAIF, Number_of_cycle, t_begin_Camp, Dur_cycle,
                                 mu_h, delta_h, bar_beta_h, beta_h, nu_hm, nu_hf, gamma_h, theta, beta_m, Wedge_m, mu_m,
                                 ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy, dtau, tau_max, Vtau, Ntau,
                                 dsigma, sigma_max, Vsigma, Nsigma, strategy, IVM_field_dependancy, t_begin_sum)
          ModelLongLasting06[[Number_of_cycle]] = ModelOutput
        }
        
        # Sauvegarde résultats LAIF
        simulation_results[[counter]] = list(
          PropIVM = PropIVM, t_begin_Camp = t_begin_Camp, p = p,
          type = "laif", ModelLongLasting06 = ModelLongLasting06,
          VectTime_between_cycles_LAIF = VectTime_between_cycles_LAIF
        )
      }
      
      if (type == "oral") {
        
        # --- BOHEMIA ---
        {
          strategy = 2; NbCycle = 3
          VectTime_between_cycles_Bohemia=c(30,30,30)
          print("Bohemia")
          print(VectTime_between_cycles_Bohemia)
          
          Saisonality = 1
          theta = Model_theta(Saisonality, Prev_mos10)
          Number_of_cycle = 3
          ModelBohemia = ModelIVM(PropIVM, VectTime_between_cycles_Bohemia, Number_of_cycle, t_begin_Camp, Dur_cycle,
                                  mu_h, delta_h, bar_beta_h, beta_h, nu_hm, nu_hf, gamma_h, theta, beta_m, Wedge_m, mu_m,
                                  ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy, dtau, tau_max, Vtau, Ntau,
                                  dsigma, sigma_max, Vsigma, Nsigma, strategy, IVM_field_dependancy, t_begin_sum)
        }
        
        # --- RIMDAMAL ---
        {
          strategy = 3; NbCycle = 4
          VectTime_between_cycles_Rimdamal=c(30,30,30)
          print("Rimdamal")
          print(VectTime_between_cycles_Rimdamal)
          
          Saisonality = 1
          theta = Model_theta(Saisonality, Prev_mos10)
          Number_of_cycle = 4
          ModelRimdamal = ModelIVM(PropIVM, VectTime_between_cycles_Rimdamal, Number_of_cycle, t_begin_Camp, Dur_cycle,
                                   mu_h, delta_h, bar_beta_h, beta_h, nu_hm, nu_hf, gamma_h, theta, beta_m, Wedge_m, mu_m,
                                   ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy, dtau, tau_max, Vtau, Ntau,
                                   dsigma, sigma_max, Vsigma, Nsigma, strategy, IVM_field_dependancy, t_begin_sum)
        }
        
        # --- KAMAU RIMDAMAL ---
        {
          strategy = 4; NbCycle = 4
          VectTime_between_cycles_KamauRimdamal=c(30,30,30)
          print("KamauRimdamal")
          print(VectTime_between_cycles_KamauRimdamal)
          
          Saisonality = 1
          theta = Model_theta(Saisonality, Prev_mos10)
          Number_of_cycle = 4
          ModelKamauRimdamal = ModelIVM(PropIVM, VectTime_between_cycles_KamauRimdamal, Number_of_cycle, t_begin_Camp, Dur_cycle,
                                        mu_h, delta_h, bar_beta_h, beta_h, nu_hm, nu_hf, gamma_h, theta, beta_m, Wedge_m, mu_m,
                                        ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy, dtau, tau_max, Vtau, Ntau,
                                        dsigma, sigma_max, Vsigma, Nsigma, strategy, IVM_field_dependancy, t_begin_sum)
        }
        
        # Sauvegarde résultats ORAL
        simulation_results[[counter]] = list(
          PropIVM = PropIVM, t_begin_Camp = t_begin_Camp, p = p,
          type = "oral",  ModelBohemia = ModelBohemia,
          ModelRimdamal = ModelRimdamal, ModelKamauRimdamal = ModelKamauRimdamal,
          VectTime_between_cycles_Bohemia = VectTime_between_cycles_Bohemia,
          VectTime_between_cycles_Rimdamal = VectTime_between_cycles_Rimdamal,
          VectTime_between_cycles_KamauRimdamal = VectTime_between_cycles_KamauRimdamal
        )
      }
      
      counter = counter + 1
    }
  }
  
  return(simulation_results)
}

plot_results = function(simulation_results, type) {
  
  if (type == "laif") {
    FigName = paste0("Scenario_SeaS_Ih_LAIF06_vk5_field_dep", IVM_field_dependancy, "_Pregnant",
                     IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".pdf")
    pdf(FigName, width = 10.5, height = 8.2)
    par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.5, 2, .1))
    par(mfrow = c(length(VectPropIVM), length(Vect_t_begin_Camp)))
    
    print(paste0("CoefDelay=", CoefDelay))
    
    GrNumber = 0
    counter = 1
    
    for (PropIVM in VectPropIVM) {
      for (p in 1:length(Vect_t_begin_Camp)) {
        
        # Récupération des données
        data_point = simulation_results[[counter]]
        ModelLongLasting06 = data_point$ModelLongLasting06
        
        MaxY = 0.155 #max(ModelZero[[p]]$PropIhTot[t_begin_plot/dt:tmax/dt]) * 1.1
        
        GrNumber = GrNumber + 1
        plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_plot, tmax), ylab = "",
             ylim = c(0, MaxY), cex.lab = 1.2, yaxt = "n", xaxt = "n")
        axis(1, at = seq(t_begin_plot, tmax, by = 60), labels = seq(0, tmaxmonths, by = 2), las = 1)
        y_vals = c(0, 0.03, 0.06, 0.09, 0.12, 0.15)
        axis(2, at = y_vals, labels = y_vals * 100, las = 1)
        par(xpd = NA)
        text(t_begin_plot - 0.05, MaxY * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.6)
        par(xpd = FALSE)
        
        lines(time, ModelZero[[p]]$PropIhTot, lwd = LC, col = ColVect[1], lty = LineVect[1])
        
        rect(t_begin_Higly_Seas, par("usr")[3], t_begin_Higly_Seas + T_theta, par("usr")[4], 
             col = adjustcolor("#D3D3D3", alpha.f = 0.3), border = NA)
        
        abline(v = data_point$t_begin_Camp, col = adjustcolor("black", alpha.f = 0.8), lty = 1, lwd = 2)
        
        Nbb = 1
        for (nc in 1:4) {
          Nbb = Nbb + 1
          lines(time, ModelLongLasting06[[nc]]$PropIhTot, lwd = LC, cex.lab = 1.2, 
                lty = LineVect[Nbb], col = ColVect[Nbb])
        }
        
        Peak_season = ifelse(p == 1, "IVM starts 30 days before peak season", 
                             ifelse(p == 2, "IVM starts at peak season", 
                                    "IVM starts 30 days after peak season"))
        
        if (is.element(GrNumber, c(1, 2, 3))) {
          mtext(Peak_season, side = 3, adj = 0.5, cex = 0.8, line = 1, font = 1)
        }
        
        if (is.element(GrNumber, c(1, 4, 7))) {
          mtext("Clinical cases (all ages; %)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8)
          mtext(bquote("Coverage of the TP: " ~ .(PropIVM * 100) ~ "%"), side = 2, adj = 0.5, cex = 0.9, line = 3, font = 1)
        }
        
        if (is.element(GrNumber, c(7, 8, 9))) {
          mtext("Time (Months)", side = 1, adj = 0.5, cex = 0.9, line = 3, font = 0.8)
        }
        
        counter = counter + 1
      }
    }
    
    # Légende LAIF
    LegendTex = c("Without IVM", expression(n[c] == 1),
                  expression(n[c] == 2), expression(n[c] == 3), expression(n[c] == 4))
    
    reset()
    legend("bottom", legend = c(bquote("mdc-STM-001-0.6-vk5")), 
           xpd = NA, horiz = TRUE, inset = c(0, 0.035), bty = "n", cex = 1)
    
    legend("bottom", legend = LegendTex, xpd = NA, horiz = TRUE, inset = c(0, 0.01), 
           bty = "n", lty = LineVect, lw = 1.75, col = ColVect, cex = 1)
    
    dev.off()
    getwd()
  }
  
  if (type == "oral") {
    FigName = paste0("Scenario_SeaS_Ih_oralF_field_dep", IVM_field_dependancy, "_Pregnant",
                     IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".pdf")
    pdf(FigName, width = 10.5, height = 8.2)
    par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.5, 2, .1))
    par(mfrow = c(length(VectPropIVM), length(Vect_t_begin_Camp)))
    
    GrNumber = 0
    counter = 1
    
    for (PropIVM in VectPropIVM) {
      for (p in 1:length(Vect_t_begin_Camp)) {
        
        # Récupération des données
        data_point = simulation_results[[counter]]
        ModelBohemia = data_point$ModelBohemia
        ModelRimdamal = data_point$ModelRimdamal
        ModelKamauRimdamal = data_point$ModelKamauRimdamal
        
        MaxY =0.152 #max(ModelZero[[p]]$PropIhTot) * 1.1
        
        GrNumber = GrNumber + 1
        plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_plot, tmax), ylab = "",
             ylim = c(0, MaxY), cex.lab = 1.2, yaxt = "n", xaxt = "n")
        axis(1, at = seq(t_begin_plot, tmax, by = 60), labels = seq(0, tmaxmonths, by = 2), las = 1)
        y_vals = c(0, 0.03, 0.06, 0.09, 0.12, 0.15)
        axis(2, at = y_vals, labels = y_vals * 100, las = 1)
        par(xpd = NA)
        text(t_begin_plot - 0.05, MaxY * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.6)
        par(xpd = FALSE)
        
        lines(time, ModelZero[[p]]$PropIhTot, lwd = LC, col = ColVect[1], lty = LineVect[1])
        
        rect(t_begin_Higly_Seas, par("usr")[3], t_begin_Higly_Seas + T_theta, par("usr")[4], 
             col = adjustcolor("#D3D3D3", alpha.f = 0.3), border = NA)
        
        abline(v = data_point$t_begin_Camp, col = adjustcolor("black", alpha.f = 0.8), lty = 1, lwd = 2)
        
        lines(time, ModelBohemia$PropIhTot, lwd = LC, cex.lab = 1.2, 
              lty = LineVect[3], col = ColVect[2])
        lines(time, ModelRimdamal$PropIhTot, lwd = LC, cex.lab = 1.2, 
              lty = LineVect[4], col = ColVect[5])
        lines(time, ModelKamauRimdamal$PropIhTot, lwd = LC, cex.lab = 1.2, 
              lty = LineVect[2], col = ColVect[4])
        
        Peak_season = ifelse(p == 1, "IVM starts 30 days before peak season", 
                             ifelse(p == 2, "IVM starts at peak season", 
                                    "IVM starts 30 days after peak season"))
        
        if (is.element(GrNumber, c(1, 2, 3))) {
          mtext(Peak_season, side = 3, adj = 0.5, cex = 0.8, line = 1, font = 1)
        }
        
        if (is.element(GrNumber, c(1, 4, 7))) {
          mtext("Clinical cases (all ages; %)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8)
          mtext(bquote("Coverage of the TP: " ~ .(PropIVM * 100) ~ "%"), side = 2, adj = 0.5, cex = 0.9, line = 3, font = 1)
        }
        
        if (is.element(GrNumber, c(7, 8, 9))) {
          mtext("Time (Months)", side = 1, adj = 0.5, cex = 0.9, line = 2, font = 0.8)
        }
        
        counter = counter + 1
      }
    }
    LegendTex = c("Without IVM", "BOHEMIA", expression("RIMDAMAL II"[S]), expression("RIMDAMAL II"[K]))
    
    reset()
    legend("bottom", legend = LegendTex,
           xpd = NA, horiz = TRUE, inset = c(0, 0.03), 
           bty = "n", lty = c(1, 2, 4, 3), lw = 1.75,
           col = c("#ff3355", "#8B4513", "#205072", "#037153", "#5f05f3", "#8B008B", "#D98E04"), cex = 1)
    
    dev.off()
  }
}

# Function pour run + save + plot 
run_and_plot = function(type) {
  # run
  results = run_simulation(type)
  
  # Sauvegarde
  filename_rdata = paste0("IVM_Seas_", type, ".RData")
  
  simulation_data = list(
    results = results,
    type = type,
    parameters = list(
      VectPropIVM = VectPropIVM,  Vect_t_begin_Camp = Vect_t_begin_Camp,
      CoefDelay = CoefDelay, IVM_field_dependancy = IVM_field_dependancy,
      IVM_Pregnancy = IVM_Pregnancy, Prev_mos10 = Prev_mos10, simulation_date = Sys.time()
    )
  )
  
  save(simulation_data, ModelZero, time, t_begin_plot, tmax, tmaxmonths, 
       t_begin_Higly_Seas, T_theta, LC, ColVect, LineVect, 
       file = filename_rdata)
  
  cat("Données sauvegardées dans:", filename_rdata, "\n")
  
  # Plot
  plot_results(results, type)
  
  return(list(results = results, filename = filename_rdata))
}
laif_data = run_and_plot("laif")
oral_data = run_and_plot("oral")

# Function reload and replot
reload_and_plot = function(rdata_filename) {
  load(rdata_filename, envir = .GlobalEnv)
  
  type = simulation_data$type
  
  plot_results(simulation_data$results, type)
}
reload_and_plot("IVM_Seas_laif.RData")
reload_and_plot("IVM_Seas_oral.RData")

sink()

### LAIF 3 months Between cycles
run_simulation = function(type) {
  
  simulation_results = list()
  counter = 1
  
  for (PropIVM in VectPropIVM) {
    for (p in 1:length(Vect_t_begin_Camp)) {
      
      t_begin_Camp = Vect_t_begin_Camp[p]
      
      print(paste0("Time begin cycle=", t_begin_Camp))
      print(paste0("PropIVM=", PropIVM))
      
      if (type == "laif") {
        print("Simulation mdc-STM-001-0.6-vk5")
        ModelLongLasting06 = list()
        
        strategy = 0; NbCycle = 4
        VectTime_between_cycles=c(90,90,90)
        print("LongLasting0.6")
        print(VectTime_between_cycles)
        VectTime_between_cycles_LAIF = VectTime_between_cycles
        # Simulation avec seasonalité
        Saisonality = 1
        theta = Model_theta(Saisonality, Prev_mos10)
        
        for (Number_of_cycle in 1:NbCycle) {
          ModelOutput = ModelIVM(PropIVM, VectTime_between_cycles_LAIF, Number_of_cycle, t_begin_Camp, Dur_cycle,
                                 mu_h, delta_h, bar_beta_h, beta_h, nu_hm, nu_hf, gamma_h, theta, beta_m, Wedge_m, mu_m,
                                 ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy, dtau, tau_max, Vtau, Ntau,
                                 dsigma, sigma_max, Vsigma, Nsigma, strategy, IVM_field_dependancy, t_begin_sum)
          ModelLongLasting06[[Number_of_cycle]] = ModelOutput
        }
        
        # Sauvegarde résultats LAIF
        simulation_results[[counter]] = list(
          PropIVM = PropIVM, t_begin_Camp = t_begin_Camp, p = p,
          type = "laif", ModelLongLasting06 = ModelLongLasting06,
          VectTime_between_cycles_LAIF = VectTime_between_cycles_LAIF
        )
      }
      
      if (type == "oral") {
        
        # --- BOHEMIA ---
        {
          strategy = 2; NbCycle = 3
          VectTime_between_cycles_Bohemia=c(30,30,30)
          print("Bohemia")
          print(VectTime_between_cycles_Bohemia)
          
          Saisonality = 1
          theta = Model_theta(Saisonality, Prev_mos10)
          Number_of_cycle = 3
          ModelBohemia = ModelIVM(PropIVM, VectTime_between_cycles_Bohemia, Number_of_cycle, t_begin_Camp, Dur_cycle,
                                  mu_h, delta_h, bar_beta_h, beta_h, nu_hm, nu_hf, gamma_h, theta, beta_m, Wedge_m, mu_m,
                                  ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy, dtau, tau_max, Vtau, Ntau,
                                  dsigma, sigma_max, Vsigma, Nsigma, strategy, IVM_field_dependancy, t_begin_sum)
        }
        
        # --- RIMDAMAL ---
        {
          strategy = 3; NbCycle = 4
          VectTime_between_cycles_Rimdamal=c(30,30,30)
          print("Rimdamal")
          print(VectTime_between_cycles_Rimdamal)
          
          Saisonality = 1
          theta = Model_theta(Saisonality, Prev_mos10)
          Number_of_cycle = 4
          ModelRimdamal = ModelIVM(PropIVM, VectTime_between_cycles_Rimdamal, Number_of_cycle, t_begin_Camp, Dur_cycle,
                                   mu_h, delta_h, bar_beta_h, beta_h, nu_hm, nu_hf, gamma_h, theta, beta_m, Wedge_m, mu_m,
                                   ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy, dtau, tau_max, Vtau, Ntau,
                                   dsigma, sigma_max, Vsigma, Nsigma, strategy, IVM_field_dependancy, t_begin_sum)
        }
        
        # --- KAMAU RIMDAMAL ---
        {
          strategy = 4; NbCycle = 4
          VectTime_between_cycles_KamauRimdamal=c(30,30,30)
          print("KamauRimdamal")
          print(VectTime_between_cycles_KamauRimdamal)
          
          Saisonality = 1
          theta = Model_theta(Saisonality, Prev_mos10)
          Number_of_cycle = 4
          ModelKamauRimdamal = ModelIVM(PropIVM, VectTime_between_cycles_KamauRimdamal, Number_of_cycle, t_begin_Camp, Dur_cycle,
                                        mu_h, delta_h, bar_beta_h, beta_h, nu_hm, nu_hf, gamma_h, theta, beta_m, Wedge_m, mu_m,
                                        ah_max, dah, Nah, tmax, Gap, dt, time, Ntime, p_f, IVM_Pregnancy, dtau, tau_max, Vtau, Ntau,
                                        dsigma, sigma_max, Vsigma, Nsigma, strategy, IVM_field_dependancy, t_begin_sum)
        }
        
        # Sauvegarde résultats ORAL
        simulation_results[[counter]] = list(
          PropIVM = PropIVM, t_begin_Camp = t_begin_Camp, p = p,
          type = "oral",  ModelBohemia = ModelBohemia,
          ModelRimdamal = ModelRimdamal, ModelKamauRimdamal = ModelKamauRimdamal,
          VectTime_between_cycles_Bohemia = VectTime_between_cycles_Bohemia,
          VectTime_between_cycles_Rimdamal = VectTime_between_cycles_Rimdamal,
          VectTime_between_cycles_KamauRimdamal = VectTime_between_cycles_KamauRimdamal
        )
      }
      
      counter = counter + 1
    }
  }
  
  return(simulation_results)
}

plot_results = function(simulation_results, type) {
  
  if (type == "laif") {
    FigName = paste0("Scenario_SeaS_Ih_LAIF06_vk5_3mois_field_dep", IVM_field_dependancy, "_Pregnant",
                     IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".pdf")
    pdf(FigName, width = 10.5, height = 8.2)
    par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.5, 2, .1))
    par(mfrow = c(length(VectPropIVM), length(Vect_t_begin_Camp)))
    
    print(paste0("CoefDelay=", CoefDelay))
    
    GrNumber = 0
    counter = 1
    
    for (PropIVM in VectPropIVM) {
      for (p in 1:length(Vect_t_begin_Camp)) {
        
        # Récupération des données
        data_point = simulation_results[[counter]]
        ModelLongLasting06 = data_point$ModelLongLasting06
        
        MaxY = 0.155 #max(ModelZero[[p]]$PropIhTot[t_begin_plot/dt:tmax/dt]) * 1.1
        
        GrNumber = GrNumber + 1
        plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_plot, tmax), ylab = "",
             ylim = c(0, MaxY), cex.lab = 1.2, yaxt = "n", xaxt = "n")
        axis(1, at = seq(t_begin_plot, tmax, by = 60), labels = seq(0, tmaxmonths, by = 2), las = 1)
        y_vals = c(0, 0.03, 0.06, 0.09, 0.12, 0.15)
        axis(2, at = y_vals, labels = y_vals * 100, las = 1)
        par(xpd = NA)
        text(t_begin_plot - 0.05, MaxY * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.6)
        par(xpd = FALSE)
        
        lines(time, ModelZero[[p]]$PropIhTot, lwd = LC, col = ColVect[1], lty = LineVect[1])
        
        rect(t_begin_Higly_Seas, par("usr")[3], t_begin_Higly_Seas + T_theta, par("usr")[4], 
             col = adjustcolor("#D3D3D3", alpha.f = 0.3), border = NA)
        
        abline(v = data_point$t_begin_Camp, col = adjustcolor("black", alpha.f = 0.8), lty = 1, lwd = 2)
        
        Nbb = 1
        for (nc in 1:4) {
          Nbb = Nbb + 1
          lines(time, ModelLongLasting06[[nc]]$PropIhTot, lwd = LC, cex.lab = 1.2, 
                lty = LineVect[Nbb], col = ColVect[Nbb])
        }
        
        Peak_season = ifelse(p == 1, "IVM starts 30 days before peak season", 
                             ifelse(p == 2, "IVM starts at peak season", 
                                    "IVM starts 30 days after peak season"))
        
        if (is.element(GrNumber, c(1, 2, 3))) {
          mtext(Peak_season, side = 3, adj = 0.5, cex = 0.8, line = 1, font = 1)
        }
        
        if (is.element(GrNumber, c(1, 4, 7))) {
          mtext("Clinical cases (all ages; %)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8)
          mtext(bquote("Coverage of the TP: " ~ .(PropIVM * 100) ~ "%"), side = 2, adj = 0.5, cex = 0.9, line = 3, font = 1)
        }
        
        if (is.element(GrNumber, c(7, 8, 9))) {
          mtext("Time (Months)", side = 1, adj = 0.5, cex = 0.9, line = 3, font = 0.8)
        }
        
        counter = counter + 1
      }
    }
    
    # Légende LAIF
    LegendTex = c("Without IVM", expression(n[c] == 1),
                  expression(n[c] == 2), expression(n[c] == 3), expression(n[c] == 4))
    
    reset()
    legend("bottom", legend = c(bquote("mdc-STM-001-0.6-vk5")), 
           xpd = NA, horiz = TRUE, inset = c(0, 0.035), bty = "n", cex = 1)
    
    legend("bottom", legend = LegendTex, xpd = NA, horiz = TRUE, inset = c(0, 0.01), 
           bty = "n", lty = LineVect, lw = 1.75, col = ColVect, cex = 1)
    
    dev.off()
    getwd()
  }
  
  if (type == "oral") {
    FigName = paste0("Scenario_SeaS_Ih_oralF_field_dep", IVM_field_dependancy, "_Pregnant",
                     IVM_Pregnancy, "Prev_mos", Name_Prev_mos(Prev_mos10), ".pdf")
    pdf(FigName, width = 10.5, height = 8.2)
    par(oma = c(5.5, 1, 2, .1), mar = c(3, 3.5, 2, .1))
    par(mfrow = c(length(VectPropIVM), length(Vect_t_begin_Camp)))
    
    GrNumber = 0
    counter = 1
    
    for (PropIVM in VectPropIVM) {
      for (p in 1:length(Vect_t_begin_Camp)) {
        
        # Récupération des données
        data_point = simulation_results[[counter]]
        ModelBohemia = data_point$ModelBohemia
        ModelRimdamal = data_point$ModelRimdamal
        ModelKamauRimdamal = data_point$ModelKamauRimdamal
        
        MaxY =0.152 #max(ModelZero[[p]]$PropIhTot) * 1.1
        
        GrNumber = GrNumber + 1
        plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_plot, tmax), ylab = "",
             ylim = c(0, MaxY), cex.lab = 1.2, yaxt = "n", xaxt = "n")
        axis(1, at = seq(t_begin_plot, tmax, by = 60), labels = seq(0, tmaxmonths, by = 2), las = 1)
        y_vals = c(0, 0.03, 0.06, 0.09, 0.12, 0.15)
        axis(2, at = y_vals, labels = y_vals * 100, las = 1)
        par(xpd = NA)
        text(t_begin_plot - 0.05, MaxY * 1.15, paste("(", LETTERS[GrNumber], ")", sep = ""), cex = 1.6)
        par(xpd = FALSE)
        
        lines(time, ModelZero[[p]]$PropIhTot, lwd = LC, col = ColVect[1], lty = LineVect[1])
        
        rect(t_begin_Higly_Seas, par("usr")[3], t_begin_Higly_Seas + T_theta, par("usr")[4], 
             col = adjustcolor("#D3D3D3", alpha.f = 0.3), border = NA)
        
        abline(v = data_point$t_begin_Camp, col = adjustcolor("black", alpha.f = 0.8), lty = 1, lwd = 2)
        
        lines(time, ModelBohemia$PropIhTot, lwd = LC, cex.lab = 1.2, 
              lty = LineVect[3], col = ColVect[2])
        lines(time, ModelRimdamal$PropIhTot, lwd = LC, cex.lab = 1.2, 
              lty = LineVect[4], col = ColVect[5])
        lines(time, ModelKamauRimdamal$PropIhTot, lwd = LC, cex.lab = 1.2, 
              lty = LineVect[2], col = ColVect[4])
        
        Peak_season = ifelse(p == 1, "IVM starts 30 days before peak season", 
                             ifelse(p == 2, "IVM starts at peak season", 
                                    "IVM starts 30 days after peak season"))
        
        if (is.element(GrNumber, c(1, 2, 3))) {
          mtext(Peak_season, side = 3, adj = 0.5, cex = 0.8, line = 1, font = 1)
        }
        
        if (is.element(GrNumber, c(1, 4, 7))) {
          mtext("Clinical cases (all ages; %)", side = 2, adj = 0.5, cex = .85, line = 2, font = .8)
          mtext(bquote("Coverage of the TP: " ~ .(PropIVM * 100) ~ "%"), side = 2, adj = 0.5, cex = 0.9, line = 3, font = 1)
        }
        
        if (is.element(GrNumber, c(7, 8, 9))) {
          mtext("Time (Months)", side = 1, adj = 0.5, cex = 0.9, line = 2, font = 0.8)
        }
        
        counter = counter + 1
      }
    }
    LegendTex = c("Without IVM", "BOHEMIA", expression("RIMDAMAL II"[S]), expression("RIMDAMAL II"[K]))
    
    reset()
    legend("bottom", legend = LegendTex,
           xpd = NA, horiz = TRUE, inset = c(0, 0.03), 
           bty = "n", lty = c(1, 2, 4, 3), lw = 1.75,
           col = c("#ff3355", "#8B4513", "#205072", "#037153", "#5f05f3", "#8B008B", "#D98E04"), cex = 1)
    
    dev.off()
  }
}

# Function pour run + save + plot 
run_and_plot = function(type) {
  # run
  results = run_simulation(type)
  
  # Sauvegarde
  filename_rdata = paste0("IVM_Seas3mois_", type, ".RData")
  
  simulation_data = list(
    results = results,
    type = type,
    parameters = list(
      VectPropIVM = VectPropIVM,  Vect_t_begin_Camp = Vect_t_begin_Camp,
      CoefDelay = CoefDelay, IVM_field_dependancy = IVM_field_dependancy,
      IVM_Pregnancy = IVM_Pregnancy, Prev_mos10 = Prev_mos10, simulation_date = Sys.time()
    )
  )
  
  save(simulation_data, ModelZero, time, t_begin_plot, tmax, tmaxmonths, 
       t_begin_Higly_Seas, T_theta, LC, ColVect, LineVect, 
       file = filename_rdata)
  
  cat("Données sauvegardées dans:", filename_rdata, "\n")
  
  # Plot
  plot_results(results, type)
  
  return(list(results = results, filename = filename_rdata))
}
laif_data = run_and_plot("laif")
oral_data = run_and_plot("oral")

# Function reload and replot
reload_and_plot = function(rdata_filename) {
  load(rdata_filename, envir = .GlobalEnv)
  
  type = simulation_data$type
  
  plot_results(simulation_data$results, type)
}
reload_and_plot("IVM_Seas3mois_laif.RData")


sink()
