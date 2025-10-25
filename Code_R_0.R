# CALCUL DU RAYON SPECTRAL DU R_0 
setwd("C:/Users/franc/Desktop/PowerPoint PhD/R studio/Project IVM4/R0")

library(dplyr)
library(readxl)
library(drc)
library(caTools)
library(stringr)
library(pracma)
library(expm)

trapz = function(x = NULL, y) {
  if (is.null(x)) {x = 1:length(y)}
  if (length(y) == 1) { z = 0}
  else {z = caTools::trapz(x, y)}
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
    sub_PK_F31 = P4D_PK %>%
      filter(strain == "vk5" & IVM_formulation == CritereIVMformulation)
    sub_PK_F31["DAI"] = as.numeric(sub_PK_F31$DAI)
    #fitting the formula
    SimulDrmLongLasting06= drm(prop_dead ~ DAI, data = sub_PK_F31,
                               fct = LL.4(fixed = c(NA, 0, 1, NA),
                                          names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                               type = "binomial")
    summary(SimulDrmLongLasting06)
    Dmax=1
    tau50=SimulDrmLongLasting06$coefficients[2]
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
    Dmax=1
    tau50=SimulDrmLongLasting1$coefficients[2]
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
                         fct = LL.4(fixed = c(NA, 0, 0.8632, NA),
                                    names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    summary(SimulDrmBohemia)
    Dmax=0.8626
    tau50=SimulDrmBohemia$coefficients[2]
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
  
  if(strategy==3){#Pour RIMDAMAL smitETAl 2018
    RIMDAMAL= Oral_Formulation %>%filter(Article == "RIMDAMAL")
    RIMDAMAL["DAI"] = as.numeric(RIMDAMAL$DAI)
    SimulDrmRimdamal= drm(prop_dead ~ DAI, data = RIMDAMAL,
                          fct = LL.4(fixed = c(NA, 0, 0.8632, NA),
                                     names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    summary(SimulDrmRimdamal)
    Dmax=0.8632
    tau50=SimulDrmRimdamal$coefficients[2]
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
  
  if(strategy==4){#Pour KamauRIMDAMAL  Kamau et Al 2024
    KamauRIMDAMAL= Oral_Formulation %>%filter(Article == "KamauRIMDAMAL")
    KamauRIMDAMAL["DAI"] = as.numeric(KamauRIMDAMAL$DAI)
    SimulDrmKamauRIMDAMAL= drm(prop_dead ~ DAI, data = KamauRIMDAMAL,
                               fct = LL.4(fixed = c(NA, 0, 0.6744, NA),
                                          names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
    summary(SimulDrmKamauRIMDAMAL)
    Dmax=0.6744
    tau50=SimulDrmKamauRIMDAMAL$coefficients[2]
    Alpha=SimulDrmKamauRIMDAMAL$coefficients[1]
    
    #plot(EDM_KamauRIMDAMAL)
    DAI = seq(0,30, 1)
    ndat = data.frame(DAI=DAI, prop_dead=NA)
    
    pred_KamauRIMDAMAL=predict(SimulDrmKamauRIMDAMAL,newdata=ndat,interval="confidence",level = 0.3)
    
    pred_KamauRIMDAMAL = data.frame(pred_KamauRIMDAMAL) %>%
      cbind(ndat) %>%
      mutate(Formulation = "KamauRIMDAMAL")
    predprop=pred_KamauRIMDAMAL
    
    sub_PK=KamauRIMDAMAL %>% mutate(Formulation = "KamauRIMDAMAL")
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
    sub_PK_F31 = P4D_PK_kis %>%
      filter(strain == "KIS" & IVM_formulation == CritereIVMformulation)
    sub_PK_F31["DAI"] = as.numeric(sub_PK_F31$DAI)
    #fitting the formula
    SimulDrmLongLasting06_Kis= drm(prop_dead ~ DAI, data = sub_PK_F31,
                                   fct = LL.4(fixed = c(NA, 0, 1, NA),
                                              names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                   type = "binomial")
    summary(SimulDrmLongLasting06_Kis)
    Dmax=1
    tau50=SimulDrmLongLasting06_Kis$coefficients[2]
    Alpha=SimulDrmLongLasting06_Kis$coefficients[1]
    
    #pas très clair ici: juste pour faire de plot visibelement?
    gam_F31 = gam_all %>%
      filter(IVM_formulation == "F31-0.6") %>%
      mutate(DAI = Hour / 24, prop_dead = NA)
    Day_seq=length(unique(gam_F31$DAI))
    pred_F31 = predict(SimulDrmLongLasting06_Kis, newdata = gam_F31, interval = "confidence")
    predprop_F31 = data.frame(pred_F31) %>%
      cbind(gam_F31) %>%
      mutate(Formulation = "F31-0.6")
    predprop_F31 = predprop_F31[, !names(predprop_F31) %in% 
                                  c("IVM_concentration", "IVM_formulation", "Hour", "prop_dead")]
    predprop=predprop_F31
    
    sub_PK_F31 = sub_PK_F31 %>%mutate(DAI = as.numeric(DAI))
    # sub_PK_F31 = sub_PK_F31 %>% mutate(Day = as.numeric(Day))
    sub_PK=sub_PK_F31
    DAI = seq(0,166, 1)
  }
  
  FunList = list("tau50" = tau50, "alpha"=Alpha, "Dmax"=Dmax, "predprop"=predprop,
                 "sub_PK"=sub_PK,"DAI"=DAI)
  return(FunList)
}

RhoFunction=function(Nah,Ntau,Vtau,Vah,IVM_field_dependancy,strategy){
  ParamStrategy=ModelEtau50Alpha(strategy)
  tau50=ParamStrategy$tau50 
  alpha=ParamStrategy$alpha
  
  
  rho=matrix(0,nrow = Nah, ncol = Ntau)
  
  if(IVM_field_dependancy){
    tau_eff=0.5*tau50
    for (a in 1:Nah) {
      for (j in 1:Ntau) {
        rho[a,j] = alpha * (Vtau[j] / tau_eff)^(alpha - 1) /
          (tau_eff * (1 + (Vtau[j] / tau_eff)^alpha)) #10*(Vtau[j]>tau_eff)
      }
    }
  } else {
    tau_eff=tau50
    for (a in 1:Nah) {
      for (j in 1:Ntau) {
        rho[a,j] = alpha * (Vtau[j] / tau_eff)^(alpha - 1) /
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

TheIvmStrategy=function(t_begin_Camp,Number_of_cycle,VectTime_between_cycles,
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

setup_parameters = function() {
  params = list()
  
  params$Gap = 1; params$dt = params$Gap; params$tmax = 400
  params$time = seq(0, params$tmax, by = params$dt); params$Ntime = length(params$time)
  
  params$dtau = 1; params$tau_max = 300; params$Vtau = seq(0, params$tau_max, by = params$dtau)
  params$Ntau = length(params$Vtau); params$dah = 1; params$ah_max = 90
  params$Vah = seq(0, params$ah_max, by = params$dah); params$Nah = length(params$Vah)
  
  params$dsigma = params$dtau; params$sigma_max = 60
  params$Vsigma = seq(0, params$sigma_max, by = params$dsigma)
  params$Nsigma = length(params$Vsigma)
  
  params$mu_m = 0.13; params$IVM_Pregnancy = 0; params$EIP = 8; params$theta = 0.5
  params$Dur_cycle = 7; params$t_begin_Camp = 10; params$q = 0.48  # Proportion des hommes
  
  params$Lambda_h = 13250
  
  # Mortalité humaine mu_h(a)
  {
    mu_h_values = c(66.8, 7.6, 1.7, 0.9, 1.3, 1.9, 2.4, 2.8, 3.6, 4.7, 6.3, 8.9,
                    13.2, 19.8, 31.1, 47.7, 71.3, 110.5, 186.7) / 1000 
    age_labels0 = c(0, 1, seq(5, 90, by = 5))
    params$mu_h = rep(NA, params$Nah)
    for (i in 1:params$Nah) {
      for (IdAgeGroup in 1:(length(age_labels0) - 1)) {
        if (age_labels0[IdAgeGroup] <= params$Vah[i] & params$Vah[i] < age_labels0[IdAgeGroup + 1]) {
          params$mu_h[i] = mu_h_values[IdAgeGroup]
        }
      }
      if (i >= params$Nah) {
        params$mu_h[i] = mu_h_values[length(mu_h_values)]
      }
    }
  }
  
  # Mortalité induite par la maladie delta_h(a)
  {
    delta_h_values = c(1.07 * 10^(-3), 7.02 * 10^(-4), 4.55 * 10^(-4), 5.73 * 10^(-5))
    params$delta_h = rep(NA, params$Nah)
    for (i in 1:params$Nah) {
      if (0 <= params$Vah[i] & params$Vah[i] <= 1) {
        params$delta_h[i] = delta_h_values[1]
      } else if (1 < params$Vah[i] & params$Vah[i] <= 5) {
        params$delta_h[i] = delta_h_values[2]
      } else if (5 < params$Vah[i] & params$Vah[i] <= 15) {
        params$delta_h[i] = delta_h_values[3]
      } else {
        params$delta_h[i] = delta_h_values[4]
      }
    }
  }
  
  {
    # Infectiosité humaine beta_h(a) et recrutement moustiques
    Prev_mos10 = 1
    if (Prev_mos10) {
      params$beta_m = 0.0416 * (params$Vsigma > params$EIP)
      Wedge_m <- 1.33 * 10^6   
      alpha_1 <- 0.122       # 0.071[0.023; 0.175]
      alpha_2 <- 0.17         # 0.302[0.16; 0.475]
    } else {
      params$beta_m = 0.063 * (params$Vsigma > params$EIP)
      Wedge_m <- 1.33 * 10^6    
      alpha_1 <- 0.071        # 0.071[0.023; 0.175]
      alpha_2 <- 0.15 
    }
    
    G_func = function(a) {
      return(22.7 * a * exp(-0.0934 * a))
    }
    params$beta_h = alpha_1 * (G_func(params$Vah)^alpha_2)
    params$bar_beta_h = 0.8 * params$beta_h
  }
  
  # Progression asymptomatique vers clinique nu_h(a)
  {
    p_f = 0.15; nu_0 = 1/11; nu_1 = 1/170; nu_2 = 1/3
    params$nu_hm = rep(NA, params$Nah)
    params$nu_hf = rep(NA, params$Nah)
    for (i in 1:params$Nah) {
      if (params$Vah[i] <= 5) { params$nu_hm[i] = nu_2;params$nu_hf[i] = nu_2
      } else if (params$Vah[i] > 5 && params$Vah[i] <= 15) {
        params$nu_hm[i] = nu_0;params$nu_hf[i] = nu_0
      } else {params$nu_hm[i] = nu_1}
      if (params$Vah[i] > 15 && params$Vah[i] <= 45) {
        params$nu_hf[i] = p_f * nu_0 + (1 - p_f) * nu_1
      } else if (params$Vah[i] > 40) {
        params$nu_hf[i] = nu_1
      }
    }
  }
  
  # Taux de guérison gamma_h(a)
  {
    recovery_ranges = list(c(163, 345), c(555, 714), c(344, 400), c(181, 204), c(82, 92),
                           c(56, 61), c(48, 55))
    age_ranges = list(c(0, 1),c(1, 5),c(5, 8),c(8, 18),c(18, 28),c(28, 43), c(43, params$ah_max))
    params$gamma_h = rep(NA, params$Nah)
    for (i in 1:params$Nah) {
      for (j in 1:length(age_ranges)) {
        age_min = age_ranges[[j]][1]; age_max = age_ranges[[j]][2]
        if (params$Vah[i] >= age_min & params$Vah[i] <= age_max) {
          params$gamma_h[i] = 1 / mean(recovery_ranges[[j]])
        }
      }
    }
    
    # Taux de perte de l'uminite k_h(a)
    params$k_h = rep(5.5*10^(-4), params$Nah)
  }
  
  return(params)
}

# MATRICES D_h, Phi, E_1, E_2, Id_2, Id_4
construct_matrices = function(params, a_idx, mphi_vec, fphi_vec) {
  # Matrices 2x2
  Id_2 = diag(2)
  
  # nu_h(a) diagonal 2x2
  nu_h_mat = diag(c(params$nu_hm[a_idx], params$nu_hf[a_idx]))
  
  # phi(a) diagonal 2x2 
  phi_mat = diag(c(mphi_vec[a_idx], fphi_vec[a_idx]))
  
  # D_h(a) 4x4
  D_h = matrix(0, 4, 4)
  D_h[1:2, 1:2] = params$mu_h[a_idx] * Id_2 + nu_h_mat
  D_h[1:2, 3:4] = matrix(0, 2, 2)
  D_h[3:4, 1:2] = -nu_h_mat
  D_h[3:4, 3:4] = (params$mu_h[a_idx] + params$gamma_h[a_idx] + params$delta_h[a_idx]) * Id_2
  
  # Phi(a) 4x4
  Phi_mat = matrix(0, 4, 4)
  Phi_mat[1:2, 1:2] = phi_mat; Phi_mat[3:4, 3:4] = phi_mat
  
  # E_1 4x4
  E_1 = matrix(0, 4, 4)
  E_1[1:2, 1:2] = Id_2
  
  # E_2 4x4
  E_2 = matrix(0, 4, 4)
  E_2[3:4, 3:4] = Id_2
  
  # Id_4
  Id_4 = diag(4)
  
  return(list(D_h = D_h, Phi = Phi_mat, E_1 = E_1, E_2 = E_2, Id_4 = Id_4, nu_h_mat = nu_h_mat, phi_mat = phi_mat))
}

# FONCTIONS GAMMA 
# Gamma_h(a) - MATRICE 4x4
compute_Gamma_h = function(params, a_idx, mphi_vec, fphi_vec) {
  if (a_idx == 1) {return(diag(4))}
  
  # Intégrale de 0 à a de [Phi(zeta)E_1 + D_h(zeta)]
  integral_mat = matrix(0, 4, 4)
  for (i in 1:(a_idx - 1)) {
    mats = construct_matrices(params, i, mphi_vec, fphi_vec)
    integral_mat = integral_mat + (mats$Phi %*% mats$E_1 + mats$D_h) * params$dah
  }
  Gamma_h = expm(-integral_mat)
  return(Gamma_h)
}

# Gamma_m(sigma) - SCALAIRE
compute_Gamma_m = function(params, sigma_idx, lambda_h_ivm_S0) {
  if (sigma_idx == 1) {return(1)}
  sigma = params$Vsigma[sigma_idx]
  integral_lambda = trapz(params$Vtau, lambda_h_ivm_S0)
  Gamma_m = exp(-(params$mu_m + integral_lambda) * sigma)
  return(Gamma_m)
}

# Gamma_{h,ivm}(a,tau) - MATRICE 4x4
compute_Gamma_h_ivm = function(params, a_idx, tau_idx, mphi_vec, fphi_vec, rho) {
  if (tau_idx == 1) { return(diag(4))}
  
  tau = params$Vtau[tau_idx]; a = params$Vah[a_idx]
  
  # Intégrale de 0 à tau de [D_h(zeta + a - tau) + rho(zeta + a - tau, zeta)Id_4]
  integral_mat = matrix(0, 4, 4)
  for (i in 1:(tau_idx - 1)) {
    zeta = params$Vtau[i]
    age_idx = max(1, min(params$Nah, which.min(abs(params$Vah - (zeta + a - tau)))))
    mats = construct_matrices(params, age_idx, mphi_vec, fphi_vec)
    integral_mat = integral_mat + (mats$D_h + rho[age_idx, i] * mats$Id_4) * params$dtau
  }
  
  Gamma_h_ivm = expm(-integral_mat)
  return(Gamma_h_ivm)
}

# tilde Gamma_{h,ivm}(s,a,tau) - MATRICE 4x4 (cas tau <= a)
compute_tilde_Gamma_h_ivm_case1 = function(params, s_idx, a_idx, tau_idx, mphi_vec, fphi_vec, rho) {
  if (s_idx >= tau_idx) { return(diag(4))}
  tau = params$Vtau[tau_idx]; a = params$Vah[a_idx]
  
  # Intégrale de s à tau de [D_h(zeta + a - tau) + rho(zeta + a - tau,zeta)Id_4]
  integral_mat = matrix(0, 4, 4)
  for (i in s_idx:(tau_idx - 1)) {
    zeta = params$Vtau[i]
    age_idx = max(1, min(params$Nah, which.min(abs(params$Vah - (zeta + a - tau)))))
    mats = construct_matrices(params, age_idx, mphi_vec, fphi_vec)
    integral_mat = integral_mat + (mats$D_h + rho[age_idx, i] * mats$Id_4) * params$dtau
  }
  
  tilde_Gamma = expm(-integral_mat)
  return(tilde_Gamma)
}

# tilde Gamma_{h,ivm}(s,tau,a) - MATRICE 4x4 (cas a <= tau)
compute_tilde_Gamma_h_ivm_case2 = function(params, s_idx, tau_idx, a_idx, mphi_vec, fphi_vec, rho) {
  if (s_idx >= a_idx) {return(diag(4))}
  tau = params$Vtau[tau_idx];a = params$Vah[a_idx]
  
  # Intégrale de s à a 
  integral_mat = matrix(0, 4, 4)
  for (i in s_idx:(a_idx - 1)) {
    age = params$Vah[i]
    tau_temp_idx = max(1, min(params$Ntau, which.min(abs(params$Vtau - tau))))
    mats = construct_matrices(params, i, mphi_vec, fphi_vec)
    integral_mat = integral_mat + (mats$D_h + rho[i, tau_temp_idx] * mats$Id_4) * params$dah
  }
  
  tilde_Gamma = expm(-integral_mat)
  return(tilde_Gamma)
}

# Gamma_{m,ivm}(sigma,tau) - SCALAIRE
compute_Gamma_m_ivm = function(params, sigma_idx, tau_idx, strategy) {
  if (sigma_idx == 1) {return(1)}
  sigma = params$Vsigma[sigma_idx]; tau = params$Vtau[tau_idx]
  
  # mu_{m,ivm}(zeta) 
  ParamIvmFormulation = ModelEtau50Alpha(strategy)
  p_vec_full = P_tau.fun(params$Vtau, tau50 = ParamIvmFormulation$tau50, 
                         alpha = ParamIvmFormulation$alpha, Dmax = ParamIvmFormulation$Dmax)
  mu_m_ivm_vec = params$mu_m - log(1 - p_vec_full)
  
  # Intégrale de 0 à sigma de mu_{m,ivm}(zeta + tau - sigma)
  integral = 0
  for (i in 1:(sigma_idx - 1)) {
    zeta = params$Vsigma[i]
    tau_idx_temp = max(1, min(params$Ntau, which.min(abs(params$Vtau - (zeta + tau - sigma)))))
    integral = integral + mu_m_ivm_vec[tau_idx_temp] * params$dsigma
  }
  
  Gamma_m_ivm = exp(-integral)
  return(Gamma_m_ivm)
}

# tilde Gamma_{m,ivm}(s,sigma,tau) - SCALAIRE (cas sigma <= tau)
compute_tilde_Gamma_m_ivm_case1 = function(params, s_idx, sigma_idx, tau_idx, strategy) {
  if (s_idx >= sigma_idx) {return(1)}
  sigma = params$Vsigma[sigma_idx]; tau = params$Vtau[tau_idx]
  
  ParamIvmFormulation = ModelEtau50Alpha(strategy)
  p_vec_full = P_tau.fun(params$Vtau, tau50 = ParamIvmFormulation$tau50,
                         alpha = ParamIvmFormulation$alpha, Dmax = ParamIvmFormulation$Dmax)
  mu_m_ivm_vec = params$mu_m - log(1 - p_vec_full)
  
  # Intégrale de s à sigma de mu_{m,ivm}(zeta + tau - sigma)
  integral = 0
  for (i in s_idx:(sigma_idx - 1)) {
    zeta = params$Vsigma[i]
    tau_idx_temp = max(1, min(params$Ntau, which.min(abs(params$Vtau - (zeta + tau - sigma)))))
    integral = integral + mu_m_ivm_vec[tau_idx_temp] * params$dsigma
  }
  
  tilde_Gamma = exp(-integral)
  return(tilde_Gamma)
}

# tilde Gamma_{m,ivm}(s,tau,sigma) - SCALAIRE (cas tau <= sigma)
compute_tilde_Gamma_m_ivm_case2 = function(params, s_idx, tau_idx, sigma_idx, strategy) {
  if (s_idx >= tau_idx) {return(1)}
  
  sigma = params$Vsigma[sigma_idx]
  tau = params$Vtau[tau_idx]
  
  ParamIvmFormulation = ModelEtau50Alpha(strategy)
  p_vec_full = P_tau.fun(params$Vtau, tau50 = ParamIvmFormulation$tau50,
                         alpha = ParamIvmFormulation$alpha,
                         Dmax = ParamIvmFormulation$Dmax)
  mu_m_ivm_vec = params$mu_m - log(1 - p_vec_full)
  
  # Intégrale de s à tau
  integral = 0
  for (i in s_idx:(tau_idx - 1)) {
    integral = integral + mu_m_ivm_vec[i] * params$dtau
  }
  
  tilde_Gamma = exp(-integral)
  return(tilde_Gamma)
}

# ÉQUILIBRE SANS MALADIE (DFE) 
compute_DFE = function(params, phi_const, strategy, IVM_field_dependancy = FALSE) {
  Nah = params$Nah; Ntau = params$Ntau; Vah = params$Vah
  Vtau = params$Vtau; dah = params$dah; dtau = params$dtau
  
  # rho dépend de la stratégie
  rho_result = RhoFunction(Nah,Ntau,Vtau,Vah,IVM_field_dependancy,strategy)
  rho = rho_result$rho
  
  #  mphi_vec et fphi_vec 
  p_f = 0.15
  mphi_vec = rep(NA, Nah); fphi_vec = rep(NA, Nah)
  
  for (i in 1:Nah) {
    if (Vah[i] <= 5) {
      mphi_vec[i] = 0
      fphi_vec[i] = 0
    } else if (Vah[i] >= 15 & Vah[i] <= 45) {
      if (params$IVM_Pregnancy == 1) {
        mphi_vec[i] = params$q * phi_const
        fphi_vec[i] = (1 - params$q) * (1 - p_f) * phi_const
      } else {
        mphi_vec[i] = params$q * phi_const
        fphi_vec[i] = 0
      }
    } else {
      mphi_vec[i] = params$q * phi_const
      fphi_vec[i] = (1 - params$q) * phi_const
    }
  }  
  
  mSh = rep(0, Nah); fSh = rep(0, Nah)
  
  mSh[1] = params$q * params$Lambda_h
  fSh[1] = (1 - params$q) * params$Lambda_h
  # Itération pour DFE
  tol = 1e-6; max_iter = 1000
  
  for (iter in 1:max_iter) {
    mSh_old = mSh; fSh_old = fSh
    for (i in 2:Nah) {
      a = Vah[i]
      
      # MALE # Terme exp(-int_0^a (mu_h + phi^male))
      integrand_surv = params$mu_h[1:(i-1)] + mphi_vec[1:(i-1)]
      integrand_surv = integrand_surv * dah
      surv_0_a_male = exp(-sum(integrand_surv))
      
      # Intégrale rho
      integral_rho_male = 0
      for (s_idx in 1:(i-1)) {
        s = Vah[s_idx]
        
        # exp(-int_s^a (mu_h + phi^male))
        integrand_s_a = params$mu_h[s_idx:(i-1)] + mphi_vec[s_idx:(i-1)]
        integrand_s_a = integrand_s_a * dah
        exp_s_a = exp(-sum(integrand_s_a))
        
        # Intégrale sur tau de 0 à s
        inner_integral = 0
        for (tau_idx in 1:s_idx) {
          tau = Vtau[tau_idx]
          if (tau <= s) {
            s_minus_tau_idx = max(1, which.min(abs(Vah - (s - tau))))
            
            # exp(-int_0^tau (mu_h(zeta + s - tau) + rho(zeta)))
            integrand_tau = rep(0, tau_idx)
            for (k in 1:tau_idx) {
              zeta = Vtau[k]
              age_temp_idx = max(1, min(Nah, which.min(abs(Vah - (zeta + s - tau)))))
              integrand_tau[k] = (params$mu_h[age_temp_idx] + rho[age_temp_idx, k]) * dtau
            }
            exp_tau = exp(-sum(integrand_tau))
            
            inner_integral = inner_integral + rho[s_idx, tau_idx] * mphi_vec[s_minus_tau_idx] * 
              mSh_old[s_minus_tau_idx] * exp_tau * dtau
          }
        }
        
        integral_rho_male = integral_rho_male + exp_s_a * inner_integral * dah
      }
      
      mSh[i] = surv_0_a_male * mSh[1] + integral_rho_male
      
      # FEMALE (même structure)
      integrand_surv_f = params$mu_h[1:(i-1)] + fphi_vec[1:(i-1)]
      integrand_surv_f = integrand_surv_f * dah
      surv_0_a_female = exp(-sum(integrand_surv_f))
      
      integral_rho_female = 0
      for (s_idx in 1:(i-1)) {
        s = Vah[s_idx]
        
        integrand_s_a_f = params$mu_h[s_idx:(i-1)] + fphi_vec[s_idx:(i-1)]
        integrand_s_a_f = integrand_s_a_f * dah
        exp_s_a_f = exp(-sum(integrand_s_a_f))
        
        inner_integral_f = 0
        for (tau_idx in 1:s_idx) {
          tau = Vtau[tau_idx]
          if (tau <= s) {
            s_minus_tau_idx = max(1, which.min(abs(Vah - (s - tau))))
            
            integrand_tau_f = rep(0, tau_idx)
            for (k in 1:tau_idx) {
              zeta = Vtau[k]
              age_temp_idx = max(1, min(Nah, which.min(abs(Vah - (zeta + s - tau)))))
              integrand_tau_f[k] = (params$mu_h[age_temp_idx] + rho[age_temp_idx, k]) * dtau
            }
            exp_tau_f = exp(-sum(integrand_tau_f))
            
            inner_integral_f = inner_integral_f + rho[s_idx, tau_idx] * fphi_vec[s_minus_tau_idx] *
              fSh_old[s_minus_tau_idx] * exp_tau_f * dtau
          }
        }
        
        integral_rho_female = integral_rho_female + exp_s_a_f * inner_integral_f * dah
      }
      
      fSh[i] = surv_0_a_female * fSh[1] + integral_rho_female
    }
    
    # Convergence check
    if (max(abs(mSh - mSh_old)) < tol && max(abs(fSh - fSh_old)) < tol) {
      break
    }
  }
  
  # S_{h,ivm}^{g0}(a,tau) = 1_{[0,a]} phi(a-tau) S^{g0}_h(a-tau) exp(-int_0^tau (mu_h(zeta+a-tau) + rho(zeta)))
  mSh_ivm = matrix(0, nrow = Nah, ncol = Ntau)
  fSh_ivm = matrix(0, nrow = Nah, ncol = Ntau)
  
  for (i in 1:Nah) {
    for (j in 1:Ntau) {
      a = Vah[i]
      tau = Vtau[j]
      
      if (tau <= a) {  # 1_{[0,a]}(tau)
        a_minus_tau_idx = max(1, which.min(abs(Vah - (a - tau))))
        
        # exp(-int_0^tau (mu_h(zeta+a-tau) + rho(zeta)))
        integrand = rep(0, j)
        for (k in 1:j) {
          zeta = Vtau[k]
          age_idx_temp = max(1, min(Nah, which.min(abs(Vah - (zeta + a - tau)))))
          integrand[k] = (params$mu_h[age_idx_temp] + rho[age_idx_temp, k]) * dtau
        }
        exp_factor = exp(-sum(integrand))
        
        mSh_ivm[i, j] = mphi_vec[a_minus_tau_idx] * mSh[a_minus_tau_idx] * exp_factor
        fSh_ivm[i, j] = fphi_vec[a_minus_tau_idx] * fSh[a_minus_tau_idx] * exp_factor
      }
    }
  }
  
  # Population totale N_h^0
  N_h_0 = sum(mSh) * dah + sum(fSh) * dah +
    sum(mSh_ivm) * dah * dtau + sum(fSh_ivm) * dah * dtau
  
  # v_h^0 et v_{h,ivm}^0 sont des VECTEURS
  # v_h^0(a) = (S_h^male(a), S_h^female(a), 0, 0)^T
  # v_{h,ivm}^0(a,tau) = (S_{h,ivm}^male(a,tau), S_{h,ivm}^female(a,tau), 0, 0)^T
  
  # lambda_{h,ivm}^{S0}(tau)
  lambda_h_ivm_S0 = rep(0, Ntau)
  e_vec = c(1, 1, 1, 1)
  
  for (j in 1:Ntau) {
    # e . v_{h,ivm}^0(a,tau)
    numerator_integrand = rep(0, Nah)
    for (i in 1:Nah) {
      v_h_ivm_0 = c(mSh_ivm[i, j], fSh_ivm[i, j], 0, 0)
      numerator_integrand[i] = sum(e_vec * v_h_ivm_0)
    }
    numerator = trapz(Vah, numerator_integrand)
    
    # Dénominateur avec intégrale sur tau aussi
    denom_integrand = rep(0, Nah)
    for (i in 1:Nah) {
      v_h_0 = c(mSh[i], fSh[i], 0, 0)
      e_vh0 = sum(e_vec * v_h_0)
      
      e_vhivm0_sum = 0
      for (k in 1:Ntau) {
        v_h_ivm_0_k = c(mSh_ivm[i, k], fSh_ivm[i, k], 0, 0)
        e_vhivm0_sum = e_vhivm0_sum + sum(e_vec * v_h_ivm_0_k) * dtau
      }
      
      denom_integrand[i] = e_vh0 + e_vhivm0_sum
    }
    denominator = trapz(Vah, denom_integrand)
    
    lambda_h_ivm_S0[j] = params$theta * numerator / denominator
  }
  
  # S_m^0
  integral_lambda = trapz(Vtau, lambda_h_ivm_S0)
  S_m_0 = params$Wedge_m / (params$mu_m + integral_lambda)
  
  # S_{m,ivm}(tau) avec mu_{m,ivm} dépendant de la stratégie
  ParamIvmFormulation = ModelEtau50Alpha(strategy)
  p_vec = P_tau.fun(Vtau, tau50 = ParamIvmFormulation$tau50,
                    alpha = ParamIvmFormulation$alpha, Dmax = ParamIvmFormulation$Dmax)
  mu_m_ivm = params$mu_m - log(1 - p_vec)
  
  S_m_ivm = rep(0, Ntau)
  S_m_ivm[1] = 0
  
  if (is.null(S_m_0) || is.na(S_m_0) || length(S_m_0) == 0) {
    warning("S_m_0 est invalide")
    return(NULL)
  }
  
  for (j in 2:Ntau) {
    tau = Vtau[j]
    
    integrand = rep(0, j)
    for (k in 1:j) {
      s = Vtau[k]
      
      # exp(-int_s^tau mu_{m,ivm}(zeta))
      exp_integral = 0
      for (l in k:j) {
        exp_integral = exp_integral + mu_m_ivm[l] * dtau
      }
      
      integrand[k] = lambda_h_ivm_S0[k] * exp(-exp_integral)
    }
    
    S_m_ivm[j] = S_m_0 * trapz(Vtau[1:j], integrand)
  }
  
  return(list( mSh = mSh, fSh = fSh, mSh_ivm = mSh_ivm,fSh_ivm = fSh_ivm, N_h_0 = N_h_0, 
               lambda_h_ivm_S0 = lambda_h_ivm_S0, S_m_0 = S_m_0, S_m_ivm = S_m_ivm, rho = rho, 
               mphi_vec = mphi_vec, fphi_vec = fphi_vec, mu_m_ivm = mu_m_ivm))
}

# OPÉRATEUR G NGO
compute_NGO = function(params, DFE, strategy) {
  
  Nah = params$Nah; Ntau = params$Ntau; Nsigma = params$Nsigma; Vah = params$Vah
  Vtau = params$Vtau; Vsigma = params$Vsigma; dah = params$dah; dtau = params$dtau
  dsigma = params$dsigma
  
  rho = DFE$rho; mSh = DFE$mSh; fSh = DFE$fSh; mSh_ivm = DFE$mSh_ivm
  fSh_ivm = DFE$fSh_ivm; N_h_0 = DFE$N_h_0; S_m_0 = DFE$S_m_0; lambda_h_ivm_S0 = DFE$lambda_h_ivm_S0
  mphi_vec = DFE$mphi_vec; fphi_vec = DFE$fphi_vec
  
  #  psi_h et psi_h_ivm sont des VECTEURS R^4
  
  dim_psi_h = Nah * 4 ; dim_psi_h_ivm = Nah * Ntau * 4 
  dim_psi_m = Nsigma; dim_psi_m_ivm = Nsigma * Ntau   
  dim_total = dim_psi_h + dim_psi_h_ivm + dim_psi_m + dim_psi_m_ivm
  G_matrix = matrix(0, nrow = dim_total, ncol = dim_total)
  
  # BLOC K : psi_{h,ivm} → psi_h 
  for (a_idx in 1:Nah) {
    a = Vah[a_idx]
    for (comp in 1:4) {
      idx_row_h = (a_idx - 1) * 4 + comp
      
      for (tau_idx in 1:Ntau) {
        tau = Vtau[tau_idx]; rho_tau = rho[a_idx, tau_idx]
        
        # CAS 1: tau <= a
        if (tau <= a) {
          # Gamma_{h,ivm}(a,tau) : matrice 4×4
          Gamma_h_ivm_mat = compute_Gamma_h_ivm(params, a_idx, tau_idx, mphi_vec, fphi_vec, rho)
          
          # Terme Pi_{2,tau} psi_{h,ivm} : trace à tau=0
          for (a2_idx in 1:Nah) {
            if (a2_idx == a_idx) {  # Condition de trace
              for (comp2 in 1:4) {
                idx_col_hivm = dim_psi_h + (a2_idx - 1) * Ntau * 4 + (1 - 1) * 4 + comp2
                contribution = rho_tau * Gamma_h_ivm_mat[comp, comp2] * dtau
                G_matrix[idx_row_h, idx_col_hivm] = G_matrix[idx_row_h, idx_col_hivm] + contribution
              }
            }
          }
          
          # Terme avec double intégrale
          for (s_idx in 1:tau_idx) {
            s = Vtau[s_idx]
            tilde_Gamma_mat = compute_tilde_Gamma_h_ivm_case1(params, s_idx, a_idx, tau_idx,
                                                              mphi_vec, fphi_vec, rho)
            
            s_plus_a_minus_tau = s + a - tau
            if (s_plus_a_minus_tau >= 0 && s_plus_a_minus_tau <= params$ah_max) {
              a_s_idx = max(1, min(Nah, which.min(abs(Vah - s_plus_a_minus_tau))))
              
              for (comp2 in 1:4) {
                idx_col_hivm = dim_psi_h + (a_s_idx - 1) * Ntau * 4 + (s_idx - 1) * 4 + comp2
                contribution = rho_tau * tilde_Gamma_mat[comp, comp2] * dtau * dtau
                G_matrix[idx_row_h, idx_col_hivm] = G_matrix[idx_row_h, idx_col_hivm] + contribution
              }
            }
          }
        }
        
        # CAS 2: a <= tau
        if (a <= tau) {
          for (s_idx in 1:min(a_idx, Ntau)) {
            s = Vah[s_idx]
            tilde_Gamma_mat = compute_tilde_Gamma_h_ivm_case2(params, s_idx, tau_idx, a_idx,
                                                              mphi_vec, fphi_vec, rho)
            
            s_plus_tau_minus_a = s + tau - a
            if (s_plus_tau_minus_a >= 0 && s_plus_tau_minus_a <= params$tau_max) {
              tau_s_idx = max(1, min(Ntau, which.min(abs(Vtau - s_plus_tau_minus_a))))
              
              for (comp2 in 1:4) {
                idx_col_hivm = dim_psi_h + (s_idx - 1) * Ntau * 4 + (tau_s_idx - 1) * 4 + comp2
                contribution = rho_tau * tilde_Gamma_mat[comp, comp2] * dtau * dah
                G_matrix[idx_row_h, idx_col_hivm] =  G_matrix[idx_row_h, idx_col_hivm] + contribution
              }
            }
          }
        }
      }
    }
  }
  
  # BLOC M : psi_m, psi_{m,ivm} → psi_h et psi_{h,ivm}
  
  # M1 : psi_m, psi_{m,ivm} → psi_h
  
  for (a_idx in 1:Nah) {
    # v_h^0(a) : vecteur R^4
    v_h_0_vec = c(mSh[a_idx], fSh[a_idx], 0, 0)
    
    # E_1 matrix
    mats = construct_matrices(params, a_idx, mphi_vec, fphi_vec)
    E_1_mat = mats$E_1
    
    # E_1 * v_h^0(a) : vecteur R^4
    E1_vh0_vec = as.vector(E_1_mat %*% v_h_0_vec)
    
    for (comp in 1:4) {
      idx_row_h = (a_idx - 1) * 4 + comp
      
      # Contribution de psi_m
      for (sigma_idx in 1:Nsigma) {
        idx_col_m = dim_psi_h + dim_psi_h_ivm + sigma_idx
        
        Gamma_m_sigma = compute_Gamma_m(params, sigma_idx, lambda_h_ivm_S0)
        
        lambda_m_contrib = params$theta / N_h_0 * params$beta_m[sigma_idx] * 
          Gamma_m_sigma * dsigma
        
        contribution = lambda_m_contrib * E1_vh0_vec[comp]
        
        G_matrix[idx_row_h, idx_col_m] = G_matrix[idx_row_h, idx_col_m] + contribution
      }
      
      # Contribution de psi_{m,ivm} 
      for (sigma_idx in 1:Nsigma) {
        for (tau_idx in 1:Ntau) {
          idx_col_mivm = dim_psi_h + dim_psi_h_ivm + dim_psi_m + 
            (sigma_idx - 1) * Ntau + tau_idx
          
          sigma = Vsigma[sigma_idx]
          tau = Vtau[tau_idx]
          
          lambda_m_ivm_contrib = 0
          
          # CAS 1: sigma <= tau
          if (sigma <= tau) {
            Gamma_m_ivm_scalar = compute_Gamma_m_ivm(params, sigma_idx, tau_idx, strategy)
            lambda_m_ivm_contrib = lambda_m_ivm_contrib + params$theta / N_h_0 * 
              params$beta_m[sigma_idx] *  Gamma_m_ivm_scalar * dtau * dsigma
            
            # Intégrale sur s
            for (s_idx in 1:sigma_idx) {
              tilde_Gamma_m_ivm_scalar = compute_tilde_Gamma_m_ivm_case1(
                params, s_idx, sigma_idx, tau_idx, strategy)
              
              lambda_m_ivm_contrib = lambda_m_ivm_contrib +  params$theta / N_h_0 * 
                params$beta_m[sigma_idx] * tilde_Gamma_m_ivm_scalar * dsigma * dtau * dsigma
            }
          }
          
          # CAS 2: tau <= sigma
          if (tau <= sigma) {
            for (s_idx in 1:tau_idx) {
              tilde_Gamma_m_ivm_scalar = compute_tilde_Gamma_m_ivm_case2(
                params, s_idx, tau_idx, sigma_idx, strategy)
              
              lambda_m_ivm_contrib = lambda_m_ivm_contrib + params$theta / N_h_0 * 
                params$beta_m[sigma_idx] * tilde_Gamma_m_ivm_scalar * dtau * dtau * dsigma
            }
          }
          
          contribution = lambda_m_ivm_contrib * E1_vh0_vec[comp]
          G_matrix[idx_row_h, idx_col_mivm] = G_matrix[idx_row_h, idx_col_mivm] + contribution
        }
      }
    }
  }
  
  # M2 : psi_m, psi_{m,ivm} → psi_{h,ivm}
  
  for (a_idx in 1:Nah) {
    for (tau_idx in 1:Ntau) {
      # v_{h,ivm}^0(a,tau) : vecteur R^4
      v_h_ivm_0_vec = c(mSh_ivm[a_idx, tau_idx], 
                        fSh_ivm[a_idx, tau_idx], 0, 0)
      
      mats = construct_matrices(params, a_idx, mphi_vec, fphi_vec)
      E_1_mat = mats$E_1
      E1_vhivm0_vec = as.vector(E_1_mat %*% v_h_ivm_0_vec)
      
      # Boucler sur les 4 composantes du vecteur
      for (comp in 1:4) {
        idx_row_hivm = dim_psi_h + (a_idx - 1) * Ntau * 4 + (tau_idx - 1) * 4 + comp
        
        # Contribution de psi_m (même structure que M1)
        for (sigma_idx in 1:Nsigma) {
          idx_col_m = dim_psi_h + dim_psi_h_ivm + sigma_idx
          
          Gamma_m_scalar = compute_Gamma_m(params, sigma_idx, lambda_h_ivm_S0)
          lambda_m_contrib = params$theta / N_h_0 * params$beta_m[sigma_idx] *
            Gamma_m_scalar * dsigma
          
          contribution = lambda_m_contrib * E1_vhivm0_vec[comp]
          G_matrix[idx_row_hivm, idx_col_m] = G_matrix[idx_row_hivm, idx_col_m] + contribution
        }
        
        # Contribution de psi_{m,ivm} (même structure)
        for (sigma_idx in 1:Nsigma) {
          for (tau_idx2 in 1:Ntau) {
            idx_col_mivm = dim_psi_h + dim_psi_h_ivm + dim_psi_m + 
              (sigma_idx - 1) * Ntau + tau_idx2
            
            sigma = Vsigma[sigma_idx]
            tau2 = Vtau[tau_idx2]
            
            lambda_m_ivm_contrib = 0
            
            if (sigma <= tau2) {
              Gamma_m_ivm_scalar = compute_Gamma_m_ivm(params, sigma_idx, tau_idx2, strategy)
              lambda_m_ivm_contrib = lambda_m_ivm_contrib +params$theta / N_h_0 * 
                params$beta_m[sigma_idx] *  Gamma_m_ivm_scalar * dtau * dsigma
              
              for (s_idx in 1:sigma_idx) {
                tilde_Gamma_m_ivm_scalar = compute_tilde_Gamma_m_ivm_case1(
                  params, s_idx, sigma_idx, tau_idx2, strategy)
                lambda_m_ivm_contrib = lambda_m_ivm_contrib +params$theta / N_h_0 * 
                  params$beta_m[sigma_idx] * tilde_Gamma_m_ivm_scalar * dsigma * dtau * dsigma
              }
            }
            
            if (tau2 <= sigma) {
              for (s_idx in 1:tau_idx2) {
                tilde_Gamma_m_ivm_scalar = compute_tilde_Gamma_m_ivm_case2(
                  params, s_idx, tau_idx2, sigma_idx, strategy)
                lambda_m_ivm_contrib = lambda_m_ivm_contrib +params$theta / N_h_0 * 
                  params$beta_m[sigma_idx]*tilde_Gamma_m_ivm_scalar * dtau * dtau * dsigma
              }
            }
            
            contribution = lambda_m_ivm_contrib * E1_vhivm0_vec[comp]
            G_matrix[idx_row_hivm,idx_col_mivm]=G_matrix[idx_row_hivm,idx_col_mivm] + contribution
          }
        }
      }
    }
  }
  
  
  # COUPLAGE bar_psi_m (psi_h → psi_h) nouveaux infected
  
  # bar_psi_m crée un couplage : psi_h → lambda_h^0 → bar_psi_m → lambda_m^0 → psi_h
  
  for (a_idx in 1:Nah) {  # Cible
    v_h_0_vec = c(mSh[a_idx], fSh[a_idx], 0, 0)
    mats = construct_matrices(params, a_idx, mphi_vec, fphi_vec)
    E_1_mat = mats$E_1
    E1_vh0_vec = as.vector(E_1_mat %*% v_h_0_vec)
    
    for (comp in 1:4) {
      idx_row_h = (a_idx - 1) * 4 + comp
      
      for (a2_idx in 1:Nah) {  # Source
        for (comp2 in 1:4) {
          idx_col_h = (a2_idx - 1) * 4 + comp2
          
          # Contribution de psi_h(a2, comp2) à lambda_h^0
          a2 = Vah[a2_idx]
          beta_h_T_a2 = c(params$bar_beta_h[a2_idx], params$bar_beta_h[a2_idx],
                          params$beta_h[a2_idx], params$beta_h[a2_idx])
          
          Gamma_h_a2 = compute_Gamma_h(params, a2_idx, mphi_vec, fphi_vec)
          
          # d(lambda_h^0) / d(psi_h[a2, comp2])
          # Approximation : contribution directe de psi_h(a2)
          dlambda_h0_dpsi_h = params$theta / N_h_0 * beta_h_T_a2[comp2] * 
            Gamma_h_a2[comp2, comp2] * dah
          
          # d(bar_psi_m) / d(lambda_h^0) = S_m^0
          dbar_psi_m_dlambda_h0 = S_m_0
          
          # d(lambda_m^0) / d(bar_psi_m)
          # = int beta_m(sigma) * Gamma_m(sigma) dsigma
          dlambda_m0_dbar_psi_m = 0
          for (sigma_idx in 1:Nsigma) {
            Gamma_m_sigma = compute_Gamma_m(params, sigma_idx, lambda_h_ivm_S0)
            dlambda_m0_dbar_psi_m = dlambda_m0_dbar_psi_m + params$theta / N_h_0 * 
              params$beta_m[sigma_idx] * Gamma_m_sigma * dsigma
          }
          
          # Contribution totale : chaîne de dérivées
          contribution = dlambda_m0_dbar_psi_m * dbar_psi_m_dlambda_h0 * dlambda_h0_dpsi_h * 
            E1_vh0_vec[comp]
          
          G_matrix[idx_row_h, idx_col_h] = G_matrix[idx_row_h, idx_col_h] + contribution
        }
      }
    }
  }
  
  # Même couplage pour psi_{h,ivm}
  for (a_idx in 1:Nah) {
    for (tau_idx in 1:Ntau) {
      v_h_ivm_0_vec = c(mSh_ivm[a_idx, tau_idx], 
                        fSh_ivm[a_idx, tau_idx], 0, 0)
      mats = construct_matrices(params, a_idx, mphi_vec, fphi_vec)
      E_1_mat = mats$E_1
      E1_vhivm0_vec = as.vector(E_1_mat %*% v_h_ivm_0_vec)
      
      for (comp in 1:4) {
        idx_row_hivm = dim_psi_h + (a_idx - 1) * Ntau * 4 + (tau_idx - 1) * 4 + comp
        
        for (a2_idx in 1:Nah) {
          for (comp2 in 1:4) {
            idx_col_h = (a2_idx - 1) * 4 + comp2
            
            beta_h_T_a2 = c(params$bar_beta_h[a2_idx], params$bar_beta_h[a2_idx],
                            params$beta_h[a2_idx], params$beta_h[a2_idx])
            Gamma_h_a2 = compute_Gamma_h(params, a2_idx, mphi_vec, fphi_vec)
            
            dlambda_h0_dpsi_h = params$theta / N_h_0 * beta_h_T_a2[comp2] * 
              Gamma_h_a2[comp2, comp2] * dah
            dbar_psi_m_dlambda_h0 = S_m_0
            
            dlambda_m0_dbar_psi_m = 0
            for (sigma_idx in 1:Nsigma) {
              Gamma_m_sigma = compute_Gamma_m(params, sigma_idx, lambda_h_ivm_S0)
              dlambda_m0_dbar_psi_m = dlambda_m0_dbar_psi_m +  params$theta / N_h_0 * 
                params$beta_m[sigma_idx] * Gamma_m_sigma * dsigma
            }
            
            contribution = dlambda_m0_dbar_psi_m * dbar_psi_m_dlambda_h0 * 
              dlambda_h0_dpsi_h * E1_vhivm0_vec[comp]
            
            G_matrix[idx_row_hivm, idx_col_h] = G_matrix[idx_row_hivm, idx_col_h] + contribution
          }
        }
      }
    }
  }
  
  # BLOC N : psi_m → psi_{m,ivm}
  
  for (sigma_idx in 1:Nsigma) {
    for (tau_idx in 1:Ntau) {
      idx_row_mivm = dim_psi_h + dim_psi_h_ivm + dim_psi_m + (sigma_idx - 1) * Ntau + tau_idx
      
      # Contribution de psi_m
      for (jsigma in 1:Nsigma) {
        idx_col_m = dim_psi_h + dim_psi_h_ivm + jsigma
        
        Gamma_m_scalar = compute_Gamma_m(params, jsigma, lambda_h_ivm_S0)
        
        contribution = lambda_h_ivm_S0[tau_idx] * Gamma_m_scalar * dsigma * dtau
        G_matrix[idx_row_mivm, idx_col_m] =  G_matrix[idx_row_mivm, idx_col_m] + contribution
      }
    }
  }
  
  if (sum(G_matrix != 0) > 0) {
    cat("  • Range (non-zero):", sprintf("[%.2e, %.2e]", 
                                         min(G_matrix[G_matrix != 0]), max(G_matrix[G_matrix != 0])), "\n")
  }
  return(G_matrix)
}

# RAYON SPECTRAL R_0
compute_R0 = function(G_matrix) {
  eigenvalues = eigen(G_matrix, only.values = TRUE)$values
  R0 = max(abs(eigenvalues))
  return(R0)
}

# ANALYSES ET PLOTS
analyze_constant_IVM = function(params, strategies = 0:5, PropIVM_range = seq(0, 0.9, by = 0.1),
                                IVM_field_dependancy = FALSE) {
  results = data.frame()
  
  for (strat in strategies) {
    cat("\n=== Stratégie", strat, "===\n")
    R0_values = c()
    
    for (prop in PropIVM_range) {
      cat("PropIVM =", prop, "...\n")
      
      DFE = compute_DFE(params, prop, strat, IVM_field_dependancy)
      G_matrix = compute_NGO(params, DFE, strat)
      R0 = compute_R0(G_matrix)
      R0_values = c(R0_values, R0)
      
      cat("R_0 =", R0, "\n")
    }
    
    results = rbind(results, data.frame(
      Strategy = strat,
      PropIVM = PropIVM_range,
      R0 = R0_values
    ))
  }
  
  return(results)
}

analyze_time_dependent_IVM = function(params, VectPropIVM = c(0.5, 0.7, 0.9), 
                                      IVM_field_dependancy = FALSE) {
  results_list = list()
  idx = 1
  
  # Stratégies 0, 1, 5: nc varie de 1 à 4 # Stratégie 2: nc fixe à 3 # Stratégies 3, 4: nc fixe à 4
  
  for (PropIVM in VectPropIVM) {
    
    for (nc in 1:4) {
      strategies_to_run = c()
      
      # Stratégies 0, 1, 5: toujours inclure pour tous les nc
      strategies_to_run = c(strategies_to_run, 0, 1, 5)
      
      # Stratégie 2 (BOHEMIA): seulement quand nc = 3
      if (nc == 3) {
        strategies_to_run = c(strategies_to_run, 2)
      }
      
      # Stratégies 3, 4 (RIMDAMAL): seulement quand nc = 4
      if (nc == 4) {
        strategies_to_run = c(strategies_to_run, 3, 4)
      }
      
      for (strat in strategies_to_run) {
        # Temps entre cycles selon stratégie
        if (strat %in% c(0, 1, 5)) {
          time_between = 60
        } else {
          time_between = 30
        }
        
        VectTime_between_cycles = rep(time_between, nc - 1)
        
        # Calculer phi(t) pour cette configuration
        phi_result = TheIvmStrategy(t_begin_Camp = params$t_begin_Camp,Number_of_cycle = nc,
                                    VectTime_between_cycles = VectTime_between_cycles,Dur_cycle = params$Dur_cycle,
                                    Gap = params$Gap, time = params$time, Ntime = params$Ntime)
        
        phi_t = phi_result$PhiIvm
        time_sampled = seq(1, params$Ntime, by = 1)
        R0_effective = c()
        
        # Calcul R0 effectif pour chaque temps
        for (t_idx in time_sampled) {
          phi_current = PropIVM * phi_t[t_idx] / params$Dur_cycle
          
          DFE = compute_DFE(params, phi_current, strat, IVM_field_dependancy)
          G_matrix = compute_NGO(params, DFE, strat)
          R0 = compute_R0(G_matrix)
          R0_effective = c(R0_effective, R0)
        }
        
        # Stocker les résultats
        results_list[[idx]] = data.frame(PropIVM = PropIVM, Strategy = strat,
                                         nc = nc, time = params$time[time_sampled], R0_effective = R0_effective,
                                         phi = phi_t[time_sampled] )
        idx = idx + 1
      }
    }
  }
  return(results_list)
}

plot_R0_constant_IVM = function(results) {
  LC = 2.5
  LEGEND = c("mdc-STM-001−0.6-vk5", "mdc-STM-001−0.6-kis", "mdc-STM-001−1.0",
             expression("RIMDAMAL II"[S]), expression("RIMDAMAL II"[K]), "BOHEMIA")
  colors = c("#ff3355", "#8B4513", "#037153", "#205072", "#8B008B", "#D98E04")
  linetypes = c(1, 4, 3, 4, 4, 2)
  
  plot(NULL, xlim = c(0, max(results$PropIVM)),
       ylim = c(min(results$R0), max(results$R0)),
       xlab = "IVM exposure rate", ylab = expression(R[0]),
       main = expression("Impact of IVM formulation on" ~ R[0]))
  
  for (i in 0:5) {
    subset_data = results[results$Strategy == i, ]
    lines(subset_data$PropIVM, subset_data$R0,
          col = colors[i + 1], lty = linetypes[i + 1], lwd = LC)
  }
  
  legend("topright", legend = LEGEND, col = colors, lty = linetypes,
         lwd = LC, bty = "n", cex = 0.8)
  abline(h = 1, lty = 2, col = "red", lwd = 2)
}

plot_R0_time_dependent = function(results_list, params, FigName = "R0_time_dependent.pdf") {
  colors = c("#ff3355", "#8B4513", "#037153", "#205072", "#8B008B", "#D98E04")
  linetypes = c(1, 4, 3, 4, 4, 2)
  LEGEND = c("mdc-STM-001−0.6-vk5", "mdc-STM-001−0.6-kis", "mdc-STM-001−1.0",
             "RIMDAMAL II (Smit)", "RIMDAMAL II (Kamau)", "BOHEMIA")
  LC = 2.5
  VectPropIVM = c(0.5, 0.7, 0.9)
  
  t_begin_Camp = params$t_begin_Camp
  tmax = max(params$time)
  tmaxmonths = (tmax - t_begin_Camp) / 30
  all_R0 = sapply(results_list, function(x) x$R0_effective)
  Max_y = max(all_R0, na.rm = TRUE) * 1.05
  
  pdf(FigName, width = 15, height = 8.7)
  par(oma = c(6, 1, 2, 0.1), mar = c(3, 3.5, 2, 0.1))
  par(mfrow = c(length(VectPropIVM), 4))
  
  GrNumber = 0
  for (PropIVM in VectPropIVM) {
    for (nc in 1:4) {
      GrNumber = GrNumber + 1
      
      plot(-1, 1, type = "l", xlab = "", xlim = c(t_begin_Camp, tmax), 
           ylim = c(0, Max_y), ylab = "", cex.lab = 1.2, yaxt = "n", xaxt = "n")
      axis(1, at = seq(t_begin_Camp, tmax, by = 60), labels = seq(0, tmaxmonths, by = 2), las = 1)
      axis(2, las = 1)
      if (GrNumber <= 4) {
        mtext(paste("n_c =", nc), side = 3, line = 0.5, cex = 0.9, font = 2)
      }
      
      if (nc == 1) {
        mtext(paste("PropIVM =", PropIVM), side = 2, line = 3, cex = 0.85, font = 2)
        mtext(expression(R[0]^{eff}), side = 2, adj = 0.5, cex = 0.85, line = 2)
      }
      if (PropIVM == 0.9) {
        mtext("Time (months)", side = 1, adj = 0.5, cex = 0.85, line = 2)
      }
      
      for (strat in c(0, 1, 5)) {
        subset_data = Filter(function(x) {
          x$PropIVM[1] == PropIVM && x$Strategy[1] == strat && x$nc[1] == nc
        }, results_list)
        if (length(subset_data) > 0) {
          data = subset_data[[1]]
          lines(data$time, data$R0_effective, lwd = LC, lty = linetypes[strat + 1], col = colors[strat + 1])
        }
      }
      
      # Stratégie 2 (BOHEMIA) 
      subset_data = Filter(function(x) {
        x$PropIVM[1] == PropIVM && x$Strategy[1] == 2 && x$nc[1] == 3
      }, results_list)
      if (length(subset_data) > 0) {
        data = subset_data[[1]]
        lines(data$time, data$R0_effective, 
              lwd = LC, lty = linetypes[2 + 1], col = colors[2 + 1])
      }
      
      
      # Stratégies 3, 4 (RIMDAMAL)
      for (strat in c(3, 4)) {
        subset_data = Filter(function(x) {
          x$PropIVM[1] == PropIVM && x$Strategy[1] == strat && x$nc[1] == 4
        }, results_list)
        
        if (length(subset_data) > 0) {
          data = subset_data[[1]]
          lines(data$time, data$R0_effective, 
                lwd = LC, lty = linetypes[strat + 1], col = colors[strat + 1])
        }
      }
      
      abline(h = 1, lty = 2, col = "black", lwd = 1)
    }
  }
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = LEGEND, xpd = NA, horiz = FALSE, bty = "n", lty = linetypes, 
         lwd = LC, col = colors, ncol = 3, cex = 1.1)
  dev.off()
}

main_analysis = function() {
  params = setup_parameters()
  
  results_constant = analyze_constant_IVM(params, strategies = 0:5)
  
  pdf("R0_vs_PropIVM_strategies.pdf", width = 10, height = 7)
  plot_R0_constant_IVM(results_constant)
  dev.off()
  
  results_time = analyze_time_dependent_IVM(params, VectPropIVM = c(0.5, 0.7, 0.9))
  
  plot_R0_time_dependent(results_time, params, "R0_effective_plot.pdf")
  
  return(list( constant = results_constant, time_dependent = results_time  ))
}

main_analysis()
