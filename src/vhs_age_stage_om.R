# vhs_age_stage_om.R
# Created by John Trochta
# Date created: 04/02/2020
# Operating model of VHS outbreaks in each year

# Things to do (05/13/2020)
# Reorder for loop for simulating outbreak
# - Consider rewriting equations like BASA so projections are not i+1 (but going from past)


vhs.asa = function(nyr=200,nage=5,nstage=3,ndays=30,
                   avg_waa,
                   vhs_trans_rate_I=1,
                   vhs_trans_rate_C=0.1,
                   vhs_mort_rate=91.32/365,vhs_rec_rate=91.32/365,
                   dep_scaling = 0,
                   fishing_mort=0,
                   survey_selA50=3,
                   survey_selA95=4,
                   maturity_A50=1,
                   maturity_A95=2,
                   disease_vulA50=3,
                   disease_vulA95=4,
                   dis_survey_selA50=2.5,
                   dis_survey_selA95=3.5,
                   nonlinear_exp_1 = 1,
                   nonlinear_exp_2 = 1,
                   sig_nat_mor = 0,
                   selA50 = 2,
                   selA95 = 3,
                   ignore.carryover.inf,
                   inf_prev_survey,
                   log_sigma_R,
                   rseed){
  
  set.seed(rseed)
  
  #### Age-structured population dynamics variables ####
  nyr <- nyr
  nage <- nage
  nstage <- nstage
  ndays <- ndays
  
  # Key population and fishing parameters
  ab_logit_midpoint <- 0.5
  ab_logit_k <- 10
  ab_logit_max <- 0.3
  ac_coef_rec <- 0.6
  
  # For time-varying option on natural mortality
  mor_devs <- rep(1,times=nyr)
  if(sig_nat_mor>0){
    # mor_devs <- rbeta(nyr,0.5,18)
    mor_devs <- exp(rnorm(nyr,0,sig_nat_mor))
    # mor_devs[1] <- rnorm(1,0.25,sig_nat_mor)
    # for(i in 2:nyr){
    #   mor_devs[i] <- rnorm(1,mor_devs[i-1],sig_nat_mor)
    # }
  }
  
  natural_mortality <- 0.25*mor_devs
  #plot(natural_mortality,type="l")

  
  log_rbar <- 5.2
  log_q_survey <- -0.5
  
  #plus_group_mortality <- 1.8
  plus_group_mortality <- natural_mortality
  
  log_sigma_R <- log_sigma_R
  selA50 <- selA50
  selA95 <- selA95
  
  survey_selA50 <- survey_selA50
  survey_selA95 <- survey_selA95
  
  maturity_A50 <- maturity_A50
  maturity_A95 <- maturity_A95
  
  disease_vulA50 <- disease_vulA50
  disease_vulA95 <- disease_vulA95

  # Transformed parameters
  sigma_R <- exp(log_sigma_R)
    
  log_rbar_devs <- rnorm(nyr,0,sigma_R)
  
  # Demographic variables - scalar OR vector values
  mature = round(1.0 / (1.0 + exp(-log(19)*((0:(nage-1))-maturity_A50)/(maturity_A95-maturity_A50))),2)
  surv_daily <- exp(-natural_mortality/365) # Daily survival probability based on an annual baseline natural mortality rate
  
  surv_annual_plus <- exp(-plus_group_mortality) # For one full year - natural mortality is an age-specific vector
  surv_daily_plus <- exp(-plus_group_mortality/365) # Daily survival probability based on an annual baseline natural mortality rate
  
  fish_slx <- round(1.0 / (1.0 + exp(-log(19)*((0:(nage-1))-selA50)/(selA95-selA50))),2)
  survey_slx <- round(1.0 / (1.0 + exp(-log(19)*((0:(nage-1))-survey_selA50)/(survey_selA95-survey_selA50))),2)
  disease_slx <- round(1.0 / (1.0 + exp(-log(19)*((0:(nage-1))-disease_vulA50)/(disease_vulA95-disease_vulA50))),2)
  dis_survey_slx <- round(1.0 / (1.0 + exp(-log(19)*((0:(nage-1))-dis_survey_selA50)/(dis_survey_selA95-dis_survey_selA50))),2)
  
  # fish_slx <- 1.0 / (1.0 + exp(-log(19)*((0:(nage-1))-selA50)/(selA95-selA50)))
  # survey_slx <- 1.0 / (1.0 + exp(-log(19)*((0:(nage-1))-survey_selA50)/(survey_selA95-survey_selA50)))
  # disease_slx <- 1.0 / (1.0 + exp(-log(19)*((0:(nage-1))-disease_vulA50)/(disease_vulA95-disease_vulA50)))
  # dis_survey_slx <- 1.0 / (1.0 + exp(-log(19)*((0:(nage-1))-dis_survey_selA50)/(dis_survey_selA95-dis_survey_selA50)))
  
  # Correct Age 0 when each age has 100% vulnerability
  disease_slx[is.nan(disease_slx) | is.infinite(disease_slx)] <- 1.0
  
  avg_waa = avg_waa
  
  # Set exploitation rate
  fish_rate <- fishing_mort #exp(-fishing_mort) # The annual exploitation rate
  # exp_daily <- c(rep(0,ndays/2-1),1,rep(0,ndays/2))# Distribution of exploitation over the season of disease epidemics
  
  # Log negative numbers starting (inserted 01/26/2021)
  log_N_init <- c(log_rbar,
                  log_rbar*exp(-natural_mortality)^(1:(nage-2)),
                  log_rbar*exp(-natural_mortality)^(nage-2)*exp(-natural_mortality)/(1-exp(-natural_mortality)))
  
  # State variables
  catches <- rep(0,nyr)
  catch_comps <- matrix(0,nyr,nage)
  Cya <- matrix(0,nyr,nage)
  N_catch <- rep(0,nyr)
  Nya <- matrix(0,nyr+1,nage)
  Post_out_Nya <- matrix(0,nyr,nage)
  predicted_comps <- matrix(0,nyr,nage)
  predicted_immune_comps <- matrix(0,nyr,nage)
  predicted_survey <- rep(0,nyr)
  obs_prop_immune <- rep(0,nyr)
  true_prop_immune <- rep(0,nyr)
  recruits <- rep(0,nyr)
  season_survival <- matrix(0,nyr,nage)
  spawners <- matrix(0,nyr,nage)
  SSB_atage <- matrix(0,nyr,nage)
  SSB <- rep(0,nyr)
  Sya <- matrix(0,nyr,nage)
  # u <- rep(0,nyr)
  # vulnerable <- rep(0,nage)
  
  #### Disease-related variables ####
  
  # Disease age and stage structure by year (beginning of each year)
  N_yas <- matrix(0,nyr+1,nage*nstage)
  
  # Disease age and stage structure by day - changes within the year
  N_das <- matrix(0,ndays,nage*nstage)
  
  # Vector of indices where # susceptible, # infected, and # recovered are placed within N_das
  infected_indices <- seq(from=2,to=nage*nstage,by=nstage)
  sus_indices <- seq(from=1,to=nage*nstage,by=nstage)
  rec_indices <- seq(from=3,to=nage*nstage,by=nstage)
  
  # Numbers of daily disease deaths by age (this matrix is reset for every year)
  N_daily_disease_deaths <- matrix(0,ndays,nage)
  
  # Numbers of NEW infections by age (this matrix is reset for every year)
  N_daily_infections <- matrix(0,ndays,nage)
  
  # Numbers of NEW infections by age (this matrix is reset for every year)
  N_daily_recoveries <- matrix(0,ndays,nage)
  
  # Daily trajectory of susceptible, infected, and recovered stages within each year (summed across all ages)
  N_y_susceptible <- matrix(0,nyr,ndays)
  N_y_infected <- matrix(0,nyr,ndays)
  N_y_recovered <- matrix(0,nyr,ndays)
  N_y_dis_deaths <- matrix(0,nyr,ndays)
  
  incidence_prop <- matrix(0,nyr,ndays)
  point_prevalence <- matrix(0,nyr,ndays)
  
  # Cumulative disease deaths and infections per age per year
  Nya_sus <- matrix(0,nyr+1,nage)
  Nya_sel_sus <- matrix(0,nyr+1,nage)
  Nya_inf <- matrix(0,nyr+1,nage)
  Nya_imm <- matrix(0,nyr+1,nage)
  Nya_dis_death <- matrix(0,nyr,nage)
  Nya_new_infect <- matrix(0,nyr,nage)
  Nya_new_recoveries <- matrix(0,nyr,nage)
  Pya_immune_postout <- matrix(0,nyr,nage)
  
  # Average rate of transmission from infected ages jth column to to susceptible ages in ith row
  mu_base_I <- vhs_trans_rate_I # actual rate days^-1 = 91.32/365 = 0.2502- This is age-invariant
  mu_base_C <- vhs_trans_rate_C # actual rate days^-1 = 91.32/365 = 0.2502- This is age-invariant
  mu_I <- matrix(mu_base_I,nrow=nage,ncol=nage,byrow=TRUE)
  mu_C <- matrix(mu_base_C,nrow=nage,ncol=nage,byrow=TRUE)
  
  # Probability of mortality given infection
  alpha <- t(sapply(vhs_mort_rate,function(x) rep(1-exp(-x),times=nage)))
  # alpha <- rep(1-exp(-91.32/365),times=nage)
  # alpha <- rep(1-exp(-0.166),times=nage)
  
  # Transition probability from infection to recovered
  rho <- t(sapply(vhs_rec_rate,function(x) rep(1-exp(-x),times=nage)))
  # rho <- rep(1-exp(-91.32/365),times=nage)
  
  # Calculate initial numbers
  #Nya[1,1] = exp(N_tot*cust_exp_fun(c(1,2:nage))[1])
  #Nya[1,2:nage] = exp(log_N_init)
  
  Nya[1,1] = exp(5.2)
  Nya[1,2:(nage-1)] = exp(-0.25)^(1:length(2:(nage-1)))*Nya[1,1]
  Nya[1,nage] = (exp(-0.25)/(1-exp(-0.25)))* Nya[1,(nage-1)]
  
  # init_stage_proportions <- c(rep(c(1,0,0),3),rep(c(0.4,0.3,0.3),nage-3)) # based on 3 stages & 10 ages
  init_stage_proportions <- rep(c(0.9,0,0.1),nage) # based on 3 stages
  N_yas[1,] <- rep(Nya[1,],times=rep(3,nage))*init_stage_proportions
  N_das[1,] <- N_yas[1,]
  Nya_sus[1,] <- N_yas[1,sus_indices]
  Nya_sel_sus[1,] <- N_yas[1,sus_indices]*disease_slx
  Nya_inf[1,] <- N_yas[1,infected_indices]
  Nya_imm[1,] <- N_yas[1,rec_indices]
  
  # Initial age-stage vulnerability to transmission
  # - describes the proportion of each age-stage class that is mixed with the population where
  #   transmission occurs
  dis_vul_vector <- c(sapply(disease_slx,FUN=function(x) c(x,1,1)))
  
  # Natural survival of non-vulnerable age-stage groups (susceptibles) during transmission season
  # season_base_survival <- t(sapply(natural_mortality,function(x) c(rep(exp(-x/(365)*(ndays)),times=nstage*(nage-1)),rep(exp(-plus_group_mortality/(365)*(ndays)),times=nstage))))
  season_base_survival <- t(sapply(natural_mortality,function(x) c(rep(exp(-x/(365)*(ndays)),times=nstage*(nage)))))
  
  ## START HERE - 09/23/2020 (Try running this bit to make sure matrix comes out alright)
  
  i <- 1
  beta_I <- matrix(0,ndays,nage)
  beta_C <- matrix(0,ndays,nage)
  
  beta_I_y <- matrix(0,nyr,ndays)
  beta_C_y <- matrix(0,nyr,ndays)
  
  for (i in 1:nyr){
     # Sample first
     # Assume survey immediately AFTER disease and fishing
     predicted_comps[i,] <- (survey_slx*Nya[i,])/sum(survey_slx*Nya[i,])
     predicted_survey[i] = exp(log_q_survey)*sum((Nya[i,])*survey_slx*avg_waa)
     
     # This bit is takes vector N_yas[i,] and sums every nstage elements (tapply command) to produce vector of length nage
     # seqs = seq_along(N_yas[i,])
     # tapply(N_yas[i,],rep(seqs,each=nstage)[seqs],FUN=sum)
     
     # This the the proportion of each age that show antibodies
     predicted_immune_comps[i,] <- N_yas[i,rec_indices]/Nya[i,]
     
     # Since antibody test may not represent all forms of immunity, we scale true immune down according to logistic curve
     true_prop_immune[i] = sum(N_yas[i,rec_indices])/sum(Nya[i,])
     obs_prop_immune[i] = sum(dis_survey_slx*N_yas[i,rec_indices])/sum(dis_survey_slx*Nya[i,])
     # predicted_immune_comps[i,]=ab_logit_max/(1+exp(-ab_logit_k*(true_prop_immune[i]-ab_logit_midpoint)))*predicted_immune_comps[i,]
     
     spawners[i,]=mature*(Nya[i,])
     SSB_atage[i,]=spawners[i,]*avg_waa
     
     # Now this could be use for stock-recruit function
     SSB[i]=sum(SSB_atage[i,])#prop_female*sum(SSB_atage[i,])
     
     # Apply fishing - calculate catch-at-age calculation
     Cya[i,] <- Nya[i,] * fish_slx * (fish_rate[i])
     N_catch[i] <- sum(Cya[i,])
     catches[i] <- sum(Cya[i,]*avg_waa)
     catch_comps[i,] <- (Cya[i,]*avg_waa)/catches[i]
     
     # Initialize numbers in each stage at the start of the next year/season
     N_das[1,] <- N_yas[i,]*rep((1-fish_slx * (fish_rate[i])),each=nstage)*dis_vul_vector
     beta_I <- matrix(0,ndays,nage)
     beta_C <- matrix(0,ndays,nage)
     
    for(d in 1:(ndays-1)){
      
      for(j in 1:nage){
        age_vector <- ((j-1)*nstage+1):((j-1)*nstage+nstage)
        
        # Frequency dependent
        beta_I[d,j] <- 1-exp(-sum(mu_I[j,]*(N_das[d,sus_indices])*N_das[d,infected_indices]^nonlinear_exp_1)*sum(N_das[d,])^dep_scaling/sum(N_das[d,]))
        beta_C[d,j] <- 1-exp(-sum(mu_C[j,]*(N_das[d,sus_indices])*N_das[d,rec_indices]^nonlinear_exp_2)*sum(N_das[d,])^dep_scaling/sum(N_das[d,]))
        # SIR is set-up so that transition probability of each column into row i
        sir <- t(matrix(c(1-beta_I[d,j]-beta_C[d,j],   0,                    0,
                          beta_I[d,j]+beta_C[d,j],      (1-alpha[i,j]-rho[i,j]),  0,
                          0,                           rho[i,j],               1),nrow=nstage,byrow = TRUE))
        
        # SIR is set-up so that transition probability of each column into row i
        N_das_temp <- ((N_das[d,age_vector]) %*% sir)
        
        # On last day, assume remaining infected either die or recover to the exact proportion from the combined rate
        # This insures no infecteds remain, as typically observed in the field
        if(d==(ndays-1) & ignore.carryover.inf){
            final_recov <- ifelse(rho[i,j]==0,0,rho[i,j]/(alpha[i,j]+rho[i,j]))
            N_das_temp <- N_das_temp + c(0,-N_das_temp[2],final_recov*N_das_temp[2])
        }
        
        # CHECK THIS HARD!!!!!
        # Getting some weird thing where N_das will go negative (numbers in certain ages and stages below 0, ONLY when incorporating some beta_C value>0)
        # N_das_temp[N_das_temp<0] <- 0
        
        # Now calculate natural survival
        if(j<nage){
          surv_temp <- surv_daily[i]
        }else{
          surv_temp <- surv_daily_plus[i]
        }
         
        # Number new infections per day - take the 1st index of the N_das[d,age_vector] vector which is the # susceptible in age j
        # N_daily_infections[d+1,j] <- (beta_I[d,j]+beta_C[d,j])*N_das[d,age_vector][1]
        N_daily_infections[d+1,j] <- -(N_das_temp[1] - N_das[d,age_vector][1])
        
        # Number daily recoveries
        # N_daily_recoveries[d+1,j] <- rho[i,j]*N_das[d,age_vector][2]
        N_daily_recoveries[d+1,j] <- (N_das_temp[3] - N_das[d,age_vector][3])
        
        # Number deaths due to disease - take the 2nd index of the N_das[d,age_vector] vector which is the # infected in age j
        # N_daily_disease_deaths[d+1,j] <- alpha[i,j]*N_das[d,age_vector][2]
        N_daily_disease_deaths[d+1,j] <- -(N_das_temp[2] - N_das[d,age_vector][2] - N_daily_infections[d+1,j] + N_daily_recoveries[d+1,j])
        
        # Now calculate natural survival
        N_das[d+1,age_vector] <- N_das_temp * surv_temp
        
        if(d==(ndays-1)){
          Post_out_Nya[i,j] <- sum((N_das[d+1,]+N_yas[i,]*rep((1-fish_slx * (fish_rate[i])),each=nstage)*(1-dis_vul_vector)*season_base_survival[i,])[age_vector])
          # age-specific (row) stage (column) proportions at the end of the season
          
          # % that are recovered within each age group at the end of the season
          Pya_immune_postout[i,] <- N_das[d+1,age_vector][3]/sum(N_das[d+1,age_vector])
        }
      }
    }
    
    # Store infection probabilities
    beta_I_y[i,] <- beta_I[,nage]
    beta_C_y[i,] <- beta_C[,nage]
     
    # Daily numbers in each stage by year
    N_y_susceptible[i,] <- apply(N_das[,sus_indices],1,sum)
    N_y_infected[i,] <- apply(N_das[,infected_indices],1,sum)
    N_y_recovered[i,] <- apply(N_das[,rec_indices],1,sum)
    
    # Daily # of disease deaths throughout the season (year x day within season)
    N_y_dis_deaths[i,] <- apply(N_daily_disease_deaths,1,sum)
    
    # Total numbers of deaths, new infections, and new recoveries (immune) of fish by age in each year i
    Nya_dis_death[i,] <- apply(N_daily_disease_deaths,2,sum)
    Nya_new_infect[i,] <- apply(N_daily_infections,2,sum)
    Nya_new_recoveries[i,] <- apply(N_daily_recoveries,2,sum)
    
    # Incidence proportion
    incidence_prop[i,] <- apply(N_daily_infections,1,sum)/apply(N_das,1,sum)
    
    # Point prevalence
    point_prevalence[i,] <- N_y_infected[i,]/apply(N_das,1,sum)
    
    # Now add in those non-vulnerable susceptibles back into the final day of the N_das matrix
    N_das[ndays,] <- N_das[ndays,] + N_yas[i,]*rep((1-fish_slx * (fish_rate[i])),each=nstage)*(1-dis_vul_vector)*season_base_survival[i,]
    
    # Cumulative survival over this period
    season_survival[i,] <- (Post_out_Nya[i,])/Nya[i,]
    
    for(j in 1:(nage-1)){
      age_vector_1 <- ((j-1)*nstage+1):((j-1)*nstage+nstage)
      age_vector_2 <- ((j)*nstage+1):((j)*nstage+nstage)
      N_yas[i+1,age_vector_2] <- N_das[ndays,age_vector_1]*exp(-natural_mortality[i]/(365)*(365-ndays+1))
      Nya[i+1,j+1] <- sum(N_yas[i+1,age_vector_2])
      
      # Realized total survival rate - I add catches back in an project their survival to not bias this estimate
      # Sya[i,j] <- (Nya[i,j] - Nya[i+1,j+1] + Cya[i,j])/Nya[i,j]
      Sya[i,j] <- (Nya[i+1,j+1] + Cya[i,j]*exp(-natural_mortality[i]/(365)*ndays))/Nya[i,j]
    }
    # Now calculate realized survival of plus group only
    # Sya[i,j+1] <- (Nya[i,j+1]- sum(N_das[ndays,age_vector_2]*exp(-natural_mortality/(365)*ndays)) )/Nya[i,j+1] - INCORRECT; used in results before 12/16/2020
    Sya[i,nage] <- (sum(N_das[ndays,age_vector_2]*exp(-plus_group_mortality[i]/(365)*(365-ndays+1))) + Cya[i,j+1]*exp(-plus_group_mortality[i]/(365)*ndays))/Nya[i,nage]
    
    # Also account for numbers in plus group (sum numbers from nage-1 moving into nage, and this year's nage)
    # N_yas[i+1,age_vector_2] <- N_yas[i+1,age_vector_2]+N_das[ndays,age_vector_2]*exp(-plus_group_mortality/(365)*ndays) - INCORRECT; used in results before 12/16/2020
    N_yas[i+1,age_vector_2] <- N_yas[i+1,age_vector_2]+N_das[ndays,age_vector_2]*exp(-plus_group_mortality[i]/(365)*(365-ndays+1))
    
    # Reassign value in final age group
    Nya[i+1,j+1] <- sum(N_yas[i+1,age_vector_2])
    
    # NEXT YEAR'S RECRUITMENT!
    #recruits(i) = alpha*SSB(i)*exp(-beta*P_indices(0,2)*SSB(i)) #  Plain Ricker - uses initial year P_index to appropriately scale beta so I don't have to change PIN values or bounds
    recruits[i] = exp(log_rbar) # Stationary recruitment - independent of biomass
    
    # Nya[i+1,1] = exp(log_rbar)
    Nya[i+1,1] = recruits[i]*exp(ac_coef_rec*(log(Nya[i,1])-log_rbar+0.5*sigma_R^2) + sqrt(1-ac_coef_rec^2)*log_rbar_devs[i+1]-0.5*sigma_R^2)
    
    # Log-normal independent recruitment
    N_yas[i+1,1] = Nya[i+1,1]
    N_yas[i+1,2:nstage] <- 0
    
    # Stage proportions in next year
    Nya_sus[i+1,] <- N_yas[i+1,sus_indices]
    Nya_sel_sus[i+1,] <- N_yas[i+1,sus_indices]*disease_slx
    Nya_inf[i+1,] <- N_yas[i+1,infected_indices]
    Nya_imm[i+1,] <- N_yas[i+1,rec_indices]
  }
  
  # Remove extra year
  N_yas <- N_yas[-(nyr+1),]
  Nya <- Nya[-(nyr+1),]
  Nya_imm <- Nya_imm[-(nyr+1),]
  Nya_inf <- Nya_inf[-(nyr+1),]
  Nya_sel_sus <- Nya_sel_sus[-(nyr+1),]
  Nya_sus <- Nya_sus[-(nyr+1),]
  
  # Epidemic characteristics:
  # Mean duration (based on incidence rate exceeding 0.01)
  duration.incidence <- apply(incidence_prop,1,function(x) sum(x>=0.1))
  duration.prevalence <- apply(point_prevalence,1,function(x) sum(x>=0.1))
  
  # Peak prevalence/incidence
  max.incidence <- apply(incidence_prop,1,max)
  max.prevalence <- apply(point_prevalence,1,max)
  
  # Timing of peak prevalence
  day.max.incidence <- apply(incidence_prop,1,function(x) which.max(x))
  day.max.prevalence <- apply(point_prevalence,1,function(x) which.max(x))
  
  # Mock VHSV prevalence survey
  N_inf_prev <- N_y_infected/(N_y_susceptible + N_y_infected + N_y_recovered)
  N_inf_prev <- cbind(1:nrow(N_inf_prev),N_inf_prev)
  true_samp_prev <- t(apply(N_inf_prev,1,function(x) x[-1][inf_prev_survey[x[1],]]))
  
  # Function for simultaneously creating names of variables/elements within a list
  listN <- function(...){
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    anonList
  }
  return(listN(beta_C_y,beta_I_y,catches,catch_comps,Cya,day.max.incidence,day.max.prevalence,disease_slx,dis_survey_slx,duration.incidence,duration.prevalence,fish_slx,
               incidence_prop,max.incidence,max.prevalence,natural_mortality,
               N_catch,Nya,Nya_dis_death,Nya_inf,Nya_imm,Nya_new_infect,
               Nya_new_recoveries,N_yas,Nya_sus,Nya_sel_sus,N_y_dis_deaths,
               N_y_infected,N_y_recovered,N_y_susceptible,obs_prop_immune,
               point_prevalence,Post_out_Nya,Pya_immune_postout,season_survival,
               predicted_comps,predicted_immune_comps,predicted_survey,recruits,true_prop_immune,
               spawners,SSB_atage,SSB,survey_slx,Sya,true_samp_prev))
}
