# fun_write_dat.R
# Created by John Trochta
# Date created: 06/20/2020
# Writes observations generated from OM to .dat file for input into TMB

fun_obs_mod = function(obs_years,Nya,comp_samp_size,catch_comp_samp_size,
                       predicted_comps,survey_slx,dis_survey_slx,
                       predicted_immune_comps,antibody_comp_samp_size,
                       predicted_survey,survey_cv,Cya,N_catch,
                       true_samp_prev,
                       rseed){
set.seed(rseed)
survey_obs = exp(rnorm(n=length(predicted_survey),mean=log(predicted_survey),sd=rep(survey_cv,length(predicted_survey))))
comp_obs = t(apply(predicted_comps,1,function(x) rmultinom(n=1,size=comp_samp_size,prob=x)))

# Fishery age compositions
catch_comps <- matrix(NaN,nrow = nrow(Cya),ncol = ncol(Cya))
if(any(Cya>0)){
  catch_comps = t(apply(Cya,1,function(x) rmultinom(n=1,size=catch_comp_samp_size,prob=x/sum(x))))
}

# predicted_immuned_comps gives age-specific proportions immune
# For EM model, have to convert these to age-stage specific proportions of entire population
immune_mat = cbind(sweep(predicted_immune_comps*Nya,MARGIN=2,dis_survey_slx,'*'),sweep((1-predicted_immune_comps)*Nya,MARGIN=2,dis_survey_slx,'*'))
imm_ind = c(rbind(1:ncol(predicted_immune_comps),(ncol(predicted_immune_comps)+1):(2*ncol(predicted_immune_comps))))
antibody = immune_mat[,imm_ind]

# The column ordering is such: Age 0 immune, Age 0 non-immune, Age 1 immune, Age 1 non-immune, ..., Age nage immune
antibody_obs = t(apply(antibody,1,function(x) rmultinom(n=1,size=antibody_comp_samp_size,prob=x)))

# Infection prevalence indices
obs_samp_prev <- apply(true_samp_prev,1,function(x) mean(rbinom(3,60,x)/60))

return(list(survey_obs=survey_obs,
            comp_obs=comp_obs,
            antibody_obs=antibody_obs,
            catch_comps=catch_comps,
            obs_samp_prev=obs_samp_prev))
}