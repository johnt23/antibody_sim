# fun_write_dat.R
# Created by John Trochta
# Date created: 06/20/2020
# Writes observations generated from OM to .dat file for input into TMB

fun_write_dat = function(obs_years,nyr,sage,nage,comp_samp_size,catch_comp_samp_size,antibody_comp_samp_size,
                         survey_obs,catches,catch_comps,
                         comp_obs,antibody_obs,avg_waa,maturity_A50,maturity_A95,
                         obs_samp_prev,
                         rseed,flag){
na.check = function(x){return(ifelse(is.na(x),-999,x))}
  dat <- list("# Simulated data from SIR age-structured hybrid model",
        paste0("# random seed = ",rseed),
        paste0("# Flag=0 if no issues, =1 if negative pop #'s, =2 if NaN, or=3 if NA. Flag:"),
        flag,"",
        "# Number of years (nyr)",
        length(obs_years),"",
        "# Starting age",
        sage,"",
        "# Ending age",
        nage,"",
        "# Sample size for comp data",
        comp_samp_size,"",
        "# Sample size for fishery age comp data",
        catch_comp_samp_size,"",
        "# Sample size for antibody prevalence data",
        antibody_comp_samp_size,"",
        "# Survey data",
        na.check(survey_obs),"",
        "# Catches",
        na.check(catches),"",
        "# Catch compositions",
        na.check(catch_comps),"",
        "# Age composition from fishery-independent survey",
        na.check(comp_obs),"",
        "# Antibody prevalence values",
        na.check(antibody_obs),"",
        "# Average weight-at-age",
        avg_waa,"",
        "# Maturity logistic function parameter age-at-50%",
        maturity_A50,"",
        "# Maturity logistic function parameter age-at-95%",
        maturity_A95,"",
        "# Infection prevalence index (avg of 3 different samples)",
        obs_samp_prev,
        "")
  lapply(dat,write.table,file="vhs_asa_em.dat",append=TRUE,sep = " ",row.names=FALSE,col.names=FALSE,quote=F)
  # write.table(dat,file = "vhs_asa_em.dat", append = F, sep = " ",
  #             row.names=FALSE,col.names=FALSE,quote=F)
}