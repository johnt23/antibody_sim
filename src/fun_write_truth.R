# fun_write_truth.R
# Created by John Trochta
# Date created: 06/20/2020
# Writes observations generated from OM to .dat file for input into TMB

fun_write_truth = function(SSB,Nya,Nya_new_infect,Nya_sel_sus,Nya_sus,Sya,
                           outbreak.dur.inci,outbreak.dur.prev,peak.inci,peak.prev,true_immune,obs_immune,day.peak.prevalence,
                           rseed,flag){
  nan.check = function(x){return(ifelse(is.nan(x),-999,x))}
  dat <- list("# Simulated data from SIR age-structured hybrid model",
        paste0("# random seed = ",rseed),
        paste0("# Flag=0 if no issues, =1 if negative pop #'s, =2 if NaN, or=3 if NA. Flag:"),
        flag,"",
        "# True SSB",
        SSB,"",
        "# True R (numbers)",
        Nya[,1],"",
        "# True infection rate given (selected) susceptible",
        nan.check(Nya_new_infect/Nya_sel_sus),"",
        "# Duration of outbreak (days) based on incidence",
        nan.check(outbreak.dur.inci),"",
        "# Duration of outbreak (days) based on prevalence",
        nan.check(outbreak.dur.prev),"",
        "# Peak (point) prevalence",
        nan.check(peak.prev),"",
        "# Day of peak prevalence",
        nan.check(day.peak.prevalence),"",
        "# Total immunity",
        nan.check(true_immune),"",
        "# Observed immunity (based on survey selectivity)",
        nan.check(obs_immune),"",
        "# True total proportion of susceptible",
        nan.check(Nya_sus/Nya),"",
        "# True vulnerable proportion of susceptible",
        nan.check(Nya_sel_sus/Nya),"",
        "# True survival by age",
        Sya,"",
        "# True Numbers-at-age",
        nan.check(Nya),
        "")
  lapply(dat,write.table,file="truth.dat",append=TRUE,sep = " ",row.names=FALSE,col.names=FALSE,quote=F)
  # write.table(dat,file = "VHS_ASA.truth", append = F, sep = " ",
  #             row.names=FALSE,col.names=FALSE,quote=F)
}