# plot_results.R
#  Created by John Trochta
#  This code loads in results from simulation runs and plots summary figures

library(dplyr)

options(scipen=9)
setwd(here::here())

# dir.ls <- c("simulations/Scenario_ignore_disease_multinom_runs.csv",
#             "simulations/Scenario_ignore_disease_runs.csv")
# scenario.names <- c("Multinomial\nlikelihood",
#                     "Binomial\nlikelihood")
 
# dir.ls <- c("simulations/Scenario_incorporate_seroprevalence_runs.csv",
#             "simulations/om_base_jan2021/Scenario_base_runs.csv"
#             )
# scenario.names <- c("Incorporate serorpevalence (Mar 2021)",
#                     "Incorporate serorpevalence (Jan 2021)")

# dir.ls <- c("simulations/Scenario_ignore_disease_runs.csv",
#             "simulations/Scenario_ignore_disease_greater survival changes_runs.csv")
# scenario.names <- c("Ignore disease (ver 1)",
#                     "Ignore disease (ver 2)")

base.om <- "simulations/base OM_age specific mixing/"

dir.ls <- c("Scenario_ignore_disease_runs.csv",
            "Scenario_incorporate_infection_prevalence_runs.csv",
            "Scenario_incorporate_seroprevalence_runs.csv",
            "Scenario_small_sample_size_runs.csv",
            "Scenario_early_mixing_of_susceptible_runs.csv",
            "Scenario_age_specific_mixing_ignored_runs.csv",
            "Scenario_time_varying_background_mortality_runs.csv",
            "Scenario_time_varying_background_mortality_ignore_disease_runs.csv",
            "Scenario_timevarying_disease_mortality_recovery_runs.csv")
dir.ls <- paste0(base.om,dir.ls)

scenario.names <- c("Ignore disease",
                    "Incorporate\ninfection\nprevalence",
                    "Incorporate\nseroprevalence",
                    "Small\nsample size",
                    "Early mixing\nof susceptibles",
                    "Age-specific\nmixing\nignored",
                    "Time-varying\nbackground\nmortality",
                    "Time-varying\nbackground\nmortality and\ndisease ignored",
                    "Time-varying\ndisease\nmortality/recovery")

# dir.ls <- c("simulations/base_age specific mixing/Scenario_time_varying_disease_mortality_recovery_greater survival changes_runs.csv")
# scenario.names <- c("Time-varying disease mortality ver2")

# Save em_runs to .csv in OM#_runs folder
sim.outputs <- data.frame(scenario=scenario.names[1],
                          read.csv(paste0(here::here("results/"),dir.ls[1]),stringsAsFactors=FALSE))


for(i in 2:length(dir.ls)){
  temp <- data.frame(scenario=scenario.names[i],
                            read.csv(paste0(here::here("results/"),dir.ls[i]),stringsAsFactors=FALSE))
  # sim.outputs <- rbind(sim.outputs,select(temp,-year))
  sim.outputs <- rbind(sim.outputs,temp)
  rm(temp)
}

# Check (and remove) non-converged models
nonconverged <- sim.outputs[sim.outputs$convergence!=0,]
sim.outputs <- sim.outputs[sim.outputs$convergence==0,]

convergence.check <- sim.outputs %>% group_by(scenario,seed) %>% summarize(con.rate = all(convergence==0)) 
convergence.rate <- convergence.check %>% group_by(scenario) %>% summarize(rate=sum(con.rate,na.rm = TRUE)/500)

##################################
# Plotting Relative Error of Key outputs (SSB, REC, Infection rates)
##################################
library(reshape2)
library(ggplot2)
library(grid)
library(gtable)
library(foreach)

# Remove NAs, but come back and check why there are NAs in the first place....
sim.outputs <- sim.outputs[!is.na(sim.outputs$scenario),]

# Factor scenario to arrange names of scenario
sim.outputs$scenario_fac <- factor(sim.outputs$scenario,levels=scenario.names)

# Scenarios where no disease is estimated (and should be ignored in calculating errors for disease-related quanities)
scen.no.dis <- c("Ignore disease","Incorporate\ninfection\nprevalence","Time-varying\nbackground\nmortality and\ndisease ignored")

# Calculate relative error for SSB and Recruitment
rel.err <- sim.outputs %>% group_by(scenario_fac,seed) %>% transmute('Year'=1:length(true_ssb),
                                                    'SB'=(est_ssb - true_ssb)/true_ssb,
                                                    'Recruitment'=(est_rec - true_rec)/true_rec,
                                                    'Infection rate'=ifelse(scenario_fac %in% scen.no.dis,NA,(est_infection-true_infection)))

# Melt the data frame
rel.err.2 <- melt(as.data.frame(rel.err),id=1:3)

# Calc 95% quantiles for each year
rel.err.3 <- rel.err.2 %>% group_by(scenario_fac,variable,Year) %>% summarize(Q.025=quantile(value,probs=0.025,na.rm=TRUE),
                                                                 Q.25=quantile(value,probs=0.25,na.rm=TRUE),
                                                                 Q.50=quantile(value,probs=0.5,na.rm=TRUE),
                                                                 Q.75=quantile(value,probs=0.75,na.rm=TRUE),
                                                                 Q.975=quantile(value,probs=0.975,na.rm=TRUE)) %>%
  filter(!is.na(Q.50))

rel.err.3$scenario <- sapply(rel.err.3$scenario_fac,FUN=function(x) paste(strwrap(x,width=25),collapse="\n"))

# Now plot
font.size <- 12

dat.2.plot <- filter(rel.err.3,scenario_fac %in% scenario.names[1:5])

re.plot <- ggplot(data=dat.2.plot,aes(x=Year,y=Q.50)) + 
  geom_ribbon(aes(ymin=Q.025,ymax=Q.975),fill="grey70")+
  geom_ribbon(aes(ymin=Q.25,ymax=Q.75),fill="grey85")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_line(size=1)+
  coord_cartesian(ylim=c(-1,1))+
  facet_grid(variable~scenario_fac,switch="y", drop=T)+
  scale_x_continuous(limits=c(0,50),breaks=c(0,50),expand=c(0.05,0.05))+
  theme_classic()+
  theme(strip.background = element_blank(),
        panel.border=element_rect(fill=NA),
        plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(size=font.size,vjust = 0),
        strip.text.x = element_text(size=font.size),
        strip.placement = "outside",
        axis.text.y = element_text(size=font.size-2),
        axis.text.x = element_text(size=font.size-2),
        axis.ticks.x= element_line(color="black"),
        # axis.title.y = element_text(size=font.size+2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=font.size),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.spacing = unit(0.25, "lines"),
        legend.position ="none")

# Remove empty facets
grob <- ggplotGrob(re.plot)
idx <- which(grob$layout$name %in% c("panel-3-1", "panel-3-2", "panel-3-8"));
for (i in idx) grob$grobs[[i]] <- nullGrob()

# Find and move x axes up of empty facets
idx <- which(grob$layout$name %in% c("axis-b-1", "axis-b-2", "axis-b-8"));
grob$layout[idx, c("t", "b")] <- grob$layout[idx, c("t", "b")] - c(2, 2, 2);

# Find and move y axes right to infection rate facets
idx <- which(grob$layout$name %in% c("axis-l-3","strip-l-3"));
grob$layout[idx, c("l","r")] <- grob$layout[idx, c("l","r")] + c(4,5);

# Remove x axis label
idx <- which(grob$layout$name %in% c("axis-b-2"));
for (i in idx) grob$grobs[[i]] <- nullGrob()

grid.newpage();
grid.draw(grob);

ggsave(plot=grob,filename=here::here(paste0("results/figures/Figure_time trajectories of relative error_set1_aug2021.png")),
       width=8, height=5, units="in",dpi=600)
# ggsave(plot=grob,filename=here::here(paste0("results/figures/Figure_time trajectories of relative error_set1.png")),
#        width=12, height=5, units="in",dpi=600)


dat.2.plot <- filter(rel.err.3,scenario_fac %in% scenario.names[c(3,6:9)])
re.plot <- ggplot(data=dat.2.plot,aes(x=Year,y=Q.50)) + 
  geom_ribbon(aes(ymin=Q.025,ymax=Q.975),fill="grey70")+
  geom_ribbon(aes(ymin=Q.25,ymax=Q.75),fill="grey85")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_line(size=1)+
  coord_cartesian(ylim=c(-1,1))+
  facet_grid(variable~scenario_fac,switch="y", drop=T)+
  scale_x_continuous(limits=c(0,50),breaks=c(0,50),expand=c(0.05,0.05))+
  theme_classic()+
  theme(strip.background = element_blank(),
        panel.border=element_rect(fill=NA),
        plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(size=font.size),
        strip.text.x = element_text(size=font.size),
        strip.placement = "outside",
        axis.text.y = element_text(size=font.size-2),
        axis.text.x = element_text(size=font.size-2),
        axis.ticks.x= element_line(color="black"),
        # axis.title.y = element_text(size=font.size+2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=font.size),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.spacing = unit(0.25, "lines"),
        legend.position ="none")

# Remove empty facets
grob <- ggplotGrob(re.plot)
idx <- which(grob$layout$name %in% c("panel-3-4"));
for (i in idx) grob$grobs[[i]] <- nullGrob()

# Find and move x axes up of empty facets
idx <- which(grob$layout$name %in% c("axis-b-4"));
grob$layout[idx, c("t", "b")] <- grob$layout[idx, c("t", "b")] - c(2);

grid.newpage();
grid.draw(grob);

ggsave(plot=grob,filename=here::here(paste0("results/figures/Figure_time trajectories of relative error_set2_aug2021.png")),
       width=8, height=5, units="in",dpi=600)

##################################
# MAREs of SSB, Recruitment, and infection
##################################
rel.err <- sim.outputs %>% group_by(scenario,seed) %>% 
  transmute(Year=1:length(true_ssb),
           ssb=(est_ssb - true_ssb)/true_ssb,
           rec=(est_rec - true_rec)/true_rec,
           inf=ifelse(scenario_fac %in% scen.no.dis,NA,(est_infection-true_infection)))

mare <- rel.err %>% group_by(scenario,Year) %>%
  summarise(mare.ssb = median(abs(ssb)),
    mare.rec = median(abs(rec)),
    mare.inf = median(abs(inf)))

mare.sum.derived <- mare %>% group_by(scenario) %>%
  summarise(finalyr.ssb = mare.ssb[Year==max(Year)],
            median.ssb = median(mare.ssb),
            median.rec = median(mare.rec),
            median.inf = median(mare.inf[Year %in% 2:(max(Year)-1)]))

##################################
# Probabilities of estimated biomass exceeding true biomass
##################################

threshold <- 0.4

# Calculate relative error for SSB and Recruitment
P.risk <- sim.outputs %>% group_by(scenario,seed) %>% 
  transmute('Year'=1:length(true_ssb),
            'Binary.above'=est_ssb > (true_ssb + threshold*true_ssb),
            'Binary.below'=est_ssb < (true_ssb - threshold*true_ssb))

P.risk <- P.risk %>% group_by(scenario,Year) %>% 
  summarise(P.above=sum(Binary.above)/length(Binary.above),
            #length.Pabove = length(Binary.above),
            P.below=sum(Binary.below)/length(Binary.below))

P.summary.0.4 <- P.risk %>% group_by(scenario) %>% 
  summarise(median.P.above.0.4=median(P.above),
            finalyr.P.above.0.4=P.above[Year==max(Year)],
            median.P.below.0.4=median(P.below),
            finalyr.P.below.0.4=P.below[Year==max(Year)])

threshold <- 0.1
P.risk <- sim.outputs %>% group_by(scenario,seed) %>% 
  transmute('Year'=1:length(true_ssb),
            'Binary.above'=est_ssb > (true_ssb + threshold*true_ssb),
            'Binary.below'=est_ssb < (true_ssb - threshold*true_ssb))
P.risk <- P.risk %>% group_by(scenario,Year) %>% 
  summarise(P.above=sum(Binary.above)/length(Binary.above),
            P.below=sum(Binary.below)/length(Binary.below))
P.summary.0.1 <- P.risk %>% group_by(scenario) %>% 
  summarise(median.P.above.0.1=median(P.above),
            finalyr.P.above.0.1=P.above[Year==max(Year)],
            median.P.below.0.1=median(P.below),
            finalyr.P.below.0.1=P.below[Year==max(Year)])

##################################
# Stats for describing bias in SB and R based on RE
##################################
threshold <- 0
P.risk <- sim.outputs %>% group_by(scenario,seed) %>% 
  transmute('Year'=1:length(true_ssb),
            'SB.above'=est_ssb > (true_ssb + threshold*true_ssb),
            'SB.below'=est_ssb < (true_ssb - threshold*true_ssb),
            'R.above'=est_rec > (true_rec + threshold*true_rec),
            'R.below'=est_rec < (true_rec - threshold*true_rec))
P.risk <- P.risk %>% group_by(scenario,Year) %>% 
  summarise(P.SB.above=sum(SB.above)/length(SB.above),
            P.SB.below=sum(SB.below)/length(SB.below),
            P.R.above=sum(R.above)/length(R.above),
            P.R.below=sum(R.below)/length(R.below))
P.summary.0.0 <- P.risk %>% group_by(scenario) %>% 
  summarise(median.SB.above.0.0=median(P.SB.above),
            finalyr.SB.above.0.0=P.SB.above[Year==max(Year)],
            median.SB.below.0.0=median(P.SB.below),
            finalyr.SB.below.0.0=P.SB.below[Year==max(Year)],
            median.R.above.0.0=median(P.R.above),
            finalyr.R.above.0.0=P.R.above[Year==max(Year)],
            median.R.below.0.0=median(P.R.below),
            finalyr.R.below.0.0=P.R.below[Year==max(Year)])
write.csv(P.risk,file=here::here(paste0("results/figures/Annual_Bias in SB and R by scenario.csv")),
          row.names=FALSE)
write.csv(P.summary.0.0,file=here::here(paste0("results/figures/Summary_Bias in SB and R by scenario.csv")),
          row.names=FALSE)
##################################
# Load and plot key parameters, particularly regarding disease
##################################
# Read in parameter estimates
# dir.ls <- c("om1_large samp_runs",
#             "om1_large samp_runs",
#             "om1_small samp_runs",
#             "om1_young sus vulnerable_runs")
# 
# file.ls <- c("vhs_asa_em1_parameters.RDS",
#              "vhs_asa_em1_no disease est_parameters.RDS",
#              "vhs_asa_em1_parameters.RDS",
#              "vhs_asa_em1_young susceptibles mixing_parameters.RDS")

dir.ls <- c("ignore_disease",
            "incorporate_infection_prevalence",
            "incorporate_seroprevalence",
            "small_sample_size",
            "early_mixing_of_susceptible",
            "age_specific_mixing_ignored",
            "time_varying_background_mortality",
            "time_varying_background_mortality_ignore_disease",
            "time_varying_disease_mortality_recovery")
dir.ls <- paste0(base.om,dir.ls)

file.ls <- c("vhs_asa_em_v1_ignore_disease_parameters.RDS",
             "vhs_asa_em_v1_incorporate_infection_prevalence_parameters.RDS",
             "vhs_asa_em_v1_incorporate_seroprevalence_parameters.RDS",
             "vhs_asa_em_v1_small_sample_size_parameters.RDS",
             "vhs_asa_em_v1_early_mixing_of_susceptible_parameters.RDS",
             "vhs_asa_em_v1_age_specific_mixing_ignored_parameters.RDS",
             "vhs_asa_em_v1_time_varying_background_mortality_parameters.RDS",
             "vhs_asa_em_v1_time_varying_background_mortality_ignore_disease_parameters.RDS",
             "vhs_asa_em_v1_timevarying_disease_mortality_recovery_parameters.RDS")

par_names <- c('Recovery rate',
               'Age 50% of susceptible mix',
               'Age 95% of susceptible mix',
               'Age 50% fish available to survey',
               'Age 95% fish available to survey',
               'Recruitment std dev',
               'Recruitment mean')
par_ID <- paste0("par_",1:7)

for(j in 1:length(dir.ls)){
  scen.dir <- paste0(here::here("results/"),dir.ls[j])
  setwd(scen.dir)
  folder.ls <- list.files()[1:500] # WATCH THIS
  est_pars <- foreach(i=1:length(folder.ls),.combine=rbind) %do% {
    # Set WD to current simulation rseed
    rseed = as.numeric(substring(folder.ls[i], 7))
    #setwd(folder.ls[i])
    if(file.exists(paste0(scen.dir,"/",folder.ls[i],"/",file.ls[j]))){
      pars <- readRDS(paste0(scen.dir,"/",folder.ls[i],"/",file.ls[j]))
      sel_pars <- pars$par[c("fix_recov_par","disease_selA50","disease_selA95","survey_selA50","survey_selA95","log_sigma_R","log_rbar")]
    }else{
      sel_pars <- rep(NA,times=7)
      names(sel_pars) <- c("fix_recov_par","disease_selA50","disease_selA95","survey_selA50","survey_selA95","log_sigma_R","log_rbar")
    }
    
    results <- data.frame(seed=rseed,t(sel_pars))
    
    results
  }
  
  colnames(est_pars) <- c("seed",par_ID)
  
  if(j==1){
    sim_pars <- data.frame(scenario=scenario.names[j],est_pars)
  }else{
    sim_pars <- rbind(sim_pars,data.frame(scenario=scenario.names[j],est_pars))
  }
  
}

# Store truth matrix
# Recovery rate, A50 disease, A95 disease, A50 survey sel, A95 survey sel, sigma R, mean R)
tru_pars <- matrix(c(26/(26+11), 3, 4, 3, 4, exp(0.1823216), exp(5.2),
                     26/(26+11), 3, 4, 3, 4, exp(0.1823216), exp(5.2),
                     26/(26+11), 3, 4, 3, 4, exp(0.1823216), exp(5.2),
                     26/(26+11), 3, 4, 3, 4, exp(0.1823216), exp(5.2),
                     26/(26+11), 1, 4, 3, 4, exp(0.1823216), exp(5.2),
                     26/(26+11), 3, 4, 3, 4, exp(0.1823216), exp(5.2),
                     26/(26+11), 3, 4, 3, 4, exp(0.1823216), exp(5.2),
                     26/(26+11), 3, 4, 3, 4, exp(0.1823216), exp(5.2),
                     NA,         3, 4, 3, 4, exp(0.1823216), exp(5.2)),ncol=7,byrow=TRUE)
colnames(tru_pars) <- par_ID
tru_pars <- data.frame(scenario=scenario.names,tru_pars)

# Melt data frames
sim_pars <- melt(as.data.frame(sim_pars),id=c("scenario","seed"))
tru_pars <- melt(tru_pars,id=c("scenario"))

# Convert logged parameter values to normal scale
sim_pars$value[sim_pars$variable=="par_6"] <- exp(sim_pars$value[sim_pars$variable=="par_6"])
sim_pars$value[sim_pars$variable=="par_7"] <- exp(sim_pars$value[sim_pars$variable=="par_7"])

# Calculate error
par_abs_err <- sim_pars %>% group_by(scenario,variable,seed) %>%
  transmute(are=abs((value-tru_pars$value[tru_pars$scenario==scenario & tru_pars$variable==variable])/tru_pars$value[tru_pars$scenario==scenario & tru_pars$variable==variable]))

par_mare <- par_abs_err %>% group_by(scenario,variable) %>% summarize(mare=median(are,na.rm=TRUE))

# Now re-insert names
par_mare$variable <- par_names[match(par_mare$variable,par_ID)]

# Extract only those parameters which I am presenting & cast to table
par_mare_2 <- filter(par_mare,variable %in% c('Recovery rate','Age 50% of susceptible mix','Age 95% of susceptible mix'))
par_mare_2 <- dcast(par_mare_2,scenario ~ variable)


##################################
# Combine errors for all quantities and scenarios
##################################
all_metrics <- left_join(left_join(par_mare_2,mare.sum.derived,by="scenario"),
                         P.summary.0.4,by="scenario")
all_metrics <- left_join(all_metrics,P.summary.0.1,by="scenario")
write.csv(all_metrics,file=here::here(paste0("results/figures/Table_performance metrics by scenario_aug2021.csv")),
          row.names=FALSE)













##################################
# EXTRA:  Check out individual replicates (for bimodal errors in estimation across simulations)
##################################
dir.ls <- "no_disease_multinom"
file.ls <- "vhs_asa_em_multinom_no_disease_multinom_ver2_parameters.RDS"

dir.ls <- "incorporate_seroprevalence"
file.ls <- "vhs_asa_em_v1_incorporate_seroprevalence_parameters.RDS"

dir.ls <- "incorporate_seroprevalence_prior"
file.ls <- "vhs_asa_em_v1_incorporate_seroprevalence_parameters.RDS"

dir.ls <- "om_base_jan2021"
file.ls <- "vhs_asa_em_v1_base_parameters.RDS"

seed <- "rseed_34029"
scen.dir <- paste0(here::here("results/simulations/"),dir.ls)

single_est_out_3 <- readRDS(paste0(scen.dir,"/",seed,"/",file.ls))


ggplot(data=test.plot.dat,aes(x=year,y=est_ssb,group=seed,color=seed)) + geom_line()
trial <- filter(test.plot.dat, scenario=='Multinomial\nlikelihood' & seed==34029)
trial.2 <- filter(sim.outputs, scenario=='Multinomial\nlikelihood' & seed==66556)
trial.3 <- filter(sim.outputs, scenario=='Multinomial\nlikelihood' & seed==113690)

plot(trial$true_rec,type="l")
lines(x=1:50,y=trial$est_rec,col="red")
plot(trial.2$true_rec,type="l")
lines(x=1:50,y=trial.2$est_rec,col="red")
plot(trial.3$true_rec,type="l")
lines(x=1:50,y=trial.3$est_rec,col="red")

plot(trial.2$true_ssb,trial.2$est_ssb)
plot(trial.2$true_ssb)
plot(trial.2$est_ssb)
plot(trial.2$true_rec)
plot(trial.2$est_rec)


##################################
# Extra plots
##################################

# Plot only disease related parameters
dat.4.plot <- filter(par_abs_err,variable %in% c('Recovery rate','Age 50% of susceptible mix',
                                                   'Age 95% of susceptible mix'))
ggplot(data=dat.4.plot,aes(x=scenario,y=error)) + 
  geom_boxplot()+
  geom_hline(yintercept=0,linetype="dashed")+
  #coord_cartesian(ylim=c(-1,1))+
  facet_wrap(variable~.,scales="free_x")+
  theme_classic()+
  theme(strip.background = element_blank(),
        panel.border=element_rect(fill=NA),
        plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(size=font.size-2),
        strip.text.x = element_text(size=font.size-2),
        axis.text.y = element_text(size=font.size-2),
        axis.text.x = element_text(size=font.size-2,angle=45,vjust=0.5),
        axis.ticks.x= element_line(color="black"),
        # axis.title.y = element_text(size=font.size+2),
        axis.title.y = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.spacing = unit(0.5, "lines"),
        legend.position ="none")+
  coord_flip()
ggsave(filename=here::here(paste0("results/figures/Figure_disease parameters absolue error.png")),
       width=10, height=5, units="in",dpi=600)






# Plot only survey selectivity parameters
dat.4.plot <- filter(par_abs_err,variable %in% c('Age 50% fish available to survey',
                                                 'Age 95% fish available to survey'))
ggplot(data=dat.4.plot,aes(x=scenario,y=error)) + 
  geom_boxplot()+
  geom_hline(yintercept=0,linetype="dashed")+
  #coord_cartesian(ylim=c(-1,1))+
  facet_wrap(variable~.,scales="free_x")+
  theme_classic()+
  theme(strip.background = element_blank(),
        panel.border=element_rect(fill=NA),
        plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(size=font.size-2),
        strip.text.x = element_text(size=font.size-2),
        axis.text.y = element_text(size=font.size-2),
        axis.text.x = element_text(size=font.size-2,angle=45,vjust=0.5),
        axis.ticks.x= element_line(color="black"),
        # axis.title.y = element_text(size=font.size+2),
        axis.title.y = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.spacing = unit(0.5, "lines"),
        legend.position ="none")+
  coord_flip()
ggsave(filename=here::here(paste0("results/figures/Figure_survey parameters absolue error.png")),
       width=10, height=5, units="in",dpi=600)
