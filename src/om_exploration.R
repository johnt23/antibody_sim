# om_exploration.R
# Created by John Trochta
# Date created: 06/20/2020
# Run simulation model from R of the disease antibodies
rm(list=ls(all=TRUE)) #clears workspace

source(here::here("src/vhs_age_stage_om.R"))
source(here::here("src/fun_obs_mod.R"))

nyr=100
nage=8
nstage=3
ndays=120
avg_waa=c(70, 94, 115, 134, 150, 160, 165, 168)

survey_selA50 <- 3
survey_selA95 <- 4
maturity_A50 <- 3
maturity_A95 <- 4
disease_vulA50 <- 3
disease_vulA95 <- 4
dis_survey_selA50 <- 3
dis_survey_selA95 <- 4

selA50 = 3
selA95 = 4
log_sigma_R = 0.1823216
nonlinear_exp_1 = 1
nonlinear_exp_2 = 1
sig_nat_mor = 0

ignore.carryover.inf <- TRUE

# VHS per capita epidemiological rates (day^-1)
# Frequency dependence value
# vhs_trans_rate_I=0
# vhs_trans_rate_C=0

vhs_trans_rate_I=0.01
vhs_trans_rate_C=0.000001

# vhs_trans_rate_I=0.0001
# vhs_trans_rate_C=0.00000001

#vhs_mort_rate=rep(30/365,times=nyr) # 11/365
#vhs_rec_rate=rep(5/365,times=nyr) # 26/365

vhs_mort_rate=rep(11/365,times=nyr) # 11/365
vhs_rec_rate=rep(26/365,times=nyr) # 26/365
# vhs_mort_rate = runif(nyr,9/365,21/365)
# vhs_rec_rate = runif(nyr,20/365,70/365) # 26/365

# Dependency scaling:
# This is an exponent in the transmission rate function that allows
# scaling from frequency dependent (dep_scaling=0) to density 
# dependent (dep_scaling=1)
dep_scaling=0

# Fishing mortality rate (year^-1)
# fishing_mort=c(rep(0,times=150),seq(0.05,0.4,length.out=20),rep(0.4,times=10),rep(0.05,times=10),rep(0.2,times=10))
fishing_mort=rep(0,times=nyr)

# Disease prevalence sampling dates (Generate random sequential dates within the first month of transmission)
inf_prev_survey <- round(runif(nyr,5,30))
inf_prev_survey <- cbind(inf_prev_survey,inf_prev_survey + round(runif(nyr,1,5)))
inf_prev_survey <- cbind(inf_prev_survey,inf_prev_survey[,2] + round(runif(nyr,1,5)))

# rseed = ceiling(runif(1,50,10000000)) # 90854094 # 2770367
rseed = 2770367
set.seed(rseed)

#rseed <- runif(1,100,10^7) #18314367
ptm <- proc.time()
run.1 = vhs.asa(nyr,nage,nstage,ndays,avg_waa,vhs_trans_rate_I,vhs_trans_rate_C,vhs_mort_rate,vhs_rec_rate,dep_scaling,
                fishing_mort,survey_selA50,survey_selA95,maturity_A50,maturity_A95,disease_vulA50,disease_vulA95,dis_survey_selA50,dis_survey_selA95,
                nonlinear_exp_1,nonlinear_exp_2,sig_nat_mor,selA50,selA95,ignore.carryover.inf,inf_prev_survey,log_sigma_R,rseed)
proc.time()-ptm
# Model checks - make sure dynamics from model are sensible

list2env(run.1,globalenv())

# All years - divide by Post_out_Nya from previous year
# Infection intensity by year
infection.intensity <- (apply(Nya_new_infect[,-nage],1,sum)/apply(Nya[,-nage],1,sum))
# Disease mortality by year
disease.mortality <- (apply(Nya_dis_death[,-nage],1,sum)/apply(Nya[,-nage],1,sum))
# Natural survival by year
# total.survival <- (rowSums(Nya[,-nage]*Sya[,-nage]) + Nya[,nage]*Sya[,nage])/(rowSums(Nya[,-nage]) + Nya[,nage])
# total.survival <- (rowSums(Nya[,-nage]*Sya[,-nage]))/(rowSums(Nya[,-nage]))
total.survival <- (rowSums(Nya*Sya))/(rowSums(Nya))

# (Nya[,nage]*Sya[,nage])/Nya[,nage]
# Immunity (Antibody prevalence)
# immunity <- (N_y_recovered[,ndays]/apply(Post_out_Nya,1,sum))
immunity <- (apply(Nya_imm[,-nage],1,sum)/apply(Nya[,-nage],1,sum))

N_y_daily <- N_y_susceptible + N_y_infected + N_y_recovered

##################################
# Daily and annual dynamics
##################################
# years = (nyr-499):nyr
years = (nyr-49):nyr
# years = (151):200

fig.name <- paste0(here::here("results/figures/"),"Figure_Example 50 year dynamics from seed ",round(rseed),".tiff")
tiff(file=fig.name,width=9, height=6, units="in", res=200)

par(xpd=NA,oma=c(4,0,0,0))
layout(matrix(1:6,3,2,byrow=FALSE))

# Plot prevalence
par(mar=c(0,4,1,2))
matplot(t(N_y_susceptible/N_y_daily)[,years],type="l",lty=1,ylab="",col=colorRampPalette(c("grey20","grey90"))(50), xaxt='n', yaxt='n', lwd=2,ylim=c(0,1))
axis(side =1, seq(0,ndays,by=ndays/3), labels = F)
axis(2, at = c(0,0.5,1),labels=c(0,0.5,1))
text(0,0.95,labels='A',cex=1.6,adj=c(0,0.5),font=2)
title(ylab="Susceptible prevalence", line=2.5, cex.lab=1.4)

matplot(t(N_y_infected/N_y_daily)[,years],type="l",lty=1,ylab="",col=colorRampPalette(c("grey20","grey90"))(50), xaxt='n', yaxt='n',lwd=2,ylim=c(0,1))
axis(side =1, seq(0,ndays,by=ndays/3), labels = F)
axis(2, at = c(0,0.5,1),labels=c(0,0.5,1))
text(0,0.95,labels='B',cex=1.6,adj=c(0,0.5),font=2)
title(ylab="Infection prevalence", line=2.5, cex.lab=1.4)

par(mar=c(0,4,1,2))
matplot(t(N_y_recovered/N_y_daily)[,years],type="l",lty=1,ylab="", xlab="",col=colorRampPalette(c("grey20","grey90"))(50),lwd=2,ylim=c(0,1), xaxt='n', yaxt='n')
axis(side =1, seq(0,ndays,by=ndays/3), labels = T)
axis(2, at = c(0,0.5,1),labels=c(0,0.5,1))
text(0,0.95,labels='C',cex=1.6,adj=c(0,0.5),font=2)
title(ylab="Seroprevalence", line=2.5, cex.lab=1.4)
title(xlab="Days", line=2.5, cex.lab=1.4)

# Plot muti-plot of spawning biomass, recruitment, and outbreak severity in the last 50 years
# par(mfrow=c(3,1))
plot(apply(SSB_atage,1,sum)[years],type="l",xlab="",ylab="",
     yaxt="n",
     xaxt = "n",
     ylim=c(0,signif(ceiling(max(apply(SSB_atage,1,sum)[years])),digits=2)),lwd=2)
points(apply(SSB_atage,1,sum)[years],pch=21,bg=colorRampPalette(c("grey20","grey90"))(50),cex=1.5)
axis(2, at = seq(0, signif(ceiling(max(apply(SSB_atage,1,sum)[years])),digits=1), length.out = 3),
     labels=seq(0, signif(ceiling(max(apply(SSB_atage,1,sum)[years])),digits=1), length.out = 3)/1000)
axis(side =1, seq(0,length(years),by=25), labels = F)
title(ylab="Spawning biomass\n(1000 tons)", line=2.5, cex.lab=1.4)
text(0,0.95*max(apply(SSB_atage,1,sum)[years]),labels='D',cex=1.6,adj=c(0,0.5),font=2)

plot(Nya[,1][years],type="l",xlab="",ylab="",ylim=c(0,signif(ceiling(max(Nya[,1][years])),digits=2)), xaxt='n',lwd=2)
points(Nya[,1][years],pch=21,bg=colorRampPalette(c("grey20","grey90"))(50),cex=1.5)
axis(side =1, seq(0,length(years),by=25), labels = F)
title(ylab="Recruits (millions)", line=2.5, cex.lab=1.4)
text(0,0.95*max(Nya[,1][years]),labels='E',cex=1.6,adj=c(0,0.5),font=2)

plot(immunity[years],type="l",xlab="",ylab="",ylim=c(0,1.1),lwd=3,lty="solid", xaxt='n',yaxt='n',col="purple4")
lines(infection.intensity[years],lwd=4.5,col="darkorange2")
lines(total.survival[years],lwd=3, col="seagreen3")
axis(side =2, c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100))
axis(side =1, seq(0,length(years),by=25), labels = c(50, 75, 100))
title(ylab="Percentages", line=2.5, cex.lab=1.4)
title(xlab="Years", line=2.5, cex.lab=1.4)
text(0,0.95*1.1,labels='F',cex=1.6,adj=c(0,0.5),font=2)

legend(28,1.15,
       c('% survived','% with antibodies','% infected'),
       #horiz = TRUE,
       #cex=1.2,
       lwd=2.5,
       inset=c(0,0),
       x.intersp=0.3,
       y.intersp = 0.9,
       # col=c('grey40','grey70','black'),
       col=c("seagreen3","purple4","darkorange2"),
       #lty = c("twodash","solid","dotted"),
       box.lty=0, bg=NA)

dev.off()

##################################
# Cross-correlations
##################################

custom.ccf <- function(x,y,lags,...){
  rho <- vector(length=2*lags+1)
  sig.level <- vector(length=2*lags+1)
  k <- 1
  for(i in -lags:lags){
    indexing <- (1:length(x))+i # reference indices
    indexing <- ifelse(indexing<=0 | is.na(indexing),NA,indexing)
    block <- data.frame(x=x[indexing],y=y)
    block <- block[!is.na(block$x),]
    if(NROW(block)<=5){
      rho[k] <- NA
      sig.level[k] <- NA
    }else{
      temporary <- cor.test(x=as.ts(block$x),y=as.ts(block$y),method="spearman",na.action=na.pass)
      rho[k] <- temporary$estimate
      sig.level[k] <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(block$x)))
    }
    k <- k+1
  }
  
  output <- data.frame(lags=-lags:lags,
                       correlations=rho,
                       significance.level=sig.level)
  
  return(output)
}

years.2 <- 100:nyr
# cc.inf.vs.rec <- ccf(Nya[,1][years.2],infection.intensity[years.2],plot=FALSE,lag.max=10)
# cc.inf.vs.rec <- custom.ccf(x=(Nya[,1]/apply(Nya,1,sum))[years.2],y=infection.intensity[years.2],lags=10)
cc.inf.vs.rec <- custom.ccf(x=(Nya[,1])[years.2],y=infection.intensity[years.2],lags=10)
cc.imm.vs.sur <- custom.ccf(x=immunity[years.2],y=total.survival[years.2],lags=10)
cc.sur.vs.inf <- custom.ccf(x=infection.intensity[years.2],y=total.survival[years.2],lags=10)
cc.imm.vs.inf <- custom.ccf(x=infection.intensity[years.2],y=immunity[years.2],lags=10)

cc.rec.vs.rec <- custom.ccf(x=(Nya[,1])[years.2],y=(Nya[,1])[years.2],lags=10)

fig.name <- paste0(here::here("results/figures/"),"Figure_Cross correlations from seed ",round(rseed),".tiff")
tiff(file=fig.name,width=8, height=6, units="in", res=200)

par(xpd=NA,mfrow=c(2,2))
par(oma=c(3,3,1,1),mar=c(0,0,0,0))
plot(x=cc.inf.vs.rec$lags,y=cc.inf.vs.rec$correlations,ylab="",
     xlab="",t="n",main=NA,ylim=c(-1,1),cex.lab=1,cex.axis=0.9,xaxt='n')
polygon(x=c(cc.inf.vs.rec$lags,rev(cc.inf.vs.rec$lags)),y=c(cc.inf.vs.rec$significance.level,-rev(cc.inf.vs.rec$significance.level)),col="grey80",border=NA)
segments(x0=-10,x1=10,y0=0)
segments(x0=cc.inf.vs.rec$lags,y0=0,y1=cc.inf.vs.rec$correlations,lwd=5,lend=1)
text(-10, 0.9, label="A. Infection incidence vs Recruits (lagged)",adj=0,cex=1,font=2)
# points(x=cc.inf.vs.rec$lags,y=cc.inf.vs.rec$correlations,pch=21,bg="white",cex=1.2)
# lines(x=cc.inf.vs.rec$lags,y=cc.inf.vs.rec$significance.level,col="blue",lty="dashed")
# lines(x=cc.inf.vs.rec$lags,y=-cc.inf.vs.rec$significance.level,col="blue",lty="dashed")

plot(x=cc.sur.vs.inf$lags,y=cc.sur.vs.inf$correlations,ylab="",xlab="",
     t="n",main=NA,ylim=c(-1,1),yaxt='n',xaxt='n', cex.lab=1,cex.axis=0.9)
polygon(x=c(cc.sur.vs.inf$lags,rev(cc.sur.vs.inf$lags)),y=c(cc.sur.vs.inf$significance.level,-rev(cc.sur.vs.inf$significance.level)),col="grey80",border=NA)
segments(x0=-10,x1=10,y0=0)
segments(x0=cc.imm.vs.sur$lags,y0=0,y1=cc.imm.vs.sur$correlations,lwd=5,lend=1)
text(-10, 0.9, label="B. Survival vs Infection incidence (lagged)",adj=0,cex=1,font=2)

plot(x=cc.imm.vs.inf$lags,y=cc.imm.vs.inf$correlations,ylab="",xlab="",
     t="n",main=NA,ylim=c(-1,1),cex.lab=1,cex.axis=0.9)
polygon(x=c(cc.imm.vs.inf$lags,rev(cc.imm.vs.inf$lags)),y=c(cc.imm.vs.inf$significance.level,-rev(cc.imm.vs.inf$significance.level)),col="grey80",border=NA)
segments(x0=-10,x1=10,y0=0)
segments(x0=cc.imm.vs.inf$lags,y0=0,y1=cc.imm.vs.inf$correlations,lwd=5,lend=1)
text(-10, 0.9, label="C. Seroprevalence vs Infection incidence (lagged)",adj=0,cex=1,font=2)

plot(x=cc.imm.vs.sur$lags,y=cc.imm.vs.sur$correlations,ylab="",xlab="",
     t="n",main=NA,ylim=c(-1,1),yaxt='n',cex.lab=1,cex.axis=0.9)
polygon(x=c(cc.imm.vs.sur$lags,rev(cc.imm.vs.sur$lags)),y=c(cc.imm.vs.sur$significance.level,-rev(cc.imm.vs.sur$significance.level)),col="grey80",border=NA)
segments(x0=-10,x1=10,y0=0)
segments(x0=cc.imm.vs.sur$lags,y0=0,y1=cc.imm.vs.sur$correlations,lwd=5,lend=1)
text(-10, 0.9, label="D. Survival vs Seroprevalence (lagged)",adj=0,cex=1,font=2)
# points(x=cc.imm.vs.sur$lags,y=cc.imm.vs.sur$correlations,pch=21,bg="white",cex=1.2)
# lines(x=cc.imm.vs.sur$lags,y=cc.imm.vs.sur$significance.level,col="blue",lty="dashed")
# lines(x=cc.imm.vs.sur$lags,y=-cc.imm.vs.sur$significance.level,col="blue",lty="dashed")

title(ylab="Cross correlation",xlab="Lags (years)",outer=TRUE,line=2)
dev.off()

#plot(lag(immunity[years.2],7),total.survival[years.2])
#plot(rank(lag(immunity[years.2],7)),rank(total.survival[years.2]))


##################################
# Age-specific mortality over time
##################################
library(dplyr)
library(ggplot2)
library(reshape2)

# Melt matrices (nyr x nage) of survival, mortality rates, and immunity
age_survival <- reshape2::melt(Sya,varnames=c("Year","Age"),value.name="Survival")
add_mort <- t(apply(Sya,1,function(x) -log(x)))
add_mort <- reshape2::melt(add_mort,varnames=c("Year","Age"),value.name="Mortality")
age_immune <- reshape2::melt(predicted_immune_comps,varnames=c("Year","Age"),value.name="Seroprevalence")
age_comp <- t(apply(Nya,1,function(x) x/sum(x)))
age_comp <- reshape2::melt(age_comp,varnames=c("Year","Age"),value.name="Age_Comp")
infection_incidence <- reshape2::melt(Nya_new_infect/Nya,varnames=c("Year","Age"),value.name="Infection incidence")

# Combine everything
mort_vs_immune <- full_join(full_join(age_survival,add_mort,by=c("Year","Age")),age_immune,by=c("Year","Age"))
mort_vs_immune <- full_join(mort_vs_immune,age_comp,by=c("Year","Age"))
mort_vs_immune <- full_join(mort_vs_immune,infection_incidence,by=c("Year","Age"))

mort_vs_immune$Age <- mort_vs_immune$Age-1 # Correct to actual ages

# Omit first 150 years
years.4.analysis <- 50:74

# Now melt with just immune, survival, and infection probabilities
# mort_vs_immune <- reshape2::melt(select(mort_vs_immune,-Age_Comp,-Mortality),id.vars=c("Year","Age"))
dat.2.plot <- reshape2::melt(select(mort_vs_immune,-Age_Comp,-Mortality,-'Seroprevalence',-'Infection incidence'),id.vars=c("Year","Age"))

# This is the immunity following this year's mortality
baseline<- data.frame(Year=years.4.analysis, value=exp(-0.25))

dat.2.plot <- filter(dat.2.plot,Year %in% years.4.analysis) %>% mutate(Year=Year,
                                                                           Age=Age,
                                                                           Age_lwd=(7-Age)*0.1+0.2)
font.size <- 14

A.plot <- ggplot(data=dat.2.plot) +
  geom_line(aes(x=Year,y=value,group=interaction(Age,variable),colour=factor(Age),size=factor(Age)))+
  #geom_line(data=baseline,aes(x=Year,y=value),color="black",linetype="dashed",lwd=2)+
  ylab("Survival") + 
  scale_size_manual(breaks=0:7,values=seq(2.25,0.5,length.out=8),
                    labels=c(0:6,"7+"))+
  scale_color_manual(values = c("olivedrab2","lightgoldenrod", "orange", "firebrick2", "darkmagenta", "blue", "cyan3","black"),
                     labels=c(0:6,"7+")) +
  theme_classic()+
  guides(colour = guide_legend(nrow = 1), size=guide_legend(nrow=1)) +
  #scale_y_continuous(limits=c(0,1))+
  #scale_y_continuous(limits=c(0.25,1),breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  scale_y_continuous(limits=c(0.5,0.9))+
  theme(#axis.title.y.left = element_text(color = "darkblue"),
    #axis.title.y.right = element_text(color = "darkorange"),
    #axis.text.y.left = element_text(color = "darkblue"),
    #axis.text.y.right = element_text(color = "orange"),
    #axis.line.y.right = element_line(color = "orange"),
    #axis.ticks.y.right = element_line(color = "orange"),
    axis.title.x = element_text(size=font.size),
    axis.title.y = element_text(size=font.size),
    axis.text.y = element_text(size=font.size-4,hjust=0.5,angle=90),
    axis.text.x = element_text(size=font.size-4),
    text=element_text(size=15),
    legend.position=c(0.05,0.1),
    legend.title = element_blank(),
    legend.justification = c(0,0.5),
    legend.text=element_text(size=font.size-3),
    legend.background = element_blank(),
    panel.border = element_rect(color="black", fill=NA))
# dev.off()

##################################
# Cohort effects
##################################

# Create indices for cohorts
cohort.indices <- data.frame(Cohort=1,Age=0:(nage-1),Year=1:nage)
for(y in 2:nyr){
  cohort.indices <- rbind(cohort.indices,data.frame(Cohort=y,Age=0:(nage-1),Year=y:(y-1+nage)))
}

# Label as cohorts
mort_vs_immune_2 <- mort_vs_immune %>% mutate(Cohort=cohort.indices$Cohort[prodlim::row.match(select(mort_vs_immune,Year,Age),select(cohort.indices,Year,Age))])

# Now melt with just immune, survival, and infection probabilities
mort_vs_immune_2 <- reshape2::melt(select(mort_vs_immune_2,-Age_Comp,-Mortality),id.vars=c("Cohort","Year","Age"))

# This is the immunity following this year's mortality
baseline <- data.frame(Age=0:(nage-1),value=c(rep(exp(-0.25),times=nage-1),exp(-1.8)))

dat.2.plot <- filter(mort_vs_immune_2,Cohort %in% years.4.analysis[1:9]) %>% mutate(Cohort=Cohort)

B.plot <- ggplot(data=dat.2.plot) +
  geom_point(aes(x=Age,y=value,group=Cohort,colour=variable,shape=variable),size=2)+
  geom_line(aes(x=Age,y=value,group=interaction(Cohort,variable),colour=variable),size=1)+
  #geom_line(data=baseline,aes(x=Age,y=value),color="black",linetype="dashed",size=1)+
  facet_wrap(.~Cohort,dir="v")+
  geom_text(aes(x=0,y=1,label=paste0("Cohort ",Cohort)),fontface="plain",hjust=0)+
  ylab("Proportions") + 
  # scale_colour_manual(values=c("orange","blue","red")) +
  scale_colour_manual(values=c("seagreen3","purple4","darkorange2")) +
  scale_y_continuous(limits=c(0,1.05),breaks=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1))+
  theme_classic()+
  theme(axis.title.x = element_text(size=font.size),
        axis.title.y = element_text(size=font.size),
        axis.text.y = element_text(size=font.size-4,hjust=0.5,angle=90),
        axis.text.x = element_text(size=font.size-4),
        legend.text = element_text(size=font.size-2),
        legend.position="bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.margin = margin(0,0,0,0,unit="mm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        panel.spacing = unit(0,"lines"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(color="black", fill=NA))
#dev.off()


##################################
# Plot antibody observations
##################################
comp_samp_size = 50
antibody_comp_samp_size = 200
survey_cv = 0.3
obs_years = 1:nyr
observations <- fun_obs_mod(obs_years,Nya,comp_samp_size,catch_comp_samp_size,
                            predicted_comps,survey_slx,dis_survey_slx,
                            predicted_immune_comps,antibody_comp_samp_size,
                            predicted_survey,survey_cv,Cya,N_catch,
                            true_samp_prev,
                            rseed)

antibodies <- observations$antibody_obs

# Calculate age-specific serpositive from overall sample
seropositive <- t(apply(antibodies,1,function(x) x[seq(1,2*nage,by=2)]/sum(x)))
colnames(seropositive) <- paste0(0:(nage-1))

# Calculate age composition of overall sample
samp_comp <- t(apply(antibodies,1,function(x) (x[seq(1,2*nage,by=2)]+x[seq(2,2*nage,by=2)])/sum(x)))
colnames(samp_comp) <- paste0(0:(nage-1))

# Melt these objects and joiong together
seropositive <- melt(data.frame(year=1:nyr,
                                seropositive),id='year',variable.name='age',value.name='seropositive')
samp_comp <- melt(data.frame(year=1:nyr,
                             samp_comp),id='year',variable.name='age',value.name='age_comp')

final_antibodies <- inner_join(samp_comp,seropositive)
final_antibodies <- filter(final_antibodies,!age%in%c('Age.0','Age.1'))

final_antibodies$age <- (0:7)[match(final_antibodies$age,unique(final_antibodies$age))]

dat.2.plot <- filter(final_antibodies,year %in% years.4.analysis) %>% mutate(year=year)
dat.2.plot <- melt(dat.2.plot,id = c('year','age'))

# Now plot
C.plot <- ggplot(data=dat.2.plot) + 
  geom_col(aes(x=age, y=value, fill=variable, group=variable),position = 'identity')+
  # geom_col(aes(x=age, y=seropositive,group=year),fill="grey75")+
  # facet_wrap(year~.,dir="v",nrow=10)+
  scale_fill_manual(values = c('age_comp'="black",'seropositive'='grey75'),
                    labels = c('- antibodies','+ antibodies')) + 
  facet_wrap(year~.,dir="v",nrow=13)+
  geom_text(aes(x=0,y=0.95,label=paste0("Year ",year)),fontface="plain",hjust=0)+
  labs(y="Sample seroprevalence (size=200)")+
  xlab("Age")+
  theme_classic()+
  scale_y_continuous(limits=c(0,1.15),breaks=c(0,0.5,1),labels=c(0,0.5,1))+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border=element_rect(fill=NA),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size=font.size),
        axis.title.x = element_text(size=font.size),
        axis.text.y = element_text(size=font.size-5,hjust=0.5,angle=90),
        axis.text.x = element_text(size=font.size-4),
        axis.ticks.x= element_line(color="black"),
        #axis.title.y = element_blank(),
        #legend.position ="none",
        legend.title = element_blank(),
        legend.position=c(0.55,0.01),
        legend.justification = c(0,0.5),
        legend.text=element_text(size=font.size-3),
        legend.background = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.spacing = unit(0, "lines"))

# ggsave(filename=here::here(paste0("results/figures/Figure_example seroprevalence seed ",rseed,".png")),width=8, height=6, units="in",dpi=600)

##################################
# Combine plots into one giant one
##################################

profile.plot.1 <- cowplot::ggdraw() +
  cowplot::draw_plot(A.plot, 0, 2/3, 2/3, 1/3) +
  cowplot::draw_plot(B.plot, 0, 0, 2/3, 2/3) +
  cowplot::draw_plot(C.plot, 2/3, 0, 1/3, 1) +
  cowplot::draw_plot_label(c("A", "B", "C"), c(0, 0, 2/3), c(1, 2/3, 1), size = 18)
ggsave(filename=here::here(paste0("results/figures/Figure_age and cohort profiles_seed ",rseed,".png")),
       plot=profile.plot.1,width=10, height=8, units="in",dpi=600)  








# EXPLORATORY
##################################
# Peak prevalance VS Max infection probability (Beta_I)
##################################

sel_years <- 100:nyr

plot((apply(beta_I_y,1,max)+apply(beta_C_y,1,max))[sel_years],max.prevalence[sel_years],ylab="Prevalence",xlab="Infection probability")
abline(v=((1-exp(-vhs_mort_rate))+(1-exp(-vhs_rec_rate))),lty="dashed")

plot(apply(beta_C_y,1,max),max.prevalence,ylab="Prevalence",xlab="Infection probability")

plot(apply(beta_I_y,1,max),apply(beta_C_y,1,max),ylab="Infection probability form Infected",xlab="Infection probability from Carrier")

# Frequency of times infections exceed mortality & recovery
outbreaks <- (apply(beta_I_y,1,max) + apply(beta_C_y,1,max)) > ((1-exp(-vhs_mort_rate))+(1-exp(-vhs_rec_rate)))
# outbreaks <- max.prevalence >=0.2
sum(outbreaks[sel_years])/length(sel_years)
median(diff(which(outbreaks[sel_years])))

# plot years following this definition
matplot(t(N_y_infected/N_y_daily)[,outbreaks],type="l",lty=1,ylab="",col=colorRampPalette(c("grey20","grey90"))(sum(outbreaks)), xaxt='n', yaxt='n',lwd=2,ylim=c(0,1))
axis(side =1, seq(0,ndays,by=50), labels = F)
axis(2, at = c(0,0.5,1),labels=c(0,0.5,1))
text(0,0.95,labels='b)',cex=1.4,adj=c(0,0.5))
title(ylab="Infection prevalence", line=2.5, cex.lab=1.4)

##################################
# R0
##################################

# Calculate R0 based on Final Epidemic Size
S0 <- (N_y_susceptible/N_y_daily)[,1]
Sinf <- (N_y_susceptible/N_y_daily)[,ndays]
R0 <- (log(S0)-log(Sinf))/(S0-Sinf)
plot(density(R0))
abline(v=1)

# Calculate R0 form SIR Flows
# Beta = mean transmission rate, gamma = recovery rate, and mu = nautral mortality rate
# R0 = beta/(gamma+mu)
vhs_trans_rate_I=0.01
vhs_trans_rate_C=0.000001

dis_mort_rate <- 21/365 # 11/365
rec_rate <- 20/365 # 26/365
nat_mort_rate <- 0.25/365

R0 <- (vhs_trans_rate_I + vhs_trans_rate_C)/(dis_mort_rate + rec_rate + nat_mort_rate)


##################################
# Age-specific mortality resulting from disease over time
##################################
library(reshape2)
library(ggplot2)
library(ggridges)
#add_mort <- t(apply(Sya,1,function(x) -log(x)-c(rep(0.25,times=nage-1),1.8)))
add_mort <- t(apply(Sya[(nyr-50):nyr,],1,function(x) -log(x)))
add_mort <- melt(add_mort,varnames=c("Year","Age"))

# One with ridge lines
ggplot(add_mort,aes(x=Age,y=Year,height=value, group=Year)) +
  geom_ridgeline()+
  theme_classic()

# One with facets
ggplot(add_mort,aes(x=Year,y=value,group=Age)) +
  geom_line()+
  facet_wrap(.~Age)+
  theme_classic()


##################################
# Age-specific mortality VS age-specific immunity
##################################
library(dplyr)
age_survival <- melt(Sya[(nyr-50):nyr,],varnames=c("Year","Age"),value.name="Survival")

add_mort <- t(apply(Sya[(nyr-50):nyr,],1,function(x) -log(x)))
add_mort <- melt(add_mort,varnames=c("Year","Age"),value.name="Mortality")

age_immune <- melt(predicted_immune_comps[(nyr-50):nyr,],varnames=c("Year","Age"),value.name="Immunity")

mort_vs_immune <- full_join(full_join(age_survival,add_mort,by=c("Year","Age")),age_immune,by=c("Year","Age"))

# Arrange by year first, then age
mort_vs_immune <- arrange(mort_vs_immune,Year,Age) %>% 
  mutate(Year_lag1=Year-1,
         Age_lag1=Age-1)

mort_vs_immune <- mort_vs_immune %>% mutate(Immunity_lead1=Immunity[prodlim::row.match(select(mort_vs_immune,Year,Age),select(mort_vs_immune,Year_lag1,Age_lag1))])

# This is the immunity preceding this year's mortality
ggplot(mort_vs_immune,aes(x=Immunity,y=Survival,group=Age)) +
  geom_point()+
  geom_line()+
  facet_wrap(.~Age)+
  theme_classic()

# This is the immunity following this year's mortality
ggplot(mort_vs_immune,aes(x=Mortality,y=Immunity_lead1,group=Age)) +
  geom_point()+
  geom_line()+
  facet_wrap(.~Age)+
  theme_classic()


##################################
# Age-structure VS age-immunity
##################################
library(dplyr)
age_comp <-t(apply(Nya[(nyr-50):nyr,],1,function(x) x/sum(x)))
age_comp <- melt(age_comp,varnames=c("Year","Age"),value.name="Age_structure")

age_immune <- melt(predicted_immune_comps[(nyr-50):nyr,],varnames=c("Year","Age"),value.name="Immunity")

structure_vs_immune <- full_join(age_comp,age_immune,by=c("Year","Age"))
trial <- structure_vs_immune %>% group_by(Age) %>% mutate(lagged_structure=lag(Age_structure))

# One with facets
ggplot(structure_vs_immune,aes(x=lagged_structure,y=Immunity,group=Age)) +
  geom_point()+
  geom_line()+
  facet_wrap(.~Age)+
  theme_classic()
