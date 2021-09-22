// vhs_asa_em.cpp
// COMMON ERRORS
// Check indexing (especially for loops): starts at 0 and ends at <N; WILL CRASH IF INDEXING IS OFF!!!
// int vs Type specification - either leads to different calculated values in some cases (e.g. int n=2 then pow(x,n) is different than Type n=2 with pow(x,n))
// When building model with dummy variable, do NOT map it to NA, it will abort (just map all other variables to NA)
// Initialize all derived variables with .setZero() - not doing so will cause GREAT agony (e.g. weird values that propagate to nonmatching results or NaNs)
#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template <class Type> 
Type posfun(Type x,Type eps,Type &pen){
   pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
   return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER(nyr)
  DATA_INTEGER(sage)  // Only use sage for calculations, NOT indexing
  DATA_INTEGER(nage)
  DATA_INTEGER(comp_samp_size)
  DATA_INTEGER(catch_comp_samp_size)
  DATA_INTEGER(antibody_comp_samp_size)
  DATA_VECTOR(survey_obs)
  DATA_VECTOR(catches)
  DATA_MATRIX(catch_comps)
  DATA_MATRIX(comp_obs)
  DATA_MATRIX(antibody_obs)
  DATA_VECTOR(avg_waa)
  DATA_SCALAR(maturity_A50)
  DATA_SCALAR(maturity_A95)
  DATA_VECTOR(obs_samp_prev)

  DATA_INTEGER(dis_mix_age_thresh)
  DATA_INTEGER(fix_rec_rate)
  DATA_INTEGER(fix_dis_survey_slx)
  DATA_SCALAR(eps)

  // std::cout << "data loaded" << "\n";

  // Parameters
  PARAMETER(dummy);
  PARAMETER(ac_coef_rec);
  PARAMETER(natural_mortality);
  PARAMETER(log_Ninit);
  PARAMETER(log_rbar);
  PARAMETER(log_q_survey);
  PARAMETER(plus_group_mortality);
  PARAMETER(log_SD_survey);
  PARAMETER(log_sigma_R);
  PARAMETER(log_sigma_Ninit);
  PARAMETER(selA50);
  PARAMETER(selA95);
  PARAMETER(survey_selA50);
  PARAMETER(survey_selA95);
  PARAMETER(disease_selA50);
  PARAMETER(disease_selA95);
  PARAMETER(dis_survey_selA50);
  PARAMETER(dis_survey_selA95);
  PARAMETER(fix_recov_par);
  PARAMETER_VECTOR(log_Ninit_devs);
  PARAMETER_VECTOR(log_rbar_devs);
  PARAMETER_VECTOR(init_immune);
  PARAMETER_VECTOR(tran_infec_rate);
  
  PARAMETER(beta_prev_index);
  

  Type nll_survey=0.0;
  Type nll_fish_comps_total=0.0;
  Type nll_comps_total=0.0;
  Type pen_nll_devs=0.0;
  Type antibody_nll=0.0;
  
  vector<Type> N_catch(nyr);
  vector<Type> predicted_survey(nyr);
  vector<Type> nll_comps_2(nyr);
  vector<Type> nll_fish_comps_2(nyr);
  vector<Type> nll_antibody_comps_2(nyr);
  vector<Type> resid_rec(nyr);
  vector<Type> mean_rec(nyr);
  vector<Type> SSB(nyr);
  vector<Type> total_prop_imm(nyr);
  vector<Type> u(nyr);
  vector<Type> vulnerable(nage);
  vector<Type> infec_rate(nyr);
  vector<Type> recov_rate(nyr);
  
  matrix<Type> mature(nyr,nage);
  matrix<Type> Sya(nyr,nage);
  matrix<Type> Vya(nyr,nage);
  matrix<Type> Cya(nyr,nage);
  matrix<Type> Nya(nyr,nage);
  matrix<Type> PostNya(nyr,nage);
  matrix<Type> spawners(nyr,nage);
  matrix<Type> predicted_sur_comps(nyr,nage);
  matrix<Type> predicted_comps(nyr,nage);
  matrix<Type> predicted_antibody(nyr,2*nage);
  matrix<Type> nll_fish_comps(nyr,nage);
  matrix<Type> nll_comps(nyr,nage);
  matrix<Type> nll_antibody_comps(nyr,2*nage);
  matrix<Type> dis_survey_slx(nyr,nage);
  matrix<Type> survey_slx(nyr,nage);
  matrix<Type> disease_vul(nyr,nage);
  matrix<Type> Slxya(nyr,nage);
  matrix<Type> SSB_atage(nyr,nage);

  matrix<Type> prop_sus(nyr,nage);
  matrix<Type> prop_imm(nyr,nage);
  matrix<Type> Nya_imm(nyr,nage);
  matrix<Type> Dis_Sya(nyr,nage);
  
  Type f_llk = 0.0;
  Type pi = 3.141593;

// Transform infection rate
  infec_rate.setZero();
  for (int i=0;i<nyr;i++){
    //infec_rate(i) = (atan(tran_infec_rate(i)*tran_infec_rate(i)*tran_infec_rate(i))+0.5*pi) / pi;
    infec_rate(i) = tran_infec_rate(i);
  }

// Fill in recovery rate using either the constant parameter (fixed rate across time) or parameter vector (varies across time)
  recov_rate.setZero();
  for (int i=0;i<nyr;i++){
    recov_rate(i) = fix_recov_par;
  }

// Calculate 'disease vulnerability'
// The proportion of each age class vulnerable to transmissions (i.e. well mixed with the infected population)  
  disease_vul.setZero();
  for (int i=0;i<nyr;i++){
    for (int j=0;j<nage;j++){
      disease_vul(i,j) = (j<dis_mix_age_thresh) ? 1.0 / (1.0 + exp(-log(19)*(sage+j-disease_selA50)/(disease_selA95-disease_selA50))) : 1.0;
      // disease_vul(i,j) = 1.0 / (1.0 + exp(-log(19)*(sage+j-disease_selA50)/(disease_selA95-disease_selA50)));
    }
  }

// Calculate Maturity
  mature.setZero();
  for (int i=0;i<nyr;i++){
    for (int j=2;j<nage;j++){
      mature(i,j) = 1.0 / (1.0 + exp(-log(19)*(sage+j-maturity_A50)/(maturity_A95-maturity_A50)));
    }
  }
  
// Calculate annual survival rate (constant across time and age)
  Sya.setZero();
  for(int i=0;i<nyr;i++){
    for (int j=0;j<(nage-1);j++){
      Sya(i,j) = (j>2 & j<5) ? exp(-(natural_mortality+beta_prev_index*obs_samp_prev[i])) : exp(-(natural_mortality));
    }
    Sya(i,nage-1) = exp(-plus_group_mortality);  
  }

// Calculate Selectivity
// Gear Selectivity
dis_survey_slx.setZero();
survey_slx.setZero();
Slxya.setZero();
for (int i=0;i<nyr;i++){
    for (int j=0;j<nage;j++){
      Slxya(i,j) = (j<5) ? 1.0 / (1.0 + exp(-log(19)*(sage+j-selA50)/(selA95-selA50))) : 1.0;
      // Slxya(i,j) = 1.0 / (1.0 + exp(-log(19)*(sage+j-selA50)/(selA95-selA50)));
    }
}
for (int i=0;i<nyr;i++){
    for (int j=0;j<nage;j++){
      survey_slx(i,j) = (j<5) ? 1.0 / (1.0 + exp(-log(19)*(sage+j-survey_selA50)/(survey_selA95-survey_selA50))) : 1.0;
      // survey_slx(i,j) = 1.0 / (1.0 + exp(-log(19)*(sage+j-survey_selA50)/(survey_selA95-survey_selA50)));
      if(fix_dis_survey_slx){
        dis_survey_slx(i,j) = survey_slx(i,j);
      }else{
        dis_survey_slx(i,j) = (j<5) ? 1.0 / (1.0 + exp(-log(19)*(sage+j-dis_survey_selA50)/(dis_survey_selA95-dis_survey_selA50))) : 1.0;
        // dis_survey_slx(i,j) = 1.0 / (1.0 + exp(-log(19)*(sage+j-dis_survey_selA50)/(dis_survey_selA95-dis_survey_selA50)));
      }
    }
}

// Calculate initial numbers
Nya.setZero();
Nya(0,0) = exp(log_rbar + log_rbar_devs(0) - 0.5 * square(exp(log_sigma_R)));
for (int j=1;j<nage;j++){
  Nya(0,j) = exp(log_Ninit + log_Ninit_devs(j-1) - 0.5 * square(exp(log_sigma_R)));
  Type pen=0.0;
  Nya(0,j)=posfun(Nya(0,j), eps, pen);
  f_llk+=1000000*pen;   
}

// Calculate initial susceptibles
prop_imm.setZero();
prop_sus.setZero();
prop_imm.row(0) = init_immune;
prop_sus.row(0) = 1.0 - init_immune;

// Now project forward
Dis_Sya.setZero();
Vya.setZero();
Cya.setZero();
PostNya.setZero();
predicted_comps.setZero();
u.setZero();
N_catch.setZero();
spawners.setZero();
SSB_atage.setZero();
SSB.setZero();
mean_rec.setZero();
predicted_survey.setZero();
predicted_sur_comps.setZero();
predicted_antibody.setZero();
total_prop_imm.setZero();

Type temp_var;
Type N_plus_temp;
int i_anti;

for (int i=0;i<nyr;i++){

  vulnerable.setZero();
  for(int j=0; j<nage; j++){
    predicted_survey(i) += Nya(i,j)*survey_slx(i,j)*avg_waa(j);
    predicted_sur_comps(i,j) = (Nya(i,j)*survey_slx(i,j));
    Nya_imm(i,j) = prop_imm(i,j)*Nya(i,j);
    vulnerable(j) = Nya(i,j)*Slxya(i,j);
  }
  predicted_survey(i)=exp(log_q_survey)*predicted_survey(i);  
  predicted_sur_comps.row(i) = predicted_sur_comps.row(i)/predicted_sur_comps.row(i).sum();
  predicted_comps.row(i)=vulnerable/vulnerable.sum();

  temp_var=0.0;
  for (int j=0;j<nage;j++){
    temp_var+=predicted_comps(i,j)*avg_waa(j);
  }
  
  // Antibody prevalence
  // First calculate total populatin immunity
  total_prop_imm(i) = Nya_imm.row(i).sum()/Nya.row(i).sum();
  i_anti=0;
  for(int j=0; j<nage; j++){
    // Calculate proportion of immune with antibodies in age j from total population 
    // This assumes the logistic 'scalar' for antibodies scales with the total proportion immune, not by age.
    // So to get the age-specfic antibodies, we first calculate total proportion of pop. with antibodies,
    // The proportion of age a with antibodies then is directly proportional to the ratio of #'s immune age a over the total #'s across all ages

    // predicted_antibody(i,i_anti) = (ab_logit_max/(1+exp(-ab_logit_k*(total_prop_imm(i)-ab_logit_midpoint))))*prop_imm(i,j)*Nya(i,j)/Nya.row(i).sum();
    predicted_antibody(i,i_anti) = prop_imm(i,j)*Nya(i,j)*dis_survey_slx(i,j);
    i_anti+=1;
    
    // Calculate proportion of fish WITHOUT antibodies in age j from total population
    // predicted_antibody(i,i_anti) = (1-ab_logit_max/(1+exp(-ab_logit_k*(total_prop_imm(i)-ab_logit_midpoint)))*prop_imm(i,j))*Nya(i,j)/Nya.row(i).sum();
    predicted_antibody(i,i_anti) = (1-prop_imm(i,j))*Nya(i,j)*dis_survey_slx(i,j);
    i_anti+=1;
  }

  predicted_antibody.row(i) = predicted_antibody.row(i)/predicted_antibody.row(i).sum();

  N_catch(i) = catches(i)/temp_var; 
  Cya.row(i) = N_catch(i)*predicted_comps.row(i);
  u(i) = Cya.row(i).sum()/vulnerable.sum();

  PostNya.row(i)=Nya.row(i)-Cya.row(i);
  for(int j=0; j<nage; j++){
    Type pen=0.0;
    PostNya(i,j)=posfun(PostNya(i,j), eps, pen);
    f_llk+=1000000*pen; 

    spawners(i,j)=mature(i,j)*PostNya(i,j);
  } 

  for(int j=0; j<nage; j++){
    SSB_atage(i,j)=spawners(i,j)*avg_waa(j);
  }
  SSB(i)=SSB_atage.row(i).sum();

  mean_rec(i) = exp(log_rbar); // Stationary recruitment - independent of biomass

  if(i<(nyr-1)){
    // Autocorrelated, density-INDEPENDENT recruitment
    Nya(i+1,0) = mean_rec(i)*exp(ac_coef_rec*(log(Nya(i,0))-log_rbar+0.5*square(exp(log_sigma_R))) + sqrt(1-square(ac_coef_rec))*log_rbar_devs(i+1)-0.5*square(exp(log_sigma_R)));
    
    // Project proportions immune and susceptible in the next year
    prop_imm(i+1,0) = 0;
    prop_sus(i+1,0) = 1-prop_imm(i+1,0);

    // Project post fishing numbers to next year undergoing disease and natural survival
    for(int j=1; j<nage; j++){

      Dis_Sya(i,j-1) = ((1-prop_sus(i,j-1)*disease_vul(i,j-1)*infec_rate(i)) + (prop_sus(i,j-1)*disease_vul(i,j-1)*infec_rate(i)*recov_rate(i)));

      Nya(i+1,j) = PostNya(i,j-1) * Sya(i,j-1) * Dis_Sya(i,j-1);

      prop_imm(i+1,j) = (prop_imm(i,j-1) + prop_sus(i,j-1)*disease_vul(i,j-1)*infec_rate(i)*recov_rate(i))/(1-prop_sus(i,j-1)*disease_vul(i,j-1)*infec_rate(i)*(1-recov_rate(i)));
      prop_sus(i+1,j) = 1-prop_imm(i+1,j);
    }

    // Those in plus group that survive to next year (stay in plus group though)
    // Dis_Sya(i,nage-1) = ((1-prop_sus(i,nage-1)*disease_vul(i,nage-1)*infec_rate(i)) + (prop_sus(i,nage-1)*disease_vul(i,nage-1)*infec_rate(i)*recov_rate(i)));

    N_plus_temp = PostNya(i,nage-1) * Sya(i,nage-1) * (1-prop_sus(i,nage-1)*disease_vul(i,nage-1)*infec_rate(i)*(1-recov_rate(i)));
    
    // Now combine the previous nage-1 group with nage (plus) group
    Nya(i+1,nage-1) += N_plus_temp;

    // Calculate plus-group proportions immune - combine those immune from preceding age group with those immune in current age group
    // prop_imm(i+1,nage-1) = (Nya(i+1,nage-1)*prop_imm(i+1,nage-1) + N_plus_temp*(prop_imm(i,nage-1) + prop_sus(i,nage-1)*disease_vul(i,nage-1)*infec_rate(i)*recov_rate(i))/((1-prop_sus(i,nage-1)*disease_vul(i,nage-1)*infec_rate(i)) + (prop_sus(i,nage-1)*disease_vul(i,nage-1)*infec_rate(i)*recov_rate(i))))/(Nya(i+1,nage-1)+N_plus_temp);
    prop_imm(i+1,nage-1) = (Sya(i,nage-2)*(prop_imm(i,nage-2)+disease_vul(i,nage-2)*prop_sus(i,nage-2)*infec_rate(i)*recov_rate(i))*PostNya(i,nage-2) + 
                            Sya(i,nage-1)*(prop_imm(i,nage-1)+disease_vul(i,nage-1)*prop_sus(i,nage-1)*infec_rate(i)*recov_rate(i))*PostNya(i,nage-1))/(Nya(i+1,nage-1));
    
    prop_sus(i+1,nage-1) = 1-prop_imm(i+1,nage-1);

    
    for(int j=0; j<nage; j++){
      Type pen=0.0;
      Nya(i+1,j)=posfun(Nya(i+1,j), eps, pen);
      f_llk+=1000000*pen; 
    }

  }
}

nll_fish_comps.setZero();
nll_fish_comps_2.setZero();
nll_comps.setZero();
nll_comps_2.setZero();
nll_antibody_comps.setZero();
nll_antibody_comps_2.setZero();
Type holder;
resid_rec.setZero();
 
for(int j=1; j<nage; j++){
  pen_nll_devs+=log_sigma_R+0.5*square(log_Ninit_devs(j-1)/exp(log_sigma_R));
} 

for(int i=0;i<nyr;i++){
  // Age comps from fishery independent survey
  for(int j=0; j<nage; j++){
  	if(predicted_sur_comps(i,j)==0 | comp_obs(i,j)==0){
  		nll_comps(i,j)=0;
  	}else{
  		nll_comps(i,j)=comp_obs.row(i).sum()*(comp_obs(i,j)/comp_obs.row(i).sum())*(log(predicted_sur_comps(i,j))-log(comp_obs(i,j)/comp_obs.row(i).sum()));
  	}
  }
  nll_comps_2(i)=nll_comps.row(i).sum();

  // Age comps from fishery 
  for(int j=0; j<nage; j++){
    if(catch_comps(i,j)<=0){
      nll_fish_comps(i,j)=0;
    }else{
      nll_fish_comps(i,j)=catch_comps.row(i).sum()*(catch_comps(i,j)/catch_comps.row(i).sum())*(log(predicted_comps(i,j))-log(catch_comps(i,j)/catch_comps.row(i).sum()));
    }
  }
  nll_fish_comps_2(i)=nll_fish_comps.row(i).sum();
}

for(int i=1;i<nyr;i++){
  // Antibody prevalence
  for(int j=0; j<(nage); j++){
   if((predicted_antibody(i,j*2)==0 & antibody_obs(i,j*2)==0) | (predicted_antibody(i,j*2+1)==0 & antibody_obs(i,j*2)==0)){
     nll_antibody_comps(i,j)=0;
   }else if((predicted_antibody(i,j*2)==0 & antibody_obs(i,j*2)>0) | (predicted_antibody(i,j*2+1)==0 & antibody_obs(i,j*2)>0)){
     nll_antibody_comps(i,j)=-1000;
   }else{
     // The Binomial
     nll_antibody_comps(i,j) = antibody_obs(i,j*2)*log( predicted_antibody(i,j*2)/(predicted_antibody(i,j*2)+predicted_antibody(i,j*2+1)) ) + 
                               (antibody_obs(i,j*2+1))*log( 1-predicted_antibody(i,j*2)/(predicted_antibody(i,j*2)+predicted_antibody(i,j*2+1)) );
     // nll_antibody_comps(i,j) = log(dbinom(antibody_obs(i,j*2),antibody_obs(i,j*2) + antibody_obs(i,j*2+1), predicted_antibody(i,j*2)/(predicted_antibody(i,j*2)+predicted_antibody(i,j*2+1))));
  
   }
  }

  // This was the former multinomial assumption
  // for(int j=0; j<(2*nage); j++){
  // 	// if(antibody_obs(i,j)==0 & predicted_antibody(i,j)==0){
  // 	//	nll_antibody_comps(i,j)=0;
  // 	//}else if(predicted_antibody(i,j)==0 & antibody_obs(i,j)>0){
  //   //  nll_antibody_comps(i,j)=-1000;
  //   // }else if(predicted_antibody(i,j)>0 & antibody_obs(i,j)==0){
  //   //  nll_antibody_comps(i,j)=0;
  //   if(predicted_antibody(i,j)==0 | antibody_obs(i,j)==0){
  //     nll_antibody_comps(i,j)=0;
  //   }else{
  //     nll_antibody_comps(i,j) = antibody_obs.row(i).sum()*(antibody_obs(i,j)/antibody_obs.row(i).sum())*(log(predicted_antibody(i,j))-log(antibody_obs(i,j)/antibody_obs.row(i).sum()));
  // 	}
  // }


  nll_antibody_comps_2(i)=nll_antibody_comps.row(i).sum();
}

for(int i=0;i<nyr;i++){
  holder=0.0;
  holder = (log(predicted_survey(i))-log(survey_obs(i)))/(exp(log_SD_survey));
  nll_survey += log_SD_survey+0.5*holder*holder;

  pen_nll_devs+=log_sigma_R+0.5*square(log_rbar_devs(i)/exp(log_sigma_R));

  //holder=0.0;
  //if(Nya(i,0)==0 | mean_rec(i)==0){
  //	resid_rec(i)=0;
  //}else{
  //	resid_rec(i) = log(Nya(i,0))-log(mean_rec(i));
  //}
}
  
  nll_fish_comps_total=-nll_fish_comps_2.sum();
  nll_comps_total=-nll_comps_2.sum();
  antibody_nll=-nll_antibody_comps_2.sum();


  // The last value is a penalty to the likelihood when Numbers go negative, or survival exceeds 100%
  f_llk = square(dummy) + nll_survey + nll_fish_comps_total + nll_comps_total + antibody_nll + pen_nll_devs;
  //f_llk = dummy*dummy;
  //obj_fun = dummy*dummy;
REPORT(f_llk);
REPORT(nll_survey);
REPORT(nll_fish_comps_total);
REPORT(nll_comps_total);
REPORT(antibody_nll);
REPORT(pen_nll_devs);
REPORT(nll_antibody_comps);

REPORT(infec_rate);
REPORT(recov_rate);
REPORT(mature);
REPORT(Dis_Sya);
REPORT(Sya);
REPORT(Nya);
REPORT(Cya);
REPORT(PostNya);
REPORT(prop_imm);
REPORT(prop_sus);
REPORT(SSB_atage);
REPORT(SSB);
REPORT(spawners)
REPORT(mean_rec);
vector<Type> recruits_obs = Nya.col(0);
REPORT(recruits_obs);
REPORT(predicted_comps);
REPORT(predicted_antibody);
REPORT(predicted_survey);
REPORT(disease_vul);
REPORT(survey_slx);
REPORT(dis_survey_slx);
REPORT(total_prop_imm);
REPORT(u);

return(f_llk);
}
