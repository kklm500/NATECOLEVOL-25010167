
library(foreach)
library(doParallel)
library(dplyr)
library(data.table)
#Details refer to Method in Article
#the net biodiversity effect on community variability (NBECV) and its multiple components as functions of species richness in LV model
#Parallel processing
detectCores()
stopImplicitCluster()
cl=makeCluster(10)
registerDoParallel(cl)

# Set repetitions
reps=10000
# scaler of environmental noise
miu_e = 0.2 
tmp_l=foreach (s = c(1:reps))%dopar%{
  fail=c()
  tmp<- data.frame()
  for(d in c(0,0.1)){
# scaler of demographic noise; O represents that demographic stochasticity is excluded
      miu_d = d
    for(i in 1:9){
      alpha_mean = 0.2
      flag_coex = 0; ## test of coexistence by deterministic model
      try=0
      while(flag_coex %in% c(0, NA) & try<1000){
        n = i+1;# species number

        alpha_ij<- runif(n*(n-1), max = 2*alpha_mean+0.01, min= 0.01)
        alpha<- diag(n)
        alpha[alpha == 0]<- alpha_ij 
  
        r = 2^runif(n,min=-1,max=0) ## random r
        K = 2^runif(n,min=-1,max=0) ## random K
       
        sigma_e = miu_e*runif(n,min=0,max=0.1); 
        sigma_d = miu_d*runif(n,min=0,max=0.1);  
        
        tsteps = 3e3;
        tsteps_r = 1e2;
        X = array(0,dim=c(tsteps,n));
        X[1,] = K/n
        for(t in c(1:(tsteps-1))){
          f = r*(1-alpha%*%X[t,]/K);
          X[t+1,] = X[t,]*exp(f);
        }
        XX = X[c((tsteps-(tsteps_r-1)):tsteps),]
        
        flag_coex = min(XX)>1e-2; ## testing coexistence
        
        range_X = apply(XX,2,max)-apply(XX,2,min)
        flag_coex = flag_coex * (max(range_X)<1e-2)
  
        try=try+1
      }
      
      if (flag_coex==1){
        ## simulation of mixture community
        tsteps = 1e3;
        tsteps_r = 1e2
        X = array(0,dim=c(tsteps,n));
        X[1,] = K/n
        xi_e = sigma_e*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n))  ## environmental noises (independent)
        xi_d = sigma_d*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## demographic noises
        for(t in c(1:(tsteps-1))){
          f = r*(1-alpha%*%X[t,]/(K)) + xi_e[t,] + xi_d[t,]*X[t,]^(-0.5);
          X[t+1,] = X[t,]*exp(f);
        }
        XX.mix = X[c((tsteps-(tsteps_r-1)):tsteps),];
        flag_coex = min(XX.mix)>1e-2; ## testing coexistence
    
        ## simulation of monoculture community
        X = array(0,dim=c(tsteps,n));
        X[1,] = K/n
        xi_e = sigma_e*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## environmental noises (independent)
        xi_d = sigma_d*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## demographic noises
        for(t in c(1:(tsteps-1))){
          f = r*(1-diag(n)%*%X[t,]/(K/n)) + xi_e[t,] + xi_d[t,]*X[t,]^(-0.5);
          X[t+1,] = X[t,]*exp(f);
        }
        XX.mono = X[c((tsteps-(tsteps_r-1)):tsteps),]
        
        print(i+1)
        print(s)
        CVpartition<- variability_partition(XX.mono, XX.mix)
        CESE<- BEF(XX.mono, XX.mix)
        tmp<- rbind(tmp, c(n,s,alpha_mean,miu_d,CVpartition, CESE))
      }else{
        fail=rbind(fail,data.frame(run=i+1,s=s,a=alpha_mean,miu_d=d,r=r,K=K))
      }
    }
  }
  return(tmp)
}

sim_results_richness_D0_D1=rbindlist(tmp_l[1:10000])
colnames(sim_results_richness_D0_D1)<- c("Richness","rep","alpha","miu_d","CV.o", "CV.e","CV.b", "CVS.o", "CVS.e", "phi.o", "phi.e",
                           "AE.SpVar", "SE.SpVar", "AE.Syn", "SE.Syn","NBE.F","SE","CE")


##Relationships between the strength of interspecific competition and the relative contributions of the multiple components behind biodiversity effects on community variability in the competition model
stopImplicitCluster()
cl=makeCluster(10)
registerDoParallel(cl)
reps = 3000
# scaler of environmental noise
miu_e = 0.2 
tmp_l=foreach (s = c(1:reps))%dopar%{
  fail=c()
  tmp<- data.frame()
  for(d in c(0,0.1)){
 # scaler of demographic noise; O represents that demographic stochasticity is excluded
    miu_d = d
    for(i in c(3, 6, 9)){
      for (j in seq(0.02,0.38,0.02) ) {
      alpha_mean = j
      flag_coex = 0; ## test of coexistence by deterministic model
      try=0
      while(flag_coex %in% c(0, NA) & try<1000){
        n = i;# species number
        
        alpha_ij<- runif(n*(n-1), max = 2*alpha_mean+0.01, min= 0.01)
        alpha<- diag(n)
        alpha[alpha == 0]<- alpha_ij 
        
        r = 2^runif(n,min=-1,max=0) ## random r
        K = 2^runif(n,min=-1,max=0) ## random K
        
        sigma_e = miu_e*runif(n,min=0,max=0.1); 
        sigma_d = miu_d*runif(n,min=0,max=0.1);  
        
        tsteps = 3e3;
        tsteps_r = 1e2;
        X = array(0,dim=c(tsteps,n));
        X[1,] = K/n
        for(t in c(1:(tsteps-1))){
          f = r*(1-alpha%*%X[t,]/K);
          X[t+1,] = X[t,]*exp(f);
        }
        XX = X[c((tsteps-(tsteps_r-1)):tsteps),]
        flag_coex = min(XX)>1e-2; ## testing coexistence
        range_X = apply(XX,2,max)-apply(XX,2,min)
        flag_coex = flag_coex * (max(range_X)<1e-2)
        try=try+1
      }
      
      if (flag_coex==1){
        ## simulation of mixture community
        tsteps = 1e3;
        tsteps_r = 1e2
        X = array(0,dim=c(tsteps,n));
        X[1,] = K/n
        xi_e = sigma_e*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n))
        xi_d = sigma_d*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## demographic noises
        for(t in c(1:(tsteps-1))){
          f = r*(1-alpha%*%X[t,]/(K)) + xi_e[t,] + xi_d[t,]*X[t,]^(-0.5);
          X[t+1,] = X[t,]*exp(f);
        }
        XX.mix = X[c((tsteps-(tsteps_r-1)):tsteps),];
        flag_coex = min(XX.mix)>1e-2; ## testing coexistence
      
        ## simulation of monoculture community
        X = array(0,dim=c(tsteps,n));
        X[1,] = K/n
        xi_e = sigma_e*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## environmental noises (independent)
        xi_d = sigma_d*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## demographic noises
        for(t in c(1:(tsteps-1))){
          f = r*(1-diag(n)%*%X[t,]/(K/n)) + xi_e[t,] + xi_d[t,]*X[t,]^(-0.5);
          X[t+1,] = X[t,]*exp(f);
        }
        XX.mono = X[c((tsteps-(tsteps_r-1)):tsteps),]
        print(i+1)
        print(s)
        CVpartition<- variability_partition(XX.mono, XX.mix)
        CESE<- BEF(XX.mono, XX.mix)
        tmp<- rbind(tmp, c(n,s,alpha_mean,miu_d,CVpartition, CESE))
      }else{
        fail=rbind(fail,data.frame(run=i+1,s=s,a=alpha_mean,miu_d=d,r=r,K=K))
      }
      }
    }
  }
  return(tmp)
}

sim_results_competition=rbindlist(tmp_l[1:3000])
colnames(sim_results_competition)<- c("Richness","rep","alpha","miu_d","CV.o", "CV.e","CV.b", "CVS.o", "CVS.e", "phi.o", "phi.e",
                                         "AE.SpVar", "SE.SpVar", "AE.Syn", "SE.Syn","NBE.F","SE","CE")

#test the effect of the correlation between r and K on SE.SpVar and SE.Syn
stopImplicitCluster()
cl=makeCluster(10)
registerDoParallel(cl)
reps=10000
# scaler of environmental noise
miu_e = 0.2 
tmp_l=foreach (s = c(1:reps))%dopar%{
  fail=c()
  tmp<- data.frame()
  for(d in c(0,0.1)){
# scaler of demographic noise; O represents that demographic stochasticity is excluded
    miu_d = d
    for(i in 1:9){
        alpha_mean = 0.2
        flag_coex = 0; ## test of coexistence by deterministic model
        try=0
        while(flag_coex %in% c(0, NA) & try<1000){
          n = i+1;# species number
          
          alpha_ij<- runif(n*(n-1), max = alpha_mean, min= alpha_mean)
          alpha<- diag(n)
          alpha[alpha == 0]<- alpha_ij 
          
          r = 2^runif(n,min=-1,max=0)## random r
          
          K = 2^runif(n,min=-1,max=1) ## random K
          cor_rK<- cor(r,K)
          sigma_e = miu_e*runif(n,min=0.05,max=0.05); 
          sigma_d = miu_d*runif(n,min=0.05,max=0.05);  
          
          tsteps = 3e3;
          tsteps_r = 1e2;
          X = array(0,dim=c(tsteps,n));
          X[1,] = K/n
          for(t in c(1:(tsteps-1))){
            f = r*(1-alpha%*%X[t,]/K);
            X[t+1,] = X[t,]*exp(f);
          }
          XX = X[c((tsteps-(tsteps_r-1)):tsteps),]

          flag_coex = min(XX)>1e-2; ## testing coexistence
          range_X = apply(XX,2,max)-apply(XX,2,min)
          flag_coex = flag_coex * (max(range_X)<1e-2)
          try=try+1
        }
        
        if (flag_coex==1){
          ## simulation of mixture community
          tsteps = 1e3;
          tsteps_r = 1e2
          X = array(0,dim=c(tsteps,n));
          X[1,] = K/n
          xi_e = sigma_e*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n))
          xi_d = sigma_d*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## demographic noises
          for(t in c(1:(tsteps-1))){
            f = r*(1-alpha%*%X[t,]/(K)) + xi_e[t,] + xi_d[t,]*X[t,]^(-0.5);
            X[t+1,] = X[t,]*exp(f);
          }
          XX.mix = X[c((tsteps-(tsteps_r-1)):tsteps),];
          flag_coex = min(XX.mix)>1e-2; ## testing coexistence
          
          ## simulation of monoculture community
          X = array(0,dim=c(tsteps,n));
          X[1,] = K/n
          xi_e = sigma_e*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## environmental noises (independent)
          xi_d = sigma_d*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## demographic noises
          for(t in c(1:(tsteps-1))){
            f = r*(1-diag(n)%*%X[t,]/(K/n)) + xi_e[t,] + xi_d[t,]*X[t,]^(-0.5);
            X[t+1,] = X[t,]*exp(f);
          }
          XX.mono = X[c((tsteps-(tsteps_r-1)):tsteps),]
          
          print(i+1)
          print(s)
          CVpartition<- variability_partition(XX.mono, XX.mix)
          CESE<- BEF(XX.mono, XX.mix)
          tmp<- rbind(tmp, c(n,s,alpha_mean,miu_d,cor_rK,CVpartition, CESE))
        }else{
          fail=rbind(fail,data.frame(run=i+1,s=s,a=alpha_mean,miu_d=d,r=r,K=K))
        }
    }
  }
  return(tmp)
}

sim_rK_correlation=rbindlist(tmp_l[1:10000])
colnames(sim_rK_correlation)<- c("Richness","rep","alpha","miu_d","cor_rK","CV.o", "CV.e","CV.b", "CVS.o", "CVS.e", "phi.o", "phi.e",
                                      "AE.SpVar", "SE.SpVar", "AE.Syn", "SE.Syn","NBE.F","SE","CE")

#test the effect of the correlation between r and a on SE.SpVar and SE.Syn
stopImplicitCluster()
cl=makeCluster(10)
registerDoParallel(cl)
reps=10000
# scaler of environmental noise
miu_e = 0.2 
tmp_l=foreach (s = c(1:reps))%dopar%{
  fail=c()
  tmp<- data.frame()
  for(d in c(0,0.1)){
    # a competitive matrix, in which a strong competitor (1) competes more with others, while (2) is less affected by others.
    miu_d = d
    for(i in 1:9){
      alpha_mean = 0.2
      flag_coex = 0; ## test of coexistence by deterministic model
      try=0
      while(flag_coex %in% c(0, NA) & try<1000){
        n = i+1;#Species number
        # set competition matrix, where strong competitor is (1) more competition to others, but (2) less affected by others
        alpha_ij<- runif(n*(n-1), max = 2*alpha_mean+0.01, min= 0.01)
        alpha<- diag(n)
        alpha_ij<- runif( n*(n-1)/2  , max = 2*alpha_mean+0.01, min= 0.01)
        alpha_ji<- 0.42-alpha_ij
        alpha[upper.tri(alpha)==1]<-alpha_ij
        alpha.t<- t(alpha)
        alpha.t[upper.tri(alpha.t)==1]<-alpha_ji
        alpha<-  alpha.t
        a_iE<- rowSums(alpha)-1

        r = 2^runif(n,min=-1,max=0)## random r
        K = 2^runif(n,min=0,max=0)## K is set to 1 for all species
        cor_ra_iE<- cor(r,a_iE)
        
        sigma_e = miu_e*runif(n,min=0.05,max=0.05); 
        sigma_d = miu_d*runif(n,min=0.05,max=0.05);  
        
        tsteps = 3e3;
        tsteps_r = 1e2;
        X = array(0,dim=c(tsteps,n));
        X[1,] = K/n
        for(t in c(1:(tsteps-1))){
          f = r*(1-alpha%*%X[t,]/K);
          X[t+1,] = X[t,]*exp(f);
        }
        XX = X[c((tsteps-(tsteps_r-1)):tsteps),]
        
        flag_coex = min(XX)>1e-2; ## testing coexistence
        range_X = apply(XX,2,max)-apply(XX,2,min)
        flag_coex = flag_coex * (max(range_X)<1e-2)
        try=try+1
      }
      
      if (flag_coex==1){
        ## simulation of mixture community
        tsteps = 1e3;
        tsteps_r = 1e2
        X = array(0,dim=c(tsteps,n));
        X[1,] = K/n
        xi_e = sigma_e*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n))
        xi_d = sigma_d*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## demographic noises
        for(t in c(1:(tsteps-1))){
          f = r*(1-alpha%*%X[t,]/(K)) + xi_e[t,] + xi_d[t,]*X[t,]^(-0.5);
          X[t+1,] = X[t,]*exp(f);
        }
        XX.mix = X[c((tsteps-(tsteps_r-1)):tsteps),];
        flag_coex = min(XX.mix)>1e-2; ## testing coexistence
      
        ## simulation of monoculture community
        X = array(0,dim=c(tsteps,n));
        X[1,] = K/n
        xi_e = sigma_e*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## environmental noises (independent)
        xi_d = sigma_d*array(rnorm(n*tsteps,0,1),dim=c(tsteps,n)); ## demographic noises
        for(t in c(1:(tsteps-1))){
          f = r*(1-diag(n)%*%X[t,]/(K/n)) + xi_e[t,] + xi_d[t,]*X[t,]^(-0.5);
          X[t+1,] = X[t,]*exp(f);
        }
        XX.mono = X[c((tsteps-(tsteps_r-1)):tsteps),]
        
        print(i+1)
        print(s)
        CVpartition<- variability_partition(XX.mono, XX.mix)
        CESE<- BEF(XX.mono, XX.mix)
        tmp<- rbind(tmp, c(n,s,alpha_mean,miu_d,cor_ra_iE,CVpartition, CESE))
      }else{
        fail=rbind(fail,data.frame(run=i+1,s=s,a=alpha_mean,miu_d=d,r=r,K=K))
      }
    }
  }
  return(tmp)
}

sim_ra_correlation=rbindlist(tmp_l[1:10000])
colnames(sim_ra_correlation)<- c("Richness","rep","alpha","miu_d","cor_ra_iE","CV.o", "CV.e","CV.b", "CVS.o", "CVS.e", "phi.o", "phi.e",
                           "AE.SpVar", "SE.SpVar", "AE.Syn", "SE.Syn","NBE.F","SE","CE")

