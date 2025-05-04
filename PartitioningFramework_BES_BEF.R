#The function used to calculate the biodiversity effect of interspecific interaction on community variability
#matrix input
#X_mix is the observed community
#X_mono is the hypothetical community; note that, such as in grassland experiments, the hypothetical species biomass should be miu_mono = 1/n*miu_mix
#columns is time series, rows is species
# Note: covariance is calculated by n-1 (estimate from a sample) in cor()
variability_partition<- function(X_mono,X_mix){
  
  X_mono<- as.matrix(X_mono)
  X_mix<- as.matrix(X_mix)
  
  CV.mix = sd(rowSums(X_mix))/mean(rowSums(X_mix))
  CV.mono = sd(rowSums(X_mono))/mean(rowSums(X_mono))
# average biomass of community 
  miu.mix = mean(rowSums(X_mix))
  miu.mono = mean(rowSums(X_mono))
  
  tmp<- data.frame()
  for (i in 1:ncol(X_mix)) {
    CV.mix.i = sd(X_mix[,i])/mean(X_mix[,i])
    CV.mono.i = sd(X_mono[,i])/mean(X_mono[,i])
    P.mix.i = mean(X_mix[,i])/miu.mix
    P.mono.i = mean(X_mono[,i])/miu.mono 
    delta.CV = CV.mix.i-CV.mono.i
    delta.P = P.mix.i-P.mono.i
    tmp<- rbind(tmp, c(CV.mix.i,CV.mono.i,P.mix.i,P.mono.i,delta.CV,delta.P))
  }
#statistic Base line
  CV.base = mean(tmp[,2])
# partitioning in species variability  
  CVS.mix = sum(tmp[,1]*tmp[,3], na.rm = T)#species variability in mixture (observed community)
  CVS.mono = sum(tmp[,2]*tmp[,4], na.rm = T)#species variability in monoculture (expected community)
  
  AE.SpVar = sum(tmp[,5]*tmp[,3], na.rm = T)/CVS.mono
  SE.SpVar = sum(tmp[,2]*tmp[,6], na.rm = T)/CVS.mono
  cov.mix <- cov(X_mix)
  cov.mono <- cov(X_mono)
  cor.mix <- cor(X_mix)
  cor.mono <- cor(X_mono)
  delta.cor<- cor.mix-cor.mono

# partitioning in species synchrony   
  phi.mix = sum(cov.mix)/(sum(sqrt(diag(cov.mix))))^2 #species synchrony in mixture (observed community); there is the φ^2
  phi.mono = sum(cov.mono)/(sum(sqrt(diag(cov.mono))))^2 #species synchrony in monoculture (expected community); there is the φ^2
  
  omega.mix<- (cov.mix/cor.mix)/(sum(sqrt(diag(cov.mix))))^2
  omega.mono<- (cov.mono/cor.mono)/(sum(sqrt(diag(cov.mono))))^2
  omega.mix[is.na(omega.mix)]<-0
  omega.mono[is.na(omega.mono)]<-0
  delta.omega<- omega.mix - omega.mono

  SE.Syn = sum(cor.mono * delta.omega, na.rm =T) / phi.mono
  AE.Syn = sum(delta.cor * omega.mix, na.rm =T) / phi.mono

  return(c(CV.o=CV.mix, CV.e=CV.mono, CV.b= CV.base, CVS.o=CVS.mix, CVS.e=CVS.mono, phi.o=sqrt(phi.mix), phi.e= sqrt(phi.mono),
           AE.SpVar = AE.SpVar, SE.SpVar = SE.SpVar, AE.Syn = AE.Syn, SE.Syn = SE.Syn))
  # CV.o is the variability of observed community
  # CV.e is the variability of expected community without interspecific interaction; ecological reference
  # CV.b or CV.null is the average monoculture variability across component species; Statistical baseline
  # CVS.o is the species variability in observed community
  # CVS.e is the species variability in expected community
  # phi.o is the species synchrony in observed community
  # phi.e is the species synchrony in expected community
  # AE.SpVar is the average effect of interspecific interaction on species variability
  # SE.SpVar is the selection effect of interspecific interaction on species variability
  # AE.Syn is the average effect of interspecific interaction on species synchrony
  # SE.Syn is the selection effect of interspecific interaction on species synchrony
  }




### The function used to calculate the biodiversity effect on community biomass
### we used an alternative definition of net biodiversity effect by the ratio of mixture and monoculture biomass, see Method
### NBE= miu_O/ miu_E= 1+ CE +SE
BEF<- function(X_mono,X_mix){
  miu_O<- colMeans(X_mix)
  miu_E<- colMeans(X_mono)
  NBE.F= sum( (miu_O/miu_E-1) * miu_E)/sum(miu_E)+1
  SE =   cov( (miu_O/miu_E-1) , miu_E)*(n-1)/sum(miu_E)
  CE =   sum(miu_O/miu_E-1) * sum(miu_E)/n/sum(miu_E)
  return(c(NBE.F=NBE.F,SE=SE,CE=CE))
}