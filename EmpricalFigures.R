library("nlme")
library("tidyverse")
library("MuMIn")
library("piecewiseSEM")
library("patchwork")
library("rgdal")
## The function of approximations
nmax=5
f_x <- function(x, nmax) {
  s <- 0
  for (n in 1:nmax) {
    s <- s + (-1)^(n+1) * (x^n) / n
  }
  return(s)
}
##  Different forms of scaling
scaleFUN <- function(x) sprintf("%.1f", x)
#reading results of partitioning component in Biodiversity experiments
#NBE.CV, NIIE, NIIE.SpVar, NIIE.Syn,  NIDE, NIDE.SpVar, NIDE.Syn are already logarithmically transformed
BiodivExps_results<- read.csv(".../Partitioning_BiodivExps_results.csv")
#BiodivExps_results<- Partitioning_BiodivExps_results

#Fig.3e
NBECV.R<-(lme(NBE.CV~Richness, random = ~1|(Study/Exp), BiodivExps_results))
predict.NBECV.R<-predict(NBECV.R, BiodivExps_results , level = 0,se.fit=T)
NIIE.R<-(lme(NIIE~Richness, random = ~1|(Study/Exp), BiodivExps_results))
predict.NIIE.R<-predict(NIIE.R, BiodivExps_results , level = 0,se.fit=T)
NIDE.R<-(lme(NIDE~Richness, random = ~1|(Study/Exp), BiodivExps_results))
predict.NIDE.R<-predict(NIDE.R, BiodivExps_results , level = 0,se.fit=T)

Fig.3e<-ggplot(BiodivExps_results)+ 
  geom_hline(aes(yintercept = 0), linetype =2)+
  geom_ribbon(aes(x=Richness, ymin=predict.NBECV.R$fit-1.96*predict.NBECV.R$se.fit ,
                  ymax=predict.NBECV.R$fit+1.96*predict.NBECV.R$se.fit), fill = "#3A3A3A",alpha = 0.15,linetype=0)+
  geom_ribbon(aes(x=Richness, ymin=predict.NIDE.R$fit-1.96*predict.NIDE.R$se.fit ,
                  ymax=predict.NIDE.R$fit+1.96*predict.NIDE.R$se.fit), fill ="#B86500",alpha = 0.15,linetype=0)+
  geom_ribbon(aes(x=Richness, ymin=predict.NIIE.R$fit-1.96*predict.NIIE.R$se.fit ,
                  ymax=predict.NIIE.R$fit+1.96*predict.NIIE.R$se.fit), fill = "#571C72",alpha = 0.15,linetype=0)+
  geom_line(aes(x = Richness, y=predict.NBECV.R$fit),color = "#3A3A3A",linewidth = 0.7)+
  geom_line(aes(x = Richness, y=predict.NIDE.R$fit),color ="#B86500" , linewidth = 0.7, linetype= 2)+
  geom_line(aes(x = Richness, y=predict.NIIE.R$fit), color ="#571C72", linewidth = 0.7)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("Richness")+
  scale_y_continuous(limits=c(-1.2, 1), labels=scaleFUN)+
  theme_bw()+
  theme(panel.border = element_rect(color="black", size=0.4),
        panel.grid.major = element_line(linewidth = 0.1),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.position="")

NIDE.SpVar.R<-(lme(NIDE.SpVar~Richness, random = ~1|(Study/Exp), BiodivExps_results))
predict.NIDE.SpVar.R<-predict(NIDE.SpVar.R, BiodivExps_results , level = 0,se.fit=T)
NIDE.Syn.R<-(lme(NIDE.Syn~Richness, random = ~1|(Study/Exp), BiodivExps_results))
predict.NIDE.Syn.R<-predict(NIDE.Syn.R, BiodivExps_results , level = 0,se.fit=T)

Fig.3f<-ggplot(BiodivExps_results)+ 
  geom_hline(aes(yintercept = 0), linetype =2)+
  geom_ribbon(aes(x=Richness, ymin=predict.NIDE.SpVar.R$fit-1.96*predict.NIDE.SpVar.R$se.fit ,
                  ymax=predict.NIDE.SpVar.R$fit+1.96*predict.NIDE.SpVar.R$se.fit), fill = "#8D0E21",alpha = 0.15,linetype=0)+
  geom_ribbon(aes(x=Richness, ymin=predict.NIDE.R$fit-1.96*predict.NIDE.R$se.fit ,
                  ymax=predict.NIDE.R$fit+1.96*predict.NIDE.R$se.fit), fill ="#D27300",alpha = 0.1,linetype=0)+
  geom_ribbon(aes(x=Richness, ymin=predict.NIDE.Syn.R$fit-1.96*predict.NIDE.Syn.R$se.fit ,
                  ymax=predict.NIDE.Syn.R$fit+1.96*predict.NIDE.Syn.R$se.fit), fill = "#1A5602",alpha = 0.15,linetype=0)+
  geom_line(aes(x = Richness, y=predict.NIDE.SpVar.R$fit),color = "#8D0E21",linewidth = 0.7)+
  geom_line(aes(x = Richness, y=predict.NIDE.R$fit),color = "#D27300", linewidth = 0.7, linetype= 2)+
  geom_line( aes(x =Richness, y=predict.NIDE.Syn.R$fit), color ="#1A5602", linewidth = 0.7)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("Richness")+
  scale_y_continuous(limits=c(-1.2, 1), labels=scaleFUN)+
  theme_bw()+
  theme(panel.border = element_rect(color="black", size=0.4),
        panel.grid.major = element_line(linewidth = 0.1),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.position="")

AE.SpVar.R<-(lme(AE.SpVar~Richness, random = ~1|(Study/Exp), BiodivExps_results))
predict.AE.SpVar.R<-predict(AE.SpVar.R, BiodivExps_results , level = 0,se.fit=T)
SE.SpVar.R<-(lme(SE.SpVar~Richness, random = ~1|(Study/Exp), BiodivExps_results))
predict.SE.SpVar.R<-predict(SE.SpVar.R, BiodivExps_results, level = 0,se.fit=T)

Fig.3g<-ggplot(BiodivExps_results)+ 
  geom_hline(aes(yintercept = 0), linetype =2)+
  geom_ribbon(aes(x=Richness, ymin=f_x(predict.AE.SpVar.R$fit,2) - f_x(1.96*predict.AE.SpVar.R$se.fit,2) ,
                  ymax=f_x(predict.AE.SpVar.R$fit,2) + f_x(1.96*predict.AE.SpVar.R$se.fit,2)), fill ="#D40808" ,alpha = 0.15,linetype=0)+
  geom_ribbon(aes(x=Richness, ymin=f_x(predict.SE.SpVar.R$fit,2)-f_x(1.96*predict.SE.SpVar.R$se.fit,2) ,
                  ymax=f_x(predict.SE.SpVar.R$fit,2)+f_x(1.96*predict.SE.SpVar.R$se.fit,2)), fill = "#F59A67",alpha = 0.15,linetype=0)+
  geom_line(aes(x = Richness, y=f_x(predict.AE.SpVar.R$fit,2)),color = "#D40808",linewidth = 0.7)+
  geom_line(aes(x = Richness, y=f_x(predict.SE.SpVar.R$fit,2)),color ="#F59A67", linewidth = 0.7, linetype =2)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("Richness")+
  scale_y_continuous(limits=c(-1.2, 1), labels=scaleFUN)+
  theme_bw()+
  theme(panel.border = element_rect(color="black", size=0.4),
        panel.grid.major = element_line(linewidth = 0.1),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.position="")

AE.Syn.R<-lme(AE.Syn~Richness, random = ~1|(Study/Exp),BiodivExps_results)
predict.AE.Syn.R<-predict(AE.Syn.R, BiodivExps_results , level = 0,se.fit=T)
SE.Syn.R<-lme(SE.Syn~Richness, random = ~1|(Study/Exp),BiodivExps_results)
predict.SE.Syn.R<-predict(SE.Syn.R, BiodivExps_results , level = 0,se.fit=T)

Fig.3h<-ggplot(BiodivExps_results)+ 
  geom_hline(aes(yintercept = 0), linetype =2)+
  geom_ribbon(aes(x=Richness, ymin=f_x(predict.AE.Syn.R$fit,2)/2 - f_x(1.96*predict.AE.Syn.R$se.fit,2)/2 ,
                  ymax=f_x(predict.AE.Syn.R$fit,2)/2 + f_x(1.96*predict.AE.Syn.R$se.fit,2)/2 ), fill ="#6D921A" ,alpha = 0.15,linetype=0)+
  geom_ribbon(aes(x=Richness, ymin=f_x(predict.SE.Syn.R$fit,2)/2 - f_x(1.96*predict.SE.Syn.R$se.fit,2)/2 ,
                  ymax=f_x(predict.SE.Syn.R$fit,2)/2 + f_x(1.96*predict.SE.Syn.R$se.fit,2)/2 ), fill ="#13C78B",alpha = 0.15,linetype=0)+
  geom_line(aes(x = Richness, y=f_x(predict.AE.Syn.R$fit,2)/2),color = "#6D921A",linewidth = 0.7)+
  geom_line(aes(x = Richness, y=f_x(predict.SE.Syn.R$fit,2)/2),color ="#13C78B", linewidth = 0.7, linetype= 2)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("Richness")+
  scale_y_continuous(limits=c(-1.2, 1), labels=scaleFUN)+
  theme_bw()+
  theme(panel.border = element_rect(color="black", size=0.4),
        panel.grid.major = element_line(linewidth = 0.1),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.position="")
Fig.3efgh<-Fig.3e+Fig.3f+Fig.3g+Fig.3h+plot_layout(ncol = 4)

#statistical results of the relationship between richness and the multiple components of biodiversity effect
#It is shown in Table S2 S3 S4
#note that: AE.SpVar, SE,SpVar, AE.Syn and SE.Syn is tested on natural scale; Others are tested on logarithmic scale

#Table S2
summary(NBECV.R)
summary(lme((NBE.CV) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((NBE.CV) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((NBE.CV) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((NBE.CV) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((NBE.CV) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
summary(NIDE.R)
summary(lme((NIDE) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((NIDE) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((NIDE) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((NIDE) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((NIDE) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
summary(NIIE.R)
summary(lme((NIIE) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((NIIE) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((NIIE) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((NIIE) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((NIIE) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))

#Table S3
summary(lme(NIIE.SpVar~Richness, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme((NIIE.SpVar) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((NIIE.SpVar) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((NIIE.SpVar) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((NIIE.SpVar) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((NIIE.SpVar) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
summary(lme(NIIE.Syn~Richness, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme((NIIE.Syn) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((NIIE.Syn) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((NIIE.Syn) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((NIIE.Syn) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((NIIE.Syn) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
summary(NIDE.SpVar.R)
summary(lme((NIDE.SpVar) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((NIDE.SpVar) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((NIDE.SpVar) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((NIDE.SpVar) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((NIDE.SpVar) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
summary(NIDE.Syn.R)
summary(lme((NIDE.Syn) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((NIDE.Syn) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((NIDE.Syn) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((NIDE.Syn) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((NIDE.Syn) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
#Table S4
summary(AE.SpVar.R)
summary(lme((AE.SpVar) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((AE.SpVar) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((AE.SpVar) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((AE.SpVar) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((AE.SpVar) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
summary(SE.SpVar.R)
summary(lme((SE.SpVar) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((SE.SpVar) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((SE.SpVar) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((SE.SpVar) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((SE.SpVar) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
summary(AE.Syn.R)
summary(lme((AE.Syn) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((AE.Syn) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((AE.Syn) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((AE.Syn) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((AE.Syn) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
summary(SE.Syn.R)
summary(lme((SE.Syn) ~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
summary(lm((SE.Syn) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
summary(lm((SE.Syn) ~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(lm((SE.Syn) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
summary(lm((SE.Syn) ~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))

#fig.4a, b, c
a.df.1<- BiodivExps_results%>%
  dplyr:: select(Study, Plot, Richness,year,NBE.CV)%>%
  filter(NBE.CV >= quantile(NBE.CV, probs = 0.005),
         NBE.CV <= quantile(NBE.CV, probs = 0.995))%>%
  gather(Item,Value,NBE.CV)
a.df.2<- BiodivExps_results%>%
  dplyr:: select(Study, Plot, Richness,year,NIIE)%>%
  filter(NIIE >= quantile(NIIE, probs = 0.005),
         NIIE <= quantile(NIIE, probs = 0.995))%>%
  gather(Item,Value,NIIE)
a.df.3<- BiodivExps_results%>%
  dplyr:: select(Study, Plot, Richness,year,NIDE)%>%
  filter(NIDE >= quantile(NIDE, probs = 0.005),
         NIDE <= quantile(NIDE, probs = 0.995))%>%
  gather(Item,Value,NIDE)
a.df.4<- BiodivExps_results%>%
  dplyr:: select(Study, Plot, Richness,year,NIIE.SpVar)%>%
  filter(NIIE.SpVar>= quantile(NIIE.SpVar, probs = 0.005),
         NIIE.SpVar <= quantile(NIIE.SpVar, probs = 0.995))%>%
  gather(Item,Value,NIIE.SpVar)
a.df.5<- BiodivExps_results%>%
  dplyr:: select(Study, Plot, Richness,year,NIIE.Syn)%>%
  filter(NIIE.Syn >= quantile(NIIE.Syn, probs = 0.005),
         NIIE.Syn <= quantile(NIIE.Syn, probs = 0.995))%>%
  gather(Item,Value,NIIE.Syn)
Fig4a.df<- rbind(a.df.1,a.df.2,a.df.3,a.df.4,a.df.5)

Fig4a.df$Order<- NA
Fig4a.df[Fig4a.df$Item=="NBE.CV",]$Order<- 1
Fig4a.df[Fig4a.df$Item%in% c("NIIE","NIDE"),]$Order<- 2
Fig4a.df[Fig4a.df$Item%in%c("NIIE.SpVar","NIIE.Syn"),]$Order<- 3
Fig4a.df$Item<- factor(Fig4a.df$Item,levels = c("NBE.CV","NIIE","NIDE","NIIE.SpVar","NIIE.Syn"))

NBE.CV_Bigbio<-(lm(NBE.CV ~1, BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NBE.CV_Wageningen<-(lm(NBE.CV ~1, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NBE.CV_Texas<-(lm(NBE.CV ~1, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NBE.CV_Agro<-(lme(NBE.CV ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Agro",]))
NBE.CV_Biodepth<-(lme(NBE.CV ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NBE.CV_Jena<-(lm(NBE.CV ~1, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

NIDE_Bigbio<-(lm(NIDE ~1,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIDE_Wageningen<-(lm(NIDE ~1, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIDE_Texas<-(lm(NIDE ~1, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIDE_Agro<-(lme(NIDE ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Agro",]))
NIDE_Biodepth<-(lme(NIDE ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIDE_Jena<-(lm(NIDE ~1, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

NIIE_Bigbio<-(lm(NIIE ~1,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIIE_Wageningen<-(lm(NIIE ~1, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIIE_Texas<-(lm(NIIE ~1,  BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIIE_Agro<-(lme(NIIE ~1, random = ~1|(Exp),BiodivExps_results[BiodivExps_results$Study == "Agro",]))
NIIE_Biodepth<-(lme(NIIE ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIIE_Jena<-(lm(NIIE ~1,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))

NIIE.SpVar_Bigbio<-(lm(NIIE.SpVar ~1,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIIE.SpVar_Wageningen<-(lm(NIIE.SpVar ~1,  BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIIE.SpVar_Texas<-(lm(NIIE.SpVar ~1, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIIE.SpVar_Agro<-(lme(NIIE.SpVar ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Agro",]))
NIIE.SpVar_Biodepth<-(lme(NIIE.SpVar ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIIE.SpVar_Jena<-(lm(NIIE.SpVar ~1,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))

NIIE.Syn_Bigbio<-(lm(NIIE.Syn ~1,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIIE.Syn_Wageningen<-(lm(NIIE.Syn ~1,  BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIIE.Syn_Texas<-(lm(NIIE.Syn ~1,  BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIIE.Syn_Agro<-(lme(NIIE.Syn ~1, random = ~1|(Exp),  BiodivExps_results[BiodivExps_results$Study == "Agro",]))
NIIE.Syn_Biodepth<-(lme(NIIE.Syn ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIIE.Syn_Jena<-(lm(NIIE.Syn ~1, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

a1.CI<- data.frame(Item = c("NBE.CV"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                  upper = c( confint(NBE.CV_Bigbio)[2], confint(NBE.CV_Wageningen)[2],confint(NBE.CV_Texas)[2],intervals(NBE.CV_Agro,which = "fixed")$fixed[3],intervals(NBE.CV_Biodepth,which = "fixed")$fixed[3],confint(NBE.CV_Jena,)[2]),
                  mean =  c(coefficients(NBE.CV_Bigbio), coefficients(NBE.CV_Wageningen),coefficients(NBE.CV_Texas),intervals(NBE.CV_Agro,which = "fixed")$fixed[2],intervals(NBE.CV_Biodepth,which = "fixed")$fixed[2],coefficients(NBE.CV_Jena)),
                  lower =  c( confint(NBE.CV_Bigbio)[1], confint(NBE.CV_Wageningen)[1],confint(NBE.CV_Texas)[1],intervals(NBE.CV_Agro,which = "fixed")$fixed[1],intervals(NBE.CV_Biodepth,which = "fixed")$fixed[1],confint(NBE.CV_Jena)[1]))
a2.CI<- data.frame(Item = c("NIDE"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                   upper = c( confint(NIDE_Bigbio)[2], confint(NIDE_Wageningen)[2],confint(NIDE_Texas)[2],intervals(NIDE_Agro,which = "fixed")$fixed[3],intervals(NIDE_Biodepth,which = "fixed")$fixed[3],confint(NIDE_Jena)[2]),
                   mean =  c( coefficients(NIDE_Bigbio), coefficients(NIDE_Wageningen),coefficients(NIDE_Texas),intervals(NIDE_Agro,which = "fixed")$fixed[2],intervals(NIDE_Biodepth,which = "fixed")$fixed[2],coefficients(NIDE_Jena)),
                   lower =  c( confint(NIDE_Bigbio)[1], confint(NIDE_Wageningen)[1],confint(NIDE_Texas)[1],intervals(NIDE_Agro,which = "fixed")$fixed[1],intervals(NIDE_Biodepth,which = "fixed")$fixed[1],confint(NIDE_Jena)[1]))
a3.CI<- data.frame(Item = c("NIIE"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                   upper = c(confint(NIIE_Bigbio)[2], confint(NIIE_Wageningen)[2],confint(NIIE_Texas)[2],intervals(NIIE_Agro,which = "fixed")$fixed[3],intervals(NIIE_Biodepth,which = "fixed")$fixed[3],confint(NIIE_Jena)[2]),
                   mean =  c(coefficients(NIIE_Bigbio), coefficients(NIIE_Wageningen),coefficients(NIIE_Texas),intervals(NIIE_Agro,which = "fixed")$fixed[2],intervals(NIIE_Biodepth,which = "fixed")$fixed[2],coefficients(NIIE_Jena)),
                   lower =  c( confint(NIIE_Bigbio)[1], confint(NIIE_Wageningen)[1],confint(NIIE_Texas)[1],intervals(NIIE_Agro,which = "fixed")$fixed[1],intervals(NIIE_Biodepth,which = "fixed")$fixed[1],confint(NIIE_Jena)[1]))

a4.CI<- data.frame(Item = c("NIIE.SpVar"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                   upper = c(confint(NIIE.SpVar_Bigbio)[2], confint(NIIE.SpVar_Wageningen)[2],confint(NIIE.SpVar_Texas)[3],intervals(NIIE.SpVar_Agro,which = "fixed")$fixed[3],intervals(NIIE.SpVar_Biodepth,which = "fixed")$fixed[3],confint(NIIE.SpVar_Jena)[3]),
                   mean =  c(coefficients(NIIE.SpVar_Bigbio), coefficients(NIIE.SpVar_Wageningen),coefficients(NIIE.SpVar_Texas)[2],intervals(NIIE.SpVar_Agro,which = "fixed")$fixed[2],intervals(NIIE.SpVar_Biodepth,which = "fixed")$fixed[2],coefficients(NIIE.SpVar_Jena)[2]),
                   lower =  c(confint(NIIE.SpVar_Bigbio)[1], confint(NIIE.SpVar_Wageningen)[1],confint(NIIE.SpVar_Texas)[1],intervals(NIIE.SpVar_Agro,which = "fixed")$fixed[1],intervals(NIIE.SpVar_Biodepth,which = "fixed")$fixed[1],confint(NIIE.SpVar_Jena)[1]))
a5.CI<- data.frame(Item = c("NIIE.Syn"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                   upper = c(confint(NIIE.Syn_Bigbio)[2], confint(NIIE.Syn_Wageningen)[2],confint(NIIE.Syn_Texas)[2],intervals(NIIE.Syn_Agro,which = "fixed")$fixed[3],intervals(NIIE.Syn_Biodepth,which = "fixed")$fixed[3],confint(NIIE.Syn_Jena)[2]),
                   mean =  c(coefficients(NIIE.Syn_Bigbio), coefficients(NIIE.Syn_Wageningen),coefficients(NIIE.Syn_Texas),intervals(NIIE.Syn_Agro,which = "fixed")$fixed[2], intervals(NIIE.Syn_Biodepth,which = "fixed")$fixed[2],coefficients(NIIE.Syn_Jena)),
                   lower =  c(confint(NIIE.Syn_Bigbio)[1], confint(NIIE.Syn_Wageningen)[1],confint(NIIE.Syn_Texas)[1],intervals(NIIE.Syn_Agro,which = "fixed")$fixed[1],intervals(NIIE.Syn_Biodepth,which = "fixed")$fixed[1],confint(NIIE.Syn_Jena)[1]))
Fig4a.CI<-a1.CI%>%rbind(a2.CI)%>%rbind(a3.CI)%>%rbind(a4.CI)%>%rbind(a5.CI)

Fig4a.CI$Order<-NA
Fig4a.CI[Fig4a.CI$Item == "NBE.CV",]$Order <- 1
Fig4a.CI[Fig4a.CI$Item %in% c("NIDE", "NIIE"),]$Order <- 2
Fig4a.CI[Fig4a.CI$Item %in% c("NIIE.SpVar","NIIE.Syn"),]$Order <- 3
Fig4a.CI$Item<- factor(Fig4a.CI$Item,levels = c("NBE.CV","NIIE","NIDE","NIIE.SpVar","NIIE.Syn"))
my6<- c("#526797","#ec9a29","#a8201a", "#0f8b8d","#cfa093","#143642")
Fig4abc<- ggplot(Fig4a.df) +
  geom_violin( aes(x = Item, y = Value),linewidth= 0.3,fill = "#F2F2F2", color = "black", trim = F, width=0.6) +
  geom_errorbar(data =Fig4a.CI, aes(x = Item, y = mean, ymax = upper, ymin = lower, color = system),linewidth= 0.3, width = 0.01,position = position_jitter(width = 0.15, height = 0,seed =2))+
  geom_point(data =Fig4a.CI, aes(x = Item, y = mean, color = system), size = 1.2, position = position_jitter(width = 0.15, height = 0,seed =2))+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values =  my6)+ 
  ylab(expression(Value))+
  xlab("")+
  theme_bw()+
  scale_y_continuous(limits=c(-2.5, 2.5), labels=scaleFUN)+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="right")+
  facet_grid(~ Order, scales = "free_x", space='free')

#Fig.4d
d.df.2<- BiodivExps_results%>%
  select(Study, Plot, Richness,year,NIDE.SpVar)%>%
  filter(NIDE.SpVar >= quantile(NIDE.SpVar, probs = 0.005),
         NIDE.SpVar <= quantile(NIDE.SpVar, probs = 0.995))%>%
  gather(Item,Value,NIDE.SpVar)
d.df.3<- BiodivExps_results%>%
  select(Study, Plot, Richness,year, NIDE.Syn)%>%
  filter( NIDE.Syn >= quantile( NIDE.Syn, probs = 0.005),
          NIDE.Syn <= quantile( NIDE.Syn, probs = 0.995))%>%
  gather(Item,Value, NIDE.Syn)
Fig4d.df<- rbind(d.df.2,d.df.3)
Fig4d.df$Item<- factor(Fig4d.df$Item,levels = c("NIDE.SpVar","NIDE.Syn" ))

NIDE.SpVar_Bigbio<-(lm(NIDE.SpVar ~1,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIDE.SpVar_Wageningen<-(lm(NIDE.SpVar ~1,  BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIDE.SpVar_Texas<-(lm(NIDE.SpVar ~1,  BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIDE.SpVar_Agro<-(lme(NIDE.SpVar ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Agro",]))
NIDE.SpVar_Biodepth<-(lme(NIDE.SpVar ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIDE.SpVar_Jena<-(lm(NIDE.SpVar ~1, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

NIDE.Syn_Bigbio<-lm( NIDE.Syn ~1,  BiodivExps_results[BiodivExps_results$Study  == "Bigbio",])
NIDE.Syn_Wageningen<-lm( NIDE.Syn ~1, BiodivExps_results[BiodivExps_results$Study  == "Wageningen",])
NIDE.Syn_Texas<-lm( NIDE.Syn ~1,  BiodivExps_results[BiodivExps_results$Study  == "Texas",])
NIDE.Syn_Agro<-lme( NIDE.Syn ~1,random = ~1|(Exp),   BiodivExps_results[BiodivExps_results$Study  == "Agro",])
NIDE.Syn_Biodepth<-lme( NIDE.Syn ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study  == "Biodepth",])
NIDE.Syn_Jena<-lm( NIDE.Syn ~1, BiodivExps_results[BiodivExps_results$Study  == "Jena",])

d1.CI<- data.frame(Item = c("NIDE.SpVar"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                   upper = c(confint(NIDE.SpVar_Bigbio)[2], confint(NIDE.SpVar_Wageningen)[2],confint(NIDE.SpVar_Texas)[2],intervals(NIDE.SpVar_Agro,which = "fixed")$fixed[3],intervals(NIDE.SpVar_Biodepth,which = "fixed")$fixed[3],confint(NIDE.SpVar_Jena)[2]),
                   mean =  c(coefficients(NIDE.SpVar_Bigbio), coefficients(NIDE.SpVar_Wageningen),coefficients(NIDE.SpVar_Texas),intervals(NIDE.SpVar_Agro,which = "fixed")$fixed[2],intervals(NIDE.SpVar_Biodepth,which = "fixed")$fixed[2], coefficients(NIDE.SpVar_Jena)),
                   lower =  c(confint(NIDE.SpVar_Bigbio)[1], confint(NIDE.SpVar_Wageningen)[1],confint(NIDE.SpVar_Texas)[1],intervals(NIDE.SpVar_Agro,which = "fixed")$fixed[1],intervals(NIDE.SpVar_Biodepth,which = "fixed")$fixed[1],confint(NIDE.SpVar_Jena)[1]))
d2.CI<- data.frame(Item = c("NIDE.Syn"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                   upper = c(confint( NIDE.Syn_Bigbio)[2], confint( NIDE.Syn_Wageningen)[2], confint( NIDE.Syn_Texas)[2], intervals(NIDE.Syn_Agro,which = "fixed")$fixed[3], intervals( NIDE.Syn_Biodepth,which = "fixed")$fixed[3], confint( NIDE.Syn_Jena)[2]),
                   mean =  c(coefficients(NIIE.SpVar_Bigbio), coefficients( NIDE.Syn_Wageningen),coefficients(NIDE.Syn_Texas),intervals(NIDE.Syn_Agro,which = "fixed")$fixed[2],intervals( NIDE.Syn_Biodepth,which = "fixed")$fixed[2], coefficients( NIDE.Syn_Jena)),
                   lower =  c(confint( NIDE.Syn_Bigbio)[1], confint( NIDE.Syn_Wageningen)[1], confint( NIDE.Syn_Texas)[1], intervals(NIDE.Syn_Agro,which = "fixed")$fixed[1],intervals( NIDE.Syn_Biodepth,which = "fixed")$fixed[1],confint( NIDE.Syn_Jena)[1]))
Fig4d.CI<- rbind(d1.CI, d2.CI)

Fig4d<- ggplot(Fig4d.df) +
  geom_violin( aes(x = Item, y = Value),linewidth = 0.3, fill = "#F2F2F2", color = "black", trim = F, width=0.6) +
  geom_errorbar(data = Fig4d.CI, aes(x = Item, y = mean, ymax = upper, ymin = lower, color = system),linewidth = 0.3, width = 0.02,position = position_jitter(width = 0.15, height = 0,seed =2))+
  geom_point(data = Fig4d.CI, aes(x = Item, y = mean, color = system), size = 1.2, position = position_jitter(width = 0.15, height = 0,seed =2))+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values =  my6)+ 
  ylab(expression(Value))+
  xlab("")+
  theme_bw()+
  scale_y_continuous(limits=c(-2.5, 2.5), labels=scaleFUN)+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")

#Fig.4e
e.df.2<- BiodivExps_results%>%
  select(Study, Plot, Richness,year,AE.SpVar)%>%
  filter(AE.SpVar >= quantile(AE.SpVar, probs = 0.005),
         AE.SpVar <= quantile(AE.SpVar, probs = 0.995))%>%
  gather(Item,Value,AE.SpVar)
e.df.3<- BiodivExps_results%>%
  select(Study, Plot, Richness,year,SE.SpVar)%>%
  filter(SE.SpVar >= quantile(SE.SpVar, probs = 0.005),
         SE.SpVar <= quantile(SE.SpVar, probs = 0.995))%>%
  gather(Item,Value,SE.SpVar)
Fig4e.df<- rbind(e.df.2,e.df.3)
Fig4e.df$Item<- factor(Fig4e.df$Item,levels = c("AE.SpVar","SE.SpVar"))

AE.SpVar_Bigbio<-(lm(AE.SpVar ~1, BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
AE.SpVar_Wageningen<-(lm(AE.SpVar ~1,  BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
AE.SpVar_Texas<-(lm(AE.SpVar ~1,  BiodivExps_results[BiodivExps_results$Study == "Texas",]))
AE.SpVar_Agro<-(lme(AE.SpVar ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Agro",]))
AE.SpVar_Biodepth<-(lme(AE.SpVar ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
AE.SpVar_Jena<-(lm(AE.SpVar ~1,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))

SE.SpVar_Bigbio<-(lm(SE.SpVar ~1,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
SE.SpVar_Wageningen<-(lm(SE.SpVar ~1, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
SE.SpVar_Texas<-(lm(SE.SpVar ~1, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
SE.SpVar_Agro<-(lme(SE.SpVar ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Agro",]))
SE.SpVar_Biodepth<-(lme(SE.SpVar ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
SE.SpVar_Jena<-(lm(SE.SpVar ~1, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

e1.CI<- data.frame(Item = c("AE.SpVar"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                   upper = c(confint(AE.SpVar_Bigbio)[2], confint(AE.SpVar_Wageningen)[2],confint(AE.SpVar_Texas)[2],intervals(AE.SpVar_Agro,which = "fixed")$fixed[3],intervals(AE.SpVar_Biodepth,which = "fixed")$fixed[3],confint(AE.SpVar_Jena)[2]),
                   mean =  c( coefficients (AE.SpVar_Bigbio), coefficients(AE.SpVar_Wageningen),coefficients(AE.SpVar_Texas),intervals(AE.SpVar_Agro,which = "fixed")$fixed[2],intervals(AE.SpVar_Biodepth,which = "fixed")$fixed[2],coefficients(AE.SpVar_Jena)),
                   lower =  c(confint(AE.SpVar_Bigbio)[1], confint(AE.SpVar_Wageningen)[1],confint(AE.SpVar_Texas)[1],intervals(AE.SpVar_Agro,which = "fixed")$fixed[1],intervals(AE.SpVar_Biodepth,which = "fixed")$fixed[1],confint(AE.SpVar_Jena)[1]))
e2.CI<- data.frame(Item = c("SE.SpVar"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                   upper = c(confint(SE.SpVar_Bigbio)[2], confint(SE.SpVar_Wageningen)[2],confint(SE.SpVar_Texas)[2],intervals(SE.SpVar_Agro,which = "fixed")$fixed[3],intervals(SE.SpVar_Biodepth,which = "fixed")$fixed[3],confint(SE.SpVar_Jena)[2]),
                   mean =  c(coefficients(SE.SpVar_Bigbio), coefficients(SE.SpVar_Wageningen),coefficients(SE.SpVar_Texas),intervals(SE.SpVar_Agro,which = "fixed")$fixed[2],intervals(SE.SpVar_Biodepth,which = "fixed")$fixed[2],coefficients(SE.SpVar_Jena)),
                   lower =  c(confint(SE.SpVar_Bigbio)[1], confint(SE.SpVar_Wageningen)[1],confint(SE.SpVar_Texas)[1],intervals(SE.SpVar_Agro,which = "fixed")$fixed[1],intervals(SE.SpVar_Biodepth,which = "fixed")$fixed[1],confint(SE.SpVar_Jena)[1]))

Fig4e.CI<- rbind(e1.CI,e2.CI)

Fig4e<- ggplot(Fig4e.df) +
  geom_violin( aes(x = Item , y = f_x(Value, 2)),linewidth = 0.3,fill = "#F2F2F2", color = "black", trim = F, width=0.6) +

  geom_errorbar(data = Fig4e.CI, aes(x = Item, y = f_x(mean, 2), ymax = f_x(upper, 2), ymin = f_x(lower, 2), color = system),linewidth = 0.3, width = 0.02, position = position_jitter(width = 0.1, height = 0,seed =2))+
  geom_point(data = Fig4e.CI, aes(x = Item, y = f_x(mean, 2), color = system), size = 1.2, position = position_jitter(width = 0.1, height = 0,seed =2))+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values =  my6)+ 
  ylab("")+
  xlab("")+
  theme_bw()+
  scale_y_continuous(limits=c(-2.5, 2.5), labels=scaleFUN)+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")
#Fig.4f
f.df.2<- BiodivExps_results%>%
  select(Study, Plot, Richness,year,AE.Syn)%>%
  filter(AE.Syn >= quantile(AE.Syn, probs = 0.005),
         AE.Syn <= quantile(AE.Syn, probs = 0.995))%>%
  gather(Item,Value,AE.Syn)
f.df.3<- BiodivExps_results%>%
  select(Study, Plot, Richness,year, SE.Syn)%>%
  filter(SE.Syn >= quantile(SE.Syn, probs = 0.005),
         SE.Syn <= quantile(SE.Syn, probs = 0.995))%>%
  gather(Item,Value,SE.Syn)
Fig4f.df<- rbind(f.df.2,f.df.3)

AE.Syn_Bigbio<-(lm((AE.Syn) ~1,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
AE.Syn_Wageningen<-(lm((AE.Syn) ~1, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
AE.Syn_Texas<-(lm((AE.Syn) ~1, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
AE.Syn_Agro<-(lme((AE.Syn) ~1, random = ~1|(Exp),BiodivExps_results[BiodivExps_results$Study == "Agro",]))
AE.Syn_Biodepth<-(lme((AE.Syn) ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
AE.Syn_Jena<-(lm((AE.Syn) ~1,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))

SE.Syn_Bigbio<-(lm((SE.Syn) ~1, BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
SE.Syn_Wageningen<-(lm((SE.Syn) ~1,  BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
SE.Syn_Texas<-(lm((SE.Syn) ~1,  BiodivExps_results[BiodivExps_results$Study == "Texas",]))
SE.Syn_Agro<-(lme((SE.Syn) ~1, random = ~1|(Exp),  BiodivExps_results[BiodivExps_results$Study == "Agro",]))
SE.Syn_Biodepth<-(lme((SE.Syn) ~1, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
SE.Syn_Jena<-(lm((SE.Syn) ~1,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))



f1.CI<- data.frame(Item = c("AE.Syn"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                   upper = c( confint(AE.Syn_Bigbio)[2], confint(AE.Syn_Wageningen)[2],confint(AE.Syn_Texas)[2],intervals(AE.Syn_Agro,which = "fixed")$fixed[3],intervals(AE.Syn_Biodepth,which = "fixed")$fixed[3],confint(AE.Syn_Jena)[2]),
                   mean =  c(coefficients(AE.Syn_Bigbio), coefficients(AE.Syn_Wageningen),coefficients(AE.Syn_Texas),intervals(AE.Syn_Agro,which = "fixed")$fixed[2],intervals(AE.Syn_Biodepth,which = "fixed")$fixed[2],coefficients(AE.Syn_Jena)),
                   lower =  c( confint(AE.Syn_Bigbio)[1], confint(AE.Syn_Wageningen)[1],confint(AE.Syn_Texas)[1],intervals(AE.Syn_Agro,which = "fixed")$fixed[1],intervals(AE.Syn_Biodepth,which = "fixed")$fixed[1],confint(AE.Syn_Jena)[1]))
f2.CI<- data.frame(Item = c("SE.Syn"),
                   system = c("Bigbio","Wageningen","Texas","Agro","Biodepth","Jena"),
                   upper = c( confint(SE.Syn_Bigbio)[2], confint(SE.Syn_Wageningen,which = "fixed")[2],confint(SE.Syn_Texas)[2],intervals(SE.Syn_Agro,which = "fixed")$fixed[3],intervals(SE.Syn_Biodepth,which = "fixed")$fixed[3],confint(SE.Syn_Jena)[2]),
                   mean =  c(coefficients(SE.Syn_Bigbio), coefficients(SE.Syn_Wageningen,which = "fixed"),coefficients(SE.Syn_Texas),intervals(SE.Syn_Agro,which = "fixed")$fixed[2],intervals(SE.Syn_Biodepth,which = "fixed")$fixed[2],coefficients(SE.Syn_Jena)),
                   lower =  c(confint(SE.Syn_Bigbio)[1], confint(SE.Syn_Wageningen,which = "fixed")[1],confint(SE.Syn_Texas)[1],intervals(SE.Syn_Agro,which = "fixed")$fixed[1],intervals(SE.Syn_Biodepth,which = "fixed")$fixed[1],confint(SE.Syn_Jena)[1]))
Fig4f.CI<- rbind(f1.CI,f2.CI)

Fig4f<- ggplot(Fig4f.df) +
  geom_violin( aes(x = Item, y = f_x(Value, 2)/2),linewidth = 0.3, fill = "#F2F2F2", color = "black", trim = F, width=0.6) +

  geom_errorbar(data = Fig4f.CI, aes(x = Item, y =f_x(mean, 2)/2 , ymax =f_x(upper, 2)/2 , ymin = f_x(lower, 2)/2, color = system),linewidth = 0.3, width = 0.02,position = position_jitter(width = 0.1, height = 0,seed =2))+
  geom_point(data = Fig4f.CI, aes(x = Item, y = f_x(mean, 2)/2, color = system), size = 1.2, position = position_jitter(width = 0.1, height = 0,seed =2))+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values =  my6)+ 
  ylab("")+
  xlab("")+
  theme_bw()+
  scale_y_continuous(limits=c(-2.5, 2.5), labels=scaleFUN)+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")
a.blank<- ggplot()
Fig.4<- (Fig4abc+a.blank+plot_layout(widths = c(1, 0.01)))/(Fig4d+Fig4e+Fig4f)

# The statistical results of Fig. 4

summary(lme(NBE.CV ~1, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme(NIDE ~1, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme(NIIE ~1, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme(NIIE.SpVar ~1, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme(NIIE.Syn ~1, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme(NIDE.SpVar ~1, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme(NIDE.Syn ~1, random = ~1|(Study/Exp), BiodivExps_results))

summary(lme(f_x(AE.SpVar,2) ~1, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme(f_x(SE.SpVar,2) ~1, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme(f_x(AE.Syn,2) ~1, random = ~1|(Study/Exp), BiodivExps_results))
summary(lme(f_x(SE.Syn,2) ~1, random = ~1|(Study/Exp), BiodivExps_results))

#Fig.5
modelO_E<- (lme(CV.o~CV.e, random = ~1|(Study/Exp/Richness), BiodivExps_results))
predict.modelO_E<- predict(modelO_E, BiodivExps_results , level = 0,se.fit=T)
modelO_E_Bigbio<-(lme(CV.o~CV.e, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
modelO_E_Wageningen<-(lme(CV.o~CV.e, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
modelO_E_Texas<-(lme(CV.o~CV.e, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Texas",]))
modelO_E_Agro<-(lme(CV.o~CV.e, random = ~1|(Exp/Richness), BiodivExps_results[BiodivExps_results$Study == "Agro",]))
modelO_E_Biodepth<-(lme(CV.o~CV.e, random = ~1|(Exp/Richness), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
modelO_E_Jena<-(lme(CV.o~CV.e, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(modelO_E)
rsquared(modelO_E)
intervals(modelO_E, which= "fixed")
Fig.5<- ggplot()+ 
  geom_point(data = BiodivExps_results, aes(x = CV.e, y = CV.o, color=Study), size=2.5, shape=21, alpha = 0.4)+
  ylim(0, 2)+
  geom_abline(intercept = 0, slope = 1, linetype = 1, linewidth = 1, color = "grey")+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = CV.e, y=predict(modelO_E_Wageningen, level = 0)), linewidth = 0.9,color= my6[6],linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = CV.e, y=predict(modelO_E_Texas, level = 0)), linewidth = 0.9,color= my6[5], linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Agro",],aes(x = CV.e, y=predict(modelO_E_Agro, level = 0)), linewidth = 0.9,color= my6[1], linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = CV.e, y=predict(modelO_E_Biodepth, level = 0)), linewidth = 0.9,color= my6[3], linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = CV.e, y=predict(modelO_E_Jena, level = 0)), linewidth = 0.9,color= my6[4], linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = CV.e, y=predict(modelO_E_Bigbio, level = 0)), linewidth = 0.9,color= my6[2], linetype = 1)+
  geom_line(data = BiodivExps_results, aes(x = CV.e, y=predict.modelO_E$fit), linewidth = 1.2, linetype = 1)+
  geom_ribbon(data = BiodivExps_results, aes(x=CV.e, ymin=predict.modelO_E$fit-1.96*predict.modelO_E$se.fit ,
                                             ymax=predict.modelO_E$fit+1.96*predict.modelO_E$se.fit),alpha = 0.18,linetype=0)+
  scale_alpha_continuous(range = c(1, 1))+
  ylab("Observed community variability")+
  xlab("Expected community variability")+
  scale_color_manual(values = my6) +
  theme_bw()+
  theme(legend.text = element_text(size = 14),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        legend.key.size = unit(1, 'cm'),
        legend.position="right")

#Fig.6
#Here, NIDE need to be transformed to natural scale
NBE.F.model<-(lme(exp(NIDE)~NBE.F, random = ~1|(Study/Exp/Richness), BiodivExps_results))
predict.NBE.F.model<- predict(NBE.F.model, BiodivExps_results , level = 0,se.fit=T)
NBE.F_Bigbio<-(lme(exp(NIDE)~NBE.F, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NBE.F_Wageningen<-(lme(exp(NIDE)~NBE.F, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NBE.F_Texas<-(lme(exp(NIDE)~NBE.F, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NBE.F_Agro<-(lme(exp(NIDE)~NBE.F, random = ~1|(Exp/Richness), BiodivExps_results[BiodivExps_results$Study == "Agro",]))
NBE.F_Biodepth<-(lme(exp(NIDE)~NBE.F, random = ~1|(Exp/Richness), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NBE.F_Jena<-(lme(exp(NIDE)~NBE.F, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Jena",]))
summary(NBE.F.model)
anova(NBE.F.model)
intervals(NBE.F.model,which = "fixed")
Fig6a<- ggplot()+ 
  geom_point(data = BiodivExps_results, aes(x = NBE.F, y =exp(NIDE), color=Study), size=1.4, shape=21, alpha = 0.4)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = NBE.F, y=predict(NBE.F_Bigbio, level = 0)), linewidth = 0.7,color= my6[2], linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = NBE.F, y=predict(NBE.F_Wageningen, level = 0)), linewidth = 0.7,color= my6[6],linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = NBE.F, y=predict(NBE.F_Texas, level = 0)), linewidth = 0.7,color= my6[5],linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Agro",],aes(x = NBE.F, y=predict(NBE.F_Agro, level = 0)), linewidth = 0.7,color= my6[1],linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = NBE.F, y=predict(NBE.F_Biodepth, level = 0)), linewidth = 0.7,color= my6[3],linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = NBE.F, y=predict(NBE.F_Jena, level = 0)), linewidth = 0.7,color= my6[4],linetype = 2)+
  geom_ribbon(data = BiodivExps_results, aes(x=NBE.F, ymin=predict.NBE.F.model$fit-1.96*predict.NBE.F.model$se.fit ,
                                             ymax=predict.NBE.F.model$fit+1.96*predict.NBE.F.model$se.fit),alpha = 0.15,linetype=0)+
  geom_line(data = BiodivExps_results, aes(x = NBE.F, y=predict.NBE.F.model$fit), linewidth = 0.8, linetype = 1)+
  ylab(expression(NIDE))+
  xlab(expression(NBE[Bio]))+
  scale_color_manual(values = my6) +
  scale_y_continuous(limits=c(0, 4), labels=scaleFUN)+
  scale_x_continuous(limits=c(0.5, 6.1), labels=scaleFUN)+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.key.size = unit(1, 'cm'),
        legend.position="")

AE.NIDE.SpVar.F<-(lme(AE.SpVar~CE, random = ~1|(Study/Exp/Richness), BiodivExps_results))
summary(AE.NIDE.SpVar.F)
rsquared(AE.NIDE.SpVar.F)
intervals(AE.NIDE.SpVar.F,which = "fixed")
predict.AE.NIDE.SpVar.F<-predict(AE.NIDE.SpVar.F, BiodivExps_results , level = 0,se.fit=T)
AE.NIDE.SpVar.F_Bigbio<-(lme(AE.SpVar~CE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
AE.NIDE.SpVar.F_Wageningen<-(lme(AE.SpVar~CE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study  == "Wageningen",]))
AE.NIDE.SpVar.F_Texas<-(lme(AE.SpVar~CE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study  == "Texas",]))
AE.NIDE.SpVar.F_Agro<-(lme(AE.SpVar~CE, random = ~1|(Exp/Richness), BiodivExps_results[BiodivExps_results$Study  == "Agro",]))
AE.NIDE.SpVar.F_Biodepth<-(lme(AE.SpVar~CE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study  == "Biodepth",]))
AE.NIDE.SpVar.F_Jena<-(lme(AE.SpVar~CE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study  == "Jena",]))

Fig6b<- ggplot()+ 
  geom_point(data = BiodivExps_results, aes(x = CE, y = AE.SpVar, color=Study ),shape = 21, size=1.4, alpha = 0.4)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = CE, y=predict(AE.NIDE.SpVar.F_Bigbio, level = 0)), linewidth = 0.7,color= my6[2], linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = CE, y=predict(AE.NIDE.SpVar.F_Wageningen, level = 0)), linewidth = 0.7,color= my6[6], linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = CE, y=predict(AE.NIDE.SpVar.F_Texas, level = 0)), linewidth = 0.7,color= my6[5], linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Agro",],aes(x = CE, y=predict(AE.NIDE.SpVar.F_Agro, level = 0)), linewidth = 0.7,color= my6[1], linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = CE, y=predict(AE.NIDE.SpVar.F_Biodepth, level = 0)), linewidth = 0.7,color= my6[3], linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = CE, y=predict(AE.NIDE.SpVar.F_Jena, level = 0)), linewidth = 0.7,color= my6[4], linetype = 2)+
  geom_ribbon(data = BiodivExps_results,aes(x=CE, ymin=predict.AE.NIDE.SpVar.F$fit-1.96*predict.AE.NIDE.SpVar.F$se.fit ,
                                            ymax=predict.AE.NIDE.SpVar.F$fit+1.96*predict.AE.NIDE.SpVar.F$se.fit), alpha = 0.15,linetype=0)+
  geom_line(data = BiodivExps_results,aes(x = CE, y=predict.AE.NIDE.SpVar.F$fit), linewidth = 0.8, linetype = 1)+
  ylab(expression(AE[SpVar]))+
  xlab(expression(CE[Bio]))+
  scale_alpha_continuous(range = c(0.2, 1))+
  scale_y_continuous(limits=c(-0.8, 4.3), labels=scaleFUN)+
  scale_x_continuous(limits=c(-0.8, 7.1), labels=scaleFUN)+
  scale_color_manual(values = my6) +
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.key.size = unit(1, 'cm'),
        legend.position="")

SE.NIDE.SpVar.F<-(lme(SE.SpVar~SE, random = ~1|(Study/Exp/Richness), BiodivExps_results))
predict.SE.NIDE.SpVar.F<-predict(SE.NIDE.SpVar.F, BiodivExps_results , level = 0,se.fit=T)
summary(SE.NIDE.SpVar.F)
rsquared(SE.NIDE.SpVar.F)
intervals(SE.NIDE.SpVar.F, which = "fixed")
SE.NIDE.SpVar.F_Bigbio<-(lme(SE.SpVar~SE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
SE.NIDE.SpVar.F_Wageningen<-(lme(SE.SpVar~SE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
SE.NIDE.SpVar.F_Texas<-(lme(SE.SpVar~SE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Texas",]))
SE.NIDE.SpVar.F_Agro<-(lme(SE.SpVar~SE, random = ~1|(Exp/Richness),BiodivExps_results[BiodivExps_results$Study == "Agro",]))
SE.NIDE.SpVar.F_Biodepth<-(lme(SE.SpVar~SE, random = ~1|(Exp/Richness), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
SE.NIDE.SpVar.F_Jena<-(lme(SE.SpVar~SE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Jena",]))

Fig6c<- ggplot()+ 
  geom_point(data = BiodivExps_results, aes(x = SE, y = SE.SpVar, color=Study),shape = 21, size=1.4, alpha = 0.4)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = SE, y=predict(SE.NIDE.SpVar.F_Biodepth, level = 0)), linewidth = 0.7,color= my6[3],linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = SE, y=predict(SE.NIDE.SpVar.F_Bigbio, level = 0)), linewidth = 0.7,color= my6[2],linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = SE, y=predict(SE.NIDE.SpVar.F_Texas, level = 0)), linewidth = 0.7,color= my6[5],linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Agro",],aes(x = SE, y=predict(SE.NIDE.SpVar.F_Agro, level = 0)), linewidth = 0.7,color= my6[1],linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = SE, y=predict(SE.NIDE.SpVar.F_Wageningen, level = 0)), linewidth = 0.7,color= my6[6], linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = SE, y=predict(SE.NIDE.SpVar.F_Jena)), linewidth = 0.7,color= my6[4], alpha = 0.8,linetype = 1)+
  geom_ribbon(data = BiodivExps_results, aes(x=SE, ymin=predict.SE.NIDE.SpVar.F$fit-1.96*predict.SE.NIDE.SpVar.F$se.fit ,
                                             ymax=predict.SE.NIDE.SpVar.F$fit+1.96*predict.SE.NIDE.SpVar.F$se.fit), alpha = 0.15,linetype=0)+
  geom_line(data = BiodivExps_results, aes(x = SE, y=predict.SE.NIDE.SpVar.F$fit), linewidth = 0.8, linetype = 1)+
  ylab(expression(SE[SpVar]))+
  xlab(expression(SE[Bio]))+
  scale_alpha_continuous(range = c(0.2, 1))+
  scale_y_continuous(limits=c(-0.6, 1), labels=scaleFUN)+
  scale_x_continuous(limits=c(-4.1, 2.1), labels=scaleFUN)+
  scale_color_manual(values = my6) +
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.key.size = unit(1, 'cm'),
        legend.position="")

AE.Syn.F<-(lme(AE.Syn~CE,  random = ~1|(Study/Exp/Richness), BiodivExps_results))
predict.AE.Syn.F<-predict(AE.Syn.F, BiodivExps_results , level = 0,se.fit=T)
summary(AE.Syn.F)
rsquared(AE.Syn.F)
intervals(AE.Syn.F,which = "fixed")
AE.Syn.F_Bigbio<-(lme(AE.Syn~CE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
AE.Syn.F_Wageningen<-(lme(AE.Syn~CE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study== "Wageningen",]))
AE.Syn.F_Texas<-(lme(AE.Syn~CE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Texas",]))
AE.Syn.F_Agro<-(lme(AE.Syn~CE, random = ~1|(Exp/Richness), BiodivExps_results[BiodivExps_results$Study == "Agro",]))
AE.Syn.F_Biodepth<-(lme(AE.Syn~CE, random = ~1|(Exp/Richness), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
AE.Syn.F_Jena<-(lme(AE.Syn~CE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Jena",]))

Fig6d<- ggplot()+ 
  geom_point(data = BiodivExps_results, aes(x = CE, y = AE.Syn, color=Study),shape = 21, size=1.4, alpha = 0.4)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = CE, y=predict(AE.Syn.F_Biodepth, level = 0)), linewidth = 0.7,color= my6[3],linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = CE, y=predict(AE.Syn.F_Bigbio, level = 0)), linewidth = 0.7,color= my6[2], linetype = 21)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = CE, y=predict(AE.Syn.F_Texas, level = 0)), linewidth = 0.7,color= my6[5],linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Agro",],aes(x =CE, y=predict(AE.Syn.F_Agro, level = 0)), linewidth = 0.7,color= my6[1],linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = CE, y=predict(AE.Syn.F_Wageningen, level = 0)), linewidth = 0.7,color= my6[6],linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = CE, y=predict(AE.Syn.F_Jena, level = 0)), linewidth = 0.7,color= my6[4],linetype = 2)+
  geom_ribbon(data = BiodivExps_results, aes(x=CE, ymin=predict.AE.Syn.F$fit-1.96*predict.AE.Syn.F$se.fit ,
                                             ymax=predict.AE.Syn.F$fit+1.96*predict.AE.Syn.F$se.fit), alpha = 0.15,linetype=0)+
  geom_line(data = BiodivExps_results, aes(x = CE, y=predict.AE.Syn.F$fit), linewidth = 0.8, linetype = 2)+
  ylab(expression(AE[Syn]))+
  xlab(expression(CE[Bio]))+
  scale_y_continuous(limits=c(-2.5, 4.5), labels=scaleFUN)+
  scale_x_continuous(limits=c(-0.8, 7.1), labels=scaleFUN)+
  scale_color_manual(values = my6) +
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.key.size = unit(1, 'cm'),
        legend.position="")


SE.Syn.F<-(lme(SE.Syn~SE, random = ~1|(Study/Exp/Richness), BiodivExps_results))
predict.SE.Syn.F<-predict(SE.Syn.F, BiodivExps_results , level = 0,se.fit=T)
summary(SE.Syn.F)
rsquared(SE.Syn.F)
intervals(SE.Syn.F,which = "fixed")
SE.Syn.F_Bigbio<-(lme(SE.Syn~SE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
SE.Syn.F_Wageningen<-(lme(SE.Syn~SE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
SE.Syn.F_Texas<-(lme(SE.Syn~SE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Texas",]))
SE.Syn.F_Agro<-(lme(SE.Syn~SE, random = ~1|(Exp/Richness), BiodivExps_results[BiodivExps_results$Study == "Agro",]))
SE.Syn.F_Biodepth<-(lme(SE.Syn~SE, random = ~1|(Exp/Richness), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
SE.Syn.F_Jena<-(lme(SE.Syn~SE, random = ~1|(Richness), BiodivExps_results[BiodivExps_results$Study == "Jena",]))

Fig6e<- ggplot()+ 
  geom_point(data = BiodivExps_results, aes(x = SE, y = SE.Syn, color=Study),shape = 21, size=1.4, alpha = 0.4)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = SE, y=predict(SE.Syn.F_Biodepth, level = 0)), linewidth = 0.7,color= my6[3],linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = SE, y=predict(SE.Syn.F_Bigbio, level = 0)), linewidth = 0.7,color= my6[2],linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = SE, y=predict(SE.Syn.F_Texas, level = 0)), linewidth = 0.7,color= my6[5],linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Agro",],aes(x =SE, y=predict(SE.Syn.F_Agro, level = 0)), linewidth = 0.7,color= my6[1],linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = SE, y=predict(SE.Syn.F_Wageningen, level = 0)), linewidth = 0.7,color= my6[6],linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = SE, y=predict(SE.Syn.F_Jena, level = 0)), linewidth = 0.7,color= my6[4], linetype = 2)+
  geom_ribbon(data = BiodivExps_results, aes(x=SE, ymin=predict.SE.Syn.F$fit-1.96*predict.SE.Syn.F$se.fit ,
                                             ymax=predict.SE.Syn.F$fit+1.96*predict.SE.Syn.F$se.fit), alpha = 0.15,linetype=0)+
  geom_line(data = BiodivExps_results, aes(x = SE, y=predict.SE.Syn.F$fit), linewidth = 0.8, linetype = 1)+
  ylab(expression(SE[Syn]))+
  xlab(expression(SE[Bio]))+
  scale_y_continuous(limits=c(-1, 5), labels=scaleFUN)+
  scale_x_continuous(limits=c(-4.1, 2.1), labels=scaleFUN)+
  scale_color_manual(values = my6) +
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.key.size = unit(1, 'cm'),
        legend.position="")
blank<-ggplot()
Fig6<-Fig6a+Fig6b+Fig6c+blank+Fig6d+Fig6e



#Extended data fig.5
world.polygon<- readOGR( ".../world_polygon_show.shp")
site<- read.csv(".../site_infro.csv")
#site<- site_infro
plot(world.polygon)
world.map<-fortify(world.polygon)
world.map2<-world.map%>% group_by(long,order,hole,piece,id,group)%>% summarise(lat = max(lat))
ED_fig.5<- ggplot(world.map)+geom_polygon(aes(x = long, y = lat, group= group),  fill = "#F7F7F7", colour="black",size = 0.2, linetype = 3) +
  geom_point(data = site, aes(x = longitude, y = latitude, color = Study),size = 2, alpha = 0.9)+
  scale_color_manual(values = my6)+
  xlim(-130,60)+
  ylim(10,85)+
  xlab("longitude")+
  ylab("latitude")+
  theme_bw()+
  theme(    panel.border = element_blank(),
            legend.text = element_text(size = 12),
            axis.title.x=element_text(vjust=0.3, size=16),
            axis.title.y=element_text(hjust=0.5, vjust=2,size=16),
            axis.text.x=element_text(vjust=0,size=14),
            axis.text.y=element_text(hjust=0,size=14),
            legend.key.size = unit(1, 'cm'),
            legend.position="")


#Fig.S1
#AgroDiversity is excluded from this analysis as it does not have a richness gradient.

NBE.R_Bigbio<-lm(NBE.CV~Richness, BiodivExps_results[BiodivExps_results$Study == "Bigbio",])
NBE.R_Wageningen<-(lm(NBE.CV~Richness,BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NBE.R_Texas<-(lm(NBE.CV~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NBE.R_Biodepth<-(lme(NBE.CV~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NBE.R_Jena<-(lm(NBE.CV~Richness, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1a<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = Richness))+ 
  geom_point(aes(x = Richness, y = NBE.CV, color=Study), size=2, shape = 21)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-3.6, 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(NBE.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(NBE.R_Wageningen)), linewidth = 0.8,color= my6[6],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(NBE.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(NBE.R_Biodepth), group =Exp), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(NBE.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 1)+
  ylab(expression(NBE[CV]~(log) ))+
  xlab("")+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values = my6[-1]) +
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="none")

NIIE.R_Bigbio<-(lm(NIIE~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIIE.R_Wageningen<-(lm(NIIE~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIIE.R_Texas<-(lm(NIIE~Richness,  BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIIE.R_Biodepth<-(lme(NIIE~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIIE.R_Jena<-(lm(NIIE~Richness, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1b<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = NIIE))+ 
  geom_point(aes(x = Richness, y = NIIE, color=Study), size=2, shape = 21, alpha = 0.5)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-1.8, 0.3)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(NIIE.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(NIIE.R_Wageningen)), linewidth = 0.8, color= my6[6],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(NIIE.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(NIIE.R_Biodepth), group =Exp ), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(NIIE.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 1)+
  ylab(expression(NIIE~(log)) )+
  xlab("")+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values = my6[-1]) +
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="none")

NIDE.R_Bigbio<-(lm(NIDE~Richness, BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIDE.R_Wageningen<-(lm(NIDE~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIDE.R_Texas<-(lm(NIDE~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIDE.R_Biodepth<-(lme(NIDE~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIDE.R_Jena<-(lm(NIDE~Richness, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1c<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = NIDE))+ 
  geom_point(aes(x = Richness, y =NIDE, color=Study), size=2, shape = 21, shape = 0.5)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-4, 2.5)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(NIDE.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(NIDE.R_Wageningen)), linewidth = 0.8,color= my6[6],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(NIDE.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(NIDE.R_Biodepth), group =Exp ), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(NIDE.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 2)+
  ylab(expression(NIDE~(log)))+
  scale_color_manual(values = my6[-1]) +
  xlab("")+
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="none")

NIIE.SpVar.R_Bigbio<-(lm(NIIE.SpVar~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIIE.SpVar.R_Wageningen<-(lm(NIIE.SpVar~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIIE.SpVar.R_Texas<-(lm(NIIE.SpVar~Richness,  BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIIE.SpVar.R_Biodepth<-(lme(NIIE.SpVar~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIIE.SpVar.R_Jena<-(lm(NIIE.SpVar~Richness, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1d<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = NIIE.SpVar))+ 
  geom_point(aes(x = Richness, y = NIIE.SpVar, color=Study), size=2, shape = 21, alpha = 0.2)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-0.8, 0.2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(NIIE.SpVar.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(NIIE.SpVar.R_Wageningen)), linewidth = 0.8,color= my6[6],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(NIIE.SpVar.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(NIIE.SpVar.R_Biodepth), group =Exp ), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(NIIE.SpVar.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 2)+
  ylab(expression(NIIE[SpVar]~(log)))+
  xlab("")+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values = my6[-1]) +
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="none")

NIIE.Syn.R_Bigbio<-(lm(NIIE.Syn~Richness, BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIIE.Syn.R_Wageningen<-(lm(NIIE.Syn~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIIE.Syn.R_Texas<-(lm(NIIE.Syn~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIIE.Syn.R_Biodepth<-(lme(NIIE.Syn~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIIE.Syn.R_Jena<-(lm(NIIE.Syn~Richness, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1e<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = NIIE.Syn))+ 
  geom_point(aes(x = Richness, y = NIIE.Syn, color=Study), size=2,shape = 21, alpha = 0.5)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-1.6, 0.1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(NIIE.Syn.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(NIIE.Syn.R_Wageningen)), linewidth = 0.8,color= my6[6],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(NIIE.Syn.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(NIIE.Syn.R_Biodepth), group =Exp), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(NIIE.Syn.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 1)+
  ylab(expression(NIIE[Syn]~(log)))+
  xlab("")+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values = my6[-1]) +
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="none")

NIDE.SpVar.R_Bigbio<-(lm(NIDE.SpVar~Richness, BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIDE.SpVar.R_Wageningen<-(lm(NIDE.SpVar~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIDE.SpVar.R_Texas<-(lm(NIDE.SpVar~Richness,BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIDE.SpVar.R_Biodepth<-(lme(NIDE.SpVar~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIDE.SpVar.R_Jena<-(lm(NIDE.SpVar~Richness, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1f<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = NIDE.SpVar))+ 
  geom_point(aes(x = Richness, y = NIDE.SpVar, color=Study), size=2,shape = 21, alpha = 0.5)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-2.2, 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(NIDE.SpVar.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(NIDE.SpVar.R_Wageningen)), linewidth = 0.8,color= my6[6],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(NIDE.SpVar.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(NIDE.SpVar.R_Biodepth), group =Exp ), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(NIDE.SpVar.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 1)+
  ylab(expression(NIDE[SpVar]~(log)))+
  xlab("")+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values = my6[-1]) +
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="")

NIDE.Syn.R_Bigbio<-(lm(NIDE.Syn~Richness, BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
NIDE.Syn.R_Wageningen<-(lm(NIDE.Syn~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
NIDE.Syn.R_Texas<-(lm(NIDE.Syn~Richness,  BiodivExps_results[BiodivExps_results$Study == "Texas",]))
NIDE.Syn.R_Biodepth<-(lme(NIDE.Syn~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
NIDE.Syn.R_Jena<-(lm(NIDE.Syn~Richness, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1g<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = NIDE.Syn))+ 
  geom_point(aes(x = Richness, y = NIDE.Syn, color=Study), size=2,shape = 21, alpha = 0.5)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-2.2, 1.2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(NIDE.Syn.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(NIDE.Syn.R_Wageningen)), linewidth = 0.8,color= my6[6],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(NIDE.Syn.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(NIDE.Syn.R_Biodepth), group =Exp ), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(NIDE.Syn.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 1)+
  ylab(expression(NIDE[Syn]~(log)))+
  xlab("")+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values = my6[-1]) +
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="none")

AE.SpVar.N<-(lme(AE.SpVar~Richness, random = ~1|(Exp), BiodivExps_results))
AE.SpVar.R_Bigbio<-(lm(AE.SpVar~Richness, BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
AE.SpVar.R_Wageningen<-(lm(AE.SpVar~Richness,  BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
AE.SpVar.R_Texas<-(lm(AE.SpVar~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
AE.SpVar.R_Biodepth<-(lme(AE.SpVar~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
AE.SpVar.R_Jena<-(lm(AE.SpVar~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1h<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = AE.SpVar))+ 
  geom_point(aes(x = Richness, y = AE.SpVar, color=Study), size=2,shape = 21, alpha = 0.5)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-1, 4)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(AE.SpVar.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(AE.SpVar.R_Wageningen)), linewidth = 0.8,color= my6[6],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(AE.SpVar.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(AE.SpVar.R_Biodepth), group =Exp ), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(AE.SpVar.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 1)+
  ylab(expression(AE[SpVar]))+
  xlab("Species richness")+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values = my6[-1]) +
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="none")

SE.SpVar.R_Bigbio<-(lm(SE.SpVar~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
SE.SpVar.R_Wageningen<-(lm(SE.SpVar~Richness,  BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
SE.SpVar.R_Texas<-(lm(SE.SpVar~Richness,  BiodivExps_results[BiodivExps_results$Study == "Texas",]))
SE.SpVar.R_Biodepth<-(lme(SE.SpVar~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
SE.SpVar.R_Jena<-(lm(SE.SpVar~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1i<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = SE.SpVar))+ 
  geom_point(aes(x = Richness, y = SE.SpVar, color=Study), size=2,shape = 21, alpha = 0.5)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-0.6, 0.8)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(SE.SpVar.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(SE.SpVar.R_Wageningen)), linewidth = 0.8,color= my6[6],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(SE.SpVar.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(SE.SpVar.R_Biodepth), group =Exp ), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(SE.SpVar.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 2)+
  ylab(expression(SE[SpVar]))+
  xlab("Species richness")+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values = my6[-1]) +
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="none")

AE.Syn.R_Bigbio<-(lm(AE.Syn~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
AE.Syn.R_Wageningen<-(lm(AE.Syn~Richness,  BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
AE.Syn.R_Texas<-(lm(AE.Syn~Richness,  BiodivExps_results[BiodivExps_results$Study == "Texas",]))
AE.Syn.R_Biodepth<-(lme(AE.Syn~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
AE.Syn.R_Jena<-(lm(AE.Syn~Richness,  BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1j<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = AE.Syn))+ 
  geom_point(aes(x = Richness, y = AE.Syn, color=Study), size=2,shape = 21, alpha = 0.5)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-2.5, 3)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(AE.Syn.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(AE.Syn.R_Wageningen)), linewidth = 0.8,color= my6[6],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(AE.Syn.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(AE.Syn.R_Biodepth), group =Exp ), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(AE.Syn.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 1)+
  ylab(expression(AE[Syn]))+
  xlab("Species richness")+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values = my6[-1]) +
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="none")


SE.Syn.R_Bigbio<-(lm(SE.Syn~Richness,  BiodivExps_results[BiodivExps_results$Study == "Bigbio",]))
SE.Syn.R_Wageningen<-(lm(SE.Syn~Richness, BiodivExps_results[BiodivExps_results$Study == "Wageningen",]))
SE.Syn.R_Texas<-(lm(SE.Syn~Richness, BiodivExps_results[BiodivExps_results$Study == "Texas",]))
SE.Syn.R_Biodepth<-(lme(SE.Syn~Richness, random = ~1|(Exp), BiodivExps_results[BiodivExps_results$Study == "Biodepth",]))
SE.Syn.R_Jena<-(lm(SE.Syn~Richness, BiodivExps_results[BiodivExps_results$Study == "Jena",]))

FigS1k<- ggplot(BiodivExps_results[BiodivExps_results$Study != "Agro",], aes(x = Richness, y = SE.Syn))+ 
  geom_point(aes(x = Richness, y = SE.Syn, color=Study), size=2,shape = 21, alpha = 0.5)+
  geom_hline(aes(yintercept = 0), linetype =2)+
  ylim(-1,7.4)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Bigbio",],aes(x = Richness, y=predict(SE.Syn.R_Bigbio)), linewidth = 0.8,color= my6[2],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Wageningen",],aes(x = Richness, y=predict(SE.Syn.R_Wageningen)), linewidth = 0.8,color= my6[6],alpha = 0.8, linetype = 1)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Texas",],aes(x = Richness, y=predict(SE.Syn.R_Texas)), linewidth = 0.8,color= my6[5],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Biodepth",],aes(x = Richness, y=predict(SE.Syn.R_Biodepth), group =Exp ), linewidth = 0.8,color= my6[3],alpha = 0.8, linetype = 2)+
  geom_line(data = BiodivExps_results[BiodivExps_results$Study == "Jena",],aes(x = Richness, y=predict(SE.Syn.R_Jena)), linewidth = 0.8,color= my6[4],alpha = 0.8, linetype = 1)+
  ylab(expression(SE[Syn]))+
  xlab("Species richness")+
  geom_hline(aes(yintercept = 0), linetype =2)+
  scale_color_manual(values = my6[-1]) +
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        legend.position="none")
blank<-ggplot()
FigS1<-FigS1a+FigS1b+FigS1c+blank+FigS1d+FigS1e+FigS1f+FigS1g+FigS1h+FigS1i+FigS1j+FigS1k+plot_layout( ncol = 4)
