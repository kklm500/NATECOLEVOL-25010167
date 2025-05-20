
library("nlme")
library("tidyverse")
library("MuMIn")
library("patchwork")
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
#Fig3a,b,c,d
sim_richness_D0<- read.csv(".../sim_results_Richness_D0.csv")
#NBE.CV, NIIE, NIIE.SpVar, NIIE.Syn,  NIDE, NIDE.SpVar, NIDE.Syn are logarithmically transformed
Fig3ab.df<- sim_results_richness_D0%>%
  na.omit()%>%
  mutate(NBE.CV= log(CV.o/CV.b),
         NIIE = log(CV.e/CV.b),
         NIDE = log(CV.o/CV.e),
         NIDE.SpVar = log(CVS.o/CVS.e),
         NIDE.Syn = log(phi.o/phi.e))%>%
  select(Richness, NBE.CV,NIIE,NIDE,NIDE.SpVar,NIDE.Syn)%>%
  gather(Item,Value,NBE.CV,NIIE,NIDE,NIDE.SpVar,NIDE.Syn)%>%
  mutate(Order= NA)
Fig3ab.df[Fig3ab.df$Item == "NBE.CV",]$Order <- "a"
Fig3ab.df[Fig3ab.df$Item %in% c("NIIE", "NIDE"),]$Order <- "b"
Fig3ab.df[Fig3ab.df$Item %in% c("NIDE.SpVar","NIDE.Syn"),]$Order <- "c"
Fig3ab.df$Item<- factor(Fig3ab.df$Item,levels = c("NBE.CV","NIIE","NIDE","NIDE.SpVar","NIDE.Syn"))
Fig3a.df_avg<- Fig3ab.df%>%filter(Item %in% c("NBE.CV", "NIIE", "NIDE"))%>% group_by(Richness,Item)%>% summarise(Value = mean(Value, na.rm = T))
Fig3a<- ggplot(Fig3a.df_avg)+geom_line(aes(x = Richness, y = Value, color = Item),linewidth= 0.7)+
  scale_color_manual(values = c("#3A3A3A","#571C72","#B86500"))+
  geom_hline(aes( yintercept= 0), linetype =2)+
  scale_y_continuous(limits=c(-1.4, 0.7), labels=scaleFUN)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("Richness")+
  theme_bw()+
  theme(panel.border = element_rect(color="black", size=0.4),
        panel.grid.major = element_line(linewidth = 0.1),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.position="")
Fig3b.df_avg<- Fig3ab.df%>%filter(Item %in% c( "NIDE","NIDE.SpVar","NIDE.Syn"))%>% group_by(Richness,Item)%>% summarise(Value = mean(Value, na.rm = T))
Fig3b<- ggplot(Fig3b.df_avg)+geom_line(aes(x = Richness, y = Value, color = Item),linewidth= 0.7)+
  scale_color_manual(values = c("#D27300","#8D0E21","#1A5602"))+
  geom_hline(aes( yintercept= 0), linetype =2)+
  scale_y_continuous(limits=c(-1.4, 0.7), labels=scaleFUN)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("Richness")+
  theme_bw()+
  theme(panel.border = element_rect(color="black", size=0.4),
        panel.grid.major = element_line(linewidth = 0.1),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.position="")

Fig3c.df<-sim_results_richness_D0%>%
  select(Richness, AE.SpVar, SE.SpVar)%>%
  gather(Item, Value, AE.SpVar, SE.SpVar)
Fig3c.df_avg<-Fig3c.df%>% group_by(Richness,Item)%>% summarise(Value = mean(Value, na.rm = T))
Fig3c<- ggplot(Fig3c.df_avg)+
  geom_hline(aes( yintercept= 0), linetype =2)+
  geom_line(aes(x = Richness, y = f_x(Value, 2) , color = Item),linewidth= 0.7)+
  scale_color_manual(values = c("#D40808","#F59A67"))+
  scale_y_continuous(limits=c(-1.4, 0.7), labels=scaleFUN)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("Richness")+
  theme_bw()+
  theme(panel.border = element_rect(color="black", size=0.4),
        panel.grid.major = element_line(linewidth = 0.1),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.position="")

Fig3d.df<-sim_results_richness_D0%>%
  select(Richness, AE.Syn, SE.Syn)%>%
  gather(Item, Value, AE.Syn, SE.Syn)
Fig3d.df_avg<- Fig3d.df%>% group_by(Richness,Item)%>% summarise(Value = mean(Value, na.rm = T))
Fig3d<- ggplot(Fig3d.df_avg)+geom_line(aes(x = Richness, y =  f_x(Value, 2)/2 , color = Item),linewidth= 0.7)+
  scale_color_manual(values =  c("#8E9E46","#3F8F80"))+
  geom_hline(aes( yintercept= 0), linetype =2)+
  scale_y_continuous(limits=c(-1.4, 0.7), labels=scaleFUN)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("Richness")+
  theme_bw()+
  theme(panel.border = element_rect(color="black", size=0.4),
        panel.grid.major = element_line(linewidth = 0.1),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.position="")
Fig3abcd<-Fig3a+Fig3b+Fig3c+Fig3d+plot_layout(ncol = 4)

#ED Fig 1
sim_results_Competition<- read.csv(".../sim_results_competition.csv")
EDfig1ab.df<-sim_results_competition%>%
  filter(Richness == 3)%>%
  na.omit()%>%
  mutate(NBE.CV= log(CV.o/CV.b),
         NIIE = log(CV.e/CV.b),
         NIDE = log(CV.o/CV.e),
         NIDE.SpVar = log(CVS.o/CVS.e),
         NIDE.Syn = log(phi.o/phi.e))%>%
  select(alpha,miu_d, NBE.CV,NIIE,NIDE,NIDE.SpVar,NIDE.Syn)%>%
  gather(Item,Value,NBE.CV,NIIE,NIDE,NIDE.SpVar,NIDE.Syn)%>%
  mutate(Order= NA)
EDfig1ab.df[EDfig1ab.df$Item == "NBE.CV",]$Order <- "a"
EDfig1ab.df[EDfig1ab.df$Item %in% c("NIIE", "NIDE"),]$Order <- "b"
EDfig1ab.df[EDfig1ab.df$Item %in% c("NIDE.SpVar","NIDE.Syn"),]$Order <- "c"
EDfig1ab.df$Item<- factor(EDfig1ab.df$Item,levels = c("NBE.CV","NIIE","NIDE","NIDE.SpVar","NIDE.Syn"))
EDfig1a.df_avg<- EDfig1ab.df%>%filter(Item %in% c("NBE.CV","NIIE","NIDE"))%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1a<- ggplot(EDfig1a.df_avg)+geom_line(aes(x = alpha, y = Value, color = Item, linetype=as.factor(miu_d) ),linewidth= 0.8)+
  scale_color_manual(values = c("#3A3A3A","#571C72","#B86500"))+
  geom_hline(aes( yintercept= 0), color = "grey")+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")

EDfig1b.df_avg<- EDfig1ab.df%>%filter(Item %in% c("NIDE", "NIDE.SpVar","NIDE.Syn"))%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1b<- ggplot(EDfig1b.df_avg)+
  geom_hline(aes( yintercept= 0), color = "grey")+
  geom_line(aes(x = alpha, y = Value, color = Item, linetype=as.factor(miu_d) ),linewidth= 0.8)+
  scale_color_manual(values =  c("#D27300","#8D0E21","#1A5602"))+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab("")+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")

EDfig1c.df<- sim_results_competition%>%
  filter(Richness == 3)%>%
  select(alpha,miu_d, AE.SpVar, SE.SpVar)%>%
  gather(Item, Value, AE.SpVar, SE.SpVar)
EDfig1c.df_avg<- EDfig1c.df%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1c<- ggplot(EDfig1c.df_avg)+
  geom_hline(aes( yintercept= 0), color = "grey")+
  geom_line(aes(x = alpha, y =f_x(Value, 2), color = Item, linetype=as.factor(miu_d)),linewidth= 0.8)+
  scale_color_manual(values = c("#D40808","#F59A67"))+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab("")+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")

EDfig1d.df<- sim_results_competition%>%
  filter(Richness == 3)%>%
  select(alpha,miu_d, AE.SpVar, SE.Syn)%>%
  gather(Item, Value, AE.SpVar, SE.Syn)
EDfig1d.df_avg<- EDfig1d.df%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1d<- ggplot(EDfig1d.df_avg)+
  geom_hline(aes( yintercept= 0), color = "grey")+
  geom_line(aes(x = alpha, y = f_x(Value, 2), color = Item, linetype=as.factor(miu_d)),linewidth= 0.8)+
  scale_color_manual(values =  c("#8E9E46","#3F8F80"))+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab("")+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")

EDfig1ef.df<- sim_results_competition%>%
  filter(Richness == 6)%>%
  na.omit()%>%
  mutate(NBE.CV= log(CV.o/CV.b),
         NIIE = log(CV.e/CV.b),
         NIDE = log(CV.o/CV.e),
         NIDE.SpVar = log(CVS.o/CVS.e),
         NIDE.Syn = log(phi.o/phi.e))%>%
  select(alpha,miu_d, NBE.CV,NIIE,NIDE,NIDE.SpVar,NIDE.Syn)%>%
  gather(Item,Value,NBE.CV,NIIE,NIDE,NIDE.SpVar,NIDE.Syn)%>%
  mutate(Order= NA)
EDfig1ef.df[EDfig1ef.df$Item == "NBE.CV",]$Order <- "a"
EDfig1ef.df[EDfig1ef.df$Item %in% c("NIIE", "NIDE"),]$Order <- "b"
EDfig1ef.df[EDfig1ef.df$Item %in% c("NIDE.SpVar","NIDE.Syn"),]$Order <- "c"
EDfig1ef.df$Item<- factor(EDfig1ef.df$Item,levels = c("NBE.CV","NIIE","NIDE","NIDE.SpVar","NIDE.Syn"))
EDfig1e.df_avg<- EDfig1ef.df%>%filter(Item %in% c("NBE.CV", "NIIE", "NIDE"))%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1e<- ggplot(EDfig1e.df_avg)+geom_line(aes(x = alpha, y = Value, color = Item, linetype=as.factor(miu_d) ),linewidth= 0.8)+
  scale_color_manual(values = c("#3A3A3A","#571C72","#B86500"))+
  geom_hline(aes( yintercept= 0), color = "grey")+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")
EDfig1f.df_avg<- EDfig1ef.df%>%filter(Item %in% c("NIDE", "NIDE.SpVar", "NIDE.Syn"))%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1f<- ggplot(EDfig1f.df_avg)+geom_line(aes(x = alpha, y = Value, color = Item, linetype=as.factor(miu_d) ),linewidth= 0.8)+
  scale_color_manual(values =  c("#D27300","#8D0E21","#1A5602"))+
  geom_hline(aes( yintercept= 0), color = "grey")+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab("")+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")
EDfig1g.df<- sim_results_competition%>%
  filter(Richness == 6)%>%
  select(alpha,miu_d, AE.SpVar, SE.SpVar)%>%
  gather(Item, Value, AE.SpVar, SE.SpVar)%>%
  mutate(scenario = "random")

EDfig1g_avg<- EDfig1g.df%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1g<- ggplot(EDfig1g_avg)+
  geom_hline(aes( yintercept= 0), color = "grey")+
  geom_line(aes(x = alpha, y =f_x(Value, 2), color = Item, linetype=as.factor(miu_d)),linewidth= 0.8)+
  scale_color_manual(values = c("#D40808","#F59A67"))+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab("")+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")

EDfig1h.df<- sim_results_competition%>%
  filter(Richness == 6)%>%
  select(alpha,miu_d, AE.Syn, SE.Syn)%>%
  gather(Item, Value, AE.Syn, SE.Syn)
EDfig1h.df_avg<-EDfig1h.df%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1h<- ggplot(EDfig1h.df_avg)+geom_line(aes(x = alpha, y = f_x(Value, 2), color = Item, linetype=as.factor(miu_d)),linewidth= 0.8)+
  scale_color_manual(values =  c("#8E9E46","#3F8F80"))+
  geom_hline(aes( yintercept= 0), color = "grey")+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab("")+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")

EDfig1ij.df<- sim_results_competition%>%
  filter(Richness == 9)%>%
  na.omit()%>%
  mutate(NBE.CV= log(CV.o/CV.b),
         NIIE = log(CV.e/CV.b),
         NIDE = log(CV.o/CV.e),
         NIDE.SpVar = log(CVS.o/CVS.e),
         NIDE.Syn = log(phi.o/phi.e))%>%
  select(alpha,miu_d, NBE.CV,NIIE,NIDE,NIDE.SpVar,NIDE.Syn)%>%
  gather(Item,Value,NBE.CV,NIIE,NIDE,NIDE.SpVar,NIDE.Syn)%>%
  mutate(Order= NA)
EDfig1ij.df[EDfig1ij.df$Item == "NBE.CV",]$Order <- "a"
EDfig1ij.df[EDfig1ij.df$Item %in% c("NIIE", "NIDE"),]$Order <- "b"
EDfig1ij.df[EDfig1ij.df$Item %in% c("NIDE.SpVar","NIDE.Syn"),]$Order <- "c"
EDfig1ij.df$Item<- factor(EDfig1ij.df$Item,levels = c("NBE.CV","NIIE","NIDE","NIDE.SpVar","NIDE.Syn"))
EDfig1i.df_avg<- EDfig1ij.df%>%filter(Item %in% c("NBE.CV", "NIIE", "NIDE"))%>% group_by(alpha, miu_d, Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1i<- ggplot(EDfig1i.df_avg)+geom_line(aes(x = alpha, y = Value, color = Item, linetype=as.factor(miu_d) ),linewidth= 0.8)+
  scale_color_manual(values = c("#3A3A3A","#571C72","#B86500"))+
  geom_hline(aes( yintercept= 0), color = "grey")+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab(expression(Biodiversity~effect[log]))+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")
EDfig1j.df_avg<- EDfig1ij.df%>%filter(Item %in% c("NIDE", "NIDE.SpVar", "NIDE.Syn"))%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1j<- ggplot(EDfig1j.df_avg)+geom_line(aes(x = alpha, y = Value, color = Item, linetype=as.factor(miu_d) ),linewidth= 0.8)+
  scale_color_manual(values =  c("#D27300","#8D0E21","#1A5602"))+
  geom_hline(aes( yintercept= 0), color = "grey")+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab("")+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")
EDfig1k.df<-  sim_results_competition%>%
  filter(Richness == 9)%>%
  select(alpha,miu_d, AE.SpVar, SE.SpVar)%>%
  gather(Item, Value, AE.SpVar, SE.SpVar)
EDfig1k_avg<- EDfig1k.df%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1k<- ggplot(EDfig1k_avg)+
  geom_hline(aes( yintercept= 0), color = "grey")+
  geom_line(aes(x = alpha, y = f_x(Value, 2), color = Item, linetype=as.factor(miu_d)),linewidth= 0.8)+
  scale_color_manual(values = c("#D40808","#F59A67"))+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab("")+
  xlab("The average intensity of interspecific competition")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")

EDfig1l.df<- sim_results_competition%>%
  filter(Richness == 9)%>%
  select(alpha,miu_d, AE.Syn, SE.Syn)%>%
  gather(Item, Value, AE.Syn,SE.Syn)
EDfig1l_avg<- EDfig1l.df%>% group_by(alpha,miu_d,Item)%>% summarise(Value = mean(Value, na.rm = T))
EDfig1l<- ggplot(EDfig1l_avg)+geom_line(aes(x = alpha, y = f_x(Value, 2), color = Item, linetype=as.factor(miu_d)),linewidth= 0.8)+
  scale_color_manual(values =  c("#8E9E46","#3F8F80"))+
  geom_hline(aes( yintercept= 0), color = "grey")+
  scale_y_continuous(limits=c(-2.1, 1.3), labels=scaleFUN)+
  ylab("")+
  xlab("")+
  xlim(0.01,0.39)+
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, vjust=2,size=12),
        axis.text.x=element_text(vjust=0,size=11),
        axis.text.y=element_text(hjust=0,size=11),
        legend.key.size = unit(1, 'cm'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="")

EDFig1<- EDfig1a+EDfig1b+EDfig1c+EDfig1d+EDfig1e+EDfig1f+EDfig1g+EDfig1h+EDfig1i+EDfig1j+EDfig1k+EDfig1l+plot_layout(ncol = 4)


#ED.Fig2a&b
sim_rK_correlation<- read.csv(".../sim_rK_correlation.csv")
ED.Fig2a<- ggplot()+ 
  geom_smooth(data = sim_rK_correlation, aes(x = cor_rK, y = SE.SpVar, color = as.factor(Richness)), linewidth = 0.9, method = "lm", fill = NA)+
  ylab(expression(SE[SpVar]))+
  xlab("The correlation between ri and Ki")+
  scale_color_manual(values = c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(vjust=0,size=14),
        axis.text.y=element_text(hjust=0,size=14),
        legend.key.size = unit(1, 'cm'),
        legend.position="")
ED.Fig2b<- ggplot()+ 
  geom_smooth(data = sim_rK_correlation, aes(x = cor_rK, y = SE.Syn, color = as.factor(Richness)), linewidth = 0.9, method = "lm", fill = NA)+
  ylab(expression(SE[Syn]))+
  xlab("")+
  scale_color_manual(values = c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(vjust=0,size=14),
        axis.text.y=element_text(hjust=0,size=14),
        legend.key.size = unit(1, 'cm'),
        legend.position="")
#ED.Fig2c&d
sim_ra_correlation<- read.csv(".../sim_ra_correlation.csv")
ED.Fig2c<- ggplot()+ 
  geom_smooth(data = sim_ra_correlation, aes(x = cor_ra_iE, y = SE.SpVar, color = as.factor(Richness)), linewidth = 0.9, method = "lm", fill = NA)+
  ylab(expression(SE[SpVar]))+
  xlab("The correlation between ri and ai_E")+
  scale_color_manual(values = c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(vjust=0,size=14),
        axis.text.y=element_text(hjust=0,size=14),
        legend.key.size = unit(1, 'cm'),
        legend.position="")

ED.Fig2d<- ggplot()+ 
  geom_smooth(data = sim_ra_correlation, aes(x = cor_ra_iE, y = SE.Syn, color = as.factor(Richness)), linewidth = 0.9, method = "lm", fill = NA)+
  ylab(expression(SE[Syn]))+
  xlab("")+
  scale_color_manual(values = c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(vjust=0,size=14),
        axis.text.y=element_text(hjust=0,size=14),
        legend.key.size = unit(1, 'cm'),
        legend.position="")
ED.Fig2<- ED.Fig2a+ED.Fig2b+ED.Fig2c+ED.Fig2d+plot_layout(ncol = 2)


#ED Fig3a-b, competition model without demographic stochasticity
pred.EDfig3a<- predict(lm(AE.SpVar~CE:as.factor(Richness)+as.factor(Richness)-1, sim_results_richness_D0), sim_results_richness_D0)
EDfig3a<- ggplot(sim_results_richness_D0)+ 
  geom_smooth( aes(x = CE, y = pred.EDfig3a, color = as.factor(Richness)),linewidth = 0.8, method = "lm", fill = NA)+
  ylab(expression(AE[SpVar]))+
  xlab(expression(CE[Bio]))+
  scale_color_manual(values = c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  scale_y_continuous(limits=c(-0.25, 0.6), labels=scaleFUN)+
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position="")
pred.EDfig3b<- predict(lm(AE.Syn~CE:as.factor(Richness)+as.factor(Richness)-1, sim_results_richness_D0), sim_results_richness_D0)
EDfig3b<- ggplot(sim_results_richness_D0)+ 
  geom_smooth( aes(x = CE, y = pred.EDfig3b, color = as.factor(Richness)), linewidth = 0.8, method = "lm", fill = NA)+
  ylab(expression(AE[Syn]))+
  xlab(expression(CE[Bio]))+
  scale_color_manual(values = c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  scale_y_continuous(limits=c(-0.8, 0.5),labels=scaleFUN)+
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position="")
pred.EDfig3c<- predict(lm(SE.SpVar~SE:as.factor(Richness)+as.factor(Richness)-1, sim_results_richness_D0), sim_results_richness_D0)
EDfig3c<- ggplot(sim_results_richness_D0)+ 
  geom_smooth( aes(x = SE, y = pred.EDfig3c, color = as.factor(Richness)), linewidth = 0.8, method = "lm", fill = NA)+
  ylab(expression(SE[SpVar]))+
  xlab(expression(SE[Bio]))+
  scale_y_continuous(limits=c(-0.15, 0.1), labels=scaleFUN)+
  scale_color_manual(values =  c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  #scale_color_tableau()+
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position="")

pred.EDfig3d<- predict(lm(SE.Syn~SE:as.factor(Richness)+as.factor(Richness)-1, sim_results_richness_D0), sim_results_richness_D0)
EDfig3d<- ggplot(sim_results_richness_D0)+ 
  geom_smooth( aes(x = SE, y = pred.EDfig3d, color = as.factor(Richness)), linewidth = 0.8, method = "lm", fill = NA)+
  ylab(expression(SE[Syn]))+
  xlab(expression(SE[Bio]))+
  scale_y_continuous(limits=c(-0.15, 0.55),labels=scaleFUN)+
  scale_color_manual(values =  c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  #scale_color_tableau()+
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position="")

#ED Fig3e-f, competition model with demographic stochasticity
sim_results_Richness_D1<- read.csv(".../sim_results_Richness_D1.csv")
pred.EDfig3e<- predict(lm(AE.SpVar~CE:as.factor(Richness)+as.factor(Richness)-1, sim_results_richness_D1), sim_results_richness_D1)
EDfig3e<- ggplot(sim_results_richness_D1)+ 
  geom_smooth(aes(x = CE, y = pred.EDfig3e, color = as.factor(Richness)),linewidth = 0.8, method = "lm", fill = NA)+
  ylab(expression(AE[SpVar]))+
  xlab(expression(CE[Bio]))+
  scale_color_manual(values = c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  scale_y_continuous(limits=c(-0.25, 0.6), labels=scaleFUN)+
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position="")

pred.EDfig3f<- predict(lm(AE.Syn~CE:as.factor(Richness)+as.factor(Richness)-1, sim_results_richness_D1), sim_results_richness_D1)
EDfig3f<- ggplot(sim_results_richness_D1)+ 
  geom_smooth( aes(x = CE, y =pred.EDfig3f, color = as.factor(Richness)), linewidth = 0.8, method = "lm", fill = NA)+
  ylab(expression(AE[Syn]))+
  xlab(expression(CE[Bio]))+
  scale_color_manual(values = c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  scale_y_continuous(limits=c(-0.8, 0.5),labels=scaleFUN)+
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position="")
pred.EDfig3g<- predict(lm(SE.SpVar~SE:as.factor(Richness)+as.factor(Richness)-1, sim_results_richness_D1), sim_results_richness_D1)
EDfig3g<- ggplot(sim_results_richness_D1)+ 
  geom_smooth( aes(x = SE, y = pred.EDfig3g, color = as.factor(Richness)), linewidth = 0.8, method = "lm", fill = NA)+
  ylab(expression(SE[SpVar]))+
  xlab(expression(SE[Bio]))+
  scale_y_continuous(limits=c(-0.15, 0.1), labels=scaleFUN)+
  scale_color_manual(values =  c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position="")

pred.EDfig3h<- predict(lm(SE.Syn~SE:as.factor(Richness)+as.factor(Richness)-1, sim_results_richness_D1), sim_results_richness_D1)
EDfig3h<- ggplot(sim_results_richness_D1)+ 
  geom_smooth( aes(x = SE, y = pred.EDfig3h, color = as.factor(Richness)), linewidth = 0.8, method = "lm", fill = NA)+
  ylab(expression(SE[Syn]))+
  xlab(expression(SE[Bio]))+
  scale_y_continuous(limits=c(-0.15, 0.55),labels=scaleFUN)+
  scale_color_manual(values =  c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  theme_bw()+
  theme(legend.text = element_text(size = 10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        legend.key.size = unit(1, 'cm'),
        legend.position="")

EDfig3<-EDfig3a+EDfig3b+EDfig3c+EDfig3d+EDfig3e+EDfig3f+EDfig3g+EDfig3h+plot_layout(ncol = 4)

#ED fig 4
#NIDE = CV.o/CV.e; natural scale
pred.EDfig4a<- predict(lm((CV.o/CV.e)~NBE.F:as.factor(Richness)+as.factor(Richness)-1, sim_results_richness_D0), sim_results_richness_D0)
EDfig4a<- ggplot(sim_results_richness_D0)+ 
  geom_smooth( aes(x = NBE.F , y = pred.EDfig4a, color = as.factor(Richness)),linewidth = 0.9, method = "lm", fill = NA)+
  ylab(expression(NIDE))+
  xlab(expression(NBE[Bio]))+
  scale_color_manual(values = c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  theme_bw()+
  theme(legend.text = element_text(size = 16),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x=element_text(vjust=0,size=16),
        axis.text.y=element_text(hjust=0,size=16),
        legend.key.size = unit(1, 'cm'),
        legend.position="")
pred.EDfig4b<- predict(lm( (CV.o/CV.e)~NBE.F:as.factor(Richness)+as.factor(Richness)-1, sim_results_richness_D1), sim_results_richness_D1)
EDfig4b<- ggplot(sim_results_richness_D1)+ 
  geom_smooth( aes(x = NBE.F, y = pred.EDfig4b, color = as.factor(Richness)), linewidth = 0.9, method = "lm", fill = NA)+
  ylab(expression(NIDE))+
  xlab(expression(NBE[Bio]))+
  scale_color_manual(values =  c("#E1E1E1","#B1B1B1","#B1B1B1","#A3A3A3","#777C7A","#575957","#404241","#222422","#0A0A0A"))+
  theme_bw()+
  theme(legend.text = element_text(size = 16),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x=element_text(vjust=0,size=16),
        axis.text.y=element_text(hjust=0,size=16),
        legend.key.size = unit(1, 'cm'),
        legend.position="")
EDfig4<- EDfig4a+EDfig4b



