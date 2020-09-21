#Time Point 1,  
#And shoot counts in June 2017
#And growth data from Fall 2016
# And Dry Mass from Winter 2017 (march/april) time point (added Oct 2018)
#final figure revision for Ecology July 2020

#Data Files Needed
#T1.txt (for count data)
#SM_T1_ALL.txt (for biomass and starch data)

rm(list=ls())
setwd("/Users/Katherine/Documents/Rworking/StressMem")

library(ggplot2)
library(doBy)
library(lme4)

####Initial look at Shoot Count Data, T1####

d<-read.delim(file="T1.txt", header=T)

d$ID <- as.character(d$ID)

Cw <- subset(d, d$Treatment=="H")
Cc <- subset(d, d$Treatment=="C")

hist(Cw$Count) #definitly not normal! 
hist(Cc$Count)

fun<-function(x,...){c(m=mean(x, na.rm=T), se=sd(x, na.rm=T)/sqrt(length(x)), n=length(x))}

sum<-summaryBy(Count ~ Treatment, data=d, FUN=fun)
sum$group <- as.character(c(1,1))

#This fig in MS, Fig 2C#
ggplot(sum, aes(x=Treatment, y=Count.m, group=group, color=Treatment)) + 
  geom_line(color="black") + 
  geom_errorbar(aes(ymin=Count.m-Count.se, ymax=Count.m+Count.se), 
                width=.1)+
  ylab(expression(paste("No. ", F[1]," Shoots"))) + theme_classic() +
  scale_color_manual(values=c("skyblue2","orangered2"))+
  geom_point(size=1) +
  theme(text = element_text (size=8)) +
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  theme(legend.position = "none")

ggsave("Count_SM1.png", dpi=600, units="in", width=1.75, height = 1.5)

d$Date <- as.Date(d$Date, format="%m/%d/%y")

m <- glmer (Count ~ Treatment + (1|Bin), data=d, family = "poisson" )
summary(m)


m2 <- glm(Count~Treatment, data=d, family="poisson")

AIC(m, m2) # m is fitting better.



#################################################################################
### Dry Mass Data Winter Time Point
################################################################################


d<- read.delim(file="SM_T1_ALL.txt", header=T)

fun<-function(x,...){c(m=mean(x, na.rm=T), se=sd(x, na.rm=T)/sqrt(length(x)))}

###Fo Term Shoot Biomass####
#This figure in MS, Fig 2A###

sum <- summaryBy(Shoot.1 ~ Treatment, data=d, FUN=fun)
sum$group <- as.character(c(1,1))

ggplot(sum, aes(x=Treatment, y=Shoot.1.m, group=group, color=Treatment)) + 
  geom_line(color="black") + geom_point(size=1) +
  geom_errorbar(aes(ymin=Shoot.1.m-Shoot.1.se, ymax=Shoot.1.m+Shoot.1.se), width=.1)+
  ylab(expression(paste(F[0], " Biomass (g)"))) + theme_classic() +
  scale_color_manual(values=c("skyblue2","orangered2"))+
  theme(text = element_text (size=8)) + 
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  theme(legend.position="none")


ggsave("ParentShoot_SM1.png", dpi=600, units="in", width=1.75, height = 1.5)

m <- lmer(Shoot.1 ~ Treatment + (1|Bin), data=d)
summary(m)
Anova(m)

RE <- as.data.frame(ranef(m, postVar=T))
Treatment <- rep(seq(1,2), 10)
Treatment <- Treatment[-c(6,15)]

RE$Treatment <- as.factor(Treatment)

m2 <- lm(Shoot.1 ~ Treatment, data=d)
summary(m2)


ggplot(RE, aes(x=condval, y=grp, color=Treatment)) + 
  geom_point() + geom_vline(xintercept=0)+
  geom_errorbarh(aes(xmax=condval+condsd, xmin=condval-condsd), height=.1) +
  xlab("Random Effect") + ylab("Bin") + 
  scale_color_manual(values = c("skyblue", "coral"),
                     labels = c("Cold", "Hot"))

#Remove Bins 2 and 16

dz <- subset(d, d$Bin != "2")
dz <- subset(dz, dz$Bin != "16")

m3 <- lmer(Shoot.1 ~ Treatment + (1|Bin), data=dz)
summary(m3)
Anova(m3)



###F1 Side Shoots Biomass###
#With All dates
#This Fig in MS, Fig 2B#
sum <- summaryBy(Side ~ Treatment, data=d, FUN=fun)
sum$group <- as.character(c(1,1))

ggplot(sum, aes(x=Treatment, y=Side.m, group=group, color=Treatment)) + 
  geom_line(color="black") + geom_point(size=1) +
  geom_errorbar(aes(ymin=Side.m-Side.se, ymax=Side.m+Side.se), width=.1)+
  ylab(expression(paste(F[1], " Total Biomass (g)"))) + theme_classic() +
  theme(text = element_text (size=8)) + 
  scale_color_manual(values=c("skyblue2","orangered2"))+
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  theme(legend.position = "none")

ggsave("OffspringShoot_SM1.png", dpi=600, units="in", width=1.75, height = 1.5)

m <- lmer(Side ~ Treatment + (1|Bin), data=d)
summary(m)
Anova(m)



######
#Per Shoot Biomass##
#This Fig in MS, Fig 2D#
d$Offspring <- d$Count -1

d$PSB <- d$Side / d$Offspring
d <- d[d$ID != 85,]

sum <- summaryBy(PSB ~ Treatment, data=d, FUN=fun)
sum$group <- as.character(c(1,1))

ggplot(sum, aes(x=Treatment, y=PSB.m, group=group, color=Treatment)) + 
  geom_line(color="black") + geom_point(size=1) +
  geom_errorbar(aes(ymin=PSB.m-PSB.se, ymax=PSB.m+PSB.se), width=.1)+
  ylab(expression(paste(F[1], " Size (g ",Shoot^-1,")"))) + theme_classic() +
  theme(text = element_text (size=8)) + 
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  scale_y_continuous(limits = c(.13,.16),
                     breaks = c(0.13, 0.14, 0.15, 0.16))+
  scale_color_manual(values=c("skyblue2","orangered2"))+
  theme(legend.position = "none")

ggsave("OffspringSize_SM1.png", dpi=600, units="in", width=1.8, height = 1.5)

m <- lmer(PSB ~ Treatment + (1|Bin), data=d)
summary(m)
Anova(m)



####Above Biomass####
#This Fig in MS, S1A##

sum <- summaryBy(Above ~ Treatment, data=d, FUN=fun)
sum$group <- as.character(c(1,1))

ggplot(sum, aes(x=Treatment, y=Above.m, group=group)) +
  geom_line() + geom_point(size=1) +
  geom_errorbar(aes(ymin=Above.m-Above.se, ymax=Above.m+Above.se), width=.1)+
  ylab("Tot. Above (g)") + theme_classic() +
  theme(text = element_text (size=9)) + 
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))

ggsave("Above_SM1.png", dpi=300, units="in", width=2, height = 2)

m <- lmer(Above ~ Treatment + (1|Bin), data=d)
summary(m)
Anova(m)


###Below Biomass####
#This fig in MS, Fig S1B#

sum <- summaryBy(Below ~ Treatment, data=d, FUN=fun)
sum$group <- as.character(c(1,1))

ggplot(sum, aes(x=Treatment, y=Below.m, group=group)) + 
  geom_line() + geom_point(size=1) +
  geom_errorbar(aes(ymin=Below.m-Below.se, ymax=Below.m+Below.se), width=.1)+
  ylab("Tot. Below (g)") + theme_classic() +
  theme(text = element_text (size=9)) + 
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))

ggsave("Below_SM1.png", dpi=300, units="in", width=2, height = 2)

m <- lmer(Below ~ Treatment + (1|Bin), data=d)
summary(m)
Anova(m)



###A:B####
#This figure in manuscript, Fig S1C

sum <- summaryBy(AB ~ Treatment, data=d, FUN=fun)
sum$group <- as.character(c(1,1))

ggplot(sum, aes(x=Treatment, y=AB.m, group=group)) + 
  geom_line() + geom_point(size=1) +
  geom_errorbar(aes(ymin=AB.m-AB.se, ymax=AB.m+AB.se), width=.1)+
  ylab("Above:Below") + theme_classic() +
  theme(text = element_text (size=9)) + 
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))

ggsave("AB_SM1.png", dpi=300, units="in", width=2, height = 2)

m <- lmer(AB ~ Treatment + (1|Bin), data=d)
summary(m)
Anova(m)


####Total Starch####
#This figure in MS, Fig S1D

d$TotStarch <- d$Starch_mg_g * d$Rhiz

sum <- summaryBy(TotStarch ~ Treatment, data=d, FUN=fun)
sum$group <- as.character(c(1,1))

ggplot(sum, aes(x=Treatment, y=TotStarch.m, group=group, color=Treatment)) + 
  geom_line(color="black") + geom_point(size=1) +
  geom_errorbar(aes(ymin=TotStarch.m-TotStarch.se, ymax=TotStarch.m+TotStarch.se), width=.1)+
  ylab("Total Starch (mg)") + theme_classic() +
  theme(text = element_text (size=8)) + 
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  theme(legend.position = "none")

ggsave("TotStarch_SM1.png", dpi=300, units="in", width=2, height = 2)

m <- lmer(TotStarch ~ Treatment + (1|Bin), data=d)
summary(m)
Anova(m)


####Starch Per F1 Shoot, Provision####
##This Fig in MS, Fig 2E##
d$Starch_Per <- d$TotStarch/d$Count

d<-d[d$ID != 52,]

sum <- summaryBy(Starch_Per ~ Treatment, data=d, FUN=fun)
sum$group <- as.character(c(1,1))

ggplot(sum, aes(x=Treatment, y=Starch_Per.m, group=group, color=Treatment)) + 
  geom_line(color="black") + geom_point(size=1) +
  geom_errorbar(aes(ymin=Starch_Per.m-Starch_Per.se, ymax=Starch_Per.m+Starch_Per.se), width=.1)+
  ylab(expression(paste(F[1]," Provision (mg Starch)"))) +
  theme_classic() +
  scale_color_manual(values=c("skyblue2", "orangered2"))+
  theme(text = element_text (size=8)) + 
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  scale_y_continuous(breaks=c(5, 5.5, 6))+
  theme(legend.position="none")

ggsave("Starch_Per_SM1.png", dpi=600, units="in", width=1.75, height = 1.5)

m <- lmer(Starch_Per ~ Treatment + (1|Bin), data=d)
summary(m)
Anova(m)

