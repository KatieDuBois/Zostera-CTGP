# R Code for manuscript "Previous exposure mediates the response of 
# eelgrass to future warming via clonal transgenerational plasticity"
# Authors K DuBois, SL Williams, JJ Stachowicz
# UC Davis, Bodega Marine Laboratory
# 8/19/19

rm(list=ls()) #Clear working environment
setwd("/Users/Katherine/Documents/Rworking/StressMem") #set working directory

#Load packages
library(ggplot2)
library(doBy)
library(lme4)
library(car)


#Read in data file for growth rates: SM2_Warm.txt
d<-read.delim(file="SM2_Warm.txt", header=T, stringsAsFactors = F) 

d$Growth_Day <- as.numeric(as.character(d$Growth_Day)) #Set class 
d$TreatmentID <- as.factor(d$TreatmentID) #Set class
d$TimePoint <- as.factor(d$TimePoint) #Set class

d<-subset(d, d$CONTRA!="1") #Remove samples with contra-indications from break down (e.g. rhizome broken)


fun <- function(x,...) {c(m=mean(x, na.rm=T), sd=sd(x,na.rm=T), n = length(x), se=sd(x,na.rm=T)/sqrt(length(x)))}


d$PairID <- as.factor(d$PairID) #set class
d$GrowthArea <- d$Growth_Day*d$Width_End #calculate growth area
d$ShootArea <- d$Length_End*d$Width_End*4
d$GrowthPercent <-round((d$GrowthArea/d$ShootArea)*100,2)


# fixing class again
d$TreatmentID<-as.factor(d$TreatmentID)
d$SM2Treatment<-as.factor(d$SM2Treatment)
d$SM1Treatment<-as.factor(d$SM1Treatment)
d$SM2_Bin<-as.factor(d$SM2_Bin)
d$SM1_Bin<-as.factor(d$SM1_Bin)


T2<-subset(d, d$TimePoint=="Post") # only analyzed data Post (1 month after) warming


##Identified as potiental outliers
##Row 167, and 226

T2 <- T2[!rownames(T2) %in% 167, ]
T2 <- T2[!rownames(T2) %in% 171, ]


###########################################
####RESPONSE VARIABLE: LEAF GROWTH RATE####

#Calculate treatment mean and variance for plot
sum <- summaryBy(Growth_Day ~ SM1Treatment + SM2Treatment, data=T2, FUN=fun) 

#Plot
ggplot(sum, aes(x=SM2Treatment, y=Growth_Day.m, group=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=Growth_Day.m-Growth_Day.se, ymax =Growth_Day.m+Growth_Day.se), width=.1) +
  theme_classic() +
  ylab(expression(paste("Growth Rate (cm ", day^-1, ")"))) +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=9))+
  theme(legend.position = "none")+
  scale_x_discrete("", labels = c("Control", "Warmed"))

ggsave(file="Growth_SM2.png", dpi=300, units="in", height=2, width=2)


#Mixed Effects Model, Fixed Effects of Warming (SM2Treatment) and Parent Exposure (SM1Treatment)
# Random effect of Parent Exposure nested within current Bin
m_G <- lmer (Growth_Day ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment),
             data=T2, na.action = na.exclude)

summary(m_G)
Anova(m_G)

shapiro.test(residuals(m_G)) # normal

#Test for interaction
m_G2 <- lmer (Growth_Day ~ SM2Treatment * SM1Treatment + (1|SM2_Bin/SM1Treatment),
              data=T2, na.action = na.exclude)


summary(m_G2)
Anova(m_G2)

shapiro.test(residuals(m_G2)) #normal

AIC (m_G, m_G2) # interaction is slightly better

####################################################
#### RESPONSE VARIABLE: Relative Growth Rate####

d$RGR <- (d$Growth_Day / d$Length_End)*100 # calculate RGR

d <- subset(d,d$RGR <= 10)

d <- d[!rownames(d) %in% 54, ] #qqplot identified one outlier


sum <- summaryBy(RGR ~ T1 + T2, data=d2, FUN=fun, na.rm=T)

ggplot(sum, aes(x=T2, y=RGR.m, group=T1))+
  geom_point() + geom_line(aes(linetype=T1)) +
  geom_errorbar(aes(ymin=RGR.m-RGR.se, 
                    ymax=RGR.m+RGR.se), width=.2)+
  theme_classic()+
  ylim(1.75,4.2)+
  scale_x_discrete("", labels=c("Control","Warmed"))+
  ylab(expression(paste("RGR (% ", d^-1,")")))+
  theme(text=element_text(size=9))+
  theme(legend.position = "none")

ggsave("RGR_RN.png", dpi=300, units="in", height=2, width=2)


m_RGR <- lmer (RGR ~ T1 * T2 + (1|SM2_Bin.x/T1),
               data=d2, na.action = na.exclude)

summary(m_RGR)
Anova(m_RGR)

shapiro.test(residuals(m_RGR)) #good
qqPlot(residuals(m_RGR))



################################################
######RESPONSE VARIABLE: SHOOT COUNT############

sum <- summaryBy(Shoot_Count ~ SM1Treatment + SM2Treatment, data=T2, FUN=fun) 

sum$Shoot_Count.m <- sum$Shoot_Count.m - 1

ggplot(sum, aes(x=SM2Treatment, y=Shoot_Count.m, group=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=Shoot_Count.m-Shoot_Count.se, ymax =Shoot_Count.m+Shoot_Count.se), width=.1) +
  theme_classic() +
  ylab("3rd Gen. Shoot Count") +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=9))+
  theme(legend.position = "none")+
  scale_x_discrete("", labels = c("Control", "Warmed"))

ggsave(file="Count_SM2.png", dpi=300, units="in", height=2, width=2)

#needed to simplify random effects because of singular fit
#picked SM1_Bin because it had a larger effect on model & improved AIC score

m_SH <- glmer (Shoot_Count ~ SM2Treatment + SM1Treatment + (1|SM2_Bin),
               data=T2, na.action = na.exclude, family="poisson")

summary(m_SH)

m_SH2 <- glmer (Shoot_Count ~ SM2Treatment * SM1Treatment + (1|SM2_Bin),
                data=T2, na.action = na.exclude, family="poisson")

summary(m_SH2)


AIC(m_SH, m_SH2) #basically the same


##########################################################################
########Estimated Specific Leaf Area)#####

d$SLA <- (d$Length_End*d$Width_End*3) / d$Total #calculate SLA

sum <- summaryBy(SLA ~ T1 + T2, data=d, FUN=fun, na.rm=T)

ggplot(sum, aes(x=T2, y=SLA.m, group=T1))+
  geom_point() + geom_line(aes(linetype=T1)) +
  geom_errorbar(aes(ymin=SLA.m-SLA.se, 
                    ymax=SLA.m+SLA.se), width=.2)+
  theme_classic()+
  ylim(150,300)+
  scale_x_discrete("", labels=c("Control","Warmed"))+
  ylab(expression(paste("SLA"[E], " (", cm^2," ", g^-1, ")")))+
  theme(text=element_text(size=9))+
  theme(legend.position = "none")

ggsave("SLA_RN.png", dpi=300, units="in", height=2, width=2)



m_SLA <- lmer (SLA ~ T1 + T2 + (1|SM2_Bin.x/T1),
               data=d2, na.action = na.exclude)

summary(m_SLA)
Anova(m_SLA)

shapiro.test(residuals(m_SLA)) #good


##########################################################
########RESPONSE VARIABLE: Leaf Area#####

d$LA <- (d$Length_End*d$Width_End*3) # calculate leaf area

sum <- summaryBy(LA ~ T1 + T2, data=d2, FUN=fun, na.rm=T)

ggplot(sum, aes(x=T2, y=LA.m, group=T1))+
  geom_point() + geom_line(aes(linetype=T1)) +
  geom_errorbar(aes(ymin=LA.m-LA.se, 
                    ymax=LA.m+LA.se), width=.2)+
  theme_classic()+
  ylim(115,200)+
  scale_x_discrete("", labels=c("Control","Warmed"))+
  ylab(expression(paste("LA"[E], " (", cm^2,")")))+
  theme(text=element_text(size=9))+
  theme(legend.position = "none")

ggsave("LA_RN.png", dpi=300, units="in", height=2, width=2)


m_LA <- lmer (LA ~ T1 + T2 + (1|SM2_Bin.x),
              data=d2, na.action = na.exclude)

summary(m_LA)
Anova(m_LA)

shapiro.test(residuals(m_LA)) #Good


################################################
########RESPONSE VARIABLE: Length#####

sum <- summaryBy(Length_End ~ T1 + T2, data=d2, FUN=fun, na.rm=T)

ggplot(sum, aes(x=T2, y=Length_End.m, group=T1))+
  geom_point() + geom_line(aes(linetype=T1)) +
  geom_errorbar(aes(ymin=Length_End.m-Length_End.se, 
                    ymax=Length_End.m+Length_End.se), width=.2)+
  theme_classic()+
  ylim(55,95)+
  scale_x_discrete("", labels=c("Control","Warmed"))+
  ylab("Length (cm)")+
  theme(text=element_text(size=9))+
  theme(legend.position = "none")

ggsave("Length_RN.png", dpi=300, units="in", height=2, width=2)


m_L <- lmer (Length_End ~ T1 + T2 + (1|SM2_Bin.x/T1),
             data=d2, na.action = na.exclude)

summary(m_L)
Anova(m_L)

shapiro.test(residuals(m_L)) #good barely




####################################################################################
####################################################################################

#####NEW DATA SET: SM_DryMassFinal.txt#####
####Dry Mass Data & MODELS####

rm(list=ls())
d <- read.delim(file="SM_DryMassFinal.txt", header=T)
d$DryMass <- as.numeric(as.character(d$DryMass))

fun <- function(x,...) {c(m=mean(x, na.rm=T), sd=sd(x,na.rm=T), n = length(x), se=sd(x,na.rm=T)/sqrt(length(x)))}


##Look at Histograms for each Tissue Type, then model
SH <- subset(d[d$TissueType=="SHOOT",])
hist(SH$DryMass) #normal, no outliers
shapiro.test(SH$DryMass) #yep

##############################################################
####RESPONSE VARIABLE: 1st Clonal Generation Shoot Dry Mass####

sum <- summaryBy(DryMass ~ SM1Treatment + SM2Treatment, data=SH, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=DryMass.m, group=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=DryMass.m-DryMass.se, ymax =DryMass.m+DryMass.se), width=.1) +
  theme_classic() +
  ylab("Clonal Offspring \n 2nd Gen. Shoot (g dry)") +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=9))+
  theme(legend.position = "none")+
  scale_x_discrete("", labels = c("Control", "Warmed"))

ggsave(file="TermShoot_SM2.png", dpi=300, units="in", height=2, width=2.125)


m_SHO <- lmer (DryMass ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment)  + (1|SM1_Bin) + (1|PairID),
               data=SH, na.action = na.exclude)

summary(m_SHO)
Anova(m_SHO)


m_SHO2 <- lmer (DryMass ~ SM2Treatment * SM1Treatment + (1|SM2_Bin/SM1Treatment)  + (1|SM1_Bin) + (1|PairID),
                data=SH, na.action = na.exclude)

summary(m_SHO2)
Anova(m_SHO2)

AIC(m_SHO, m_SHO2) #model without interaction is better supported


m_SHO3 <- lmer (DryMass ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment),
                data=SH, na.action = na.exclude)

summary(m_SHO3) # AIC favors this model... 
Anova(m_SHO3)

AIC(m_SHO, m_SHO3)

shapiro.test(residuals(m_SHO3)) #Good


###################################################################
###RESPONSE VARIABLE: 2nd Generation Clonal Offspring Dry Mass ####

SI <- subset(d[d$TissueType=="SIDE",])
hist(SI$DryMass) #not normal, heavily skewed positive
shapiro.test(SI$DryMass) #nope
qqPlot(SI$DryMass)

SI$Log <- log(SI$DryMass)

sum <- summaryBy(DryMass ~ SM1Treatment + SM2Treatment, data=SI, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=DryMass.m, group=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=DryMass.m-DryMass.se, ymax =DryMass.m+DryMass.se), width=.1) +
  theme_classic() +
  ylab("Clonal Offspring \n 3rd Gen. Shoots (g dry)") +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=9))+
  theme(legend.position = "none")+
  scale_x_discrete("", labels = c("Control", "Warmed"))


ggsave(file="SideShoot_SM2.png", dpi=300, units="in", height=2, width=2.125)


m_SI <- lmer (Log ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment),
              data=SI, na.action = na.exclude)

summary(m_SI)
Anova(m_SI)

shapiro.test(residuals(m_SI)) #Good



#########################################
####RESPONSE VARIABLE: Root Dry Mass####

RO <- subset(d[d$TissueType=="ROOT",])
hist(RO$DryMass) #not normal, heavily skewed positive
shapiro.test(RO$DryMass) #not normals

RO$LogRo <- log(RO$DryMass)
qqPlot(RO$LogRo)
shapiro.test(RO$LogRo) # normal.


sum <- summaryBy (DryMass ~ SM1Treatment + SM2Treatment, data=RO, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=DryMass.m, group=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=DryMass.m-DryMass.se, ymax =DryMass.m+DryMass.se), width=.1) +
  theme_classic() +
  ylab("Root (g dry)") +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=9))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x=element_blank())

ggsave(file="Root_SM2.png", dpi=300, units="in", height=1.58, width=2)

m_Ro <- lmer(LogRo ~ SM2Treatment + SM1Treatment + (1|SM2_Bin),
             data=RO, na.action = na.exclude) 

summary(m_Ro)
Anova(m_Ro)

shapiro.test(residuals(m_Ro)) #not normal...


m_Ro2 <- lmer(LogRo ~ SM2Treatment + SM1Treatment + (1|SM1_Bin),
              data=RO, na.action = na.exclude) 


summary(m_Ro2)
Anova(m_Ro2)

shapiro.test(residuals(m_Ro2))  # normal, had to simplify random effects to get model to fit


#############################################
####RESPONSE VARIABLE: Rhizome####

RH <- subset(d[d$TissueType=="RHIZ",])
hist(RH$DryMass) #not normal, heavily skewed positive, one outlier?
qqPlot(RH$DryMass)

RH <- RH[-68,] # outlier

shapiro.test(RH$DryMass) #not normal

RH$LogR <- log(RH$DryMass)

shapiro.test(RH$LogR) #now normal
qqPlot(RH$LogR)
RH <- RH[-121,]


sum <- summaryBy (DryMass ~ SM1Treatment + SM2Treatment, data=RH, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=DryMass.m, group=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=DryMass.m-DryMass.se, ymax =DryMass.m+DryMass.se), width=.1) +
  theme_classic() +
  ylab("Rhizome (g dry)") +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=9))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x=element_blank())

ggsave(file="Rhiz_SM2.png", dpi=300, units="in", height=1.58, width=2)


m_Rh <- lmer (LogR ~ SM2Treatment + SM1Treatment + (1|SM1_Bin),
              data=RH, na.action = na.exclude)


summary(m_Rh) # needed to remove random effects or singular fit.
Anova(m_Rh)
shapiro.test(residuals(m_Rh)) #not normal


m_Rh1 <- lmer (LogR ~ SM2Treatment + SM1Treatment + (1|SM2_Bin),
               data=RH, na.action = na.exclude)


summary(m_Rh1) # needed to remove random effects or singular fit.
Anova(m_Rh1)
shapiro.test(residuals(m_Rh1)) #normal


###############################################
####RESPONSE VARIABLE: Total Above Dry Mass####

SH <- subset(d[d$TissueType=="SHOOT", c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment", "DryMass")])
SI <- subset(d[d$TissueType=="SIDE", c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment", "DryMass")])

Above <- left_join(SH, SI, by=c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment"))
Above$ABOV <- rowSums(Above[,c("DryMass.x", "DryMass.y")], na.rm=TRUE)

sum <- summaryBy (ABOV ~ SM1Treatment + SM2Treatment, data=Above, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=ABOV.m, group=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=ABOV.m-ABOV.se, ymax =ABOV.m+ABOV.se), width=.1) +
  theme_classic() +
  ylab("Tot. Above (g dry)") +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=9))+
  theme(legend.position = "none")+
  scale_x_discrete("", labels = c("Control", "Warmed"))

ggsave(file="Above_SM2.png", dpi=300, units="in", height=2, width=2)


m_ABOV <- lmer (ABOV ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment),
                data=Above, na.action = na.exclude)

summary(m_ABOV)
Anova(m_ABOV)

shapiro.test(residuals(m_ABOV)) #Good


#############################################
####RESPONSE VARIABLE: Total Below Dry Mass####

RH <- subset(d[d$TissueType=="RHIZ", c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment", "DryMass")])
RO <- subset(d[d$TissueType=="ROOT", c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment", "DryMass")])

Below <- left_join(RH, RO, by=c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment"))
Below$BELO <- rowSums(Below[,c("DryMass.x", "DryMass.y")], na.rm=TRUE)


hist(Below$BELO)
qqPlot(Below$BELO) 
shapiro.test(Below$BELO) # not normal


Below$LogB <- log(Below$BELO) # transform because not normal, and model wont run
hist(Below$LogB)
qqPlot(Below$LogB) 
shapiro.test(Below$LogB) #normal...
Below <- Below[-7,] # extreme outlier removed

sum <- summaryBy (BELO ~ SM1Treatment + SM2Treatment, data=Below, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=BELO.m, group=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=BELO.m-BELO.se, ymax =BELO.m+BELO.se), width=.1) +
  theme_classic() +
  ylab("Tot. Below (g dry)") +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=9))+
  theme(legend.position = "none")+
  scale_x_discrete("", labels = c("Control", "Warmed"))

ggsave(file="Below_SM2.png", dpi=300, units="in", height=2, width=2)


m_Be <- lmer (LogB ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment),
              data=Below, na.action = na.exclude)


summary(m_Be)
Anova(m_Be)

shapiro.test(residuals(m_Be)) #Just barely... but passes


#########################################################
#RESPONSE VARIABLE: A:B Biomass Ratio################

AB <- left_join (Above, Below, by = c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment"))
AB$AB <- AB$ABOV/AB$BELO

hist(AB$AB)
qqPlot(AB$AB)
shapiro.test(AB$AB) #Almost normal

AB <- AB[-64,] #outlier

sum <- summaryBy (AB ~ SM1Treatment + SM2Treatment, data=AB, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=AB.m, group=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=AB.m-AB.se, ymax =AB.m+AB.se), width=.1) +
  theme_classic() +
  ylab("Above:Below") +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=9))+
  theme(legend.position = "none")+
  scale_x_discrete("", labels = c("Control", "Warmed"))

ggsave(file="AB_SM2.png", dpi=300, units="in", height=2, width=2)


m_AB <- lmer (AB ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment),
              data=AB, na.action = na.exclude)


summary(m_AB)
Anova(m_AB)

shapiro.test(residuals(m_AB)) #Good

##########################################################################
####RESPONSE VARIABLE: Biomass per shoot 2nd clonal generation shoots#####
T2$SM2_Bin <- as.factor(T2$SM2_Bin)
SI$SM2_Bin <- as.factor(SI$SM2_Bin)

z <- left_join(T2, SI, by = c("Shoot_ID", "SM2_Bin"))

z$Shoot_Count <- z$Shoot_Count -1

z$PerSI <- z$DryMass/z$Shoot_Count

z1<-do.call(data.frame,lapply(z, function(x) replace(x, is.infinite(x),NA)))

sum <- summaryBy(PerSI ~ SM1Treatment.x + SM2Treatment.x, data=z1, FUN=fun) 

ggplot(sum, aes(x=SM2Treatment.x, y=PerSI.m, group=SM1Treatment.x)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment.x)) +
  geom_errorbar(aes(ymin=PerSI.m-PerSI.se, ymax =PerSI.m+PerSI.se), width=.1) +
  theme_classic() +
  ylab(expression(paste("Biomass/Shoot (g ", Shoot^-1, ")"))) +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=9))+
  theme(legend.position = "none")+
  scale_x_discrete("", labels = c("Control", "Warmed"))

ggsave(file="PerShoot_SM2.png", dpi=300, units="in", height=2, width=2)


m_PS <- lmer (PerSI ~ SM2Treatment.x + SM1Treatment.x + (1|SM2_Bin/SM1Treatment.x),
              data=z1, na.action = na.exclude)


summary(m_PS)
Anova(m_PS)


###################################################################################
#################################################################
#### NEW DATA SET: Rhizome Starch####

#Read in file
d<- read.delim(file="STARCH_Final.txt", header=T)

###Decided to use replicates from APRIL RUNS only because these runs were balanced
## Meaning they had equal numbers from each treatment in a run
## Run effect was large and balanced runs were essential

d <- subset(d,d$REP!= 1) #leaves only data from April runs

ggplot(d, aes(x=RUN_DATE, y=STARCH_mgg, color=T1))+geom_point()

me<-lmer(log(STARCH_mgg) ~ T1 + T2  + (1|RUN_DATE), data=d) # log transformed for normality

summary(me)
Anova(me)

shapiro.test(residuals(me)) # normal
qqPlot(residuals(me)) # normal no outliers

sum <- summaryBy(STARCH_mgg ~ T1 + T2, data=d, FUN=fun) #calculate means

ggplot(sum, aes(x=T2, y=STARCH_mgg.m, group=T1))+
  geom_point() + geom_line(aes(linetype=T1)) +
  geom_errorbar(aes(ymin=STARCH_mgg.m-STARCH_mgg.se, 
                    ymax=STARCH_mgg.m+STARCH_mgg.se), width=.2)+
  theme_classic()+
  ylim(20,40)+
  scale_x_discrete("", labels=c("Control","Warmed"))+
  ylab(expression(paste("Starch (mg ", g^-1, ")")))+
  theme(text=element_text(size=9))+
  theme(legend.position = "none") 

ggsave("StarchSM2RN.png", dpi=300, units="in", height=2, width=2)

