#StressMem2 Growth and Count at end of second warming treatment Nov 2017
#1/9/15 KD

# added in dry mass data from final break down SM2 Dec 2017

# added in starch data April 2018

#Correct labels on figures for revision, Jan 2020
#Looked at WGP and TGP with model contrasts
#recorrected figure labels again

#Final figure revision for Ecology July 2020 KD
#Removed extraneous code and added notes for GitHub

#data files:
#SM_DryMassFinal.txt
#STARCH_Final.txt
#SM2_Warm.txt (this is growth and count data)

rm(list=ls())
setwd("/Users/Katherine/Documents/Rworking/StressMem")

install.packages("ggplot2")

library(ggplot2)
library(doBy)
library(lme4)
library(dplyr)
library(car)
library(emmeans)

fun <- function(x,...) {c(m=mean(x, na.rm=T), sd=sd(x,na.rm=T), n = length(x), se=sd(x,na.rm=T)/sqrt(length(x)))}

#################################################################################
################################################################################


##### Dry Mass Data: Figures & Models####
rm(list=ls())
d <- read.delim(file="SM_DryMassFinal.txt", header=T)
d$DryMass <- as.numeric(as.character(d$DryMass))


##Look at Histograms for each Tissue Type, then model
SH <- subset(d[d$TissueType=="SHOOT",])
hist(SH$DryMass) #normal, no outliers
shapiro.test(SH$DryMass) #yep


####################################################
#F1, terminal shoot biomass, Fig 3A
##This figure in MS###
SH <- subset(d[d$TissueType=="SHOOT",])

hist(SH$DryMass)

sum <- summaryBy(DryMass ~ SM1Treatment + SM2Treatment, data=SH, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=DryMass.m, group=SM1Treatment, color=SM1Treatment)) +
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=DryMass.m-DryMass.se, ymax =DryMass.m+DryMass.se), width=.1) +
  theme_classic() +
  ylab(expression(paste(F[1], " Biomass (g)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("skyblue2", "orangered2"))+
  theme(text=element_text(size=8))+
  theme(legend.position = "none")+
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))

ggsave(file="TermShoot_SM2.png", dpi=600, units="in", height=1.5, width=1.75)


m_SHO <- lmer (DryMass ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment),
               data=SH, na.action = na.exclude)

summary(m_SHO)
Anova(m_SHO)


m_SHO2 <- lmer (DryMass ~ SM2Treatment * SM1Treatment + (1|SM2_Bin),
                data=SH, na.action = na.exclude)

summary(m_SHO2)
Anova(m_SHO2)

anova(m_SHO, m_SHO2) #model without interaction is better supported

shapiro.test(residuals(m_SHO))

emmF1 <- emmeans(m_SHO, ~SM2Treatment|SM1Treatment)
contrast(emmF1, method="tukey")

emmF1 <- emmeans(m_SHO, ~SM1Treatment|SM2Treatment)
contrast(emmF1, method="tukey")



###F2 Side Shoots Biomass####
#This figure in MS, Fig 4A#
SI <- subset(d[d$TissueType=="SIDE",])
hist(SI$DryMass) #not normal, heavily skewed positive
shapiro.test(SI$DryMass) #nope
qqPlot(SI$DryMass)

SI$Log <- log(SI$DryMass)

sum <- summaryBy(DryMass ~ SM1Treatment + SM2Treatment, data=SI, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=DryMass.m, group=SM1Treatment, color=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=DryMass.m-DryMass.se, ymax =DryMass.m+DryMass.se), width=.1) +
  theme_classic() +
  ylab(expression(paste("Total ",F[2]," Biomass (g)"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("skyblue2","orangered2"))+
  theme(text=element_text(size=8))+
  theme(legend.position = "none")+
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  annotate("text", x=1.5, y=.20, label="*")+
  annotate("text", x=1.5, y=.35, label="*")+
  geom_segment(aes(x=.75, xend=.75, y=.3, yend=.39), size = 1, color="gray80",
               arrow = arrow(length = unit(0.25, "cm")))+
  geom_segment(aes(x=2.25, xend=2.25, y=.17, yend=.24), size = 1, color="gray80",
               arrow = arrow(length = unit(0.25, "cm")))


ggsave(file="SideShoot_SM2.png", dpi=600, units="in", height=1.5, width=1.75)


m_SI <- lmer (Log ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment),
              data=SI, na.action = na.exclude)

summary(m_SI)
Anova(m_SI)

shapiro.test(residuals(m_SI)) #Good


###############################
####Total Above####
#This Figure in MS, Fig 5A##
SH <- subset(d[d$TissueType=="SHOOT", c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment", "DryMass")])
SI <- subset(d[d$TissueType=="SIDE", c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment", "DryMass")])

Above <- left_join(SH, SI, by=c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment"))
Above$ABOV <- rowSums(Above[,c("DryMass.x", "DryMass.y")], na.rm=TRUE)

sum <- summaryBy (ABOV ~ SM1Treatment + SM2Treatment, data=Above, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=ABOV.m, group=SM1Treatment, color=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=ABOV.m-ABOV.se, ymax =ABOV.m+ABOV.se), width=.1) +
  theme_classic() +
  ylab("Tot. Above (g)") +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=8))+
  theme(legend.position = "none")+
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  scale_color_manual(values=c("skyblue2","orangered2"))+
  annotate("text", x=1.5, y=.86, label="*")+
  annotate("text", x=1.5, y=.59, label="*")+
  geom_segment(aes(x=.75, xend=.75, y=.74, yend=.95), size = 1, color="gray80",
               arrow = arrow(length = unit(0.25, "cm")))+
  geom_segment(aes(x=2.25, xend=2.25, y=.52, yend=.67), size = 1, color="gray80",
               arrow = arrow(length = unit(0.25, "cm")))


ggsave(file="Above_SM2.png", dpi=600, units="in", height=1.5, width=1.75)


m_ABOV <- lmer (ABOV ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment),
                data=Above, na.action = na.exclude)

summary(m_ABOV)
Anova(m_ABOV)

shapiro.test(residuals(m_ABOV)) #Good

m_ABOV2 <- lmer (ABOV ~ SM2Treatment * SM1Treatment + (1|SM2_Bin/SM1Treatment),
                 data=Above, na.action = na.exclude)

anova(m_ABOV, m_ABOV2) #no interaction


emmA <- emmeans(m_ABOV, ~SM2Treatment|SM1Treatment)
contrast(emmA, method="tukey")

emmA <- emmeans(m_ABOV, ~SM1Treatment|SM2Treatment)
contrast(emmA, method="tukey")



####Total Below biomass####
#this Figure in MS, Fig 5B#
RH <- subset(d[d$TissueType=="RHIZ", c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment", "DryMass")])
RO <- subset(d[d$TissueType=="ROOT", c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment", "DryMass")])

Below <- left_join(RH, RO, by=c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment"))
Below$BELO <- rowSums(Below[,c("DryMass.x", "DryMass.y")], na.rm=TRUE)


hist(Below$BELO)
qqPlot(Below$BELO) 
shapiro.test(Below$BELO) #super not normal


Below$LogB <- log(Below$BELO) # transform because not normal, and model wont run
hist(Below$LogB)
qqPlot(Below$LogB) 
shapiro.test(Below$LogB) #normal...
Below <- Below[-7,]

sum <- summaryBy (BELO ~ SM1Treatment + SM2Treatment, data=Below, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=BELO.m, group=SM1Treatment, color=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=BELO.m-BELO.se, ymax =BELO.m+BELO.se), width=.1) +
  theme_classic() +
  ylab("Tot. Below (g)") +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("skyblue2","orangered2"))+
  theme(text=element_text(size=8))+
  theme(legend.position = "none")+
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  annotate("text", x=1.5, y=.59, label="*")+
  annotate("text", x=1.5, y=.51, label="*")

ggsave(file="Below_SM2.png", dpi=600, units="in", height=1.5, width=1.75)


m_Be <- lmer (LogB ~ SM2Treatment + SM1Treatment + (1|SM2_Bin),
              data=Below, na.action = na.exclude)

summary(m_Be)
Anova(m_Be)

shapiro.test(residuals(m_Be)) #Just barely... but good

m_Be2 <- lmer (LogB ~ SM2Treatment * SM1Treatment + (1|SM2_Bin),
               data=Below, na.action = na.exclude)

anova(m_Be, m_Be2) #no interaction supported

emmB <- emmeans(m_Be, ~SM2Treatment|SM1Treatment)
contrast(emmB, method="tukey")

emmB <- emmeans(m_Be, ~SM1Treatment|SM2Treatment)
contrast(emmB, method="tukey")


##########################
#A:B
#This figure in MS, Fig 5C##
AB <- left_join (Above, Below, by = c("TreatmentID", "PairID", "SM2_Bin", "SM1_Bin", "SM1Treatment", "SM2Treatment"))
AB$AB <- AB$ABOV/AB$BELO

hist(AB$AB)
qqPlot(AB$AB)
shapiro.test(AB$AB) #close

AB <- AB[-64,]

sum <- summaryBy (AB ~ SM1Treatment + SM2Treatment, data=AB, FUN=fun)

ggplot(sum, aes(x=SM2Treatment, y=AB.m, group=SM1Treatment, color=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=AB.m-AB.se, ymax =AB.m+AB.se), width=.1) +
  theme_classic() +
  ylab("Above:Below") +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("skyblue2","orangered2"))+
  theme(text=element_text(size=8))+
  theme(legend.position = "none")+
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  geom_segment(aes(x=.75, xend=.75, y=1.32, yend=1.57), size = 1, color="gray80",
               arrow = arrow(length = unit(0.25, "cm")))+
  geom_segment(aes(x=2.25, xend=2.25, y=1.26, yend=1.52), size = 1, color="gray80",
               arrow = arrow(length = unit(0.25, "cm")))

ggsave(file="AB_SM2.png", dpi=600, units="in", height=1.5, width=1.75)


#AB Model
m_AB <- lmer (AB ~ SM2Treatment + SM1Treatment + (1|SM2_Bin/SM1Treatment),
              data=AB, na.action = na.exclude)


summary(m_AB)
Anova(m_AB)

shapiro.test(residuals(m_AB)) #Good

emmB <- emmeans(m_AB, ~SM2Treatment|SM1Treatment)
contrast(emmB, method="tukey")

emmB <- emmeans(m_AB, ~SM1Treatment|SM2Treatment)
contrast(emmB, method="tukey")


#############################################################
#################################################################
###Final Starch####

rm(list=ls())
fun <- function(x,...) {c(m=mean(x, na.rm=T), sd=sd(x,na.rm=T), n = length(x), se=sd(x,na.rm=T)/sqrt(length(x)))}

d<- read.delim(file="STARCH_Final.txt", header=T)

Bio <- read.delim(file="SM_DryMassFinal.txt", header=T)
Bio$DryMass <- as.numeric(as.character(Bio$DryMass))
Bio$SAMPLEID <- lapply (1:489, function(x) {paste(Bio$SM1_Bin[x],"-",Bio$Shoot_ID[x],Bio$SM2Treatment[x], sep="")}) # add sample ID colum
Bio$SAMPLEID <- as.character(Bio$SAMPLEID)

d$SAMPLEID <- as.character(d$SAMPLEID)

Rhiz <- subset(Bio, Bio$TissueType=="RHIZ")

d2 <- left_join(d, Rhiz, by="SAMPLEID")


###Total Starch###
#Figure in MS, Fig 5D#
d2$TotStarch <- d2$STARCH_mgg * d2$DryMass
d2 <- subset(d2, d2$RUN_DATE != "NOV19_ABC") #remove unbalanced run, human error
d2 <- subset(d2, d2$RUN_DATE != "13AUG19_AB") #remove unbalanced run, human error
d2 <- subset(d2, d2$RUN_DATE != "FEB") #remove unbalanced run, human error

sum <- summaryBy(TotStarch ~ T1 + T2, data=d2, FUN=fun)

ggplot(sum, aes(x=T2, y=TotStarch.m, group=T1, color=T1))+
  geom_point(size=1) + geom_line(aes(linetype=T1)) +
  geom_errorbar(aes(ymin=TotStarch.m-TotStarch.se, 
                    ymax=TotStarch.m+TotStarch.se), width=.1)+
  theme_classic()+
  ylim(5,14)+
  scale_x_discrete(NULL, labels=c("Control","Warmed"))+
  scale_color_manual(values=c("skyblue2", "orangered2"))+
  ylab("Total Starch (mg)")+
  theme(text=element_text(size=8))+
  theme(legend.position = "none") 

ggsave("TotStarch_SM2.png", dpi=600, width=1.75, height=1.5, units="in")


m_S <- lmer (TotStarch ~ T2 + T1 + (1|SM2_Bin/T1),
             data=d2, na.action = na.exclude)


summary(m_S)
Anova(m_S) #no effect of treatments



##Adding in Leaf growth rate data
GR <- read.delim(file="SM2_Warm.txt", header=T, stringsAsFactors = F)
GR <- subset(GR, GR$TimePoint=="Post") #only final time point

GR$SAMPLEID <- lapply (1:135, function(x) {paste(GR$SM1_Bin[x],"-",GR$Shoot_ID[x],GR$SM2Treatment[x], sep="")}) # add sample ID colum
GR$SAMPLEID <- as.character(GR$SAMPLEID)
GR$Growth_Day <- as.numeric(GR$Growth_Day)


d2 <- left_join(d2, GR, by="SAMPLEID")

d2 <- subset(d2, d2$Growth_Day > 0) #remove dead shoots


###Starch Per Shoot#### F2 provisioning###
#This fig in MS, Fig. 4D#

d2$Prov <-  d2$TotStarch/d2$Shoot_Count

sum <- summaryBy(Prov ~ T1 + T2, data=d2, FUN=fun)

ggplot(sum, aes(x=T2, y=Prov.m, group=T1, color=T1))+
  geom_point(size=1) + geom_line(aes(linetype=T1)) +
  geom_errorbar(aes(ymin=Prov.m-Prov.se, 
                    ymax=Prov.m+Prov.se), width=.1)+
  theme_classic()+
  ylim(0,4)+
  scale_x_discrete(NULL, labels=c("Control","Warmed"))+
  scale_color_manual(values=c("skyblue2", "orangered2"))+
  ylab(expression(paste(F[2]," Provision (mg Starch)")))+
  theme(text=element_text(size=8))+
  theme(legend.position = "none") 

ggsave("F2Prov_SM2.png", dpi=600, width=1.75, height=1.5, units="in")

m_S <- lmer (Prov ~ T2 + T1 + (1|SM2_Bin.x/T1),
             data=d2, na.action = na.exclude)

summary(m_S)
Anova(m_S) #no effect of treatments


########(Length)######################
#This Fig in MS, Figure 3B#
sum <- summaryBy(Length_End ~ T1 + T2, data=d2, FUN=fun, na.rm=T)

ggplot(sum, aes(x=T2, y=Length_End.m, group=T1,color=T1))+
  geom_point(size=1) + geom_line(aes(linetype=T1)) +
  geom_errorbar(aes(ymin=Length_End.m-Length_End.se, 
                    ymax=Length_End.m+Length_End.se), width=.1)+
  theme_classic()+
  ylim(55,95)+
  scale_x_discrete(NULL, labels=c("Control","Warmed"))+
  scale_color_manual(values=c("skyblue2","orangered2"))+
  ylab(expression(paste(F[1]," Length (cm)")))+
  theme(text=element_text(size=8))+
  theme(legend.position = "none")+
  annotate("text", x=1.5, y=71, label="*")+
  annotate("text", x=1.5, y=82, label="*")

ggsave("Length_RN.png", dpi=600, units="in", height=1.5, width=1.75)


m_L <- lmer (Length_End ~ T1 + T2 + (1|SM2_Bin.x/T1) + (1|SM1_Bin.x),
             data=d2, na.action = na.exclude)

summary(m_L)
Anova(m_L)

shapiro.test(residuals(m_L)) #good barely


#########################################################
### Relative Growth Rate####
d2$RGR <- (d2$Growth_Day / d2$Length_End)*100

d2 <- subset(d2,d2$RGR <= 10) #huge outliers

#RGR Figure for ms, Fig 3C####
sum <- summaryBy(RGR ~ T1 + T2, data=d2, FUN=fun, na.rm=T)

ggplot(sum, aes(x=T2, y=RGR.m, group=T1, color=T1))+
  geom_point(size=1) + geom_line(aes(linetype=T1)) +
  geom_errorbar(aes(ymin=RGR.m-RGR.se, 
                    ymax=RGR.m+RGR.se), width=.1)+
  theme_classic()+
  scale_x_discrete(NULL, labels=c("Control","Warmed"))+
  scale_color_manual(values=c("skyblue2","orangered2"))+
  ylab(expression(paste(F[1]," RGR (% ", d^-1,")")))+
  theme(text=element_text(size=8))+
  theme(legend.position = "none")+
  annotate("text", x=1.5, y=3.47, label="*")+
  geom_segment(aes(x=2.25, xend=2.25, y=2.85, yend=2.05), size = 1, color="gray80",
               arrow = arrow(length = unit(0.25, "cm")))

ggsave("RGR_RN.png", dpi=600, units="in", height=1.5, width=1.75)


m_RGR <- lmer (RGR ~ T1 * T2 + (1|SM2_Bin.x),
             data=d2, na.action = na.exclude)

summary(m_RGR)
Anova(m_RGR)

shapiro.test(residuals(m_RGR)) #ok
hist(residuals(m_RGR))



#Growth Data mixed model, 
#Count and per shoot biomass analysis below

d<-read.delim(file="SM2_Warm.txt", header=T, stringsAsFactors = F) 

d$Growth_Day <- as.numeric(as.character(d$Growth_Day))
d$TreatmentID <- as.factor(d$TreatmentID)
d$TimePoint <- as.factor(d$TimePoint)
d$PairID <- as.factor(d$PairID)

d<-subset(d, d$CONTRA!="1") #Contra = from notes during experiment when shoot was damaged but not killed by human error (i.e. Rhizome crushed)


# fixing class
d$TreatmentID<-as.factor(d$TreatmentID)
d$SM2Treatment<-as.factor(d$SM2Treatment)
d$SM1Treatment<-as.factor(d$SM1Treatment)
d$SM2_Bin<-as.factor(d$SM2_Bin)
d$SM1_Bin<-as.factor(d$SM1_Bin)

d$GrowthArea <- d$Growth_Day*d$Width_End
d$ShootArea <- d$Length_End*d$Width_End*4
d$GrowthPercent <-round((d$GrowthArea/d$ShootArea)*100,2)

#Subset out different timepoints for final breakdown
#only data from T2 included in manuscript
T2<-subset(d, d$TimePoint=="Post")

#Huge outliers
T2 <- T2[!rownames(T2) %in% 167, ]
T2 <- T2[!rownames(T2) %in% 171, ]

####RGR Post####
T2 <- subset(T2, T2$GrowthPercent < 3)
T2 <- subset(T2, T2$GrowthPercent > 0)

m_RGR <- lmer (GrowthPercent ~ SM1Treatment * SM2Treatment +
                 (1|SM2_Bin) + (1|SM1_Bin),
               data=T2, na.action = na.exclude) 

summary(m_RGR)
Anova(m_RGR)

shapiro.test(residuals(m_RGR))
hist(residuals(m_RGR))
qqPlot(T2$GrowthPercent) 


##########Shoot Count Post###########
##Figure in MS, Figure 4B##

qqPlot(T2$Shoot_Count)

sum <- summaryBy(Shoot_Count ~ SM1Treatment + SM2Treatment, data=T2, FUN=fun) 

sum$Shoot_Count.m <- sum$Shoot_Count.m - 1

ggplot(sum, aes(x=SM2Treatment, y=Shoot_Count.m, group=SM1Treatment, color=SM1Treatment)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment)) +
  geom_errorbar(aes(ymin=Shoot_Count.m-Shoot_Count.se, ymax =Shoot_Count.m+Shoot_Count.se), width=.1) +
  theme_classic() +
  ylab(expression(paste("No. ",F[2], " Shoots"))) +
  scale_linetype_manual(values=c(1,2)) +
  theme(text=element_text(size=8))+
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  scale_color_manual(values=c("skyblue2", "orangered2"))+
  annotate("text", x=1.5, y=2.5, label="*")+
  theme(legend.position = "none")+
  
  ggsave(file="Count_SM2.png", dpi=600, units="in", height=1.5, width=1.75)


#had to simplify random effects because of singular fit
#picked SM1_Bin because it had a larger effect on model & improved AIC score

m_SH <- glmer (Shoot_Count ~ SM2Treatment + SM1Treatment + (1|SM2_Bin),
               data=T2, na.action = na.exclude, family="poisson")

summary(m_SH)

m_SH2 <- glmer (Shoot_Count ~ SM2Treatment * SM1Treatment + (1|SM2_Bin),
                data=T2, na.action = na.exclude, family="poisson")

summary(m_SH2)


anova(m_SH, m_SH2) #no interaction supported


emmC<-emmeans(m_SH2, ~SM2Treatment|SM1Treatment)
contrast(emmC, method= "tukey")

emmC<-emmeans(m_SH2, ~SM1Treatment|SM2Treatment)
contrast(emmC, method= "tukey")


################################################
####Biomass per shoot 2nd generation shoots#####
##Figure in MS, Fig 4C##

T2$PairID <- as.character(T2$PairID)
d <- read.delim(file="SM_DryMassFinal.txt", header=T)
d$DryMass <- as.numeric(as.character(d$DryMass))
SI <- subset(d[d$TissueType=="SIDE",])
SI$PairID <- as.character(SI$PairID) #SI from biomass data below

z <- left_join(T2, SI, by = c("PairID", "TreatmentID"))

z$Shoot_Count <- z$Shoot_Count -1

z$PerSI <- z$DryMass/z$Shoot_Count


sum <- summaryBy(PerSI ~ SM1Treatment.x + SM2Treatment.x, data=z, FUN=fun) 

ggplot(sum, aes(x=SM2Treatment.x, y=PerSI.m, group=SM1Treatment.x, color=SM1Treatment.x)) + 
  geom_point(size=1)+ geom_line(aes(linetype=SM1Treatment.x)) +
  geom_errorbar(aes(ymin=PerSI.m-PerSI.se, ymax =PerSI.m+PerSI.se), width=.1) +
  theme_classic() +
  ylab(expression(paste(F[2], " Size (g ", Shoot^-1, ")"))) +
  scale_linetype_manual(values=c(1,2)) +
  scale_color_manual(values=c("skyblue2", "orangered2"))+
  theme(text=element_text(size=8))+
  theme(legend.position = "none")+
  scale_x_discrete(NULL, labels = c("Control", "Warmed"))+
  annotate("text", x=1.5, y=0.089, label="*")+
  geom_segment(aes(x=0.75, xend=0.75, y=0.079, yend=.1), size = 1, color="gray80",
               arrow = arrow(length = unit(0.25, "cm")))

ggsave(file="PerShoot_SM2.png", dpi=600, units="in", height=1.5, width=1.75)


qqPlot(residuals(m_PS))
shapiro.test(residuals(m_PS))
z <- z[!rownames(z) %in% 80, ]
z <- z[!rownames(z) %in% 79, ]


m_PS2 <- lmer (PerSI ~ SM2Treatment.x + SM1Treatment.x + (1|SM2_Bin.x),
               data=z, na.action = na.exclude)

Anova(m_PS2)


emmPS<-emmeans(m_PS2, ~SM2Treatment.x|SM1Treatment.x)
contrast(emmPS, method="tukey")

emmPS<-emmeans(m_PS2, ~SM1Treatment.x|SM2Treatment.x)
contrast(emmPS, method="tukey")

