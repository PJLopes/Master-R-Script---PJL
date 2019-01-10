#-----------------------------MASTER GENERAL SCRIPT----------------------------------------

setwd("C:/Users/patri/Desktop/UniL/Master/Data/Analyses R/")

library(car);library(lme4);library(nlme);library(knitr);library(MASS);library(survival);library(ranger)
library(ggplot2);library(dplyr);library(ggfortify);library(ggplot2);library(arm);library(lme4);library(car)
library(LMERConvenienceFunctions);library(languageR);library(nlme);library(plyr);library(ggplot2)
library(gplots);library(lattice);library(Rcmdr);library(scatterplot3d);library(sciplot);library(OutlierDC)
library(extremevalues);library(plot3D);library(multcomp);library(outliers);library(survival);library(coxme)
library(rms);library(lmerTest);library(survival);library(survminer);library(readr);library(dplyr)
library(foreign);library(nnet);library(ggplot2);library(reshape2);library(caret);library(survBootOutliers)
library(multcomp);library(lsmeans);library(tidyr);library(data.table);library(plyr);library(plotrix);library(survminer)
library(fields)

#------------------SURVIVAL ANALYSIS-----------------------------------------------

AllDataSurvival <- read.delim("All Data - Survival Analysis - Final File.csv", header = TRUE, sep = ";")
AllData <- read.delim("All Data - Final.csv", header = TRUE, sep = ";")
AllDataSeparated <- read.delim("All Data - Final - Separated.csv", header = TRUE, sep = ";")
AllDataSurvivalWorkers <- read.delim("Observations - Workers Survival Analysis.csv", header = TRUE, sep = ";")

summary(AllDataSeparated$hib_treatment[AllDataSeparated$treatment=="MM"])
summary(AllDataSeparated$hib_treatment[AllDataSeparated$treatment.split=="PM.M"])
summary(AllDataSeparated$hib_treatment[AllDataSeparated$treatment.split=="PM.P"])
summary(AllDataSeparated$hib_treatment[AllDataSeparated$treatment=="PP"])

AllDataSeparated$replicate_number<-as.factor(AllDataSeparated$replicate_number)
AllDataSurvival$replicate_number<-as.factor(AllDataSurvival$replicate_number)


#Taking out replicates with colony sisters 
AllDataSurvival_NoSis <- filter(AllDataSurvival, replicate_number != '22')
AllDataSurvival_NoSis <- filter(AllDataSurvival_NoSis, replicate_number != '40')

AllDataSeparated_NoSis <- filter(AllDataSeparated, replicate_number != '22')
AllDataSeparated_NoSis <- filter(AllDataSeparated_NoSis, replicate_number != '40')

AllData_NoSis <- filter(AllData, replicate_number != '22')
AllData_NoSis <- filter(AllData_NoSis, replicate_number != '40')

AllDataSurvivalWorkers_NoSis <- filter(AllDataSurvivalWorkers, replicate_number != '22')
AllDataSurvivalWorkers_NoSis <- filter(AllDataSurvivalWorkers_NoSis, replicate_number != '40')

###QUEEN SURVIVAL###
#Plots
km <- with(AllDataSurvival_NoSis, Surv(survival_time_queen, queen_status))
km_fit <- survfit(Surv(survival_time_queen, queen_status) ~ 1, data=AllDataSurvival_NoSis)
summary(km_fit, times = c(1,10,20,30,50,70,90,100))
ggsurvplot(km_fit, ylab = "survival rate", conf.int=TRUE)

km_trt_fit <-  survfit(Surv(survival_time_queen, queen_status) ~ treatment.split, data=AllDataSurvival_NoSis)
autoplot(km_trt_fit, type="single", main="Survival of Queens across treatments", ylab = "Survival rate", xlab = "Time (days)")
?autoplot
ggsurvplot(km_trt_fit, ylab = "Survival rate (Total nb queens alive/Total nb queens)", xlab = "Time (days)",legend.title = "Treatment", legend=c("right"), palette = c("indianred1","palegreen3","palegreen","lightskyblue"), legend.labs=c("MxM","M in PxM","P in PxM","PxP"), conf.int=FALSE)
survgraph1 + scale_fill_discrete(name="Treatment")
?ggsurvplot


km_origin_fit <-  survfit(Surv(survival_time_queen, queen_status) ~ colony_social_origin, data=AllDataSurvival_NoSis)
autoplot(km_origin_fit, main="Survival of Queens according to social origin (in general)", ylab = "survival rate")

AllDataSurvival_NoSisPM <- filter(AllDataSurvival_NoSis, treatment=="PM")
km_trtPM_fit <-  survfit(Surv(survival_time_queen, queen_status) ~ 1, data=AllDataSurvival_NoSisPM)
autoplot(km_trtPM_fit, main="Survival of queens among PM treatments")
km_trtPM_Origin_fit <-  survfit(Surv(survival_time_queen, queen_status) ~ colony_social_origin, data=AllDataSurvival_NoSisPM)
autoplot(km_trtPM_Origin_fit, main="Survival of Queens according to social origins among PM treatments", ylab = "survival rate")
summary(AllDataSurvival_NoSisPM)
AllDataSurvival_NoSisPM.M <- filter(AllDataSurvival_NoSisPM, colony_social_origin=="M")
AllDataSurvival_NoSisPM.P <- filter(AllDataSurvival_NoSisPM, colony_social_origin=="P")
summary(AllDataSurvival_NoSisPM.M)
summary(AllDataSurvival_NoSisPM.P)

AllDataSurvival_NoSisMM <- filter(AllDataSurvival_NoSis, treatment=="MM")
km_trtMM_fit <-  survfit(Surv(survival_time_queen, queen_status) ~ 1, data=AllDataSurvival_NoSisMM)
autoplot(km_trtMM_fit, main="Survival of queens among MM treatments")

AllDataSurvival_NoSisPP <- filter(AllDataSurvival_NoSis, treatment=="PP")
km_trtPP_fit <-  survfit(Surv(survival_time_queen, queen_status) ~ 1, data=AllDataSurvival_NoSisPP)
autoplot(km_trtPP_fit, main="Survival of queens among PP treatments")

#Models - Coxph and GLMER
model_survQ.coxph <- coxph(Surv(survival_time_queen, queen_status) ~ treatment.split + worker_num + frailty(replicate_number), data=AllDataSurvival_NoSis)
summary(model_survQ.coxph)
model_survQ.coxph
anova(model_survQ.coxph)
Anova(model_survQ.coxph, type= 2)
summary(Anova(model_survQ.coxph, type=2))
lsmeans(model_survQ.coxph, list(pairwise ~ treatment.split), adjust = 'fdr')
lsmeans(model_survQ.coxph, list(pairwise ~ worker_num), adjust = 'fdr')

cox.zph(model_survQ.coxph)
ggcoxzph(cox.zph(model_survQ.coxph))
ggcoxdiagnostics(model_survQ.coxph, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(model_survQ.coxph, type = "deviance", linear.predictions = FALSE, ggtheme = theme_bw())

model_survQ_PM.coxph <- coxph(Surv(survival_time_queen, queen_status) ~ colony_social_origin + worker_num + frailty(replicate_number), data=AllDataSurvival_NoSisPM)
lsmeans(model_survQ_PM.coxph, list(pairwise ~ colony_social_origin), adjust = 'fdr')
summary(model_survQ_PM.coxph)
Anova(model_survQ_PM.coxph, type = 2)
cox.zph(model_survQ_PM.coxph)
ggcoxzph(cox.zph(model_survQ_PM.coxph))
ggcoxdiagnostics(model_survQ_PM.coxph, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(model_survQ_PM.coxph, type = "deviance", linear.predictions = FALSE, ggtheme = theme_bw())

ifelse(AllDataSeparated_NoSis$queen_alive == 'yes', 1, 0)->queen_alive2
AllDataSeparated2<-cbind(AllDataSeparated_NoSis, queen_alive2)
model_survQ_01<- glmer(queen_alive2 ~ treatment.split + worker_num_col_2nd_count + (1|replicate_number) + (1|colony_fem_fam), family=binomial,  data=AllDataSeparated2)
summary(model_survQ_01)
Anova(model_survQ_01, type = 2)
lsmeans(model_survQ_01, list(pairwise ~ treatment.split), adjust = 'fdr')

AllDataSeparated2PM <- filter(AllDataSeparated2, treatment=="PM")
model_survQ_01PM<- glmer(queen_alive2 ~ colony_social_origin + worker_num_col_2nd_count + (1|replicate_number) + (1|colony_fem_fam), family=binomial,  data=AllDataSeparated2PM)
summary(model_survQ_01PM)
Anova(model_survQ_01PM, type = 2)
lsmeans(model_survQ_01PM, list(pairwise ~ colony_social_origin), adjust = 'fdr')

AIC(model_survQ.coxph,model_survQ_01)
AIC(model_survQ_PM.coxph, model_survQ_01PM)

###WOKER SURVIVAL###
#Plots
kmW <- with(AllDataSurvivalWorkers_NoSis, Surv(survival_time_worker, worker_status))
kmW_fit <- survfit(Surv(survival_time_worker, worker_status) ~ 1, data=AllDataSurvivalWorkers_NoSis)
summary(kmW_fit, times = c(1,10,20,30,50,70,90,100))
autoplot(kmW_fit,main = "Survival of Workers in general", ylab = "survival rate")

kmW_trt_fit <-  survfit(Surv(survival_time_worker, worker_status) ~ treatment, data=AllDataSurvivalWorkers_NoSis)
autoplot(kmW_trt_fit, main="Survival of Workers across treatments", ylab = "survival rate")
ggsurvplot(kmW_trt_fit, title = "Survival of Workers across treatments", ylab = "survival rate", conf.int=FALSE, legend=c("right"), legend.labs=c("MM","PM","PP"))


kmW_origin_fit <-  survfit(Surv(survival_time_worker, worker_status) ~ colony_social_origin, data=AllDataSurvivalWorkers_NoSis)
autoplot(kmW_origin_fit, main="Survival of Workers according to social origin in general", ylab = "survival rate")

AllDataSurvivalWorkersPM <- filter(AllDataSurvivalWorkers_NoSis, treatment=="PM")
kmW_trtPM_fit <-  survfit(Surv(survival_time_worker, worker_status) ~ 1, data=AllDataSurvivalWorkersPM)
autoplot(kmW_trtPM_fit, main="Survival of workers among PM treatments")
kmW_trtPM_Origin_fit <-  survfit(Surv(survival_time_worker, worker_status) ~ colony_social_origin, data=AllDataSurvivalWorkersPM)
autoplot(kmW_trtPM_Origin_fit, main="Survival of workers according to social origins among PM treatments", ylab = "survival rate")
summary(AllDataSurvivalWorkersPM$colony_social_origin)
summary(AllDataSurvivalWorkersPM$worker_status)
summary(AllDataSurvivalWorkersPP$worker_status)
summary(AllDataSurvivalWorkersMM$worker_status)

AllDataSurvivalWorkersMM <- filter(AllDataSurvivalWorkers_NoSis, treatment=="MM")
kmW_trtMM_fit <-  survfit(Surv(survival_time_worker, worker_status) ~ 1, data=AllDataSurvivalWorkersMM)
autoplot(kmW_trtMM_fit, main="Survival of workers among MM treatments")

AllDataSurvivalWorkersPP <- filter(AllDataSurvivalWorkers_NoSis, treatment=="PP")
kmW_trtPP_fit <-  survfit(Surv(survival_time_worker, worker_status) ~ 1, data=AllDataSurvivalWorkersPP)
autoplot(kmW_trtPP_fit, main="Survival of workers among PP treatments")

#Model
model_survW.coxph <- coxph(Surv(survival_time_worker, worker_status) ~ treatment + frailty(replicate_number), data=AllDataSurvivalWorkers_NoSis)
model_survW.coxph
anova(model_survW.coxph)
lsmeans(model_survW.coxph, list(pairwise ~ treatment), adjust = 'fdr')

cox.zph(model_survW.coxph)
ggcoxzph(cox.zph(model_survW.coxph))
ggcoxdiagnostics(model_survW.coxph, type = "dfbeta",linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(model_survW.coxph, type = "deviance", linear.predictions = FALSE, ggtheme = theme_bw())


model_survW_PM.coxph <- coxph(Surv(survival_time_worker, worker_status) ~ colony_social_origin + frailty(replicate_number), data=AllDataSurvivalWorkersPM)
lsmeans(model_survW_PM.coxph, list(pairwise ~ colony_social_origin), adjust = 'fdr')
Anova(model_survW_PM.coxph, type=2)
anova(model_survW_PM.coxph)
model_survW_PM.coxph
cox.zph(model_survW_PM.coxph)

#------------------------------FUSION ANALYSIS--------------------------------------------------

AllDataSurvival <- read.delim("All Data - Survival Analysis - Final File.csv", header = TRUE, sep = ";")
AllData <- read.delim("All Data - Final.csv", header = TRUE, sep = ";")
AllDataSeparated <- read.delim("All Data - Final - Separated.csv", header = TRUE, sep = ";")
AllDataSeparated$replicate_number<-as.factor(AllDataSeparated$replicate_number)
AllDataSurvival$replicate_number<-as.factor(AllDataSurvival$replicate_number)
AllDataSeparated_NoSis <- filter(AllDataSeparated, replicate_number != '22')
AllDataSeparated_NoSis <- filter(AllDataSeparated_NoSis, replicate_number != '40')
AllData_NoSis <- filter(AllData, replicate_number != '22')
AllData_NoSis <- filter(AllData_NoSis, replicate_number != '40')
AllData_NoSis

#Plots
library(Rmisc)

ifelse(AllData$association == 1, 'yes', 'no')->association2
tgc <- summarySE(AllData, measurevar="association", groupvars="treatment", na.rm = TRUE)
tgc

?ggplot
try <- ggplot(tgc, col = c("indianred1","palegreen2","lightskyblue"), aes(x=treatment, y=association)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=association-ci, ymax=association+ci),
                width=.2,                    
                position=position_dodge(.9)) +
  labs(x="", y="Proportion of colony fusion (± CI)") +
  scale_fill_manual(values=c("#FFFFFF", "#111214", "#111214")) +
  scale_x_discrete(labels = c('MM',
                              'PM',
                              'PP')) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.title.x =element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y =element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
try

#Model
model_assoc2<- glm(association~ treatment, family=binomial,  data=AllData_NoSis)
summary(model_assoc2)
model_assoc2
summary(anova(model_assoc2))
Anova(model_assoc2, type=2)
lsmeans(model_assoc2, list(pairwise ~ treatment), adjust = 'fdr')
summary(AllData_NoSis$treatment[AllData_NoSis$association==1])

#--------------------------------WEIGHT ANALYSIS-----------------------------------------------------------

AllDataSurvival <- read.delim("All Data - Survival Analysis - Final File.csv", header = TRUE, sep = ";")
AllData <- read.delim("All Data - Final.csv", header = TRUE, sep = ";")
AllDataSeparated <- read.delim("All Data - Final - Separated.csv", header = TRUE, sep = ";")
AllDataSeparated$replicate_number<-as.factor(AllDataSeparated$replicate_number)
AllDataSurvival$replicate_number<-as.factor(AllDataSurvival$replicate_number)
AllDataSeparated_NoSis <- filter(AllDataSeparated, replicate_number != '22')
AllDataSeparated_NoSis <- filter(AllDataSeparated_NoSis, replicate_number != '40')
AllData_NoSis <- filter(AllData, replicate_number != '22')
AllData_NoSis <- filter(AllData_NoSis, replicate_number != '40')
AllDataSeparated_no_na<-filter(AllDataSeparated_NoSis, weight_queen_post != 'NA')


AlldataPM_2Alive <- filter(AllData_NoSis, survival_time_queen1 =="108" )
AlldataPM_2Alive <- filter(AlldataPM_2Alive, survival_time_queen2 =="108" )
AlldataPM_2Alive <- filter(AlldataPM_2Alive, treatment =="PM" )
AllDataPM_2Alive_no_na <- filter(AlldataPM_2Alive, weight_queen1_post != 'NA')

AllDataPM_SurvEnd <- read.delim("All Data PM SurvEnd.csv", header = TRUE, sep = ";")

#Plots
set.panel(1,1)
summary(AllDataSeparated_NoSis)
boxplot(AllDataSeparated_NoSis$weight_queen)
AllDataSeparatedMM <- filter(AllDataSeparated_no_na, treatment=="MM")
boxplot(AllDataSeparatedMM$weight_queen, AllDataSeparatedMM$weight_queen_post, col = c("indianred1","indianred1"), names=c("MM Before", "MM After"), main="MM treatment", xlab="Stage of the experiment", ylab = "Weight of the queens (mg)")


AllDataSeparatedPM <- filter(AllDataSeparated_no_na, treatment=="PM")
AllDataSeparatedPM.P <- filter(AllDataSeparatedPM, colony_social_origin=="P")
AllDataSeparatedPM.M <- filter(AllDataSeparatedPM, colony_social_origin=="M")
boxplot(AllDataSeparatedPM.P$weight_queen, AllDataSeparatedPM.P$weight_queen_post, AllDataSeparatedPM.M$weight_queen, AllDataSeparatedPM.M$weight_queen_post, col = c("palegreen","palegreen","palegreen3","palegreen3"), names=c("PM.P Before", "PM.P After", "PM.M Before", "PM.M After"), main="PM treatment", xlab="Stage of the experiment", ylab = "Weight of the queens (mg)")

AllDataSeparatedPP <- filter(AllDataSeparated_no_na, treatment=="PP")
boxplot(AllDataSeparatedPP$weight_queen, AllDataSeparatedPP$weight_queen_post, col = c("lightskyblue","lightskyblue"), names=c("PP Before", "PP After"), main="PP treatment", xlab="Stage of the experiment", ylab = "Weight of the queens (mg)")

boxplot(AllDataSeparatedMM$weight_queen, AllDataSeparatedMM$weight_queen_post,AllDataSeparatedPM.P$weight_queen, AllDataSeparatedPM.P$weight_queen_post, AllDataSeparatedPM.M$weight_queen, AllDataSeparatedPM.M$weight_queen_post, AllDataSeparatedPP$weight_queen, AllDataSeparatedPP$weight_queen_post, col = c("indianred1","indianred1","palegreen","palegreen","palegreen3","palegreen3","lightskyblue","lightskyblue"), names=c("MM Start", "MM End", "PM.P Start", "PM.P End", "PM.M Start", "PM.M End", "PP Start", "PP End"),  xlab="Stage of the experiment", ylab = "Weight of the queens (mg)")
boxplot(AllDataSeparatedMM$weight_queen_diff, AllDataSeparatedPM.M$weight_queen_diff, AllDataSeparatedPM.P$weight_queen_diff, AllDataSeparatedPP$weight_queen_diff, col = c("indianred1","palegreen3","palegreen","lightskyblue"), names=c("MM", "PM.M", "PM.P", "PP"), main="Difference in weight of the Queens according to treatment", xlab="Treatment", ylab = "Difference in weight of the queens")


#Wilcox.tests
wilcox.test(AllDataSeparated_NoSis$weight_queen, AllDataSeparated_NoSis$weight_queen_post, paired=TRUE)
wilcox.test(AllDataSeparated_no_na$weight_queen, AllDataSeparated_no_na$weight_queen_post, paired=TRUE)

wilcox.test(AllDataSeparatedPP$weight_queen, AllDataSeparatedPP$weight_queen_post, paired=TRUE)
wilcox.test(AllDataSeparatedMM$weight_queen, AllDataSeparatedMM$weight_queen_post, paired=TRUE)
wilcox.test(AllDataSeparatedPM$weight_queen, AllDataSeparatedPM$weight_queen_post, paired=TRUE)
wilcox.test(AllDataSeparatedPM.P$weight_queen, AllDataSeparatedPM.P$weight_queen_post, paired=TRUE)
wilcox.test(AllDataSeparatedPM.M$weight_queen, AllDataSeparatedPM.M$weight_queen_post, paired=TRUE)

summary(AllDataSeparatedPM$weight_queen)
std.error(AllDataSeparatedPM$weight_queen)
summary(AllDataSeparatedPM$weight_queen_post)
std.error(AllDataSeparatedPM$weight_queen_post)

#Models
model.weight.split <- lmer(weight_queen_diff ~ treatment.split + association + weight_queen + (1|replicate_number) + (1|colony_fem_fam), data= AllDataSeparated_no_na)
summary(model.weight.split)
par(mfrow=c(1,3))
resids = resid(model.weight.split)
plot(resids)
fits = fitted(model.weight.split)
plot(fits, resids)
qqnorm(resids); qqline(resids, col = 2)
sort(resids)
Anova(model.weight.split, type=2)
lsmeans(model.weight.split, list(pairwise ~ treatment.split), adjust = 'fdr')
summary(AllDataSeparated_no_na$treatment.split)
AllDataSeparated_no_outlier<-AllDataSeparated_no_na[-c(28), ]
AllDataSeparated_no_na[28, ] #Replicate Nb 59

#2
model.weight.split2 <- lmer(weight_queen_diff ~ treatment.split + weight_queen + association + (1|replicate_number) + (1|colony_fem_fam), data= AllDataSeparated_no_outlier)
summary(model.weight.split2)
par(mfrow=c(1,3))
resids = resid(model.weight.split2)
plot(resids)
fits = fitted(model.weight.split2)
plot(fits, resids)
qqnorm(resids); qqline(resids, col = 2)
sort(resids)
Anova(model.weight.split2, type=2)
lsmeans(model.weight.split2, list(pairwise ~ treatment.split), adjust = 'fdr')
summary(AllDataSeparated_no_outlier$treatment.split)


AIC(model.weight.split, model.weight.split2)

#AllDataSeparatedPM TEST
wilcox.test(AllDataPM_SurvEnd$weight_queen1_diff, AllDataPM_SurvEnd$weight_queen2_diff, paired=TRUE)
summary(AllDataSeparatedPM$replicate_number)

#-----------------------------------DIFFERENCE IN WORKFORCE--------------------------------

AllDataSurvival <- read.delim("All Data - Survival Analysis - Final File.csv", header = TRUE, sep = ";")
AllData <- read.delim("All Data - Final.csv", header = TRUE, sep = ";")
AllDataSeparated <- read.delim("All Data - Final - Separated.csv", header = TRUE, sep = ";")
AllDataSeparated$replicate_number<-as.factor(AllDataSeparated$replicate_number)
AllDataSurvival$replicate_number<-as.factor(AllDataSurvival$replicate_number)
AllDataSeparated_NoSis <- filter(AllDataSeparated, replicate_number != '22')
AllDataSeparated_NoSis <- filter(AllDataSeparated_NoSis, replicate_number != '40')
AllData_NoSis <- filter(AllData, replicate_number != '22')
AllData_NoSis <- filter(AllData_NoSis, replicate_number != '40')

#Preparing Data
AllDataSeparated_no_na_W <-filter(AllDataSeparated_NoSis, nb_workers_colony_end != 'NA')
#To have only the queens that are alive at the end
AllDataSeparated_no_na_W_no_assoc <-filter(AllDataSeparated_no_na_W, association != '1')
#To have only the queens that are alive at the end and not associated
#Double the number of workers instead OR take out PP from the contrasts
summary(AllDataSeparated_no_na_W)
summary(AllDataSeparated_no_na_W_no_assoc)

AllDataSeparated_no_na_WDiff <- AllDataSeparated_no_na_W$nb_workers_colony_end_corr - AllDataSeparated_no_na_W$worker_num_col_2nd_count
AllDataSeparated_no_na_W$Diff.W <- AllDataSeparated_no_na_WDiff
AllDataSeparated_no_na_W$Diff.W

AllDataSeparated_no_na_W_no_assocDiff <- AllDataSeparated_no_na_W_no_assoc$nb_workers_colony_end_corr - AllDataSeparated_no_na_W_no_assoc$worker_num_col_2nd_count
AllDataSeparated_no_na_W_no_assoc$Diff.W <- AllDataSeparated_no_na_W_no_assocDiff

#Plots
plot(AllDataSeparated_no_na_W$treatment.split, AllDataSeparated_no_na_W$nb_workers_colony_end, col = c("indianred1","palegreen3","palegreen","lightskyblue"), main="Number of workers by colony at the end of experiment (by treatment)", xlab="Treatment", ylab = "Number of workers")
plot(AllDataSeparated_no_na_W$treatment.split, AllDataSeparated_no_na_W$nb_workers_colony_end,  col = c("lightcoral", "lightgreen", "lightblue"), main="Number of workers by colony at the end of experiment (by treatment, inclu.assoc)", xlab="Treatment", ylab = "Number of workers")

plot(AllDataSeparated_no_na_W$colony_social_origin, AllDataSeparated_no_na_W$nb_workers_colony_end, col = c("grey50","grey90"), main="Number of workers by colony at the end of experiment (by social origin)", xlab="Social Origin", ylab = "Number of workers")

AllDataSeparated_no_na_WPM <- filter(AllDataSeparated_no_na_W, treatment=="PM")
plot(AllDataSeparated_no_na_WPM$colony_social_origin, AllDataSeparated_no_na_WPM$nb_workers_colony_end, col = c("palegreen3","palegreen"), main="Nb of W. by colony at the end of experiment for, PM treatments (by social origin)", xlab="Social origin", ylab = "Number of workers")


par(mfrow=c(1,1))
plot(AllDataSeparated_no_na_W$treatment.split, AllDataSeparated_no_na_W$worker_num_col_2nd_count)
plot(AllDataSeparated_no_na_W$treatment.split, AllDataSeparated_no_na_W$nb_workers_colony_end)
plot(AllDataSeparated_no_na_W$treatment.split, AllDataSeparated_no_na_W$nb_workers_colony_end_corr)

plot(AllDataSeparated_no_na_W$colony_social_origin, AllDataSeparated_no_na_W$Diff.W)
plot(AllDataSeparated_no_na_W_no_assoc$treatment.split, AllDataSeparated_no_na_W_no_assoc$Diff.W, col = c("indianred1","palegreen3","palegreen","lightskyblue"), xlab="Treatment", ylab = "Number of workers")

#Models

model.workers.diff.trt <- lmer(Diff.W ~ treatment.split + weight_queen + cuck_num_col + number_dead_workers_col + (1|replicate_number) + (1|colony_fem_fam), AllDataSeparated_no_na_W)
summary(model.workers.diff.trt)
Anova(model.workers.diff.trt, type = 2)
par(mfrow=c(1,3))
resids = resid(model.workers.diff.trt)
plot(resids)
fits = fitted(model.workers.diff.trt)
plot(fits, resids)
qqnorm(resids); qqline(resids, col = 2)
sort(resids)
lsmeans(model.workers.diff.trt, list(pairwise ~ treatment.split), adjust = 'tukey')


model.workers.diff.trt.bis <- lmer(Diff.W ~ treatment.split + weight_queen + cuck_num_col + number_dead_workers_col + (1|replicate_number) + (1|colony_fem_fam), AllDataSeparated_no_na_W_no_assoc)
summary(model.workers.diff.trt.bis)
Anova(model.workers.diff.trt.bis, type = 2)
par(mfrow=c(1,3))
resids = resid(model.workers.diff.trt.bis)
plot(resids)
fits = fitted(model.workers.diff.trt.bis)
plot(fits, resids)
qqnorm(resids); qqline(resids, col = 2)
sort(resids)
lsmeans(model.workers.diff.trt.bis, list(pairwise ~ treatment.split), adjust = 'tukey')
AllDataSeparated_no_na_W_no_assoc[c(14), ] #replicate_number 37
AllDataSeparated_no_outlier_no_na_W<-AllDataSeparated_no_na_W_no_assoc[-c(14), ] #14 possible outlier

model.workers.diff.trt2 <- lmer(Diff.W ~ treatment.split + weight_queen + cuck_num_col + number_dead_workers_col + (1|replicate_number) + (1|colony_fem_fam), AllDataSeparated_no_outlier_no_na_W)
summary(model.workers.diff.trt2)
Anova(model.workers.diff.trt2, type = 2)
par(mfrow=c(1,3))
resids = resid(model.workers.diff.trt2)
plot(resids)
fits = fitted(model.workers.diff.trt2)
plot(fits, resids)
qqnorm(resids); qqline(resids, col = 2)
sort(resids)
lsmeans(model.workers.diff.trt2, list(pairwise ~ treatment.split), adjust = 'fdr')
summary(AllDataSeparated_no_outlier_no_na_W$treatment.split)

AIC(model.workers.diff.trt.bis,model.workers.diff.trt2)

cor.test(AllDataSeparated_no_outlier_no_na_W$weight_queen, AllDataSeparated_no_outlier_no_na_W$cuck_num_col)
cor.test(AllDataSeparated_no_outlier_no_na_W$weight_queen, AllDataSeparated_no_outlier_no_na_W$number_dead_workers_col)
cor.test(AllDataSeparated_no_outlier_no_na_W$cuck_num_col, AllDataSeparated_no_outlier_no_na_W$number_dead_workers_col)
?cor
?cor.test

# check residuals by variable in the model  #
boxplot(resids ~ AllDataSeparated_no_outlier_no_na_W$treatment.split) 
leveneTest(resids ~ AllDataSeparated_no_outlier_no_na_W$treatment.split)######## this is ok i guess
leveneTest(resids ~ AllDataSeparated_no_outlier_no_na_W$weight_queen)
par(mfrow=c(1,1))
plot(resids ~ AllDataSeparated_no_outlier_no_na_W$treatment)

qqp(AllDataSeparated_no_outlier_no_na_W$weight_queen) ### this one seems okey too

qqp(AllDataSeparated_no_outlier_no_na_W$worker_num_col_2nd_count) ### 
qqp(AllDataSeparated_no_outlier_no_na_W$cuck_num_col) ###
qqp(AllDataSeparated_no_outlier_no_na_W$number_dead_workers_col) ###

#---------------------------------------------------AGGRESSION ANALYSIS------------------------------------

AllDataAgr <- read.delim("Observations - Aggressivness3bis.csv", header = TRUE, sep = ";")
AllDataAgr$replicate_number<-as.factor(AllDataAgr$replicate_number)
AllDataAgr$treatment<-as.factor(AllDataAgr$treatment)
AllDataAgr$treatment.split<-as.factor(AllDataAgr$treatment.split)
AllDataAgr <- filter(AllDataAgr, replicate_number != '22')
AllDataAgr <- filter(AllDataAgr, replicate_number != '40')
?Anova


plot(AllDataAgr$treatment,AllDataAgr$Indice_5, col = c("indianred1","palegreen2","lightskyblue"), xlab="Treatment", ylab = "Aggression Index")


#Taking out the colonies with no interactions toward the other colony
AllDataAgr <- filter(AllDataAgr, Tot_num_interact_twd != '0')
summary(AllDataAgr$treatment)

summary(AllDataAgr$treatment)
n1= lmer(Indice_5 ~ treatment +  (1|replicate_number) + (1|colony_fem_fam),AllDataAgr)
Anova(n1,type="2")
summary(n1)
lsmeans(n1, list(pairwise ~ treatment), adjust = 'fdr')

par(mfrow=c(1,3))
resids = resid(n1)
plot(resids)
fits = fitted(n1)
plot(fits, resids)
qqnorm(resids); qqline(resids, col = 2)
sort(resids)

par(mfrow=c(1,1))
boxplot(resids ~ AllDataAgr$treatment) 
leveneTest(resids ~ AllDataAgr$treatment)
plot(resids ~ AllDataAgr$treatment)
plot(resids ~ AllDataAgr$Indice_5)

#total number of observations positively correlated with the total number of interactions?
cor.test(AllDataAgr$Nb_obs_events,AllDataAgr$Tot_num_interact_twd)


#total number of aggressive interactions correlated to the total number of interactions in general?
cor.test(AllDataAgr$Tot_num_agr_twd,AllDataAgr$Tot_num_interact_twd)


#------------------------------------END-----------------------------------------------------------