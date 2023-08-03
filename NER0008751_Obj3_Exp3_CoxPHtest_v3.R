#DavidHCollins
#Created: 25 February 2020
#Modified 02 Ferbruary 2021

#Aim

#The aim of this script is to analyse data from our main Drosophila experiment (please see NER0008751_Obj3_Exp3_Datasheet_v7_2018-10-18.xlsx for the original raw data and Protocol_NER0008751_Objective 3_Experiment 3_v1_2018-05-05 for the experimental protocol). This particular script details the Cox Proportional Hazards regression analysis on age specific mortality data.

#Method

#In the experiment we reared Drosophila larvae in 3 larval dietary treatments (21 vials per treatment)
#These are labelled treatment 1 (20% SYA), 2 (100% SYA) and 3 (150% SYA)
# For the life history females (see the referenced protocol) we then collected 56 females from each treatment, mated them in a mating chamber (10 females and 10 males per chamber), then placed them on individually on 110% SYA. Every day for 7 days we then measured whether they were alive or dead (scored twice daily), their egg counts and their offspring counts (11 days later). After 7 days we measured their survival twice daily and their egg and offspring counts every two days. 

#We have collected the following data:

# 1) Longevity (days of life as an adult, i.e. their date of death - date of eclosion in days)
# 2) Fecundity (egg count) on each of the measurement days (Every day for 7 days, then every two days for the rest of their lives)
# 3) Total fecundity and Mean daily fecundity across all of the measurement days
# 4) Adult offspring produced on each of the measurement days (Every day for 7 days, then every two days for the rest of their lives)
# 5) Total offspring produced and Mean daily offspring production across all of the measurement days

########################################################################################################################

#First you must always reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script
library(dplyr)
library(survival)
library("survminer") #For survival figures, also loads ggpubr which is needed for this package
library(bbmle) #produces AIC values which allows you to compare multiple models against each other rather than laboriously doing pairwise comparisons between each model
library(readr) #read in files
library(coxme) #for mixed effects coxph models

#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor

#Tell R where to look
setwd("~/Documents/UEA Work documents/Post Doc_LHT/Objective 3 Experiment 3/Analysis/")

#The file contains individual level data for fecundity and longevity
NER0008751_Obj3_Exp3_LH_data_v2 <- read_csv("NER0008751_Obj3_Exp3_LH_data_v2.csv")
View(NER0008751_Obj3_Exp3_LH_data_v2)
df<-NER0008751_Obj3_Exp3_LH_data_v2

#1st Convert the treatment variable into a factor
df$treat<-as.factor(df$treat) 
is.factor(df$treat)

#create a new dataframe where T1 is first for figures (we will continue to use the original dataframe for models)
df2 <- mutate(df, treat = relevel (treat, ref = "20% SYA", "100% SYA"))

#Check that the order is now correct
levels(df$treat)
levels(df2$treat)

#then make the survival curve, with censors
fit1 <- survfit(Surv(age,censor) ~ treat, data = df2)

(plot1 <- ggsurvplot(fit1, 
                    data=df2,
                    censor=FALSE, 
                    legend="right",
                    legend.title = "Treatment",
                    conf.int = FALSE,
                    linetype = c( "twodash", "longdash", "solid"),
                    legend.labs= c("L (20% SYA)", "M (100% SYA)", "H (120% SYA)"),
                    ggtheme = theme_classic2(base_size=12),
                    font.family = "Arial",
                    palette = c("#D55E00", "#0072B2", "#56B4E9"))+
  labs(Title = FALSE, x = "Time (days since eclosion)", y = "Proportion alive"))


ggsave('plot_survival.png', plot1$plot, width = 16, height = 10, units = "cm") #Note. important when using ggave for ggsurvplot objects to specify the '$plot', as ggsurvplot produces a plot and a list, and ggsave is not able to interpret both together

#To make figure panel
ggsave('plot_survival.svg', plot1$plot, width = 16, height = 8, units = "cm") #Note. important when using ggave for ggsurvplot objects to specify the '$plot', as ggsurvplot produces a plot and a list, and ggsave is not able to interpret both together



##############################2.0 test Cox P assumptions###########################################

##Get summary statistics
tapply(df$age[df$censor==1],df$treat_no[df$censor==1],median)
tapply(df$age[df$censor==1],df$treat_no[df$censor==1],quantile)

##Testing for proportional hazards: graphical test for all treatment categories on same plot##
fit1$strata #tells us the number of UNIQUE TIME observations in the individual strata and the order they occurred, note strata are the individual survival curves for each treatment
lines<-rep(1:3,fit1$strata) #indicator variable to identify the strata
print(lines)
plot(fit1$time,log(-log(fit1$surv)),type='n',xlab="t",ylab=expression(log(-log(S(t)))))
lines(fit1$time[lines==1],log(-log(fit1$surv[lines==1])),col=1,type='s')
lines(fit1$time[lines==2],log(-log(fit1$surv[lines==2])),col=2,type='s')
lines(fit1$time[lines==3],log(-log(fit1$surv[lines==3])),col=3,type='s')
legend(locator(1),c("20% SYA","100% SYA","120% SYA"),col=c(1:3),lty=1)

##Testing for proportional hazards: analytical test on coxph model##
coxmodel1 <- coxph(Surv(age,censor) ~ treat, data = df) #Builds the Coxph model where age is modelled against the larval treatment
#Alternative syntax: coxmodel1<-coxph(Surv(df$age,df$censor)~factor(df$treat)) 
Test1<-cox.zph(coxmodel1)
Test1  ##tells us if individual regression coefficients are constant over time, H0: regression coefficient is constant, but don't take too literally
##we can accept that the regression coefficients are constant over time (p>0.05) and that the slope is not significantly different from zero, for all treatment levels
plot(Test1)

summary(coxmodel1) #No effect on survival between different treatments
#This test shows that there was no effect of treatment on survival, the 20% and 120% treatments did not differ significantly from the 100% treatment (coxph: z = -1.234, p = 0.217; z = -1.799, p =  0.072 respectively)


####test effect of size on the cox model

#First test if there is an effect of thorax on age
thoraxmodel <- lm(age ~ treat*thorax, data = df)
summary(lm.age.th) #No significant effect, but nearly significant (df=136, t=-1.956, p=0.0526) so try model with and without size as a factor

#Create a df where only individuals with thorax observations are included

df3 <- df %>% 
  filter(!is.na(thorax))

#Without size, using 100% as the intercept
coxmodel2 <- coxph(Surv(age,censor) ~ treat, data = df3) #Builds the Coxph model where age is modeled against the larval treatment

coxmodel3 <- coxph(Surv(age,censor) ~ treat*thorax, data = df3) #Builds the Coxph model where age is modeled against the larval treatment

coxmodel4 <- coxph(Surv(age,censor) ~ treat + thorax, data = df3) 

fit1 <- coxph(Surv(age,censor) ~ thorax, data = df3) 
AIC(coxmodel2,coxmodel3,coxmodel4,coxmodel5)

#Best model that contained body size was model 5
summary (coxmodel5)

#Test coxmodel 5 assumptions
#Plot figure
(plot1 <- ggsurvplot(fit1, 
                     data=df3))

#Are regression coefficients constant with time?
Test2<-cox.zph(coxmodel2)
Test2  ##tells us if individual regression coefficients are constant over time, H0: regression coefficient is constant, but don't take too literally
##we can accept that the regression coefficients are constant over time (p>0.05) and that the slope is not significantly different from zero, for all treatment levels
plot(Test2)


# Overall Coxmodel 2 is better as the AIC shows that there is a better model fit when size is not included as a factor. Cox model 1 is better than cox model 2 as it includes a greater number of observations (individuals that were excluded from model 2 as they did not have any thorax measurements). Therefore, the best model comparing treatments does not include body size as an independent fixed effect.


###Population - build a model where population is included as a random effect

#ensure source_pop is a factor?
is.factor(df$source_pop) #No so make it a factor
df$source_pop <- as.factor(df$source_pop)
is.factor(df$source_pop) #Yes

# test if there is an effect of population on age
population.model <- lm(age ~ treat*source_pop, data = df)
summary(population.model)#no sign of an effect on age

#Build model when source population is included as a random effect
coxmodel5 <- coxme(Surv(age,censor) ~ treat + (1|source_pop), data = df) #Builds the Coxph model where age is modeled against the larval treatment
AIC(coxmodel1,coxmodel5)#no difference when source population is included as a random effect so go with simple model

#From all of these tests it looks like model 1 is the best model

###Effects when data is partitioned into early and late deaths

#We can also partition the model into two parts, part 1 where most of the individuals stay alive, with a few deaths, and part 2 the main death phase. This was done because it looked like (from the figure) there was a treatment specific difference in deaths in the first part of the experiment

#Just the first part of the model, use censor2 column which will only include the individuals that died before day 55
#Using 100% as the intercept
coxmodel6 <- coxph(Surv(age,censor2) ~ treat, data = df) #Builds the Coxph model where age is modelled against the larval treatment
summary(coxmodel6)
#This test shows that there was no effect of treatment on survival, the 20% and 120% treatments did not differ signficantly from the 100% treatment (coxph: z = -1.406, p = 0.160; z = -1.246, p =  0.213 respectively)

#Just the second part of the model, use censor3 column which will only include the individuals that died after day 55
#Using 100% as the intercept
coxmodel7 <- coxph(Surv(age,censor3) ~ treat, data = df) #Builds the Coxph model where age is modelled against the larval treatment
summary(coxmodel7)
#This test shows that there was no effect of treatment on survival, the 20% and 120% treatments did not differ signficantly from the 100% treatment (coxph: z = -0.718, p = 0.473; z = -1.367, p =  0.172 respectively)

###Was asked by reviewer to calculate mortality rate over time:

#For this load in a new df with mortality and survivorship data per day
df3 <- read_csv("NER0008751_Obj3_Exp3_mortality_rate.csv")

#Filter all values where the number alive at the end of the day is 0 (i.e., the day the last fly died in each treatment)
df4 <-   df3 %>% 
  filter (number_alive > 0)


# Calculate the mortality rate for each day
df4$mortality_rate <- df4$number_dead / (df4$number_alive + df4$number_dead)

# Calculate mux using the natural log of the mortality rate
df4$mux <- -log(1 - df4$mortality_rate)

#Make a new natural log of mux
df4$logmux <- log(df4$mux)


#Plot mortality rate against age
mort.rate.figure <- ggplot(df4, aes(x = day, y = logmux, colour = treatment, group = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Age", y = bquote('Mortality rate ('*mu[x]~')'), colour = 'Larval diet')+
  scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
  theme_classic()

ggsave("Mortality rate.png", mort.rate.figure, width = 16, height = 10, units = "cm")


#Mortality analysis
mort.model <- lm(mux ~ treatment*day, data = df4)


autoplot(mort.model) # Assumptions violated

mort.model2 <- glm(mux ~ treatment*day, data = df4)

autoplot(mort.model2)


mort.model3 <- glm(mux ~ treatment*day + day^2, data = df4)

autoplot(mort.model3)

# Does treatment affect mortality rate?
anova(mort.model)

#How do the different components of the model affect mortality rate?
summary (mort.model)

#These data, and the figure, show that mortality rate increases with age (which you would expect as flies age), but there was no effect of treatment on mortality rate, and no interaction between treatment and age in their combined effects on mortality rate.
