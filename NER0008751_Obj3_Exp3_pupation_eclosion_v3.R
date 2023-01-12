#DavidHCollins
#Created: 26 January 2021

#Aim

#The aim of this script is to analyse data from our main Drosophila experiment (please see NER0008751_Obj3_Exp3_Datasheet_v7_2018-10-18.xlsx for the original raw data and Protocol_NER0008751_Objective 3_Experiment 3_v1_2018-05-05 for the experimental protocol).

#Method

#In the experiment we reared Drosophila larvae in 3 larval dietary treatments (21 vials per treatment)
#These are labelled treatment 1 (20% SYA), 2 (100% SYA) and 3 (150% SYA)
# For the life history females (see the referenced protocol) we then collected 56 females from each treatment, mated them in a mating chamber (10 females and 10 males per chamber), then placed them on individually on 110% SYA. Every day for 7 days we then measured whether they were alive or dead (scored twice daily), their egg counts and their offspring counts (11 days later). After 7 days we measured their survival twice daily and their egg and offspring counts every two days. 

#We have collected the following data:

# 1) Pupation success for each vial
# 2) Mean pupation time (in hours) for each vial
# 3) Eclosion success for each vial
# 4) Eclosion time (in hours) for each vial
# 5) The researcher that set up the vials
# 6) The source population for each vial
# 7) Pupation duration (the eclosion time minus the pupation time)
# 8) The sex of each individual that eclosed in each vial

#We will use these data to establish whether there is an effect of treatment on pupation/eclosion success and development time. Originally I used ANOVA to analyse these data, however ANOVA is not really appropriate for count data so I will use either glms or glmms (depending on whether I need to incorportate mixed effects).


########################################################################################################################

#First you must always reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(lme4) ## for glmer
library(tidyverse)
library(AER) # for other glms on count data and dispersison test
library(jtools) # summ function which is a much nicer easier to read version of summary with some better info
library(effects) #allEffects function that allows us to visualise the model outputs
library(ggbeeswarm) ##add quasirandom function so points can be added to barplots as a 'beeswarm' pattern rather than jitter
library(bbmle) #ALlows use of AICctab
#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor

#Find file using one of two methods (actually four methods in the book but two are beyond me)

library(readr)
setwd("~/Documents/UEA Work documents/Post Doc_LHT/Objective 3 Experiment 3/Analysis/")
#Use this for old computer
#setwd("U:/Documents/Post Doc_LHT/Objective 3 Experiment 3/Analysis/")

NER0008751_Obj3_Exp3_pupation_eclosion.csv<-read_csv("NER0008751_Obj3_Exp3_pupation_eclosion.csv")

#make into a tibble

df <- as.tibble(NER0008751_Obj3_Exp3_pupation_eclosion.csv)

#First make a duplicate df (for use with figures)

dff <- df

# To look at the relationships of individual treatments, we need to ensure that r is treating treatment as a factor

is.factor(df$Treatment) 
df$Treatment <- as.factor(df$Treatment)
dff$Treatment <- as.factor(dff$Treatment)

df <- mutate(df, Treatment = relevel (Treatment, ref = "20% SYA"))
df <-  mutate(df, Treatment = relevel (Treatment, ref = "100% SYA"))
dff <- mutate(dff, Treatment = relevel (Treatment, ref = "20% SYA"))

#Check that the order is now correct

levels(df$Treatment)
levels(dff$Treatment)

glimpse(df)

###############Were there any noticable differences in the pupation and eclosion success between researchers#######

# Box plot to see if there was an effect of researcher on the survival of pupae and adults

(Pupation.Res.bp <- ggplot(dff,aes(x=Researcher, y=Pupae)) +
  geom_boxplot(inherit.aes=TRUE))

(Eclosion.Res.bp <- ggplot(df,aes(x=Researcher, y=Adults)) +
  geom_boxplot())


### GLM to test whether researcher had an effect on Pupation and Eclosion ###

# 1) Pupae

#Build model
Pupation.Res.model1 <- glm(Pupae ~ Researcher,family="poisson",data = df)

#Check for overdispersion
dispersiontest(Pupation.Res.model1,trafo=1) #data are not overdispersed (c = -0.352, when it should be zero or less)

#What is the chisq critical value
qchisq(0.95, df.residual(Pupation.Res.model1))
deviance(Pupation.Res.model1) #below the critical value, implying a good fit to the data

plot(Pupation.Res.model1) #besides the one outlier the plot looks okay (sample 27)

#Use summ function to investigate if researcher had an effect and what the effect was
summ(Pupation.Res.model1, digits=3) #shows that researcher did have an effect which may need to be accounted for in later pupation models

# 2) eclosion

#Build model
Eclosion.Res.model1 <- glm(Adults ~ Researcher,family="poisson",data = df)

dispersiontest(Eclosion.Res.model1,trafo=1) #data are not overdispersed (c = -0.387, when it should be zero or less)

#What is the chisq critical value
qchisq(0.95, df.residual(Eclosion.Res.model1))
deviance(Eclosion.Res.model1) #below the critical value, implying a good fit to the data

plot(Eclosion.Res.model1) # the plot looks less good than pupation but will do for out purposes (sample 27)

#Use summ function to investigate if researcher had an effect and what the effect was
summ(Eclosion.Res.model1, digits=3) #shows that researcher did have an effect which may need to be accounted for in later pupation models

###################################################################################################################
##################Was there an effect of the original oviposition plate on the pupation and eclosion success######

# To look at the relationships of individual oviposition plates, we need to ensure that r is treating oviposition plate as a categorical variable
is.character(df$Oviposition.Plate)
df$Oviposition.Plate <- as.character(df$Oviposition.Plate)
is.character(df$Oviposition.Plate)

# First make a plot to see if there are likely to be any differences between different oviposition plates
Pupation.Ov.bp <- ggplot(df,aes(factor(Oviposition.Plate),x=Oviposition.Plate, y=Pupae)) +
  geom_boxplot()

Eclosion.Ov.bp <- ggplot(df,aes(x=Oviposition.Plate, y=Adults)) +
  geom_boxplot()


Pupation.Ov.bp
Eclosion.Ov.bp


geom_boxplot(aes(factor(group), x))


### glms to test if there was an effect of ovposition plate on pupation and eclosion ###

# 1) Pupation

#Build model
pupation.ov.model1 <- glm(Pupae ~ Oviposition.Plate,family="poisson",data = df)

dispersiontest(pupation.ov.model1,trafo=1) #data are not overdispersed (c = -0.041, when it should be zero or less)
## For paper: Graphical and analytical tests showed that the data were not overdispersed (overdispersion test: c = -0.041, p = 0.638)


#What is the chisq critical value
qchisq(0.95, df.residual(pupation.ov.model1))
deviance(pupation.ov.model1) #below the critical value, implying a good fit to the data

plot(pupation.ov.model1) # the plot looks good (sample 27 is an outlier)

#Use summ function to investigate if researcher had an effect and what the effect was
summ(pupation.ov.model1, digits=3) #shows that oviposition plate did not have a significant effect on number of pupae counted

summary(pupation.ov.model1)

# 2) Eclosion

#Build model
Eclosion.Ov.model1 <- glm(Adults ~ Oviposition.Plate,family="poisson",data = df)

dispersiontest(Eclosion.Ov.model1,trafo=1) 
## For paper: Graphical and analytical tests showed that the data were not overdispersed (overdispersion test: c = -0.226, p = 0.878)

#What is the chisq critical value
qchisq(0.95, df.residual(pupation.ov.model1))
deviance(pupation.ov.model1) #below the critical value, implying a good fit to the data

#test assumptions
plot(Eclosion.Ov.model1, smooth.colour = NA) #these data aren't great (again sample 27 is an outlier) but will do for now

summ(Eclosion.Ov.model1)

#These tests show that there is no significant effect of the oviposition plate (with the possible exception of oviposition plate 3 having slightly fewer adults eclose than the other plates) and therefore Drosophila population used. S ovipostion plate does not need to be accounted for

###################################################################################################################
############################# Was there an effect of treatment on the pupation and eclosion success? ##############

#First use a figure to see if there is likely to be a significant effect of the treatment on the pupation and eclosion success


(Pupation.success.bp <- ggplot(dff,aes(x=Treatment, y=Pupae, fill = Treatment)) +
  geom_boxplot()+ 
  theme_classic()+
  theme(legend.position = "none")+
  ggbeeswarm::geom_quasirandom()+
  labs(x="Treatment", y = "Pupation success")+ 
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
  scale_x_discrete(breaks=c("20% SYA","100% SYA","120% SYA"), labels=c("L (20% SYA)", "M (100% SYA)", "H (120 % SYA)")))

(Eclosion.success.bp <- ggplot(dff,aes(x=Treatment, y=Adults, fill = Treatment)) +
  geom_boxplot()+ 
  theme_classic()+
  theme(legend.position = "none")+
  ggbeeswarm::geom_quasirandom()+
  labs(x="Treatment", y = "Eclosion success")+ 
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
  scale_x_discrete(breaks=c("20% SYA","100% SYA","120% SYA"), labels=c("L (20% SYA)", "M (100% SYA)", "H (120 % SYA)")))

#Save for the paper
ggsave("Pupation.success.bp.png", Pupation.success.bp, width = 16, height = 10, units = "cm")
ggsave("Eclosion.success.bp.png", Eclosion.success.bp, width = 16, height = 10, units = "cm")


### Use GLMMs to test whether the was a significant effect of treatment on pupation and eclosion success

#1) Pupation

#Note, originally I treated values as percentages, but I could also look at the raw numbers with poisson distribution, or as fractions with the binomial distribution. To make the fractions for the binomial, Ill find out what the largest value is (should be around 100 max), and then divide by that value

df %>% 
  arrange(desc(Pupae)) %>% 
  dplyr::slice(1) %>% 
  ungroup()

#102 is the highest value for pupae

df %>% 
  arrange(desc(Adults)) %>% 
  dplyr::slice(1) %>% 
  ungroup() %>% 
  dplyr::select(TubeID, Adults)

# 86 is the highest value of adults

## Each number in the pupae and adults column should be the number of developing pupae and adults out of an initial input of 100 larvae. However due to possible misscounts, the highest pupae value is 104. Therefore to make these into proportions for the binomial glms we shall divide the pupae values by 104 (as nearly 100), and the adult values by 100

df2 <- df %>% 
  mutate(prop.pupae = Pupae/104) %>% 
  mutate(prop.adults = Adults/100)
  
#standard glm with poisson distribution
Pupation.num.model1 <- glm(Pupae ~ Treatment,family="poisson",data = df)

summ(Pupation.num.model1)

plot(allEffects(Pupation.num.model1)) #100% SYA has significantly lower pupation rate than 20%

#standard glm with binomial distribution
Pupation.num.model2 <- glm(prop.pupae ~ Treatment,family="binomial",data = df2)

summ(Pupation.num.model2)

plot(allEffects(Pupation.num.model2)) #100% SYA has significantly lower pupation rate than 20%

#When we look at the plots, the models don't look all that different from each other. I will stick with the poisson as I don't need to mess with the data as much.

#Next make the pupae model into a glmm where researcher is accounted for as a fixed effect:
Pupation.num.model3 <- glm(Pupae ~ Treatment+ Researcher,family="poisson",data = df)

#Is this better than when the random effect is excluded
AIC(Pupation.num.model1,Pupation.num.model3) #model3 appears to be better

#Next make the pupae model into a glmm where researcher and plate are accounted for as additional effects:
Pupation.num.model4 <- glm(Pupae ~ Treatment+ Researcher+ Oviposition.Plate,family="poisson",data = df)


#Next make the pupae model into a glmm where just plate is accounted for as an additional effect:
Pupation.num.model5 <- glm(Pupae ~ Treatment + Oviposition.Plate,family="poisson",data = df)

bbmle::AICctab(Pupation.num.model1, Pupation.num.model3,Pupation.num.model4,Pupation.num.model5)#model 3 appears to be better with researcher as an additional effect rather than oviposition plate
rdf <- as.data.frame(rdf)
## For paper: The best model included researcher but not population as a random effect.

Anova(Pupation.num.model5, test ='Chisq')
#Make a null version of model 3 to test the significance of the fixed effects
Pupation.num.model6 <- glm(Pupae ~ 1 + Researcher,family="poisson",data = df)

AICctab(Pupation.num.model3,Pupation.num.model6) #model 6 (the null model) appears to be a better fit to the data

anova(Pupation.num.model3,Pupation.num.model6) #no difference with the null model
##For paper: the best model containing fixed effects was not significantly different from the null model (pupation success: LRT: Chisq = 0.931, p = 0.628) 

Anova(Pupation.num.model3, test ='LR')
summ(Pupation.num.model3, digits = 3) #There was no effect of treatment on number of pupae with the full model
summary(Pupation.num.model3, digits = 3)
#Plot predictions of the best model
plot(allEffects(Pupation.num.model3)) 


#Plot fitted vs residuals
 plot(fitted(Pupation.num.model3),residuals (Pupation.num.model3)) #One outlier but otherwise the model looks fine

#Overall there appears to be no effect of treatment on pupae number

###### 2) Eclosion

#standard glm with poisson distribution
Eclosion.num.model1 <- glm(Adults ~ Treatment,family="poisson",data = df)

summ(Eclosion.num.model1)

plot(allEffects(Eclosion.num.model1)) #120% SYA has significantly higher pupation rate than 20%


#glmmm with researcher as a fixed effect
Eclosion.num.model3 <- glm(Adults ~ Treatment+ Researcher,family="poisson",data = df)

#Is this better than when the random effect is excluded
AIC(Eclosion.num.model1,Eclosion.num.model3) #model3 appears to be slightly better

#glmmm with researcher and plate as a random effect
Eclosion.num.model4<- glm(Adults ~ Treatment + Researcher+ Oviposition.Plate,family="poisson",data = df)

#glmmm with researcher and plate as a random effect
Eclosion.num.model5<- glm(Adults ~ Treatment+ Oviposition.Plate,family="poisson",data = df)

AIC(Eclosion.num.model3,Eclosion.num.model5) #model 3 appears to be better, therefore researcher is improves model but oviposition plate does not

#Make a null version of model 3 to test the liklihood
Eclosion.num.model6 <- glm(Adults ~ 1 + Researcher,family="poisson",data = df)

AICctab(Eclosion.num.model3,Eclosion.num.model6) #model 3 performs worse than the null model

#Do  a liklihood ratio test to compare the models
Anova(Eclosion.num.model3, test ='LR')

summ(Eclosion.num.model3, digits = 3) #There was no effect of treatment on number of adults with the full model

plot(allEffects(Eclosion.num.model3)) 


#these data show that there was no effect of treatment on eclosion success

################## Was there an effect of treatment on the mean pupation and eclosion times for each vial? #####################################################################################################

#Visualise as boxplots

(Pupation.time.bp <- ggplot(dff,aes(x=Treatment, y=Pup_time, fill = Treatment)) +
   geom_boxplot()+ 
   theme_classic()+
   geom_quasirandom()+
   theme(legend.position = "none")+
   labs(x="Treatment diet", y = "Pupation time (hours)")+ 
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
  scale_x_discrete(breaks=c("20% SYA","100% SYA","120% SYA"), labels=c("L (20% SYA)", "M (100% SYA)", "H (120 % SYA)")))

(Eclosion.time.bp <- ggplot(dff,aes(x=Treatment, y=Ecl_time, fill = Treatment)) +
    geom_boxplot()+ 
    theme_classic()+
    theme(legend.position = "none")+
    geom_quasirandom()+
    labs(x="Treatment diet", y = "Eclosion time (hours)")+ 
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
  scale_x_discrete(breaks=c("20% SYA","100% SYA","120% SYA"), labels=c("L (20% SYA)", "M (100% SYA)", "H (120 % SYA)")))

#Save for the paper
ggsave("Pupation.time.bp.png", Pupation.time.bp, width = 16, height = 10, units = "cm")
ggsave("Eclosion.time.bp.png", Eclosion.time.bp, width = 16, height = 10, units = "cm")


#Pupation models

#standard glm with gamma distribution
Pupation.time.model1 <- glm(Pup_time ~ Treatment,family="Gamma",data = df)
Pupation.time.model2 <- glm(Pup_time ~ Treatment+Oviposition.Plate,family="Gamma",data = df)
Pupation.time.model3 <- glm(Pup_time ~ 1,family="Gamma",data = df)

AICctab(Pupation.time.model1,Pupation.time.model2,Pupation.time.model3) #model 2 performs best but not much better than simplified model, so go with simplified model

#Better than null?
Anova(Pupation.time.model1, test ='LR')
#Model summary?
summ(Pupation.time.model2, digits = 3)


#Eclosion models

#standard glm with gamma distribution
Eclosion.time.model1 <- glm(Ecl_time ~ Treatment,family="Gamma",data = df)
Eclosion.time.model2 <- glm(Ecl_time ~ Treatment+Oviposition.Plate,family="Gamma",data = df)
Eclosion.time.model3 <- glm(Ecl_time ~ 1,family="Gamma",data = df)

#Best model
AICctab(Eclosion.time.model1,Eclosion.time.model2,Eclosion.time.model3) #model 2 performs better than simplified model, so go with model2

#Better than null?
Anova(Eclosion.time.model2, test ='LR')

#Model summary
summ(Eclosion.time.model2, digits = 3)


#Summarise the mean, sd, and sample size of pupation times (for the paper)
df %>%
  group_by(Treatment) %>%
  summarise(mean_pup_time = mean(Pup_time), 
            sd_pup_time = sd (Pup_time),
            mean_ecl_time = mean(Ecl_time), 
            sd_ecl_time = sd (Ecl_time),n = n())