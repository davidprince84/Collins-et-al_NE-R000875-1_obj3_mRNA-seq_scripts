#DavidHCollins
#Created: 16 February 2021

#Aim

#The aim of this script is to analyse the size data from our main Drosophila experiment (please see NER0008751_Obj3_Exp3_Datasheet_v7_2018-10-18.xlsx for the original raw data and Protocol_NER0008751_Objective 3_Experiment 3_v1_2018-05-05 for the experimental protocol).

#The reason we are doing this is to establish if size varied according to treatment (which will then affect how we code our models)


#First you must always reset R
rm(list=ls())

#make these packages (and associated functions) available to use in the script

library(tidyverse)
library(readr)
library(ggfortify) #Allows autoplot to accept objects of type lm

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

#Omit the na values
df<-data.table(df)
df2<-na.omit(df, cols=c("thorax"))

# Note this excel file is a csv version of the raw LH data in NER0008751_Obj3_Exp3_Datasheet_v6_2018-10-09.xlsx

#do you have the right data?
glimpse(df2)

#due to a high death rate at the start of the experiment and the occasional accidental death throughout the experiment, some of the data had to be censored. Therefore, we need to create a new data frame looking at just the individuals that were NOT censored (i.e. were classified as 'd' in the censor column)
df3 <- filter (df2, censor == "1")

#also need to re-order the treatment names so they follow a logical order

#1st Convert the treatment variable into a factor
df3$treat<-as.factor(df3$treat) 
is.factor(df3$treat)

#create a new dataframe where T1 is first for figures (we will continue to use the original dataframe for models)
df3b <- mutate(df3, treat = relevel (treat, ref = "20% SYA", "100% SYA"))

#Check that the order is now correct
levels(df3$treat)
levels(df3b$treat)

# 1.6.1 boxplot: effect of treatment on average thorax length
(thoraxplot <- ggplot(df3b,aes(x=treat, y=thorax, fill = treat)) +
  geom_boxplot()+
  ggbeeswarm::geom_quasirandom()+
  labs(x="Treatment diet", y = "Thorax size (mm)")+
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_discrete(limits = c("20% SYA", "100% SYA", "120% SYA"),
                     labels = c("L (20% SYA)", "M (100% SYA)", "H (120% SYA)")))


ggsave("thoraxplot.png", thoraxplot, width = 16, height = 10, units = "cm")


# 1.6.2 ANOVA effect of treatment on average thorax length

#Build model
thoraxmodel <- lm(thorax ~ treat, data = df3)

#test assumptions
autoplot(thoraxmodel , smooth.colour = NA) #Plots look okay

# use anova() to test if there is an effect of Treatment
anova(thoraxmodel)

#use summary() to test what the nature of the effect of treatment is, and see if it makes sense with the boxplot itself
summary(thoraxmodel)

# to summarize there was a significant effect of treatment on thorax size in Drosophila melanogaster, adults from treatment 1 were smaller than adults from treatment 2 and 3. Adults from treatments 2 and 3 were not significantly different from each other.

  # 1.7.1 boxplot: effect of treatment on average wing cell length
  
  
  #boxplot
 ( wingplot <- ggplot(df3b,aes(x=treat, y=wing)) +
  geom_boxplot()+
  labs(title="The effect of larval dietary regime on female marginal cell length in Drosophila melanogaster females",x="Larval diet", y = "Marginal cell length (mm)")+
  theme_classic())


# 1.7.2 ANOVA effect of treatment on average wing length

#Build model
wingmodel <- lm(wing ~ treat, data = df3)

#test assumptions
autoplot(wingmodel , smooth.colour = NA) 

# use anova() to test if there is an effect of Treatment
anova(wingmodel)

#use summary() to test what the nature of the effect of treatment is, and see if it makes sense with the boxplot itself
summary(wingmodel)

# to summarize there was a significant effect of treatment on wing size in Drosophila melanogaster, adults from treatment 1 were smaller than adults from treatment 2 and 3. Adults from treatments 2 and 3 were not significantly different from each other.

#Correlation between wing size and thorax length:

(correlationplot<- ggplot(df3b,aes(x=wing, y=thorax)) +
  geom_point()+
  labs(title="The relationship between thorax and wing",x="Wing length (mm)", y = "Thorax length (mm)")+
  theme_classic())

correlationmodel <- lm(thorax~wing, data=df3)
autoplot(correlationmodel) #assumptions not violated, although be aware of a couple of influential datapoints
summary(correlationmodel)

# To summarise, the two size measurements are highly correlated. As wings deteriorated throughout the experiment then we will use thorax as it was a more reliable measure of size which stayed constant and for which we had more samples.


#Test this using a glm

thoraxmodel2 <- glm(thorax~treat, data = df3, family = gaussian (link = log))

summ(thoraxmodel2, digits = 3)
autoplot(thoraxmodel2)
plot(allEffects(thoraxmodel2))
anova(thoraxmodel2, test = 'Chisq')

#Make a model so we can compare 120% and 20% to each other

thoraxmodel2b <- glm(thorax~treat, data = df3b, family = gaussian (link = log))
summ(thoraxmodel2b, digits = 3)


#Note that the gaussian glm is effectively the same as the linear model

#Try a gamma glm instead
thoraxmodel3 <- glm(thorax~treat, data = df3, family = Gamma(link = 'log'))
plot(allEffects(thoraxmodel3))
AICtab(thoraxmodel,thoraxmodel2,thoraxmodel3)
#The gaussian glm (which is equivalent to a linear model) is better

#Find out the mean, medians and sds of the thorax sizes for each treatment
df3 %>% 
  group_by(treat) %>% 
  dplyr::summarise(mean=mean(thorax),
                   sd = sd(thorax),
                   median = median(thorax),
                   n = n())


