#Aim

#The aim of this script is to analyse data from our main Drosophila experiment (please see NER0008751_Obj3_Exp3_Datasheet_v7_2018-10-18.xlsx for the original raw data and Protocol_NER0008751_Objective 3_Experiment 3_v1_2018-05-05 for the experimental protocol). This particular script details methods to analyse age specific reproduction data. V3 of the script told us the following information:
  #1) A glmm was required to analyse the numerical count data
  #2) Negative binomial performed better than Poisson
  #3) With the standard Negative Binomial model there was still significant zero inflation - so a Zero-inflated negative binomial was fitted using glmmTMB package. This performed better than the standard negative binomial which was fitted using the same package.
  #4) As the models included time series data and the effect of time had a non-linear effect on egg/offspring count then a quadratic function was included. This improved the fit of the model to the data.
  #5) Two random effects were tested (using AICc). A random effect with individual ID included improved the model fit, but a random effect with source population did not.

#V4 uses this information from V3 and tests whether the following fixed effects improve the model fit: Treatment, Age, body size, and their interactions. We will test this for four different hypotheses:

      ## 1) There is an effect of treatment on egg number in the first week
      ## 2) There is an effect of treatment on egg number across the whole experiment
      ## 3) There is an effect of treatment on offspring number in the first week
      ## 4) There is an effect of treatment on offspring number across the whole experiment

#First you must always reset R
rm(list=ls())

#Load in these packages
library(lme4) ## for glmer
library(MASS) # for negative binomial glmms
library(bbmle) # for AICctab 
library(plyr)
library(tidyverse)
sessionInfo()
library(glmm)
library(readr)
library(data.table) #data.table function
library(effects) ## All effects function which gives confidence intervals of actual values when using a glm
library(jtools) ## summ function which is a much nicer easier to read version of summary with some better info
library(DHARMa) #Use for dispersion and zero-inflation test function
library(glmmTMB) ##for zero inflated glmms
library(beeswarm) ##add quasirandom function so points can be added to barplots as a 'beeswarm' pattern rather than jitter

#Tell R where to look
setwd("~/Documents/UEA Work documents/Post Doc_LHT/Objective 3 Experiment 3/Analysis/")


#Load in age dependent fertility data
NER0008751_Obj3_Exp3_dailyoff_3 <- read_csv("NER0008751_Obj3_Exp3_dailyoff_3.csv")
df<-NER0008751_Obj3_Exp3_dailyoff_3

View(df)
names(df) ##to check column titles 
str(df) ##to check all in correct format


# Remove missing values by turning df from a dataframe to a datatable, and then omit all the missing rows in the egg_tot and off_tot column. Do it using the data.table function as we do not want to lose rows which have missing size data 
df<-data.table(df)
df2<-na.omit(df, cols=c("egg_tot", "off_tot"))

#The following individuals died before day 4, these should not be counted so remove them:
#'L17', 'L20', 'L42', 'L50', 'L113', 'L142', 'L149', 'L154', 'L160', 'L162', 'L163'
df2a <- df2 %>% 
  filter(!id %in% c('L17', 'L20', 'L42', 'L50', 'L113', 'L142', 'L149', 'L154', 'L160', 'L162', 'L163'))


#ensure you still have data for every individual (there should be 168 unique IDs for df2 and 157 for df2a)
df %>% summarise(Unique_Elements = n_distinct(id)) #correct number of ids for df
df2 %>% summarise(Unique_Elements = n_distinct(id)) #correct number of ids for df2
df2a %>% summarise(Unique_Elements = n_distinct(id)) #correct number of ids for df2


#Make dataframes with week instead of day to smooth out daily trends - please note that 'week' is actually every 8 day period between each mating. Each 8 day period was used because it was easier to measure fertility every two days than alternating between 2 and 3 days.
df3<- df2a %>%
  group_by (week,id,id_no,treat,size) %>%
  dplyr::summarise(
    egg_tot = sum(egg_tot),
    off_tot = sum(off_tot),
    viability = off_tot/egg_tot 
  )

#Do we still have all of the data
df3 %>% summarise(Unique_Elements = n_distinct(id)) #correct number of ids for df3


# Make a series of new dfs that will aid in the visual representation of the data on a by-week basis rather than a two daily basis;


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = sd(x[[col]], na.rm=TRUE)/sqrt(length(x))
    )
  }
  data_sum<-plyr::ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#the new dataframe summarising mean and sd for each week and treatment for egg total
df4a <- data_summary(df3, varname="egg_tot",
                     groupnames=c("treat", "week"))

#the new dataframe summarising mean and sd for each week and treatment for offspring total
df4b <- data_summary(df3, varname="off_tot",
                     groupnames=c("treat", "week"))

#Merge the two dataframes into a new table
df4 <- merge(df4a,df4b,by=c("treat","week"))
setnames(df4, old=c("se.x","se.y"), new=c("egg_se","off_se"))

#change treat to factors so the order can be changed
df4$treat<-as.factor(df4$treat) 

# Use df5 for graphics as the order is more sensible to look at, use df 4 for modeling as the order makes more sense as a control
df5 <- mutate(df4, treat = relevel (treat, ref = "20% SYA", "100% SYA"))



#### EGG COUNTS ####

### 1) EGG COUNTS FOR THE FIRST WEEK

## Make a plot of daily egg counts for the first week

#Filter Data for just the first week
df6 <- dplyr::filter(df2a, week == 1)

#Make plot
(early_egg_plot<- df6 %>%
   group_by (day,id,treat,size) %>%
   dplyr::summarise(
     egg_tot = sum(egg_tot),
     off_tot = sum(off_tot),
     viability = off_tot/egg_tot 
   ) %>% 
   ggplot(aes(x=as.factor(df6$day), y=egg_tot, colour=treat))+
   geom_boxplot()+
   labs(title="Average number of eggs produced per treatment per week in Drosophila melanogaster females",x="Day", y = "Eggs per day")+
   theme_classic() + 
   scale_color_brewer(type = 'div', palette = 'Set1', direction = 1))

## Test out the possible models

#The full model
egg.early.model1 <- glmmTMB(egg_tot ~ treat*poly(day,2) + (1|id), data=df6, ziformula=~1, family = nbinom1)

#No interaction
egg.early.model2 <- glmmTMB(egg_tot ~ treat+poly(day,2) + (1|id), data=df6, ziformula=~1,family = nbinom1)

#Just treatment
egg.early.model3 <- glmmTMB(egg_tot ~ treat + (1|id), data=df6, ziformula=~1,family = nbinom1)

#Just day
egg.early.model4 <- glmmTMB(egg_tot ~ poly(day,2) + (1|id), data=df6, ziformula=~1,family = nbinom1)

#Null
egg.early.model5 <- glmmTMB(egg_tot ~ 1 + (1|id), data=df6, ziformula=~1,family = nbinom1)

(egg.early.selection.table <- as.data.frame(AICtab(egg.early.model1,egg.early.model2,egg.early.model3,egg.early.model4,egg.early.model5, weights = T, base = T, logLik = T)))

plot(allEffects(egg.early.model1))
plot(allEffects(egg.early.model2))
(egg.early.model1.coef<- summary(egg.early.model1)$coefficients[["cond"]])
#treatments are borderline significant, but the interaction is not so try model without the interaction
summary(egg.early.model1)
summary(egg.early.model2)
summary(allEffects(egg.early.model1))
#Best to stick with egg model 1
  
#Figure representing the raw data  
  df6%>% group_by(id,treat) %>% 
    dplyr::summarise(egg_mean=mean(egg_tot)) %>% 
    ggplot(aes(x=treat, y=egg_mean, colour=treat))+
    geom_boxplot()+
    geom_quasirandom()+
    labs(title="mean eggs first week",x="Treatment", y = "Eggs per day")+
    theme_classic() + 
    scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)

  #Determine what colours R selected:
  RColorBrewer::display.brewer.pal(n = 8, name = 'Set1')
  RColorBrewer::brewer.pal(n = 8, name = "Set1") #Colours correspond to: "#E41A1C" "#377EB8" "#4DAF4A"
  
#Find mean value per individual per treatment to quote in the paper
  df6 %>% 
    group_by(treat) %>% 
    dplyr::summarise(mean=mean(egg_tot),
                     sd = sd(egg_tot),
                     median = median(egg_tot))

### 2) EGG COUNTS FOR WHOLE EXPERIMENT

## Make a plot of egg count per treatment

  (egg_plot <- df3  %>% 
      group_by(id,treat) %>% 
    dplyr::summarise(egg_mean=median(egg_tot)) %>% 
    ggplot(aes(x=treat, y=egg_mean))+
    geom_boxplot(outlier.shape = NA)+
    labs(title="egg plot A",x="Treatment", y = "Mean eggs per period")+
    theme_classic() + 
    scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)+
    scale_y_continuous(breaks=seq(0,140,20))+
    scale_x_discrete(limits = c("20% SYA", "100% SYA", "120% SYA"),
                     labels = c("L (20% SYA)", "M (100% SYA)", "H (120% SYA)")))
  
  ## Make a plot of weekly egg counts for the entire experiment (Note that 'week' is actually an 8 day period between each mating not week as is normally meant)
  
  
  
  (egg_plot2<- ggplot(df5, aes(x=week, y=egg_tot, group = treat, shape=treat, colour=treat))+
      geom_errorbar(aes(ymin=egg_tot-egg_se, ymax=egg_tot+egg_se), width=.1, position=position_dodge(0.25)) +
      geom_line() +
      geom_point()+
      labs(x="Experiment period", y = "Whole-life egg production", shape = "Treatment", colour = "Treatment")+
      theme_classic() + 
      theme(legend.position = c(0.9,0.5))+
      scale_shape_manual(values=c(17, 16, 15), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
      scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
      scale_x_continuous(breaks=seq(0,11,1)))

  #TC prefers these data presented as time (i.e., including day number) so make a new collumn in df5 which converts 'week' into the last day number of the last day for each experimental period Easiest way to do this is to multiply 'week' by eight.
  
 df5b <- df5 %>% 
   mutate(last_day = week*8)
 
 (egg_plot3<- ggplot(df5b, aes(x=last_day, y=egg_tot, group = treat, shape=treat, colour=treat))+
     geom_errorbar(aes(ymin=egg_tot-egg_se, ymax=egg_tot+egg_se), width=.1, position=position_dodge(2)) +
     geom_line() +
     geom_point()+
     labs(x="Time (days)", y = "Whole-life egg production", shape = "Treatment", colour = "Treatment")+
     theme_classic(base_size = 18) + 
     theme(legend.position = c(0.9,0.5))+
     scale_shape_manual(values=c(17, 16, 15), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
     scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
     scale_x_continuous(breaks=seq(0,88,8)))
 
 
  
ggsave("Egg plot A.png", egg_plot3, width = 8, height = 6, units = "cm")

#To make figure panel
ggsave("Egg plot A.svg", egg_plot3, width = 21, height = 15, units = "cm")

#Find mean value per individual per treatment to quote in the paper
df3%>% group_by(id,treat) %>% 
  dplyr::summarise(egg_mean=mean(egg_tot)) %>% 
  group_by(treat) %>% 
  dplyr::summarise(egg_mean=mean(egg_mean))



## Test out the possible models

#Full model
egg.model1 <- glmmTMB(egg_tot ~ treat*poly(week,2) + (1|id),data=df3,ziformula=~1,family = nbinom1) 

#No interaction
egg.model2 <- glmmTMB(egg_tot ~ treat+poly(week,2) + (1|id),data=df3,ziformula=~1, family = nbinom1)

#Just treatment
egg.model3 <- glmmTMB(egg_tot ~ treat + (1|id),data=df3,ziformula=~1,family = nbinom1) 

#Just week
egg.model4 <- glmmTMB(egg_tot ~ poly(week,2) + (1|id),data=df3,ziformula=~1,family = nbinom1) 

#Null
egg.model5 <- glmmTMB(egg_tot ~ 1 + (1|id),data=df3,ziformula=~1,family = nbinom1) 

(egg.selection.table <- as.data.frame(AICtab(egg.model1,egg.model2,egg.model3,egg.model4,egg.model5, weights = T, base = T, logLik = T)))

#Egg model 4 is best but egg models 2 and 1 are the best that include treatment
print(summary(egg.model2), digits = 5)
summary(allEffects(egg.model2))
(egg.model2.coef<- summary(egg.model2)$coefficients[["cond"]])

plot(allEffects(egg.model2))
# The plots correspond reasonably well with the egg plot and what we know about the data

#There was no significant effect of treatment on mean egg count across the experiment


###Check if body size effects egg number

# Filter all na values in the size column for df3
df3b <- df3 %>% 
  filter(!is.na(size))

egg.size.model <- glmmTMB(egg_tot ~ size + (1|id),data=df3b,ziformula=~1,family = nbinom1) 
summary(egg.size.model) # size on egg number
plot(allEffects(egg.size.model))

# Make the two best models (including one that included treatment) that did not include size as a fixed effect, and compare these with models that do include size as a fixed effect

egg.size.model <- glmmTMB(egg_tot ~ size + (1|id),data=df3b,ziformula=~1,family = nbinom1) 
egg.size.model2 <- glmmTMB(egg_tot ~ treat+poly(week,2) + (1|id),data=df3b,ziformula=~1, family = nbinom1)
egg.size.model3 <- glmmTMB(egg_tot ~ treat+poly(week,2) + size + (1|id),data=df3b,ziformula=~1, family = nbinom1)
egg.size.model4 <- glmmTMB(egg_tot ~ poly(week,2) + (1|id),data=df3b,ziformula=~1,family = nbinom1) 
egg.size.model5 <- glmmTMB(egg_tot ~ poly(week,2) + size + (1|id),data=df3b,ziformula=~1,family = nbinom1) 
egg.size.model6 <- glmmTMB(egg_tot ~ 1 + (1|id),data=df3b,ziformula=~1,family = nbinom1) 

(egg.selection.table2 <- as.data.frame(
  AICtab(
  egg.size.model,egg.size.model2,egg.size.model3,egg.size.model4,egg.size.model5,egg.size.model6,
  weights = T, base = T, logLik = T))) 
#From this we can see that there is very little difference between each of the models, whether or not size or treatment is included. The best performing model with this dataset is model 2. The best model that included adult female body size was model3
summary (egg.size.model3)


### 3) Lifetime egg production

#Make a new dataframe that calculates the entire egg production for each fly
#dataframe <- df3 %>% 
#  group_by(id,treat) %>% 
#  dplyr::summarise (LRS_egg_tot = sum (egg_tot),
#                    LRS_offspring_tot = sum (off_tot))

#This has already been done and the data used to make the following dataframe

LRS.df <- read_csv("NER0008751_Obj3_Exp3_LH_data_v2.csv")

LRS.df <- filter (LRS.df, exclude == "n")


#Ensure the correct number of ids are being calculated
LRS.df %>% summarise(Unique_Elements = n_distinct(no)) #correct number of ids for df

#change treat to factors so the order can be changed
LRS.df$treat<-as.factor(LRS.df$treat) 



# Use LRS.df2 for graphics as the order is more sensible to look at, use LRS.df for modeling as the order makes more sense as a control
LRS.df2 <- mutate(LRS.df, treat = relevel (treat, ref = "20% SYA", "100% SYA"))


#Make a boxplot

(egg_plot4 <- LRS.df2 %>% 
    ggplot(aes(x=treat, y=tot_eggs, fill = treat))+
    geom_boxplot(outlier.shape = NA)+
    ggbeeswarm::geom_quasirandom()+
    labs(x="Treatment", y = "Total eggs produced")+
    #scale_shape_manual(values=c(17, 16, 15), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_fill_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_y_continuous(breaks=seq(0,7000,50))+
    scale_x_discrete(limits = c("20% SYA", "100% SYA", "120% SYA"),
                     labels = c("L (20% SYA)", "M (100% SYA)", "H (120% SYA)"))+
    theme_classic(base_size = 18)+
    theme(legend.position = "none"))


ggsave("Egg plot B.png", egg_plot4, width = 21, height = 15, units = "cm")

ggsave("Egg plot B.svg", egg_plot4, width = 21, height = 15, units = "cm")


LRS.egg.model1 <- lm(tot_eggs~treat, data = LRS.df)
plot(LRS.egg.model1) #Plot does not look good, try using a glm instead

LRS.egg.model2 <- glm(tot_eggs~treat, family = poisson, data = LRS.df)
sim <- simulateResiduals(LRS.egg.model2, refit=T)
testDispersion(sim)
plot(LRS.egg.model2) #PLot looks better than above, but still not great, and the data are clearly overdispersed and zero-inflated
plot(allEffects(LRS.egg.model2))


#Use quasipoisson
LRS.egg.model3 <- glm(tot_eggs~treat, family = quasipoisson, data = LRS.df)
plot(allEffects(LRS.egg.model3)) #Shows a plot that is much closer to the actual data, but as it's a quasi-model it cannot be tested against the other models

#Use negative binomial
LRS.egg.model4 <-glm.nb(tot_eggs~treat, data = LRS.df)
sim2 <- simulateResiduals(LRS.egg.model4, refit=T)
testDispersion(sim2) #The negative binomial model is not overdispersed
plot(allEffects(LRS.egg.model4))



#Add an observation level random effect (id number will suffice for this test as each id is unique)
LRS.egg.model5 <-glmer(tot_eggs~treat + (1|no), family = poisson, data = LRS.df)
plot(LRS.egg.model5)
sim3 <- simulateResiduals(LRS.egg.model5, refit=T)
testOverdispersion(sim3) # This model is overdispersed
plot(allEffects(LRS.egg.model5))


LRS.egg.model6 <-glmmTMB(tot_eggs~treat + (1|no), family = poisson,ziformula=~1, data = LRS.df)
sim4 <- simulateResiduals(LRS.egg.model6, refit=T)
testDispersion(sim4) # Will not work as model does not implement calculations
plot(LRS.egg.model6)
plot(allEffects(LRS.egg.model6))

AICctab(LRS.egg.model1,LRS.egg.model2,LRS.egg.model3, LRS.egg.model4,LRS.egg.model5, LRS.egg.model6) #Besides model 1 (which we know does not fit the assumptions of the parametric test), model 4 seems to be the best, but model 3 cannot be evaluated. Model 3 has the closest predictions to the data.

LRS.df %>% 
  group_by(treat) %>% 
  dplyr::summarise(mean=mean(tot_eggs),
                   sd = sd(tot_eggs),
                   median = median(tot_eggs))

#The main predictions of model1, model 3 and model4 are quite similar
summary(LRS.egg.model1)

summary (LRS.egg.model3)

summary(LRS.egg.model4)
#Therefore stick to model 4 in reporting the results

#Make null model and compare with that:
LRS.egg.model4.null <-glm.nb(tot_eggs~1, data = LRS.df)

anova(LRS.egg.model4, LRS.egg.model4.null) #Test shows no effect of treatment on Lifetime egg production

#### OFFSPRING COUNTS ####

### 1) OFFSPRING COUNTS FOR THE FIRST WEEK

#Make plot
(early_off_plot<- df6 %>%
   group_by (day,id,treat,size) %>%
   dplyr::summarise(
     egg_tot = sum(egg_tot),
     off_tot = sum(off_tot),
     viability = off_tot/egg_tot 
   ) %>% 
   ggplot(aes(x=as.factor(df6$day), y=off_tot, colour=treat))+
   geom_boxplot()+
   labs(title="Early Offspring Plot",x="Treatment", y = "Offspring per day")+
   theme_classic() + 
   scale_color_brewer(type = 'div', palette = 'Set1', direction = 1))


## Test out the possible models

#The full model - note  the zero inflation parameter has been dropped. Not the data are only slightly zero-inflated so this might be alright for these models
off.early.model1 <- glmmTMB(off_tot ~ treat*poly(day,2) + (1|id), data=df6, family = nbinom1)

#No interaction
off.early.model2 <- glmmTMB(off_tot ~ treat+poly(day,2) + (1|id), data=df6, family = nbinom1)

#Just treatment
off.early.model3 <- glmmTMB(off_tot ~ treat + (1|id), data=df6,family = nbinom1)

#Just day
off.early.model4 <- glmmTMB(off_tot ~ poly(day,2) + (1|id), data=df6, family = nbinom1)

#Null
off.early.model5 <- glmmTMB(off_tot ~ 1 + (1|id), data=df6, ziformula=~1,family = nbinom1)

(off.early.selection.table <- as.data.frame(AICtab(off.early.model1,off.early.model2,off.early.model3,off.early.model4,off.early.model5, weights = T, base = T, logLik = T)))

plot(allEffects(off.early.model1))
(off.early.model1.coef<- summary(off.early.model1)$coefficients[["cond"]]) #treatments are not significant and the best model did not contain treatment as a fixed effect (however there was very little to choose between each model, except model 3)


### 2) OFFSPRING COUNTS FOR WHOLE EXPERIMENT


## Make a plot of offspring count per treatment

(off_plot <- df3%>% 
    group_by(id,treat) %>% 
    dplyr::summarise(off_mean=median(off_tot)) %>% 
    ggplot(aes(x=treat, y=off_mean))+
    geom_boxplot(outlier.shape = NA)+
    ggbeeswarm::geom_quasirandom()+
    labs(title="offspring plot",x="Treatment", y = "Mean offspring per individual")+
    theme_classic() + 
    scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)+
    scale_y_continuous(breaks=seq(0,120,20))+
    scale_x_discrete(limits = c("20% SYA", "100% SYA", "120% SYA"),
                     labels = c("L (20% SYA)", "M (100% SYA)", "H (120% SYA")))

## Make a plot of weekly offspring counts for the entire experiment (Note that 'week' is actually an 8 day period between each mating not week as is normally meant)

(off_plot2<- ggplot(df5, aes(x=week, y=off_tot, group = treat, shape=treat, colour=treat))+
    geom_errorbar(aes(ymin=off_tot-off_se, ymax=off_tot+off_se), width=.1, 
                  position=position_dodge(0.25)) +
    geom_line() +
    geom_point()+
    labs(x="Experiment period", y = "Whole-life offspring production", shape = "Treatment", colour = "Treatment")+
    theme_classic() + 
    theme(legend.position = c(0.9,0.5))+
    scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_shape_manual(values=c(17, 16, 15),labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_x_continuous(breaks=seq(0,11,1)))

#Put time instead of 'week number' to accomodate TCs request

(off_plot3<- ggplot(df5b, aes(x=last_day, y=off_tot, group = treat, shape=treat, colour=treat))+
    geom_errorbar(aes(ymin=off_tot-off_se, ymax=off_tot+off_se), width=.1, 
                  position=position_dodge(2)) +
    geom_line() +
    geom_point()+
    labs(x="Time (days)", y = "Whole-life offspring production", shape = "Treatment", colour = "Treatment")+
    theme_classic(base_size = 18) + 
    theme(legend.position = c(0.9,0.5))+
    scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_shape_manual(values=c(17, 16, 15),labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_x_continuous(breaks=seq(0,88,8)))

#Find mean value per individual per treatment to quote in the paper
df3%>% group_by(id,treat) %>% 
  dplyr::summarise(off_mean=mean(off_tot)) %>% 
  group_by(treat) %>% 
  dplyr::summarise(off_mean=mean(off_mean))


ggsave("Offspring plot A.png", off_plot3, width = 16, height = 10, units = "cm")


ggsave("Offspring plot A.svg", off_plot3, width = 21, height = 15, units = "cm")

## Test out the possible models

#Full model
off.model1 <- glmmTMB(off_tot ~ treat*poly(week,2) + (1|id),data=df3,ziformula=~1,family = nbinom1) 

#No interaction
off.model2 <- glmmTMB(off_tot ~ treat+poly(week,2) + (1|id),data=df3,ziformula=~1,family = nbinom1) 

#Just treatment
off.model3 <- glmmTMB(off_tot ~ treat + (1|id),data=df3,ziformula=~1,family = nbinom1) 

#Just week
off.model4 <- glmmTMB(off_tot ~ poly(week,2) + (1|id),data=df3,ziformula=~1,family = nbinom1) 

#Null
off.model5 <- glmmTMB(off_tot ~ 1 + (1|id),data=df3,ziformula=~1,family = nbinom1) 

(off.selection.table <- as.data.frame(AICtab(off.model1,off.model2,off.model3,off.model4,off.model5, weights = T, base = T, logLik = T)))

#off model 4 is best - but does not include treatment:
summary(off.model4)
plot(allEffects(off.model4))

#Look at the best model that does include treatment
(off.model2.coef<- summary(off.model2)$coefficients[["cond"]])
allEffects(off.model2)
plot(allEffects(off.model2))
summary(off.model2)
summary(allEffects(off.model2))

#How does this compare when thorax size is included?


###Check if body size effects offspring number

off.size.model <- glmmTMB(off_tot ~ size + (1|id),data=df3b,ziformula=~1,family = nbinom1) 
summary(off.size.model) # as with eggs, no apparent effect of body size on egg number
plot(allEffects(off.size.model))

# Make the two best models (including one that included treatment) that did not include size as a fixed effect, and compare these with models that do include size as a fixed effect

off.size.model <- glmmTMB(off_tot ~ size + (1|id),data=df3b,ziformula=~1,family = nbinom1) 
off.size.model2 <- glmmTMB(off_tot ~ treat+poly(week,2) + (1|id),data=df3b,ziformula=~1, family = nbinom1)
off.size.model3 <- glmmTMB(off_tot ~ treat+poly(week,2) + size + (1|id),data=df3b,ziformula=~1, family = nbinom1)
off.size.model4 <- glmmTMB(off_tot ~ poly(week,2) + (1|id),data=df3b,ziformula=~1,family = nbinom1) 
off.size.model5 <- glmmTMB(off_tot ~ poly(week,2) + size + (1|id),data=df3b,ziformula=~1,family = nbinom1) 
off.size.model6 <- glmmTMB(off_tot ~ 1 + (1|id),data=df3b,ziformula=~1,family = nbinom1) 

(off.selection.table2 <- as.data.frame(
  AICtab(
    off.size.model,off.size.model2,off.size.model3,off.size.model4,off.size.model5,off.size.model6,
    weights = T, base = T, logLik = T))) 
#From this we can see that there is very little difference between each of the models, whether or not size or treatment is included. The best performing model with this dataset is model 2.

#best performing model that contains body size is model 3
summary (off.size.model3)


### 3) Lifetime offspring production

#Make a new dataframe that calculates the entire egg production for each fly
#dataframe <- df3 %>% 
#  group_by(id,treat) %>% 
#  dplyr::summarise (LRS_egg_tot = sum (egg_tot),
#                    LRS_offspring_tot = sum (off_tot))

#Or you can use the previously calculated LRS.df

#Ensure the correct number of ids are being calculated
LRS.df %>% summarise(Unique_Elements = n_distinct(no)) #correct number of ids for df


#Make a boxplot

(off_plot4 <- LRS.df2 %>% 
    ggplot(aes(x=treat, y=tot_offspring ,fill = treat))+
    geom_boxplot(outlier.shape = NA)+
    ggbeeswarm::geom_quasirandom()+
    labs(x="Treatment", y = "Total offspring produced")+
    scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)+
    scale_fill_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_y_continuous(breaks=seq(0,7000,50))+
    scale_x_discrete(limits = c("20% SYA", "100% SYA", "120% SYA"),
                     labels = c("L (20% SYA)", "M (100% SYA)", "H (120% SYA)"))+
    theme_classic(base_size = 18) + 
    theme(legend.position = "none"))


ggsave("Offspring plot B.png", off_plot4, width = 16, height = 10, units = "cm")

ggsave("Offspring plot B.svg", off_plot4, width = 21, height = 15, units = "cm")

LRS.off.model1 <- lm(tot_offspring~treat, data = LRS.df)
plot(LRS.off.model1) #Plot does not look good, try using a glm instead

LRS.off.model2 <- glm(tot_offspring~treat, family = poisson, data = LRS.df)
sim <- simulateResiduals(LRS.off.model2, refit=T)
testDispersion(sim)
plot(LRS.off.model2) #PLot looks better than above, but still not great, and the data are clearly overdispersed
plot(allEffects(LRS.off.model2))


#Use quasipoisson
LRS.off.model3 <- glm(tot_offspring~treat, family = quasipoisson, data = LRS.df)
plot(allEffects(LRS.off.model3)) #Shows a plot that is much closer to the actual data, but as it's a quasi-model it cannot be tested against the other models

#Use negative binomial
LRS.off.model4 <-glm.nb(tot_offspring~treat, data = LRS.df)
sim2 <- simulateResiduals(LRS.off.model4, refit=T)
testDispersion(sim2) #The negative binomial model is not overdispersed
plot(allEffects(LRS.off.model4))



#Add an observation level random effect (id number will suffice for this test as each id is unique)
LRS.off.model5 <-glmer(tot_offspring~treat + (1|no), family = poisson, data = LRS.df)
plot(LRS.off.model5)
sim3 <- simulateResiduals(LRS.off.model5, refit=T)
testOverdispersion(sim3) # This model is overdispersed
plot(allEffects(LRS.off.model5))


LRS.off.model6 <-glmmTMB(tot_offspring~treat + (1|no), family = poisson,ziformula=~1, data = LRS.df)
sim4 <- simulateResiduals(LRS.off.model6, refit=T)
testDispersion(sim4) # Will not work as model does not implement calculations
plot(LRS.off.model6)
plot(allEffects(LRS.off.model6))

AICctab(LRS.off.model1,LRS.off.model2,LRS.off.model3, LRS.off.model4,LRS.off.model5, LRS.off.model6) #Besides model 1 (which we know does not fit the assumptions of the parametric test), model 4 seems to be the best, but model 3 cannot be evaluated. Model 3 has the closest predictions to the data.

LRS.df %>% 
  group_by(treat) %>% 
  dplyr::summarise(mean=mean(tot_offspring),
                   sd = sd(tot_offspring),
                   median = median(tot_offspring))

#The main predictions of model1, model 3 and model4 are quite similar
summary(LRS.off.model1)

summary (LRS.off.model3)

summary(LRS.off.model4)
#Therefore stick to model 4 in reporting the results

#Make null model and compare with that:
LRS.off.model4.null <-glm.nb(tot_offspring~1, data = LRS.df)

anova(LRS.off.model4, LRS.off.model4.null, test = 'Chisq') #Test shows no effect of treatment on Lifetime egg production
anova(LRS.off.model4, test = 'Chisq')
Anova(LRS.off.model4, test = 'LR')

### 3) VIABILITY FOR WHOLE EXPERIMENT

#Make a new viability column calculated from the data (ignore the one that exists in the import data as it is not accurate)
df7 <- df3 %>% 
  mutate (viability = off_tot/egg_tot)

#Check all the situations when egg count is zero to see if offspring were observed
print((df7 %>% filter(egg_tot == 0)),n=nrow(df7 %>% filter(egg_tot == 0)))
# There were 86 situations where offspring = zero was divided by eggs = zero resulting in a NaN
# There was 1 situation where offspring = 1 was divided by eggs = zero resulting in INF

#Replace NA values with zero (as zero eggs hatching into zero offspring should arguably be treated as zero viability)
df7 <- df7 %>% mutate(viability = replace_na(viability, 0))

#Check all the sitations where viability is higher than 1 (presumably because the eggs were miscounted)
print((df7 %>% filter(viability>1)),n=nrow(df7 %>% filter(viability > 1)))
#55 situations where viability is greater than 1. We can either exclude these values, or simply change them to equal 1. The latter might be justified as the egg count is usually very close, and so the viability would have been 1 or very close to it in each situation.

#Make new df where viability greater than 1 situations are removed
df8 <- df7 %>% 
  filter(viability<1)

#Make new df where viability greater than 1 situations are turned into viability = 1
df9 <- df7 %>% 
  mutate(viability = replace(viability, viability > 1, 1))

#Remake df4 and df5 to include new viability data

#the new dataframe summarising mean and sd for each week and treatment for egg total
df4c <- data_summary(df8, varname="viability",
                     groupnames=c("treat", "week"))

#Merge the two dataframes into a new table
df10 <- merge(df4,df4c,by=c("treat","week"))
setnames(df10, old="se", new= "viability_se")

#change treat to factors so the order can be changed
df10$treat<-as.factor(df10$treat) 

# Use df11 for graphics as the order is more sensible to look at
df11 <- mutate(df10, treat = relevel (treat, ref = "20% SYA", "100% SYA"))

# Add a new column that converts period to day as per TCs request
df11b <- df11 %>% 
  mutate(last_day = week*8)

## Make a viability plot showing how viability changes for each treatment over time
(via_plot<- ggplot(df11b, aes(x=last_day, y=viability, group = treat, shape=treat, colour=treat))+
    geom_errorbar(aes(ymin=viability-viability_se, ymax=viability+viability_se), width=.1, 
                  position=position_dodge(2)) +
    geom_line() +
    geom_point()+
    theme_classic() + 
    labs(x="Time (days)", y = "Egg viability", colour = "Treatment", shape = "Treatment")+
    scale_shape_manual(values=c(17, 16, 15), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_y_continuous(breaks=seq(0,1.0,0.25))+
    scale_x_continuous(breaks=seq(0,88,8)))


ggsave("Viability plot.png", via_plot, width = 16, height = 10, units = "cm")


##Viability models - once again, check v3 of this script for further model selection criteria and tests to show the final model is correct. Note df8 is used instead of df9, as when df9 is used the first model fails to converge, and it is probably better to remove confirmed miscounts rather than modifying the dataset. Note that using df8 vs df9 had no effect on the order preference of models by AIC values and no effect on the significance of the fixed effects
via.model1 <- glmer(viability ~ treat*poly(week,2) + (1|id), family=binomial, data=df8)
#Note model does not converge so should not be included as the final model
#Ignore the 'Warning message: 'In eval(family$initialize, rho) : non-integer #successes in a binomial glm!' ' as it's not a problem in this model (we are modeling the proportions themselves rather than the number of successes vs failures)

#Remove interaction
via.model2 <- glmer(viability ~ treat+poly(week,2) + (1|id), family=binomial, data=df8)

#Just treatment
via.model3 <- glmer(viability ~ treat+ (1|id), family=binomial, data=df8)

#Just time
via.model4 <- glmer(viability ~ poly(week,2) + (1|id), family=binomial, data=df8)

#Interaction only
via.model5 <- glmer(viability ~ 1 + (1|id), family=binomial, data=df8)

(via.selection.table <- as.data.frame(AICtab(via.model1,via.model2,via.model3,via.model4,via.model5, weights = T, base = T, logLik = T)))

#Note. Model 4 is considered the 'best' model. Model 2 is the best model that still includes treatment. Need to make a note of this in the results

#Best model that contains treatment
summary(via.model2)

(via.model2.coef <- summary(via.model2)$coefficients)

#No effect of treatment on egg viability

