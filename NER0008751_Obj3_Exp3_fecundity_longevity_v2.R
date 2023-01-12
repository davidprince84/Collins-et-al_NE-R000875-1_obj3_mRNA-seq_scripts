#DavidHCollins
#Created: 25 February 2020

#Aim

#The aim of this script is to analyse data from our main Drosophila experiment (please see NER0008751_Obj3_Exp3_Datasheet_v7_2018-10-18.xlsx for the original raw data and Protocol_NER0008751_Objective 3_Experiment 3_v1_2018-05-05 for the experimental protocol). This is a shortened version of 'NER0008751_Obj3_Exp3_fecundity_longevity.R' which contains many more analyses. This version should be used for the paper, and the figures and analysis have been re-ordered to accomodate that.

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
library(ggplot2)

#Might be usefult to draw a survival curve for the longevity data later on, install the ggplot package fortify which has suvival curve functions
#(if on a new computer use:
#install.packages('ggfortify')

library(ggfortify)
library(survival)
library("survminer")
library(readr)

#here are some useful shortcuts to get to know
#ctrl+a highlights the whole script
#ctrl+enter runs the highlighted script or just the selected line of code if nothing is highlighted
#ctrl+1 moves cursor to the console
#ctrl+2 moves it back to the script editor


#Tell R where to look
setwd("~/Documents/UEA Work documents/Post Doc_LHT/Objective 3 Experiment 3/Analysis/")


#Load in LH data

df <- read_csv("NER0008751_Obj3_Exp3_LH_data_v2.csv")
View(df)

# Note this excel file is a csv version of the raw LH data in NER0008751_Obj3_Exp3_Datasheet_v6_2018-10-09.xlsx



# Later on we will use a different excel file that contains the daily offspring and egg counts for each D. melanogaster female (we don't need this file just yet!)
library(readr)
NER0008751_Obj3_Exp3_dailyoff <- read_csv("NER0008751_Obj3_Exp3_dailyoff_3.csv")
View(NER0008751_Obj3_Exp3_dailyoff)


#due to a high death rate at the start of the experiment and the occasional accidental death throughout the experiment, some of the data had to be censored. Therefore, we need to create a new data frame looking at just the individuals that were NOT censored (i.e. were classified as 'd' in the censor column)
df2 <- filter (df, censor == "1")

#also need to re-order the treatment names so they follow a logical order

#1st Convert the treatment variable into a factor
df2$treat<-as.factor(df2$treat) 
is.factor(df2$treat)

#create a new dataframe where T1 is first for figures (we will continue to use the original dataframe for models)
df2b <- mutate(df2, treat = relevel (treat, ref = "20% SYA", "100% SYA"))

#Check that the order is now correct
levels(df2$treat)
levels(df2b$treat)


###################################################################################################################

###################################################################################################################


                                  ###########1.0 Effect of Treatment on Longevity and Fecundity############


# 1.2.1 boxplot: effect of treatment on total fecundity 
bp.tot_fec <- ggplot(df2b,aes(x=treat, y=tot_eggs)) +
  geom_boxplot()+
  labs(title="The effect of larval dietary regime on total egg production in Drosophila melanogaster females",x="Larval diet", y = "Total number of eggs")+
  theme_classic()

bp.tot_fec

# 1.2.2 ANOVA: effect of treatment on total fecundity

#Build model
lm.tot_fec <- lm(tot_eggs ~ treat, data = df2)

#test assumptions
autoplot(lm.tot_fec , smooth.colour = NA) #these data aren't perfect but will do for now

# use anova() to test if there is an effect of Treatment
anova(lm.tot_fec )

#use summary() to test what the nature of the effect of treatment is, and see if it makes sense with the boxplot itself
summary(lm.tot_fec )

# to summarize there was no difference in the amount of eggs produced between the three treatments

----------------------------------------------------------------------------------------------------------------

  # 1.3.1 boxplot: effect of treatment on total offspring number 

(bp.tot_off <- ggplot(df2b,aes(x=treat, y=tot_offspring)) +
  geom_boxplot()+
  labs(title="The effect of larval dietary regime on total offspring production in Drosophila melanogaster females",x="Larval diet", y = "Total number of offspring")+
  theme_classic())



# 1.3.2 ANOVA: effect of treatment on offspring number 

#Build model
lm.tot_Off <- lm(tot_offspring ~ treat, data = df2)

#test assumptions
autoplot(lm.tot_Off , smooth.colour = NA) #these data aren't perfect but will do for now

# use anova() to test if there is an effect of Treatment
anova(lm.tot_Off)

#use summary() to test what the nature of the effect of treatment is, and see if it makes sense with the boxplot itself
summary(lm.tot_Off )

# to summarize there was no effect of treatment on total number of offspring produced

  
# 1.6.1 boxplot: effect of treatment on average thorax length
bp.avg_th <- ggplot(df2b,aes(x=treat, y=thorax)) +
  geom_boxplot()+
  labs(title="The effect of larval dietary regime on female thorax size in Drosophila melanogaster 
females",x="Larval diet", y = "Thorax length (mm)")+
  theme_classic()
bp.avg_th

# 1.6.2 ANOVA effect of treatment on average thorax length

#Build model
lm.avg_th <- lm(thorax ~ treat, data = df2)

#test assumptions
autoplot(lm.avg_th , smooth.colour = NA) 

# use anova() to test if there is an effect of Treatment
anova(lm.avg_th)

#use summary() to test what the nature of the effect of treatment is, and see if it makes sense with the boxplot itself
summary(lm.avg_th)

# to summarize there was a significant effect of treatment on thorax size in Drosophila melanogaster, adults from treatment 1 were smaller than adults from treatment 2 and 3. Adults from treatments 2 and 3 were not significantly different from each other.

# 1.6.3 ANOVA effect of treatment on average thorax length

#Build model
lm.avg_th2 <- lm(thorax ~ treat, data = df2)

#test assumptions
autoplot(lm.avg_th2 , smooth.colour = NA) 

# use anova() to test if there is an effect of Treatment
anova(lm.avg_th2)

#use summary() to test what the nature of the effect of treatment is, and see if it makes sense with the boxplot itself
summary(lm.avg_th2)


----------------------------------------------------------------------------------------------------------------
  
# 1.7.1 boxplot: effect of treatment on average wing cell length


#boxplot
bp.avg_wi <- ggplot(df2b,aes(x=treat, y=wing)) +
  geom_boxplot()+
  labs(title="The effect of larval dietary regime on female marginal cell length in Drosophila melanogaster 
       females",x="Larval diet", y = "Marginal cell length (mm)")+
  theme_classic()

bp.avg_wi


# 1.7.2 ANOVA effect of treatment on average wing length

#Build model
lm.avg_wi <- lm(wing ~ treat, data = df2)

#test assumptions
autoplot(lm.avg_wi , smooth.colour = NA) 

# use anova() to test if there is an effect of Treatment
anova(lm.avg_wi)

#use summary() to test what the nature of the effect of treatment is, and see if it makes sense with the boxplot itself
summary(lm.avg_wi)

# to summarize there was a significant effect of treatment on wing size in Drosophila melanogaster, adults from treatment 1 were smaller than adults from treatment 2 and 3. Adults from treatments 2 and 3 were not significantly different from each other.


##################################################################################################################

                                      ###########2. 0 Longevity-Fecundity relationships############



# First a Scatterplot to see the general pattern of longevity and total/mean offspring across all treatments

(sp.all.tot_fec<- ggplot(df2b,aes(x=age, y=tot_eggs)) +
  geom_point()+
  labs(title="The relationship between longevity and total egg production in Drosophila melanogaster females",x="Longevity (days)", y = "Total number of eggs")+
  theme_classic())

(sp.all.tot_off<- ggplot(df2b,aes(x=age, y=tot_offspring)) +
  geom_point() +
  labs(title="The relationship between longevity and total offpring production in Drosophila melanogaster females",x="Longevity (days)", y = "Total number of offspring")+
  theme_classic())


#Or plot each treatment seperately

(sp.ind.tot_fec<- ggplot(df2b,aes(x=age, y=tot_eggs, colour = treat)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="The relationship between longevity and total egg production in Drosophila melanogaster females reared
under three larval dietary regimes",x="Longevity (days)", y = "Total number of eggs", colour = 'Larval diet') +
  theme_classic())

(sp.ind.tot_off<- ggplot(df2b,aes(x=age, y=tot_offspring, colour = treat)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="The relationship between longevity and total offspring production in Drosophila melanogaster females
reared under three larval dietary regimes",x="Longevity (days)", y = "Total number of offspring", colour='Larval diet') +
  theme_classic())




# Second a scatterplot to see the general pattern of longevity and average offspring across all treatments


#To make the plots - decide a colour scheme
library(colorBlindness)

#Make cvdPlot for known colourblind friendly colours
cvdPlot(displayColors(c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))) #My preferred colours are: lightblue - "#56B4E9, darkblue - "#0072B2", and redorange - "#D55E00". Add these to the next plot


##### THESE ARE THE MOST IMPORTANT PLOTS AND ANALYSES FOR THESE DATA ########


#For colourblind testing - see plot below
# (cvdPlot(sp.ind.avg_fec<- ggplot(df2b,aes(x=age, y=avg_eggs, colour = treat)) +
#  geom_point() +
#  geom_smooth(method="lm") +
#  labs(title="The relationship between longevity and average egg offspring production in Drosophila melanogaster
# females reared under three larval dietary regimes",x="Longevity (days)", y = "Mean egg production", colour='Larval diet')+
#   scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
#  theme_classic()))


(sp.ind.avg_fec<- ggplot(df2b,aes(x=age, y=avg_eggs, colour = treat)) +
    geom_point() +
    geom_smooth(method="lm") +
    labs(title="The relationship between longevity and average egg offspring production in Drosophila melanogaster
females reared under three larval dietary regimes",x="Longevity (days)", y = "Mean egg production", colour='Larval diet')+
   scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    theme_classic())

(sp.ind.avg_off<- ggplot(df2b,aes(x=age, y=avg_offspring, colour = treat)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="The relationship between longevity and average offspring production in Drosophila 
melanogaster females reared under three larval dietary regimes",x="Longevity (days)", y = "Mean offspring production", colour='Larval diet')+
    scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    theme_classic())


# Test relationship between each linear model

egg.model1 <- lm(avg_eggs ~ treat*age, data = df2)

# Test if thorax should be included as an additional fixed effect

#Filter out all the females that did not have the have their thorax measured
df3 <- df2 %>% 
  filter(!is.na(thorax))

#Makes models with and without thorax
egg.model2 <- lm(avg_eggs ~ treat*age, data = df3)
egg.model3 <- lm(avg_eggs ~ treat*age*thorax, data = df3)
egg.model4 <- lm(avg_eggs ~ treat*age + thorax, data = df3)

#Compare AIC values
(egg.selection.table <- as.data.frame(
  AICtab(
    egg.model2,egg.model3, egg.model4,
    weights = T, base = T, logLik = T))) 

anova(egg.model2, egg.model4) #No difference in likelihood ratio between models 2 and 4

#As model 2 was better than models 3 and 4 (or at least performed no worse), stick to model 1 (which has more data as it includes individuals with no thorax measurements)

autoplot(egg.model1) # This plot shows 1) TL no humps or valleys so the line is an appropriate fit to the data; 2) TR most of the residuals fall on the expected pattern if there is normal distribution; 3) BL variance is mostly constant over all predicted values of the response variable, so assumption of equal variance is probably not violated, there is a possibility of a negative realtionship between the fitted values and the standardized residuals but very weak if there at all; 4) BR shows that a single point might be having a disproportionate influence on the patterns... perhaps this is worth highlighting in the main figure. That point has a high residual value (i.e. it is a big outlier). In addition it has a large leverage value (i.e. it has a large affect on the model that was fitted). Either of these values alone would not be a major cause for concern, but both together are a problem. The final figure on plot shows that that value has a cook's distance of greater than 1 which gives us some justification for excluding it, so long as we tailor our discussion and conclusions to both scenarios where that point is removed and not removed.

anova(egg.model1) # no overall effect of age or treatment on the average egg production, but there is a considerable interaction between the two variables

summary(egg.model1) # treatment 3 has a reversed longevity-fecundity relationship compared to the other treatments

# Test effects for offspring data
offspring.model1 <- lm(avg_offspring ~ treat*age, data = df2)

autoplot(offspring.model1) # This plot shows 1)TL no humps or valleys so the line is an appropriate fit to the data; 2) TR most of the residuals fall on the expected pattern if there is normal distribution; 3) BL variance is mostly constant over all predicted values of the response variable, so assumption of equal variance is probably not violated, there is a possibility of a negative realtionship between the fitted values and the standardized residuals but very weak if there at all; 4) BR shows that a single point might be having a disproportionate influence on the patterns... perhaps this is worth highlighting in the main figure

anova(offspring.model1) # no overall effect of age or treatment on the average offspring production, but there is a considerable interaction between the two variables

summary(offspring.model1) # treatment 3 has a reversed longevity-fecundity relationship compared to the other treatments


#Add the egg model to the figure

new.x<-expand.grid(age= seq (from = 7, to = 88, length.out = 10),
                    treat = levels(df2b$treat))

#use predict() to generate new y values
new.y <- predict(egg.model1,newdata=new.x,interval = 'confidence')

#Collect the new x and new y for plotting
addThese<-data.frame(new.x,new.y)
#change name from fit to match the original data (in this case EGGS)
addThese<-dplyr::rename(addThese, avg_eggs = fit)
#check it worked
head(addThese)

(egg.relationship<- ggplot(df2b,aes(x=age, y=avg_eggs, colour = treat)) +
    geom_point(aes(shape = treat)) +
    geom_smooth(data = addThese,
                aes(ymin = lwr, ymax = upr, fill = treat, colour = treat),
                stat = 'identity',show.legend = FALSE)+
                labs(x="Longevity (days)", y = "Mean egg production (per day)", colour = "Treatment", shape = "Treatment", fill = "Treatment")+
    scale_shape_manual(values=c(17, 16, 15), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_fill_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    theme_classic())


#Add the offspring model to the figure

#use predict() to generate new y values
new.y2 <- predict(offspring.model1,newdata=new.x,interval = 'confidence')

#Collect the new x and new y for plotting
addThese2<-data.frame(new.x,new.y2)
#change name from fit to match the original data (in this case EGGS)
addThese2<-dplyr::rename(addThese2, avg_offspring = fit)
#check it worked
head(addThese)

(offspring.relationship<- ggplot(df2b,aes(x=age, y=avg_offspring, colour = treat)) +
    geom_point(aes(shape = treat)) +
    geom_smooth(data = addThese2,
                aes(ymin = lwr, ymax = upr, fill = treat, colour = treat),
                stat = 'identity',show.legend = FALSE)+
    labs(x="Longevity (days)", y = "Mean offspring production (per day)", colour='Treatment', shape = "Treatment", fill = "Treatment")+
    scale_shape_manual(values=c(17, 16, 15), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_color_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    scale_fill_manual(values = c("#D55E00", "#0072B2", "#56B4E9"), labels = c("L (20% SYA)","M (100% SYA)", "H (120% SYA)"))+
    theme_classic())

ggsave("Egg relationship.png", egg.relationship, width = 16, height = 10, units = "cm")
ggsave("Offspring relationship.png", offspring.relationship, width = 16, height = 10, units = "cm")

#Determine the r2 values for each linear relationship

##Treatment 1
Treat1 <- df2 %>% 
  filter(treat_no == 1)

#Eggs
egg.T1 <- lm(avg_eggs ~ age, data = Treat1)
summary(egg.T1)

#Offspring
offspring.T1<- lm(avg_offspring ~age, data = Treat1)
summary(offspring.T1)

##Treatment 2
Treat2 <- df2 %>% 
  filter(treat_no == 2)

#Eggs
egg.T2 <- lm(avg_eggs ~age, data = Treat2)
summary(egg.T2,digits = 5)

#Offspring
offspring.T2<- lm(avg_offspring ~age, data = Treat2)
summary(offspring.T2)

##Treatment 3

Treat3 <- df2 %>% 
  filter(treat_no == 3)

#Eggs

egg.T3 <- lm(avg_eggs ~age, data = Treat3)
summary(egg.T3,digits = 5)


#Offspring
offspring.T3<- lm(avg_offspring ~age, data = Treat3)
summary(offspring.T3)

####### Now test once again when the outliers are removed (i.e. day 55 onwards
#Make new dataframe (with a simple name)

df2r <- filter(df2, age>= 55)

#Simplify names

(egg.relationship2 <- ggplot(df2r,aes(x=age, y=avg_eggs, colour = treat)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="The relationship between longevity and average egg offspring production in Drosophila melanogaster
       females reared under two larval dietary regimes after individuals below day 50 are removed",x="Longevity (days)", y = "Average number of eggs laid per day alive", colour='Larval diet')+
  theme_classic())

offspring.relationship2 <- ggplot(df2r,aes(x=age, y=avg_offspring, colour = treat)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(title="The relationship between longevity and average offspring production in Drosophila 
       melanogaster females reared under three larval dietary regimes after individuals below day
       50 are removed",x="Longevity (days)", y = "Average number of offspring produced per day alive", colour='Larval diet')+
  theme_classic()

sp2.3r.egg
sp2.3r.off


#linear model for egg production
egg.model2 <- lm(avg_eggs ~ treat*age, data = df2r)
autoplot(egg.model2) # outliers are now removed
anova(egg.model2) # Still a (much statistically weaker) effect of treatment interacting with age to affect fecundity
summary(egg.model2) # signficiant difference between treatment 2 and treatment 3 (but not treatment 1 and treatment 2)

#linear model for offspring production
offspring.model2 <- lm(avg_offspring ~ treat*age, data = df2r)
autoplot(offspring.model2) # outliers are now removed
anova(offspring.model2) # no significant effect of age or treatment on the offspring egg production, or on the interaction between them (but it is close to significance)
summary(offspring.model2) # no effect of longevity on average daily offspring production, but some possible evidence of an effect of treatment on longevity-offspring relationships when individuals that died before day 50 are excluded.


#Check individual relationships
Treat2r <- df2r %>% 
  filter(treat_no == 2)

egg.T2r <- lm(avg_eggs ~age, data = Treat2r)
summary(egg.T2r,digits = 5)
offspring.T2r<- lm(avg_offspring ~age, data = Treat2r)
summary(offspring.T2r)

#Check individual relationships
Treat1r <- df2r %>% 
  filter(treat_no == 1)

egg.T1r <- lm(avg_eggs ~age, data = Treat1r)
summary(egg.T1r,digits = 5)
offspring.T1r<- lm(avg_offspring ~age, data = Treat1r)
summary(offspring.T1r)


#Check individual relationships
Treat2r <- df2r %>% 
  filter(treat_no == 2)

egg.T2r <- lm(avg_eggs ~age, data = Treat2r)
summary(egg.T2r,digits = 5)
offspring.T2r<- lm(avg_offspring ~age, data = Treat2r)
summary(offspring.T2r)


#Check individual relationships
Treat3r <- df2r %>% 
  filter(treat_no == 3)

egg.T3r <- lm(avg_eggs ~age, data = Treat3r)
summary(egg.T3r,digits = 5)
offspring.T3r<- lm(avg_offspring ~age, data = Treat3r)
summary(offspring.T3r)


#Overall the relationships still hold (though are quite a bit weaker) when only the individuals that reached day 50 are included in the models

