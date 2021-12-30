################################################################
#Title: Alloparental Care in Stickleback 
#Purpose: Determine the rate of alloparental care in male stickleback in natural populations
#Created by: T. I. Ingram, D.I. Bolnick & R. E. Rangel
################################################################
##### Packages #####
library(tidyverse)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(lme4)
library(lmerTest)
library(pastecs)
library(car)
library(MCMCglmm)
library(arm)
library(beanplot)
library(PropCIs)

###Working Directory###--------------------------------------------------------
setwd("/Users/racinerangel/Desktop/UT Austin/Assort Mismatch")

#----------------------------------------------------------------------------
#RATE OF ALLOPARENTAL CARE---------------------------------------------------
male_dat<-read.csv("male_egg_mismatch_v3.csv", header=T)
str(male_dat)
#-------------------------------------------------------

#CALCULATING FREQUENCIES-------------------------------
t1 <- table(male_dat$Lake, male_dat$any_mismatch)
t1<-data.frame(match = t1[,1], mismatch = t1[,2], mismatch_freq = t1[,2] / rowSums(t1)) #original

stat.desc(t1$mismatch_freq)

#TRANSFORMING MALE MASS + NEST DEPTH--------------------
male_dat$logMass <- log(male_dat$Male_Mass_g)
male_dat$logDepth <- log(male_dat$NestDepth_cm)

hist(male_dat$Male_Zd13C)

#####RESCALE BY / 2 SD - Best practices for multiple regression especially with a mix of continuous and binary predictors------------
male_dat$Z_Male_Zd13C <- rescale(male_dat$Male_Zd13C)
male_dat$Z_logMass <- rescale(male_dat$logMass)
male_dat$Z_logDepth <- rescale(male_dat$logDepth)

#Look at the frequency in each category for the substrate/vegetation variables-------------------------
round(summary(male_dat$Substrate) / sum(summary(male_dat$Substrate)), 2)
round(summary(male_dat$VegStruct) / sum(summary(male_dat$VegStruct)), 2)
round(summary(factor(male_dat$Veg_Present)) / sum(summary(factor(male_dat$Veg_Present))), 2)
round(summary(factor(male_dat$AnyStructure))/sum(summary(factor(male_dat$AnyStructure))), 2)
table(male_dat$Lake, male_dat$Veg_Present) 

round(summary(factor(male_dat$num_mismatch))/sum(summary(factor(male_dat$num_mismatch))), 2)
table(male_dat$Lake, male_dat$num_mismatch) 


#REDUCED DATA SET FOR ANALYSIS-------------
male_dat_full <- male_dat
male_dat <- male_dat_full[complete.cases(male_dat_full[,c("any_mismatch","logDepth", "logMass", "Veg_Present", "Male_Zd13C")]),]
 apply(male_dat, 2, function(x) sum(is.na(x)))

round(summary(factor(male_dat$num_mismatch))/sum(summary(factor(male_dat$num_mismatch))), 2)
table(male_dat$Lake, male_dat$num_mismatch) 

#MODEL--------------------------
bigmodel <- glm(any_mismatch ~ Lake * Z_logMass + Lake * Z_logDepth + Lake * Z_Male_Zd13C + Lake * Veg_Present, male_dat, family = "binomial")
anova(bigmodel, test = "Chisq") #this fits terms sequentially (equivalent to Type I SS)
Anova(bigmodel) #from 'car' package = calculates marginal change in deviance for each effect, so the significance level is independent of variable order
summary(bigmodel)
AIC(bigmodel)
stepAIC(bigmodel)
#AIC = 955 
exp(coef(bigmodel))

smallmodel<-glm(any_mismatch ~ Lake + Z_Male_Zd13C + Veg_Present + Lake:Veg_Present, male_dat, family="binomial")
summary(smallmodel)
Anova(smallmodel)
summary(smallmodel)
exp(coef(smallmodel))


hist(predict.glm(bigmodel, type = "response"))
plot(predict.glm(bigmodel, type = "response") ~ male_dat$Z_Male_Zd13C)
abline(coef(lm(predict.glm(bigmodel, type = "response") ~ male_dat$Z_Male_Zd13C)))

bigmodelQ <- glm(any_mismatch ~ Lake * logMass + Lake * logDepth + Lake * Male_Zd13C + Lake * Veg_Present, male_dat, family = "quasibinomial")
Anova(bigmodelQ)
summary(bigmodelQ)#dispersion parameter is close to 1.0 indicating that the binomial family is reasonable to assume (quasi allows the variance-mean relationship to scale differently)


############################################
#################Lake Level#################
############################################
lake_dat <- read.csv("Lake_Data_v3.csv")
names(lake_dat)

hist(lake_dat$PropMismatch)
mean(lake_dat$PropMismatch)


hist(lake_dat$Surface_Area_ha)
hist(log(lake_dat$Surface_Area_ha))
hist(lake_dat$Nest_Density_invMNND)
hist(1/lake_dat$Nest_Density_invMNND) #just the mean neighbor distance in meters
hist(lake_dat$Fish_per_Transect)

lake_dat$logArea <- log(lake_dat$Surface_Area_ha)

#sculpin/trout density and lake area are correlated, since we have few observations and the transect data wasn't the strongest, probably just stick with lake area and nest density (including fish wouldn't change anything anyways)
plot(Fish_per_Transect ~ logArea, data = lake_dat)

#this is what we had before
glmLake <- glm(cbind(NumMismatch, NumMatch) ~ logArea + Nest_Density_invMNND, family = "binomial", data = lake_dat)
Anova(glmLake) 
summary(glmLake)

#however, the dispersion parameter is far from 1, indicating a poor fit for the binomial family (and the near-significance disappears in the quasi model that accounts for the dispersion)
summary(glm(cbind(NumMismatch, NumMatch) ~ logArea + Nest_Density_invMNND, family = "quasibinomial", data = lake_dat)) #dispersion doesn't fit binomial very well

#the relationship isn't very convincing at any rate
plot(PropMismatch ~ logArea, data = lake_dat)
abline(coef(lm(PropMismatch ~ logArea, data = lake_dat)))

#since the proportion mismatched is ~normal and not bounded at 0 or 1, justifiable just to use a lm
lmLake <- lm(PropMismatch ~ logArea + Nest_Density_invMNND, data = lake_dat)
Anova(lmLake)
hist(resid(lmLake))



#hard to say if watersheds differ - AdeC looks like it clusters in a narrow range compared to Campbell, while BB, M, P, and VB are all on the higher side but only 1-2 lakes each
plot(PropMismatch ~ Watershed, data = lake_dat)

lmW <- lm(PropMismatch ~ logArea + Nest_Density_invMNND + Watershed, data = lake_dat)
Anova(lmW)


#------------------------------------------------------------
#FIGURES------------------------------------------------------------
#------------------------------------------------------------

#FIGURE 1-------------------------------------------------
#Formatting Data for Figure
N_tested <- c(51,62,54,57,53,53,24,24,51,44,61,44,35,47,58)
N_mismatched <- c(14,14,27,27,35,24,12,8,18,20,34,16,18,24,32)
HA <- c(98.7, 20.6, 62.5, 52.9, 22.9, 65.6, 37.5, 4.4, 5.6, 160, 620.9, 5.7, 7.9, 369.8, 76)
name <- c("Boot", "Echo", "Gosling", "Gray", "Lawson", "Merrill", "Blackwater", "Little Mud", "Ormund", "Roberts", "Mohun", "Cranberry", "Brown's Bay", "Pye", "Village Bay")
names <- paste(name, " Lake", sep ="")


lowerCI <- c()
upperCI <- c()
for(i in 1:nrow(dat)){
  CIs <- exactci( x = N_mismatched[i], n = N_tested[i], conf.level = 0.95)
  lowerCI <- c(lowerCI, CIs[[1]][1])
  upperCI <- c(upperCI, CIs[[1]][2])
}

dat <- data.frame(names, HA, N_tested, N_mismatched, Mismatchrate = N_mismatched/N_tested, lowerCI, upperCI, logHA = log(HA))
dat <- dat[order(dat$logHA),]
dat$x <-  seq(1, nrow(dat))

#pdf("Fig1.pdf",height=6,width=7)
par(mar = c(10,5,1,1))
plot(Mismatchrate ~ x, data = dat, axes = F, xlab = "", ylab = "Frequency of alloparenting", cex.lab = 1.4, pch = 16, cex = 2, ylim = c(0,0.8))
box()
arrows(x0 = dat$x, y0 = dat$lowerCI, y1 = dat$upperCI, code = 3, length = 0.04, angle = 90)
axis(2)
axis(1, at = dat$x, labels = dat$names, las = 2)
text(x = dat$x, y = 0.0, dat$N_tested, cex = 0.8)

#dev.off()

Y <- cbind(dat$N_mismatched, dat$N_tested - dat$N_mismatched)
summary(glm(Y ~ dat$logHA, family = "binomial"))


#------------------------------------------------------------
#Figure 2-----------------------------------------------------
#------------------------------------------------------------
#fitting a model with just d13C for plotting purposes
modelC <- glm(any_mismatch ~ Male_Zd13C, male_dat, family = "binomial")
(coef(modelC))

xvalues <- seq(-2.5, 2.5, length.out = 100)
predP <- predict.glm(modelC, newdata = data.frame(Male_Zd13C = xvalues), type = "link", se.fit = TRUE)
upperP <- predP$fit + 1.96 * predP$se.fit
lowerP <- predP$fit - 1.96 * predP$se.fit

plot(predP$fit ~ xvalues, type = "l", xlab = "Male Z-d13C", ylab = "logit(Prob Mismatch)", col = "red", ylim = c(-2,2))
lines(upperP ~ xvalues, lty = "dashed", col = "red")
lines(lowerP ~ xvalues, lty = "dashed", col = "red")

pdf("Fig2.pdf",height=6,width=7)
par(mar = c(5,5,1,1))
plot(any_mismatch ~ Male_Zd13C, data = male_dat, ylim = c(0, 1), pch = "|", cex = 0.5, xlim = c(-4,4), ylab = "Mismatch Probability", xlab = expression(paste("Male Z-",delta^13,"C",sep="")))
lines(modelC$family$linkinv(predP$fit) ~ xvalues, col = "black")
lines(modelC$family$linkinv(upperP) ~ xvalues, col = "black", lty = "dashed")
lines(modelC$family$linkinv(lowerP) ~ xvalues, col = "black", lty = "dashed")

beanplot(male_dat$Male_Zd13C[male_dat$any_mismatch == FALSE], add = TRUE, horizontal = TRUE, side = "second", log = FALSE, at = 0.05, boxwex = 0.3, what = c(0,1,0,0))
beanplot(male_dat$Male_Zd13C[male_dat$any_mismatch == TRUE], add = TRUE, horizontal = TRUE, side = "first", log = FALSE, at = 0.95, boxwex = 0.3, what = c(0,1,0,0))

dev.off()

#------------------------------------------------------------
#FIGURE 3----------------------------------------------------
#------------------------------------------------------------
#simple plot of mismatch vs vegetation in each lake (also displays the among-lake variation)
probs <- data.frame(all = tapply(male_dat$any_mismatch, male_dat$Lake, mean), veg = tapply(male_dat$any_mismatch[male_dat$Veg_Present == 1], male_dat$Lake[male_dat$Veg_Present == 1], mean), noveg = tapply(male_dat$any_mismatch[male_dat$Veg_Present == 0], male_dat$Lake[male_dat$Veg_Present == 0], mean))
probs <- probs[order(probs[,2], decreasing = TRUE),]

pdf("Fig3.pdf",height=6,width=7)
plot(NA, xlim = c(-0.5,1.5), xaxt = "n", xlab = "Vegetation", ylim = c(0,1), ylab = "Mismatch Probability")
axis(1, at = c(0,1), labels = c("Absent", "Present"))
segments(rep(0,15), probs[,3], rep(1,15), probs[,2])
text(rep(1.1,15), probs[,2] + c(0,0.01,0,0.01,0.03,0.03,0.01,0,0,0,0,0,0,-0.01,-0.02), rownames(probs), adj = 0, cex = 0.5)

dev.off()


#Clutch descriptive data------------------------------------------------------------------------------------
#Male count
# After filtering above for full cases
length(male_dat$FishID)
#706 Males

male_dat %>%
  count(Num_clutches)

male_dat %>%
  count(num_mismatch)
 

#Clutch 1 Data-------------
male_dat %>%
  count(clutch1_egg1)
#clutch1_mismatch   n
#1            FALSE 399
#2             TRUE 307

male_dat %>%
  count(clutch1_egg2)


male_dat %>%
  count(clutch2_egg1)

male_dat %>%
  count(clutch3_egg2)

#Total = 1418 --> 31.7% proportion eggs are mistmatched to male + 68.3% are matched to male

