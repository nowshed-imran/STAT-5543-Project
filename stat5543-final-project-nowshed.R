
#####################################################################
######## What explains systolic blood pressure in US adults? ########
#####################################################################

# STAT 5543 Applied Regression Analysis, final project
# Spring 2025, revised October 2026
# Nowshed Hossain Imran

# Data: NHANES 2011 to 2012, adults (cardioStatsUSA nhanes_data)
# Response: systolic blood pressure (mean of up to three readings)

# Parts 1 to 6 follow the steps used in class: read the data,
# descriptive stats, initial fit and transformation, model,
# diagnostics, testing and prediction.
# Part 7 is the extension: assumption checks and fixes from
# econometrics (robust standard errors, FGLS, omitted variable bias,
# median regression).

# Code credit and AI assistance: Parts 1 to 6 (the main study) follow
# the code and conventions from the STAT 5543 class scripts by
# Dr. Pratyaydipta Rudra, Oklahoma State University. For Part 7 (the
# extension), I took help from Claude (Anthropic) to write the code. The
# study design, hypotheses, model choices and interpretation are my own.
# I have reviewed the code and would appreciate any suggestions for
# improvement.

# Reference: Wooldridge, J. M. (2020). Introductory econometrics:
# A modern approach (7th ed.). Cengage.

rm(list=ls())

#### Setting a working directory
setwd("path/to/STAT-5543-Project")
#This will be different for you depending on where you have kept the project folder

library("car") #vif, linearHypothesis, outlierTest
library("MASS") #studres
library("lmtest") #extension: resettest, bptest, coeftest
library("sandwich") #extension: robust standard errors
library("quantreg") #extension: median regression

dir.create("output", showWarnings = F)
sink("output/results.txt", split = T)
#everything printed below is also saved in output/results.txt




####################### 1. Reading the data #######################
####################################################################

dat <- read.csv("data/nhanes.csv", header = T)
dat <- dat[which(dat$svy_year == "2011-2012"), ]
nrow(dat) #5,615 adults

#Is there any missing value?
sum(is.na(dat$chol_ldl))
sum(is.na(dat$chol_trig))
#LDL and triglycerides are only measured in the morning fasting subsample


## LDL is calculated with the Friedewald formula, not measured
ok <- which(!is.na(dat$chol_ldl) & !is.na(dat$chol_total))
friedewald <- dat$chol_total[ok] - dat$chol_hdl[ok] - dat$chol_trig[ok]/5
cor(friedewald, dat$chol_ldl[ok])
#0.998, so the formula is used. It is not valid when triglycerides >= 400


## Choose the required columns and remove the missing values (NAs)
subdat <- dat[, c("bp_sys_mean", "demo_gender", "demo_pregnant", "demo_race",
                  "bp_med_use", "cc_smoke", "cc_diabetes", "chol_hdl",
                  "chol_ldl", "chol_trig", "demo_age_years")]
keep <- which(complete.cases(subdat)==T)
compdat <- cbind(subdat[keep, ], cc_bmi = dat$cc_bmi[keep])
nrow(compdat) #2,187

sum(compdat$demo_pregnant == "Yes") #19 pregnant women
sum(compdat$demo_pregnant == "No" & compdat$chol_trig >= 400) #42
compdat <- subset(compdat, demo_pregnant == "No" & chol_trig < 400)
n <- nrow(compdat)
n #2,126 adults in the analysis


## Variables for the model
compdat$sbp <- compdat$bp_sys_mean
compdat$women <- as.numeric(compdat$demo_gender == "Women")
compdat$race <- relevel(factor(compdat$demo_race), ref = "Non-Hispanic White")
compdat$smoke <- relevel(factor(compdat$cc_smoke), ref = "Never")
compdat$diabetes <- as.numeric(compdat$cc_diabetes == "Yes")
compdat$bpmed <- as.numeric(compdat$bp_med_use == "Yes")
compdat$logtrig <- log(compdat$chol_trig)
compdat$age <- compdat$demo_age_years - 50 #centered at 50
compdat$hdl <- (compdat$chol_hdl - 50)/10 #10 mg/dL units, centered at 50
compdat$ldl <- (compdat$chol_ldl - 110)/10 #10 mg/dL units, centered at 110
#Centering makes the women coefficient the gap at age 50, not at age 0




####################### 2. Descriptive stats #######################
####################################################################

skew <- function(x) mean((x - mean(x))^3)/sd(x)^3 #sample skewness

summary(compdat$sbp)
sd(compdat$sbp)
hist(compdat$sbp)
boxplot(compdat$sbp, horizontal = T)
skew(compdat$sbp) #1.07, response is positively skewed
skew(log(compdat$sbp)) #0.48, much better after log

## Mean BP by group
tapply(compdat$sbp, compdat$demo_gender, mean)
tapply(compdat$sbp, compdat$race, mean)
tapply(compdat$sbp, compdat$smoke, mean)
tapply(compdat$sbp, compdat$cc_diabetes, mean)
tapply(compdat$sbp, compdat$bp_med_use, mean)
table(compdat$demo_gender)

boxplot(sbp ~ demo_gender, data = compdat)
boxplot(sbp ~ race, data = compdat)
plot(compdat$demo_age_years, compdat$sbp, pch = 16, cex = 0.5)
#BP rises with age, and faster at older ages




#################### 3. Initial fit and transformation ####################
###########################################################################

## Blood pressure medication is left out of the model on purpose.
## People take it because their BP is high, and it then lowers BP.
## Holding it fixed would bias the other coefficients (a bad control).

fitlevel <- lm(sbp ~ age + I(age^2) + women + age*women + race + smoke + diabetes
               + hdl + ldl + ldl*women + logtrig, data = compdat)

fit <- lm(log(sbp) ~ age + I(age^2) + women + age*women + race + smoke + diabetes
          + hdl + ldl + ldl*women + logtrig, data = compdat)

par(mfrow = c(1,2))
qqnorm(fitlevel$residuals, main = "SBP")
qqline(fitlevel$residuals)
qqnorm(fit$residuals, main = "log(SBP)")
qqline(fit$residuals)
par(mfrow = c(1,1))

skew(fitlevel$residuals) #0.82
skew(fit$residuals) #0.26, the log fixes most of the skew
#Also with log(SBP), each coefficient reads as a percent change
#So we use log(SBP) as response




###################### 4. The model ######################
##########################################################

## Should men and women share one model? (Chow test)
## Same idea as testing whether two lines are identical
fit_pooled <- lm(log(sbp) ~ age + I(age^2) + race + smoke + diabetes
                 + hdl + ldl + logtrig, data = compdat)
fit_full <- lm(log(sbp) ~ (age + I(age^2) + race + smoke + diabetes
                           + hdl + ldl + logtrig)*women, data = compdat)
anova(fit_pooled, fit_full)
#F = 9.50, p < 0.001: men and women do not share one model

## Is it only the age slope that differs?
fit_age <- lm(log(sbp) ~ (age + I(age^2))*women + race + smoke + diabetes
              + hdl + ldl + logtrig, data = compdat)
anova(fit_age, fit_full)
#p = 0.002: something else differs too

## Which interactions matter? (robust p-values)
out <- coeftest(fit_full, vcov. = vcovHC(fit_full, type = "HC1"))
out[grep(":women", rownames(out)), ]
#Only age:women and ldl:women. ldl:women p = 0.0001 < 0.05/12 (Bonferroni)
summary(fit_age)$coefficients["I(age^2):women", ]
#p = 0.97, the curvature of age does not differ by gender

fit_age_ldl <- lm(log(sbp) ~ (age + I(age^2))*women + race + smoke + diabetes
                  + hdl + ldl + ldl*women + logtrig, data = compdat)
anova(fit_age_ldl, fit_full)
#p = 0.29: the other ten interactions add nothing


## Final model (same as fit above)
summary(fit)
#R2 = 0.276




###################### 5. Diagnostics ######################
############################################################

#Residuals vs fitted values
plot(fit$fitted.values, studres(fit), pch = 16, cex = 0.5, col = "darkgreen")
abline(h = 0)
#no clear pattern

#QQ-plot
qqnorm(fit$residuals, pch = 16, cex = 0.5)
qqline(fit$residuals, col = "red")
#slightly heavy right tail, but n = 2,126 is large, so normality is not a big concern

#Multicollinearity
vif(fit, type = "predictor")
#generalized VIFs; the adjusted values (last column) are 1.02 to 1.22
#nothing to worry about

#Outliers and influential observations
outlierTest(fit)
#largest studentized residual 4.57, Bonferroni p = 0.011
max(cooks.distance(fit)) #0.016, far below 1
plot(cooks.distance(fit), typ = "b")
which.max(cooks.distance(fit))

#Drop the 20 most influential observations and refit
top <- order(cooks.distance(fit), decreasing = T)[1:20]
fit_drop <- update(fit, data = compdat[-top, ])




################# 6. Testing and prediction #################
#############################################################

summary(fit)$coefficients

## Percent effects: 100*(exp(beta) - 1)
round(100*(exp(fit$coefficients) - 1), 2)

## Effect of 10% higher triglycerides
100*(exp(fit$coefficients["logtrig"]*log(1.1)) - 1) #0.25%

## Age effect per 10 years at a given age
b <- fit$coefficients
100*(exp(10*(b["age"] + 2*b["I(age^2)"]*(40 - 50))) - 1) #men at 40
100*(exp(10*(b["age"] + 2*b["I(age^2)"]*(60 - 50))) - 1) #men at 60
100*(exp(10*(b["age"] + 2*b["I(age^2)"]*(40 - 50) + b["age:women"])) - 1) #women at 40
100*(exp(10*(b["age"] + 2*b["I(age^2)"]*(60 - 50) + b["age:women"])) - 1) #women at 60

## Age where the gender gap closes (LDL at 110)
50 - b["women"]/b["age:women"] #about 69


## Testing groups of terms with the extra SS F-test
fit_norace <- update(fit, . ~ . - race)
anova(fit_norace, fit)

fit_nosmoke <- update(fit, . ~ . - smoke)
anova(fit_nosmoke, fit)

fit_nolipid <- lm(log(sbp) ~ age + I(age^2) + women + age*women + race + smoke
                  + diabetes, data = compdat)
anova(fit_nolipid, fit)
#These use the usual standard errors. Part 7 shows the errors are
#heteroskedastic, so the robust versions of these tests are in Part 7


## Prediction
newdat <- expand.grid(demo_age_years = 20:80, women = c(0, 1))
newdat$age <- newdat$demo_age_years - 50
newdat$race <- factor("Non-Hispanic White", levels = levels(compdat$race))
newdat$smoke <- factor("Never", levels = levels(compdat$smoke))
newdat$diabetes <- 0
newdat$hdl <- (median(compdat$chol_hdl) - 50)/10
newdat$ldl <- (median(compdat$chol_ldl) - 110)/10
newdat$logtrig <- median(compdat$logtrig)
newdat$pred <- exp(predict(fit, newdat))
newdat[which(newdat$demo_age_years %in% c(20, 40, 60, 80)), c("demo_age_years", "women", "pred")]

#Predicted BP and 95% prediction interval for a 60 year old woman
#with median values for the other variables
exp(predict(fit, newdat[which(newdat$demo_age_years == 60 & newdat$women == 1), ],
            interval = "prediction"))




#################### 7. Extension: econometric checks ####################
##########################################################################

#### MLR.1 Functional form ####
resettest(fit, power = 2:3)
#p = 0.34, no sign of missing nonlinear terms

fit_noage2 <- update(fit, . ~ . - I(age^2))
resettest(fit_noage2, power = 2:3)
#p = 0.02 without age squared, so the quadratic is needed


#### MLR.2 Random sampling: is the fasting sample different? ####
dat$in_sample <- rownames(dat) %in% rownames(compdat)
tapply(dat$bp_sys_mean, dat$in_sample, mean, na.rm = T)
tapply(dat$demo_age_years, dat$in_sample, mean)
tapply(dat$demo_gender == "Women", dat$in_sample, mean)
tapply(dat$bp_med_use == "Yes", dat$in_sample, mean, na.rm = T)
table(dat$in_sample)
#similar BP and gender mix, sample is 2.5 years older


#### MLR.4 Zero conditional mean ####

## (a) Bad control: what happens if medication is added?
fit_med <- update(fit, . ~ . + bpmed)
summary(fit_med)$coefficients["bpmed", ]
fit_untreated <- update(fit, data = compdat[which(compdat$bpmed == 0), ])
#diabetes effect: 3.75% (main), 3.10% (with medication), 4.73% (untreated)

## (b) Omitted variable: obesity (BMI 30 or more)
## Omitted variable bias: short coefficient = long coefficient + gamma*delta
## gamma: obesity coefficient in the long regression
## delta: coefficient from regressing obesity on the included variables
bmidat <- compdat[which(!is.na(compdat$cc_bmi)), ]
bmidat$obese <- as.numeric(bmidat$cc_bmi %in% c("30 to <35", "35+"))
nrow(bmidat) #2,088
mean(bmidat$obese)

fit_short <- update(fit, data = bmidat)
fit_long <- update(fit, . ~ . + obese, data = bmidat)
fit_aux <- lm(obese ~ age + I(age^2) + women + age*women + race + smoke + diabetes
              + hdl + ldl + ldl*women + logtrig, data = bmidat)

k <- c("diabetes", "logtrig", "hdl", "raceNon-Hispanic Black")
gamma <- fit_long$coefficients["obese"]
ovb <- cbind(short = fit_short$coefficients[k], long = fit_long$coefficients[k],
             delta = fit_aux$coefficients[k],
             predicted_bias = gamma*fit_aux$coefficients[k],
             actual_change = fit_short$coefficients[k] - fit_long$coefficients[k])
round(ovb, 5)
#predicted bias = actual change, exactly
coeftest(fit_long, vcov. = vcovHC(fit_long, type = "HC1"))["obese", ]


#### MLR.5 Homoskedasticity ####
bptest(fit) #Breusch-Pagan

#White test (special form): regress u^2 on fitted values and their square
u2 <- fit$residuals^2
yhat <- fit$fitted.values
white <- lm(u2 ~ yhat + I(yhat^2))
LM <- n*summary(white)$r.squared
LM
pchisq(LM, df = 2, lower.tail = F)
#both reject: the errors are heteroskedastic

## Fix 1: heteroskedasticity-robust standard errors (HC1)
coeftest(fit, vcov. = vcovHC(fit, type = "HC1"))
range(sqrt(diag(vcovHC(fit, type = "HC1")))/sqrt(diag(vcov(fit))))
#robust SEs are 0.92 to 1.14 times the usual ones

## Fix 2: feasible GLS
## regress log(u^2) on the regressors, then weight by 1/exp(fitted)
X <- as.data.frame(model.matrix(fit)[, -1])
fit_var <- lm(log(fit$residuals^2) ~ ., data = X)
h <- exp(fit_var$fitted.values)
fit_fgls <- update(fit, weights = 1/h)
round(cbind(ols = fit$coefficients, fgls = fit_fgls$coefficients), 5)
#FGLS estimates are close to OLS, so we report OLS with robust SEs


## Robust versions of the group tests
linearHypothesis(fit, grep("^race", names(fit$coefficients), value = T),
                 vcov. = vcovHC(fit, type = "HC1"))
linearHypothesis(fit, grep("^smoke", names(fit$coefficients), value = T),
                 vcov. = vcovHC(fit, type = "HC1"))
linearHypothesis(fit, c("hdl", "ldl", "logtrig"),
                 vcov. = vcovHC(fit, type = "HC1"))
#race p < 0.001, smoking p = 0.27, lipids p < 0.001

#LDL slope for women
linearHypothesis(fit, "ldl + women:ldl = 0", vcov. = vcovHC(fit, type = "HC1"))
100*(exp(b["ldl"] + b["women:ldl"]) - 1)
#about zero for women


#### Outliers: median regression (LAD) ####
set.seed(2026)
fit_lad <- rq(log(sbp) ~ age + I(age^2) + women + age*women + race + smoke + diabetes
              + hdl + ldl + ldl*women + logtrig, tau = 0.5, data = compdat)
lad <- summary(fit_lad, se = "boot", R = 500)$coefficients
lad


#### Comparison across models (% effects) ####
pick <- c("age", "women", "raceNon-Hispanic Black", "diabetes", "hdl", "ldl", "logtrig")
mult <- c(10, 1, 1, 1, 1, 1, 1) #age per 10 years
compare <- cbind(main = fit$coefficients[pick], with_med = fit_med$coefficients[pick],
                 untreated = fit_untreated$coefficients[pick],
                 with_obese = fit_long$coefficients[pick],
                 fgls = fit_fgls$coefficients[pick], drop20 = fit_drop$coefficients[pick],
                 lad = lad[pick, 1])
round(100*(exp(compare*mult) - 1), 2)
100*(exp(fit_med$coefficients["bpmed"]) - 1) #medication itself, 4.16%
c(nobs(fit), nobs(fit_untreated), nobs(fit_long), nobs(fit_drop))
round(lad[pick, 4], 4) #LAD p-values
round(coeftest(fit_untreated, vcov. = vcovHC(fit_untreated, type = "HC1"))[pick, 4], 4)




########################## Figure ##########################
############################################################

png("output/figure1.png", width = 2400, height = 860, res = 220, family = "serif")
par(mfrow = c(1,2), mar = c(4.4, 4.4, 3, 1), las = 1)

plot(fit$fitted.values, fit$residuals, pch = 19, cex = 0.4,
     col = adjustcolor("#2a78d6", 0.5),
     xlab = "Fitted log(BP)", ylab = "Residual",
     main = "(a) Residuals vs fitted values", font.main = 1)
abline(h = 0)

plot(20:80, newdat$pred[which(newdat$women == 0)], typ = "l", lwd = 2, col = "#2a78d6",
     ylim = range(newdat$pred), xlab = "Age (years)",
     ylab = "Predicted systolic BP (mmHg)",
     main = "(b) Predicted blood pressure by age", font.main = 1)
lines(20:80, newdat$pred[which(newdat$women == 1)], lwd = 2, lty = 2, col = "#eb6834")
legend("topleft", legend = c("Men", "Women"), lwd = 2, lty = c(1, 2),
       col = c("#2a78d6", "#eb6834"), bty = "n")
dev.off()

sink()

