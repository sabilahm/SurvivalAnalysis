library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(survminer)
library(party)
library("coin")
library(tidyverse)
library(JM)
library(lattice)
library(MASS)
library(Rcmdr)
library(car)
library(lmtest)
library(stats)
library(StepReg)


data3 = read.csv("C:/Users/Hp/OneDrive - Institut Teknologi Sepuluh Nopember/Semester 8/Analisis Survival/FP Ansur/data3.csv")
head(data3)

data3$cancer <- as.factor(data3$cancer)
data3$treatment <- as.factor(data3$treatment)
data3$status <- as.numeric(data3$status)
data3$time <- as.numeric(data3$time)
data3$age <- as.numeric(data3$age)
data3$score <- as.numeric(data3$score)
data3$month <- as.numeric(data3$month)

surv_diff <- survdiff(Surv(time, status) ~ ., data = data3)
surv_diff

surv <- survfit(Surv(data3$time, data3$status==1)~1)
summary(surv)

Y <- Surv(data3$time, data3$status==1)
Y~1
Y~data3$treatment
kmfit1=survfit(Y~1)
plot(kmfit1,lty =c("solid","dashed"),col=c("black","grey"),xlab="Survival Time in Days", ylab="Survival Probabilities")+
  legend("topright",c("Treatment 1","Treatment 2"),lty=c("solid","dashed"),col=c("black","grey"))

treatment <- survdiff(Surv(data3$time, data3$status) ~ treatment, data=data3)
treatment

cancer <- survdiff(Surv(data3$time, data3$status) ~ cancer, data=data3)
cancer

modpar1= survreg(Surv(time, status) ~ ., data3, dist="exponential")
summary(modpar1)

modpar2= survreg(Surv(time, status) ~ ., data3, dist="weibull")
summary(modpar2)

modpar3= survreg(Surv(time, status) ~ ., data3, dist="lognormal")
summary(modpar3)

extractAIC(modpar1)
extractAIC(modpar2)
extractAIC(modpar3)

pattern1 = data.frame(cancer=1,age=50,score=50,month=2,treatment=1)
pct=c(.25,.50,.75)
days=predict(modpar1, newdata=pattern1, type="quantile",p=pct)
days

require(flexsurv)
s <- with(data3,Surv(time, status))
flexsurvreg(s ~ treatment,dist='exponential',data=data3)

coxph <- coxph(formula = Surv(time, status) ~ age+score+month+cancer+treatment, data = data3)
summary(coxph)

ci <- confint(coxph)
ci

exp(cbind(coef(coxph), ci))["horThyes",]

coxzph <- cox.zph(coxph)
coxzph

ctree <- ctree(Surv(time, status) ~ ., data = data3)

GBSG2_ctree <- ctree(Surv(time, status) ~ ., data = data3)
GBSG2_ctree

layout(matrix(1:3, ncol = 3))
res <- residuals(coxph)
plot(res ~ age, data = data3, ylim = c(-2.5, 1.5),
     + pch = ".", ylab = "Martingale Residuals")
abline(h = 0, lty = 3)
plot(res ~ score, data = data3, ylim = c(-2.5, 1.5),
     + pch = ".", ylab = "")
abline(h = 0, lty = 3)
plot(res ~ month, data = data3, ylim = c(-2.5, 1.5),
     + pch = ".", ylab = "")
abline(h = 0, lty = 3)
plot(res ~ cancer, data = data3, ylim = c(-2.5, 1.5),
     + pch = ".", ylab = "")
abline(h = 0, lty = 3)
plot(res ~ treatment, data = data3, ylim = c(-2.5, 1.5),
     + pch = ".", ylab = "")
abline(h = 0, lty = 3)

ageplot <- plot(coxzph, var = 'age')

age.ph <- coxph(Surv(time, status) ~ .,  data3, method="breslow")
data3$resid<- residuals(age.ph, type="martingale", data=data3)
plot(data3$age, data3$resid, xlab="Age",ylab="Martingale Residuals")
lines(lowess(data3$age, data3$resid))

score.ph <- coxph(Surv(time, status) ~ .,  data3, method="breslow")
data3$resid<- residuals(score.ph, type="martingale", data=data3)
plot(data3$score, data3$resid, xlab="Score",ylab="Martingale Residuals")
lines(lowess(data3$score, data3$resid))

month.ph <- coxph(Surv(time, status) ~ .,  data3, method="breslow")
data3$resid<- residuals(month.ph, type="martingale", data=data3)
plot(data3$month, data3$resid, xlab="Month",ylab="Martingale Residuals")
lines(lowess(data3$month, data3$resid))

cancer.ph <- coxph(Surv(time, status) ~ .,  data3, method="breslow")
data3$resid<- residuals(cancer.ph, type="martingale", data=data3)
plot(data3$cancer, data3$resid, xlab="Cancer",ylab="Martingale Residuals")
lines(lowess(data3$cancer, data3$resid))

treatment.ph <- coxph(Surv(time, status) ~ .,  data3, method="breslow")
data3$resid<- residuals(treatment.ph, type="martingale", data=data3)
plot(data3$treatment, data3$resid, xlab="Treatment",ylab="Martingale Residuals")
lines(lowess(data3$treatment, data3$resid))

resid_mart <- residuals(coxph, type = "martingale")

ggplot(data = data3, mapping = aes(x = age, y = resid_mart)) +
  geom_point() +
  geom_smooth() +
  labs(title = "age") +
  theme_bw() + theme(legend.key = element_blank())

ggplot(data = data3, mapping = aes(x = month, y = resid_mart)) +
  geom_point() +
  geom_smooth() +
  labs(title = "month") +
  theme_bw() + theme(legend.key = element_blank())

ggplot(data = data3, mapping = aes(x = cancer, y = resid_mart)) +
  geom_point() +
  geom_smooth() +
  labs(title = "cancer") +
  theme_bw() + theme(legend.key = element_blank())

ggplot(data = data3, mapping = aes(x = treatment, y = resid_mart)) +
  geom_point() +
  geom_smooth() +
  labs(title = "treatment") +
  theme_bw() + theme(legend.key = element_blank())

## Cox-Snell residuals
resid_coxsnell <- -(data3$resid - data3$status)


## Fit model on Cox-Snell residuals (Approximately Expo(1) distributed under correct model)
fit_coxsnell <- coxph(formula = Surv(resid_coxsnell, status) ~ 1,
                      data    = data3,
                      ties    = c("efron","breslow","exact")[1])

## Nelson-Aalen estimator for baseline hazard (all covariates zero)
df_base_haz <- basehaz(fit_coxsnell, centered = FALSE)
head(df_base_haz)

ggplot(data = df_base_haz, mapping = aes(x = time, y = hazard)) +
  geom_point() +
  scale_x_continuous(limit = c(0,3.5)) +
  scale_y_continuous(limit = c(0,3.5)) +
  labs(x = "Cox-Snell residuals as pseudo observed times",
       y = "Estimated cumulative hazard at pseudo observed times") +
  theme_bw() + theme(legend.key = element_blank())


gg_coxsnell(coxph) +
  geom_abline(intercept=0, slope=1, col=2)
gg_coxsnell(coxph, type="cdf") +
  geom_line(aes(y=F), col=2)

lmeFit.aids <- lme(CD4 ~ time + time:treatment,
                   random = ~ time, data = data3)

coxFit.aids <- coxph(Surv(Time, death) ~ drug, data = aids.id,
                     x = TRUE)

jointFit.aids <- jointModel(lmeFit.aids, coxFit.aids,
                            timeVar = "obstime", method = "piecewise-PH-aGH")

resCST <- residuals(data3, process = "Event",
                    type = "CoxSnell")
sfit <- survfit(Surv(resCST, status) ~ 1, data = data3)

data3 = read.csv("C:/Users/Hp/OneDrive - Institut Teknologi Sepuluh Nopember/Semester 8/Analisis Survival/FP Ansur/data3_time.csv")
head(data3)

data3$ID <- as.factor(data3$ID)
data3$cancer <- as.factor(data3$cancer)
data3$treatment <- as.factor(data3$treatment)
data3$status <- as.numeric(data3$status)
data3$time <- as.numeric(data3$time)
data3$age <- as.numeric(data3$age)
data3$score <- as.numeric(data3$score)
data3$month <- as.numeric(data3$month)

#TIME DEPENDENT
cutpoints=unique(data3$time[data3$status==1])
datatime=survSplit(data=data3, cut=cutpoints,
                  end="time",start="time_0", event="status")
head(datatime)
datatime2=datatime[order(datatime$ID),]
datatime3=datatime2[,c(2:9)]
datatime3$time <- as.numeric(datatime3$time)
model1<-coxph(Surv(time,status)~age+score+month+cancer+treatment, data=datatime3)
summary(model1)

MT <- stepAIC(model1, direction = "backward", trace = FALSE)
MT <- stepwise(model1, direction="backward")
MT <- step(coxph(Surv(time,status)~age+score+month+cancer+treatment, data=datatime3),direction="backward")

model2<-coxph(Surv(time,status)~age+month+cancer+treatment,data=datatime3)
summary(model2)

model3 <- coxph(formula = Surv(time, status) ~ age+month+cancer+treatment, data = datatime3)
model3

coxzph2 <- cox.zph(model3)
coxzph2
