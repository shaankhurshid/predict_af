# Script to perform main analyses underlying Predict-AF project
# Script assumes the presence of a suitable data.table with specific variables obtained from UKBB source data
# Please see readme for glossary of variables

# Dependencies
library(data.table)
library(survival)
library(rms)
library(prodlim)
library(timeROC)
library(nricens)
library(epiR)
library(dca)

# Source helper functions
source(file='functions.R')

############################################# STEP 1: Create scores (these coefficients have been derived from UKBB sample and are hard-coded here for replication)
## Create scores
af_set[,predict_af := charge*1.117871 + prs_norm*185.757558]
af_set[,prs_score := prs_norm*160.2864]
af_set[,age_score := agevisit0*0.1163611]

############################################# STEP 2: Calculate predicted risks (with centering upon baseline hazard of sample)
### CHARGE
af_set[,charge_avgbeta := mean(charge)]
res <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge, data=af_set)
km <- survfit(res, data=data.frame(x1=mean(charge)),type="kaplan-meier")
charge_s0 <- summary(km, times=c(5),extend=TRUE)$surv
af_set[,charge_pred5 := (1-(charge_s0)^exp(charge - charge_avgbeta))*100]

### PRS Linear Predictor
af_set[,prs_avgbeta := mean(prs_score)]
res <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_score, data=af_set)
km <- survfit(res, data=data.frame(x1=mean(prs_score)),type="kaplan-meier")
prs_s0 <- summary(km, times=c(5),extend=TRUE)$surv
af_set[,prs_pred5 := (1-(prs_s0)^exp(prs_score - (prs_avgbeta)))*100]

### PREDICT-AF
af_set[,predict_af_avgbeta := mean(predict_af)]
res <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af, data=af_set)
km <- survfit(res, data=data.frame(x1=mean(predict_af)),type="kaplan-meier")
predict_af_s0 <- summary(km, times=c(5),extend=TRUE)$surv
af_set[,predict_af_pred5 := (1-(predict_af_s0)^exp(predict_af - predict_af_avgbeta))*100]

### AGE Linear Predictor
# Average beta for age
af_set[,age_avgbeta := mean(age_score)]
res <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ age_score, data=af_set)
km <- survfit(res, data=af_set.frame(x1=mean(age_score)),type="kaplan-meier")
age_s0 <- summary(km, times=c(5),extend=TRUE)$surv
af_set[,age_pred5 := (1-(age_s0)^exp(age_score - age_avgbeta))*100]

# Cutoff variables
af_set[,':=' (charge_pred5_above5 = ifelse(charge_pred5 >= 5,1,0),
                            predict_af_pred5_above5 = ifelse(predict_af_pred5 >= 5,1,0))]

## Add values to stroke dataset
setkey(af_set,ID); setkey(stroke_set,ID)
stroke_set[af_set,':='(predict_af_pred5 = i.predict_af_pred5,
                                                     predict_af = i.predict_af,
                                                     charge_pred5 = i.charge_pred5,
                                                     prs_pred5 = i.prs_pred5)]

############################################# STEP 3: Compare model fits and concordance
###### Create standardized scores to facilitate comparison
af_set[,':='(charge_std = (charge - mean(charge))/sd(charge),
                           prs_std = (prs_norm - mean(prs_norm))/sd(prs_norm),
                           predict_af_std = (predict_af - mean(predict_af))/sd(predict_af))]
stroke_set[,':='(charge_std = (charge - mean(charge))/sd(charge),
                                 prs_std = (prs_norm - mean(prs_norm))/sd(prs_norm),
                                 predict_af_std = (predict_af - mean(predict_af))/sd(predict_af))]

###### AF comparisons
### Age
boot_age <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='age_score',data=af_set,runs=200)
age_cox <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ age_std,data=af_set)
age_cph <- cph(Surv(incd_af_5y.t,incd_af_5y) ~ age_score,data=af_set)

## Values
cstat_age <- c(summary(age_cox)$concordance[1],
               summary(age_cox)$concordance[1]-1.96*summary(age_cox)$concordance[2],
               summary(age_cox)$concordance[1]+1.96*summary(age_cox)$concordance[2])
r2_age <- c(age_cph$stats['R2'],age_cph$stats['R2']-1.96*sd(boot_age),age_cph$stats['R2']+1.96*sd(boot_age))
hr_age <- c(exp(age_cox$coefficients[1]),exp(confint(age_cox)[1]),exp(confint(age_cox)[2]))
aic_age <- AIC(age_cox)

### CHARGE
boot_charge <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='charge_std',data=af_set,runs=200)
charge_cox <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge_std,data=af_set)
charge_cph <- cph(Surv(incd_af_5y.t,incd_af_5y) ~ charge_std,data=af_set)

## Values
cstat_charge <- c(summary(charge_cox)$concordance[1],
                  summary(charge_cox)$concordance[1]-1.96*summary(charge_cox)$concordance[2],
                  summary(charge_cox)$concordance[1]+1.96*summary(charge_cox)$concordance[2])
r2_charge <- c(charge_cph$stats['R2'],charge_cph$stats['R2']-1.96*sd(boot_charge),charge_cph$stats['R2']+1.96*sd(boot_charge))
hr_charge <- c(exp(charge_cox$coefficients[1]),exp(confint(charge_cox)[1]),exp(confint(charge_cox)[2]))
aic_charge <- AIC(charge_cox)

### PRS Linear Predictor
boot_prs <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='prs_std',data=af_set,runs=200)
prs_cox <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_std,data=af_set)
prs_cph <- cph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_std,data=af_set)

## Values
cstat_prs <- c(summary(prs_cox)$concordance[1],
               summary(prs_cox)$concordance[1]-1.96*summary(prs_cox)$concordance[2],
               summary(prs_cox)$concordance[1]+1.96*summary(prs_cox)$concordance[2])
r2_prs <- c(prs_cph$stats['R2'],prs_cph$stats['R2']-1.96*sd(boot_prs),prs_cph$stats['R2']+1.96*sd(boot_prs))
hr_prs <- c(exp(prs_cox$coefficients[1]),exp(confint(prs_cox)[1]),exp(confint(prs_cox)[2]))
aic_prs <- AIC(prs_cox)

# Predict-AF
boot_predict <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='predict_af_std',data=af_set,runs=200)
predict_cox <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_std,data=af_set)
predict_cph <- cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_std,data=af_set)

## Values
cstat_predict <- c(summary(predict_cox)$concordance[1],
                   summary(predict_cox)$concordance[1]-1.96*summary(predict_cox)$concordance[2],
                   summary(predict_cox)$concordance[1]+1.96*summary(predict_cox)$concordance[2])
r2_predict <- c(predict_cph$stats['R2'],predict_cph$stats['R2']-1.96*sd(boot_predict),predict_cph$stats['R2']+1.96*sd(boot_predict))
hr_predict <- c(exp(predict_cox$coefficients[1]),exp(confint(predict_cox)[1]),exp(confint(predict_cox)[2]))
aic_predict <- AIC(predict_cox)

### Compare Predict-AF to CHARGE-AF using bootstrapping
predict_charge <- boot_compare(time="incd_af_5y.t",status='incd_af_5y',data=af_set,
                               response1='charge',response2='predict_af',runs=200)
c_diff_predict_charge <- c(cstat_predict[1]-cstat_charge[1],(cstat_predict[1]-cstat_charge[1])-1.96*sd(predict_charge),(cstat_predict[1]-cstat_charge[1])+1.96*sd(predict_charge))
p_predict_charge <- 2*(1-pnorm(abs(c_diff_predict_charge[1]/sd(predict_charge))))

### Discrimination in age < 65 subset
charge_cox_less65 <- coxph(Surv(incd_af_before65_5y.t,incd_af_before65_5y) ~ charge_std,data=af_set[!is.na(incd_af_before65_5y.t)])
predict_cox_less65 <- coxph(Surv(incd_af_before65_5y.t,incd_af_before65_5y) ~ predict_af_std,data=af_set[!is.na(incd_af_before65_5y.t)])

## Values
cstat_charge_less65 <- c(summary(charge_cox_less65)$concordance[1],
                         summary(charge_cox_less65)$concordance[1]-1.96*summary(charge_cox_less65)$concordance[2],
                         summary(charge_cox_less65)$concordance[1]+1.96*summary(charge_cox_less65)$concordance[2])

cstat_predict_less65 <- c(summary(predict_cox_less65)$concordance[1],
                          summary(predict_cox_less65)$concordance[1]-1.96*summary(predict_cox_less65)$concordance[2],
                          summary(predict_cox_less65)$concordance[1]+1.96*summary(predict_cox_less65)$concordance[2])

###### Stroke comparisons
### Age
boot_age <- boot_r2(time='stroke_5y.t',status='stroke_5y',response='age_std',data=stroke_set,runs=200)
age_cox <- coxph(Surv(stroke_5y.t,stroke_5y) ~ age_std,data=stroke_set)
age_cph <- cph(Surv(stroke_5y.t,stroke_5y) ~ age_std,data=stroke_set)

## Values
cstat_age <- c(summary(age_cox)$concordance[1],
               summary(age_cox)$concordance[1]-1.96*summary(age_cox)$concordance[2],
               summary(age_cox)$concordance[1]+1.96*summary(age_cox)$concordance[2])
r2_age <- c(age_cph$stats['R2'],age_cph$stats['R2']-1.96*sd(boot_age),age_cph$stats['R2']+1.96*sd(boot_age))
hr_age <- c(exp(age_cox$coefficients[1]),exp(confint(age_cox)[1]),exp(confint(age_cox)[2]))
aic_age <- AIC(age_cox)

### CHARGE
boot_charge <- boot_r2(time='stroke_5y.t',status='stroke_5y',response='charge_std',data=stroke_set,runs=200)
charge_cox <- coxph(Surv(stroke_5y.t,stroke_5y) ~ charge_std,data=stroke_set)
charge_cph <- cph(Surv(stroke_5y.t,stroke_5y) ~ charge_std,data=stroke_set)

## Values
cstat_charge <- c(summary(charge_cox)$concordance[1],
                  summary(charge_cox)$concordance[1]-1.96*summary(charge_cox)$concordance[2],
                  summary(charge_cox)$concordance[1]+1.96*summary(charge_cox)$concordance[2])
r2_charge <- c(charge_cph$stats['R2'],charge_cph$stats['R2']-1.96*sd(boot_charge),charge_cph$stats['R2']+1.96*sd(boot_charge))
hr_charge <- c(exp(charge_cox$coefficients[1]),exp(confint(charge_cox)[1]),exp(confint(charge_cox)[2]))
aic_charge <- AIC(charge_cox)

### PRS Linear Predictor
boot_prs <- boot_r2(time='stroke_5y.t',status='stroke_5y',response='prs_std',data=stroke_set,runs=200)
prs_cox <- coxph(Surv(stroke_5y.t,stroke_5y) ~ prs_std,data=stroke_set)
prs_cph <- cph(Surv(stroke_5y.t,stroke_5y) ~ prs_std,data=stroke_set)

## Values
cstat_prs <- c(summary(prs_cox)$concordance[1],
               summary(prs_cox)$concordance[1]-1.96*summary(prs_cox)$concordance[2],
               summary(prs_cox)$concordance[1]+1.96*summary(prs_cox)$concordance[2])
r2_prs <- c(prs_cph$stats['R2'],prs_cph$stats['R2']-1.96*sd(boot_prs),prs_cph$stats['R2']+1.96*sd(boot_prs))
hr_prs <- c(exp(prs_cox$coefficients[1]),exp(confint(prs_cox)[1]),exp(confint(prs_cox)[2]))
aic_prs <- AIC(prs_cox)

# Predict-AF
boot_predict <- boot_r2(time='stroke_5y.t',status='stroke_5y',response='predict_af_std',data=stroke_set,runs=200)
predict_cox <- coxph(Surv(stroke_5y.t,stroke_5y) ~ predict_af_std,data=stroke_set)
predict_cph <- cph(Surv(stroke_5y.t,stroke_5y) ~ predict_af_std,data=stroke_set)

## Values
cstat_predict <- c(summary(predict_cox)$concordance[1],
                   summary(predict_cox)$concordance[1]-1.96*summary(predict_cox)$concordance[2],
                   summary(predict_cox)$concordance[1]+1.96*summary(predict_cox)$concordance[2])
r2_predict <- c(predict_cph$stats['R2'],predict_cph$stats['R2']-1.96*sd(boot_predict),predict_cph$stats['R2']+1.96*sd(boot_predict))
hr_predict <- c(exp(predict_cox$coefficients[1]),exp(confint(predict_cox)[1]),exp(confint(predict_cox)[2]))
aic_predict <- AIC(predict_cox)

############################################# STEP 4: Obtain, plot, and compare 5 year event CIs based on risk cutoffs
### Apply risk cutoffs
af_set[,':='(charge_pred5_above5 = ifelse(charge_pred5 >= 5,1,0),
                           predict_af_pred5_above5 = ifelse(predict_af_pred5 >= 5,1,0))]
stroke_set[,':='(charge_pred5_above5 = ifelse(charge_pred5 >= 5,1,0),
                                 predict_af_pred5_above5 = ifelse(predict_af_pred5 >= 5,1,0))]

# Variables for graphical display
stroke_set[,':='(charge_pred_above5_graphical = factor(ifelse(charge_pred5_above5==1,"Yes","No"),levels=c("Yes","No")))]
stroke_set[,':='(predict_af_pred_above5_graphical = factor(ifelse(predict_af_pred5_above5==1,"Yes","No"),levels=c("Yes","No")))]
stroke_set[,':='(age65_visit0_graphical = factor(ifelse(age65_visit0==1,"Yes","No"),levels=c("Yes","No")))]

af_set[,':='(charge_pred_above5_graphical = factor(ifelse(charge_pred5_above5==1,"Yes","No"),levels=c("Yes","No")))]
af_set[,':='(predict_af_pred_above5_graphical = factor(ifelse(predict_af_pred5_above5==1,"Yes","No"),levels=c("Yes","No")))]
af_set[,':='(age65_visit0_graphical = factor(ifelse(age65_visit0==1,"Yes","No"),levels=c("Yes","No")))]

###### Plot cumulative event rates
############# PLOTS FOR AF
# CHARGE
af_charge <- prodlim(Hist(incd_af_5y.t,incd_af_5y)~charge_pred_above5_graphical,data=af_set)

CairoPDF(file='af_charge.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(af_charge,"cuminc",ylim=c(0,0.10),
     axis2.at=seq(0,0.10,0.025),axis2.las=2,lwd=1.5,background=F,
     atrisk.times=c(0,1,2,3,4,5),col=c("#253494","#41b6c4"),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.10,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("CHARGE-AF risk \u2265 5%","CHARGE-AF risk < 5%"),
     atrisk.title=("                     "),atrisk.pos=-0.35,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative incidence of AF (%)",side=2,line=-1.2,at=0.05,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('AF risk \u2265 5%',side=1, line=-0.5,cex=1.8,at=-0.35)
dev.off()

## Predict-AF
af_predict <- prodlim(Hist(incd_af_5y.t,incd_af_5y)~predict_af_pred_above5_graphical,data=af_set)

CairoPDF(file='af_predict.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(af_predict,"cuminc",ylim=c(0,0.10),
     axis2.at=seq(0,0.10,0.025),axis2.las=2,lwd=1.5,background=F,
     atrisk.times=c(0,1,2,3,4,5),col=c("#253494","#41b6c4"),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.10,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Predict-AF risk \u2265 5%","Predict-AF risk < 5%"),
     atrisk.title=("                     "),atrisk.pos=-0.35,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative incidence of AF (%)",side=2,line=-1.2,at=0.05,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('AF risk \u2265 5%',side=1, line=-0.5,cex=1.8,at=-0.35)
dev.off()

## Age >= 65
af_65 <- prodlim(Hist(incd_af_5y.t,incd_af_5y)~age65_visit0_graphical,data=af_set)

CairoPDF(file='af_65.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(af_65,"cuminc",ylim=c(0,0.10),
     axis2.at=seq(0,0.10,0.025),axis2.las=2,lwd=1.5,background=F,
     atrisk.times=c(0,1,2,3,4,5),col=c("#253494","#41b6c4"),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.10,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Age \u2265 65","Age < 65"),
     atrisk.title=("                     "),atrisk.pos=-0.35,atrisk.line=c(1.2,2.5),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative incidence of AF (%)",side=2,line=-1.2,at=0.05,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Age \u2265 65',side=1, line=-0.5,cex=1.8,at=-0.35)
dev.off()

############# PLOTS FOR STROKE
# CHARGE-AF
stroke_charge5 <- prodlim(Hist(stroke_5y.t,stroke_5y)~charge_pred_above5_graphical,data=stroke_set)

CairoPDF(file='stroke_charge.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(stroke_charge5,"cuminc",ylim=c(0,0.025),xlim=c(0,5),background=F,
     axis2.at=c(0,0.005,0.01,0.015,0.02,0.025),lwd=1.5,axis2.las=2,
     atrisk.times=c(0,1,2,3,4,5),col=c("#253494","#41b6c4"),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.025,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.legend=c("CHARGE-AF risk \u2265 5%","CHARGE-AF risk < 5%"),legend.cex=2.2,
     atrisk.title=("                     "),atrisk.pos=-0.3,atrisk.line=c(1.2,2.5),
     atrisk.dist=0.6,atrisk.adj=1,
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative incidence of stroke (%)",side=2, line=-1.2,at=0.01,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('AF risk \u2265 5%',side=1, line=-0.5,cex=1.8,at=-0.35)
dev.off()

## Predict-AF 
stroke_predict5 <- prodlim(Hist(stroke_5y.t,stroke_5y)~predict_af_pred_above5_graphical,data=stroke_set)

CairoPDF(file='stroke_predict.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(stroke_predict5,"cuminc",ylim=c(0,0.025),xlim=c(0,5),background=F,
     axis2.at=c(0,0.005,0.01,0.015,0.02,0.025),lwd=1.5,axis2.las=2,
     atrisk.times=c(0,1,2,3,4,5),col=c("#253494","#41b6c4"),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.025,axis1.cex.axis=2.5,axis2.cex.axis=2.5, legend.cex=2.2,
     legend.legend=c("Predict-AF risk \u2265 5%","Predict-AF risk < 5%"),axis1.padj=0.5,
     atrisk.title=("                     "),atrisk.pos=-0.3,atrisk.line=c(1.2,2.5),
     atrisk.dist=0.6,atrisk.adj=1,
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative incidence of stroke (%)",side=2, line=-1.2,at=0.01,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('AF risk \u2265 5%',side=1, line=-0.5,cex=1.8,at=-0.35)
dev.off()

## Age >= 65
stroke_65 <- prodlim(Hist(stroke_5y.t,stroke_5y)~age65_visit0_graphical,data=stroke_set)

CairoPDF(file='stroke_65.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(stroke_65,"cuminc",ylim=c(0,0.025),xlim=c(0,5),background=F,
     axis2.at=c(0,0.005,0.01,0.015,0.02,0.025),lwd=1.5,axis2.las=2,
     atrisk.times=c(0,1,2,3,4,5),col=c("#253494","#41b6c4"),atrisk.col='black',confint=FALSE,
     legend.x=0,legend.y=0.025,axis1.cex.axis=2.5,axis2.cex.axis=2.5,
     legend.legend=c("Age \u2265 65","Age < 65"), legend.cex=2.2,atrisk.line=c(1.2,2.5),
     atrisk.title=("                     "),atrisk.pos=-0.3,axis1.padj=0.5,
     atrisk.dist=0.6,atrisk.adj=1,
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative incidence of stroke (%)",side=2, line=-1.2,at=0.01,cex=2.5)
mtext("Years",side=1, line=-0.6,cex=2.5)
mtext('Age \u2265 65',side=1, line=-0.5,cex=1.8,at=-0.35)
dev.off()

############# Now quantify the event rates at 5 years
###### AF
#### All population
af_rate_age <- cuminc(data=af_set[age65_visit0==1],time='incd_af_5y.t',status='incd_af_5y')
af_rate_charge <- cuminc(data=af_set[charge_pred5_above5==1],time='incd_af_5y.t',status='incd_af_5y')
af_rate_predict <- cuminc(data=af_set[predict_af_pred5_above5==1],time='incd_af_5y.t',status='incd_af_5y')

#### Age < 65
af_rate_charge_less65 <- cuminc(data=af_set[charge_pred5_above5==1 & !is.na(incd_af_before65_5y.t)],time='incd_af_before65_5y.t',status='incd_af_before65_5y')
af_rate_predict_less65 <- cuminc(data=af_set[predict_af_pred5_above5==1 & !is.na(incd_af_before65_5y.t)],time='incd_af_before65_5y.t',status='incd_af_before65_5y')

#### Age >= 65
af_rate_charge_gr65 <- cuminc(data=af_set[charge_pred5_above5==1 & age65_visit0==1],time='incd_af_after65_5y.t',status='incd_af_after65_5y')
af_rate_predict_gr65 <- cuminc(data=af_set[predict_af_pred5_above5==1 & age65_visit0==1],time='incd_af_after65_5y.t',status='incd_af_after65_5y')

###### Stroke
#### All population
stroke_rate_age <- cuminc(data=stroke_set[age65_visit0==1],time='stroke_5y.t',status='stroke_5y')
stroke_rate_charge <- cuminc(data=stroke_set[charge_pred5_above5==1],time='stroke_5y.t',status='stroke_5y')
stroke_rate_predict <- cuminc(data=stroke_set[predict_af_pred5_above5==1],time='stroke_5y.t',status='stroke_5y')

#### Age < 65
stroke_rate_charge_less65 <- cuminc(data=stroke_set[charge_pred5_above5==1 & !is.na(stroke_before65_5y.t)],time='stroke_before65_5y.t',status='stroke_before65_5y')
stroke_rate_predict_less65 <- cuminc(data=stroke_set[predict_af_pred5_above5==1 & !is.na(stroke_before65_5y.t)],time='stroke_before65_5y.t',status='stroke_before65_5y')

#### Age >= 65
stroke_rate_charge_gr65 <- cuminc(data=stroke_set[charge_pred5_above5==1 & age65_visit0==1],time='stroke_after65_5y.t',status='stroke_after65_5y')
stroke_rate_predict_gr65 <- cuminc(data=stroke_set[predict_af_pred5_above5==1 & age65_visit0==1],time='stroke_after65_5y.t',status='stroke_after65_5y')

############################################# STEP 6: Density Plots
### Density Plots
# Create separate data.tables for cases/controls
gp_incident_af <- af_set[incd_af_5y==1,]
gp_noaf <- af_set[incd_af_5y==0,]

# CHARGE VS Predict-AF (controls)
x <- list(v1=gp_noaf$charge_pred5,v2=gp_noaf$predict_af_pred5)
data <- melt(x)

# Plot distribution
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
  scale_x_continuous(breaks=seq(0,10,2),expand=c(0,0),limits=c(0,10)) +
  scale_y_continuous(breaks=seq(0,0.90,0.10),expand=c(0,0),limits=c(0,0.90)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('CHARGE-AF','Predict-AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.8,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted 5-year AF risk (%)',y='density')
ggsave(filename='density_predict_ctrl.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

# CHARGE vs Predict-AF (cases)
x <- list(v1=gp_incident_af$charge_pred5,v2=gp_incident_af$predict_af_pred5)
data <- melt(x)

# Plot CHARGE PRED 5 distribution
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
  scale_x_continuous(breaks=seq(0,20,1),expand=c(0,0),limits=c(0,20)) +
  scale_y_continuous(breaks=seq(0,0.25,0.05),expand=c(0,0),limits=c(0,0.25)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('CHARGE-AF','Predict-AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.8,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predicted 5-year AF risk (%)',y='density')
ggsave(filename='density_predict_case.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

## Plotting separately among scores
### Density Plots
## Raw scores
# Generate CHARGE distribution
x <- list(v1=gp_incident_af$charge,v2=gp_noaf$charge)
data <- melt(x)

# Plot CHARGE distribution
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
  scale_x_continuous(breaks=seq(8,16,1),expand=c(0,0.1),limits=c(8,16)) +
  scale_y_continuous(breaks=seq(0,0.60,0.05),expand=c(0,0),limits=c(0,0.60)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('incident AF','no incident AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='CHARGE-AF Score',y='density')
ggsave(filename='density_charge_raw.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

# Generate PRS distribution
x <- list(v1=gp_incident_af$prs,v2=gp_noaf$prs)
data <- melt(x)

# Plot PRS distribution
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
  scale_x_continuous(breaks=seq(80,120,5),expand=c(0,0.5),limits=c(80,120)) +
  scale_y_continuous(breaks=seq(0,0.12,0.01),expand=c(0,0),limits=c(0,0.12)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('incident AF','no incident AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='AF PRS Score',y='density')
ggsave(filename='density_prs_raw.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

# Generate CHARGE PRS distribution
x <- list(v1=gp_incident_af$predict_af,v2=gp_noaf$predict_af)
data <- melt(x)

# Plot Predict-AF distribution
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) +
  scale_x_continuous(breaks=seq(16,26,2),expand=c(0,0.1),limits=c(16,26)) +
  scale_y_continuous(breaks=seq(0,0.50,0.05),expand=c(0,0),limits=c(0,0.50)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('incident AF','no incident AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.20,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='Predict-AF Score',y='density')
ggsave(filename='density_predict_raw.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

############################################# STEP 7: Calibration
### Raw slopes
## Models
age <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ age_score,data=af_set)
charge <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge,data=af_set)
prs <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_score,data=af_set)
predict <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=af_set)

## Slopes
cal_age <- c(age$coefficients[1],confint(age)[1],confint(age)[2])
cal_charge <- c(charge$coefficients[1],confint(charge)[1],confint(charge)[2])
cal_prs <- c(prs$coefficients[1],confint(prs)[1],confint(prs)[2])
cal_predict <- c(predict$coefficients[1],confint(predict)[1],confint(predict)[2])

### Correct internal estimates
## Estimate bias
boot_age <- boot_cal(time='incd_af_5y.t',status='incd_af_5y',response='age_score',data=af_set,runs=200)
boot_prs <- boot_cal(time='incd_af_5y.t',status='incd_af_5y',response='prs_score',data=af_set,runs=200)
boot_predict <- boot_cal(time='incd_af_5y.t',status='incd_af_5y',response='predict_af',data=af_set,runs=200)

## Correct
corrected_age <- cal_age[1]*2 - mean(boot_age)
corrected_prs <- cal_prs[1]*2 - mean(boot_prs)
corrected_predict_pt  <- cal_predict[1]*2 - mean(boot_predict)

############################################# STEP 8: NRI
## All population
charge_v_predict <- nricens(p.std=af_set$charge_pred5_above5, p.new=af_set$predict_af_pred5_above5,
                            time=af_set$incd_af_5y.t, event=af_set$incd_af_5y, cut=0.5,
                            niter = 10, t0=5,updown='category')

age_v_predict <- nricens(p.std=af_set$age65_visit0, p.new=af_set$predict_af_pred5_above5,
                         time=af_set$incd_af_5y.t, event=af_set$incd_af_5y, cut=0.5,
                         niter = 10, t0=5,updown='category')

age_v_charge <- nricens(p.std=af_set$age65_visit0, p.new=af_set$charge_pred5_above5,
                        time=af_set$incd_af_5y.t, event=af_set$incd_af_5y, cut=0.5,
                        niter = 10, t0=5,updown='category')

## Age < 65 subset
youth <- af_set[!is.na(incd_af_before65_5y.t)]
charge_v_chargeprs <- nricens(p.std=youth$charge_pred5_above5, p.new=youth$predict_af_pred5_above5,
                              time=youth$incd_af_before65_5y.t, event=youth$incd_af_before65_5y, cut=0.5,
                              niter = 10, t0=5,updown='category')

############################################# STEP 9: Sensitivity Analysis: Age > 70
######## AF
## Get 5-year AF incidence rates in subset
af_rate_age70 <- cuminc(data=af_set[age70_visit0==1],time='incd_af_5y.t',status='incd_af_5y')

## Bootstrap to obtain event rates
# AGE 65 versus AGE 70
af_65_70 <- boot_cuminc(time='incd_af_5y.t',status='incd_af_5y',response1='age65_visit0',response2='age70_visit0',data=af_set,runs=200)

# CHARGE-AF versus AGE 70
af_charge_70 <- boot_cuminc(time='incd_af_5y.t',status='incd_af_5y',response1='charge_pred5_above5',response2='age70_visit0',data=af_set,runs=200)

# Predict-AF versus AGE 70
af_predict_70 <- boot_cuminc(time='incd_af_5y.t',status='incd_af_5y',response1='predict_af_pred5_above5',response2='age70_visit0',data=af_set,runs=200)

### Z-tests using SE of difference in event rates obtained via bootstrapping above
z_af_65_70 <- as.numeric(((af_rate_age[,3]/100  - af_rate_age70[,3]/100) / sd(af_65_70[,3])))
p_af_65_70 <- 2*(1-pnorm(abs(z_af_65_70)))

z_af_charge_70 <- as.numeric(((af_rate_charge[,3]/100  - af_rate_age70[,3]/100) / sd(af_charge_70[,3])))
p_af_charge_70 <- 2*(1-pnorm(abs(z_af_charge_70)))

z_af_predict_70 <- as.numeric(((af_rate_predict[,3]/100  - af_rate_age70[,3]/100) / sd(af_predict_70[,3])))
p_af_predict_70 <- 2*(1-pnorm(abs(z_af_predict_70)))

######## Stroke
## 5-year stroke
stroke_rate_age70 <- cuminc(data=stroke_set[age70_visit0==1],time='stroke_5y.t',status='stroke_5y')

## Bootstrap to obtain event rates
# AGE 65 versus AGE 70
stroke_65_70 <- boot_cuminc(time='stroke_5y.t',status='stroke_5y',response1='age65_visit0',response2='age70_visit0',data=stroke_set,runs=200)

# CHARGE-AF versus AGE 70
stroke_charge_70 <- boot_cuminc(time='stroke_5y.t',status='stroke_5y',response1='charge_pred5_above5',response2='age70_visit0',data=stroke_set,runs=200)

# Predict-AF versus AGE 70
stroke_predict_70 <- boot_cuminc(time='stroke_5y.t',status='stroke_5y',response1='predict_af_pred5_above5',response2='age70_visit0',data=stroke_set,runs=200)

### Z-test using SE of difference in event rates obtained via bootstrapping above
z_stroke_65_70 <- as.numeric(((stroke_rate_age[,3]/100  - stroke_rate_age70[,3]/100) / sd(stroke_65_70[,3])))
p_stroke_65_70 <- 2*(1-pnorm(abs(z_stroke_65_70)))

z_stroke_charge_70 <- as.numeric(((stroke_rate_charge[,3]/100  - stroke_rate_age70[,3]/100) / sd(stroke_charge_70[,3])))
p_stroke_charge_70 <- 2*(1-pnorm(abs(z_stroke_charge_70)))

z_stroke_predict_70 <- as.numeric(((stroke_rate_predict[,3]/100  - stroke_rate_age70[,3]/100) / sd(stroke_predict_70[,3])))
p_stroke_predict_70 <- 2*(1-pnorm(abs(z_stroke_predict_70)))

############################################# STEP 10: EM by sex analysis
mod <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af + sex + predict_af:sex, data=af_set)

m <- af_set[sex=='male']
f <- af_set[sex=='female']

# Re-standardize within groups
f[,predict_af_std := (predict_af-mean(predict_af))/sd(predict_af)]
m[,predict_af_std := (predict_af-mean(predict_af))/sd(predict_af)]

mod_m <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_std, data=m)
print(AIC(mod_m))
print(c(summary(mod_m)$concordance[1],
        summary(mod_m)$concordance[1]-1.96*summary(mod_m)$concordance[2],
        summary(mod_m)$concordance[1]+1.96*summary(mod_m)$concordance[2]))

mod_f <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_std, data=f)
print(AIC(mod_f))
print(c(summary(mod_f)$concordance[1],
        summary(mod_f)$concordance[1]-1.96*summary(mod_f)$concordance[2],
        summary(mod_f)$concordance[1]+1.96*summary(mod_f)$concordance[2]))

# M
boot_predict_m <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='predict_af',data=m,runs=200)
## Values
r2_predict_m <- c(cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=m)$stats['R2'],
                  cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=m)$stats['R2']-1.96*sd(boot_predict_m),
                  cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=m)$stats['R2']+1.96*sd(boot_predict_m))

# F
boot_predict_f <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='predict_af',data=f,runs=200)
## Values
r2_predict_f <- c(cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=f)$stats['R2'],
                  cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=f)$stats['R2']-1.96*sd(boot_predict_f),
                  cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=f)$stats['R2']+1.96*sd(boot_predict_f))

############################################# STEP 11: Summed predict
## Define simple sum
af_set[,simple_prs := charge + prs]
stroke_set[,simple_prs := charge + prs]

## Standardize simple PRS
af_set[,simple_prs_std := (simple_prs - mean(simple_prs))/sd(simple_prs)]
stroke_set[,simple_prs_std := (simple_prs - mean(simple_prs))/sd(simple_prs)]

## Associations
af_simple_prs <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ simple_prs_std,data=af_set)
stroke_simple_prs <- coxph(Surv(stroke_5y.t,stroke_5y) ~ simple_prs_std,data=stroke_set)

## Values
boot_simple <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='simple_prs',data=af_set,runs=200)
cstat_simple <- c(summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ simple_prs,data=af_set))$concordance[1],
                  summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ simple_prs,data=af_set))$concordance[1]-1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ simple_prs,data=af_set))$concordance[2],
                  summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ simple_prs,data=af_set))$concordance[1]+1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ simple_prs,data=af_set))$concordance[2])
r2_simple <- c(cph(Surv(incd_af_5y.t,incd_af_5y) ~ simple_prs,data=af_set)$stats['R2'],
               cph(Surv(incd_af_5y.t,incd_af_5y) ~ simple_prs,data=af_set)$stats['R2']-1.96*sd(boot_simple),
               cph(Surv(incd_af_5y.t,incd_af_5y) ~ simple_prs,data=af_set)$stats['R2']+1.96*sd(boot_simple))
AIC(af_simple_prs)

## Values
boot_simple_stroke <- boot_r2(time='stroke_5y.t',status='stroke_5y',response='simple_prs',data=stroke_set,runs=200)
cstat_simple <- c(summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ simple_prs,data=stroke_set))$concordance[1],
                  summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ simple_prs,data=stroke_set))$concordance[1]-1.96*summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ simple_prs,data=stroke_set))$concordance[2],
                  summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ simple_prs,data=stroke_set))$concordance[1]+1.96*summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ simple_prs,data=stroke_set))$concordance[2])
r2_simple_stroke <- c(cph(Surv(stroke_5y.t,stroke_5y) ~ simple_prs,data=stroke_set)$stats['R2'],
                      cph(Surv(stroke_5y.t,stroke_5y) ~ simple_prs,data=stroke_set)$stats['R2']-1.96*sd(boot_simple_stroke),
                      cph(Surv(stroke_5y.t,stroke_5y) ~ simple_prs,data=stroke_set)$stats['R2']+1.96*sd(boot_simple_stroke))
AIC(stroke_simple_prs)

############################################# STEP 12: Quartiles
# Generate quartiles of risk
af_set[,':='(charge_pred5_quartile = ifelse(charge_pred5 >= quantile(charge_pred5,probs=0.75),4,
                                                          ifelse(charge_pred5 >= quantile(charge_pred5,probs=0.50),3,
                                                                 ifelse(charge_pred5 >= quantile(charge_pred5,probs=0.25),2,1))),
                           prs_pred5_quartile = ifelse(prs_pred5 >= quantile(prs_pred5,probs=0.75),4,
                                                       ifelse(prs_pred5 >= quantile(prs_pred5,probs=0.50),3,
                                                              ifelse(prs_pred5 >= quantile(prs_pred5,probs=0.25),2,1))),
                           predict_af_pred5_quartile = ifelse(predict_af_pred5 >= quantile(predict_af_pred5,probs=0.75),4,
                                                              ifelse(predict_af_pred5 >= quantile(predict_af_pred5,probs=0.50),3,
                                                                     ifelse(predict_af_pred5 >= quantile(predict_af_pred5,probs=0.25),2,1))))]

stroke_set[,':='(charge_pred5_quartile = ifelse(charge_pred5 >= quantile(charge_pred5,probs=0.75),4,
                                                                ifelse(charge_pred5 >= quantile(charge_pred5,probs=0.50),3,
                                                                       ifelse(charge_pred5 >= quantile(charge_pred5,probs=0.25),2,1))),
                                 prs_pred5_quartile = ifelse(prs_pred5 >= quantile(prs_pred5,probs=0.75),4,
                                                             ifelse(prs_pred5 >= quantile(prs_pred5,probs=0.50),3,
                                                                    ifelse(prs_pred5 >= quantile(prs_pred5,probs=0.25),2,1))),
                                 predict_af_pred5_quartile = ifelse(predict_af_pred5 >= quantile(predict_af_pred5,probs=0.75),4,
                                                                    ifelse(predict_af_pred5 >= quantile(predict_af_pred5,probs=0.50),3,
                                                                           ifelse(predict_af_pred5 >= quantile(predict_af_pred5,probs=0.25),2,1))))]

# Generate quartiles of age
af_set[,':='(age_quartile = ifelse(agevisit0 >= quantile(agevisit0,probs=0.75),4,
                                                 ifelse(agevisit0 >= quantile(agevisit0,probs=0.50),3,
                                                        ifelse(agevisit0 >= quantile(agevisit0,probs=0.25),2,1))))]

stroke_set[,':='(age_quartile = ifelse(agevisit0 >= quantile(agevisit0,probs=0.75),4,
                                                       ifelse(agevisit0 >= quantile(agevisit0,probs=0.50),3,
                                                              ifelse(agevisit0 >= quantile(agevisit0,probs=0.25),2,1))))]

# Models
## AF
mod_predict <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ factor(predict_af_pred5_quartile),data=af_set)
mod_age <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ factor(age_quartile),data=af_set)

## Stroke
mod_predict_stroke <- coxph(Surv(stroke_5y.t,stroke_5y) ~ factor(predict_af_pred5_quartile),data=stroke_set)
mod_age_stroke <- coxph(Surv(stroke_5y.t,stroke_5y) ~ factor(age_quartile),data=stroke_set)

## AF
# AGE
boot_age <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='age_quartile',data=af_set,runs=200,size=nrow(af_set))
## Values
cstat_age <- c(summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ age_quartile,data=af_set))$concordance[1],
               summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ age_quartile,data=af_set))$concordance[1]-1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ age_quartile,data=af_set))$concordance[2],
               summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ age_quartile,data=af_set))$concordance[1]+1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ age_quartile,data=af_set))$concordance[2])
r2_age <- c(cph(Surv(incd_af_5y.t,incd_af_5y) ~ age_quartile,data=af_set)$stats['R2'],
            cph(Surv(incd_af_5y.t,incd_af_5y) ~ age_quartile,data=af_set)$stats['R2']-1.96*sd(boot_age),
            cph(Surv(incd_af_5y.t,incd_af_5y) ~ age_quartile,data=af_set)$stats['R2']+1.96*sd(boot_age))
AIC(mod_age)

# Predict-AF
boot_predict <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='predict_af_pred5_quartile',data=af_set,runs=200,size=nrow(af_set))
## Values
cstat_predict <- c(summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_pred5_quartile,data=af_set))$concordance[1],
                   summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_pred5_quartile,data=af_set))$concordance[1]-1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_pred5_quartile,data=af_set))$concordance[2],
                   summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_pred5_quartile,data=af_set))$concordance[1]+1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_pred5_quartile,data=af_set))$concordance[2])
r2_predict <- c(cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_pred5_quartile,data=af_set)$stats['R2'],
                cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_pred5_quartile,data=af_set)$stats['R2']-1.96*sd(boot_predict),
                cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_pred5_quartile,data=af_set)$stats['R2']+1.96*sd(boot_predict))
AIC(mod_predict)

## Stroke
# AGE
boot_age <- boot_r2(time='stroke_5y.t',status='stroke_5y',response='age_quartile',data=stroke_set,runs=200,size=nrow(stroke_set))
## Values
cstat_age <- c(summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ age_quartile,data=stroke_set))$concordance[1],
               summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ age_quartile,data=stroke_set))$concordance[1]-1.96*summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ age_quartile,data=stroke_set))$concordance[2],
               summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ age_quartile,data=stroke_set))$concordance[1]+1.96*summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ age_quartile,data=stroke_set))$concordance[2])
r2_age <- c(cph(Surv(stroke_5y.t,stroke_5y) ~ age_quartile,data=stroke_set)$stats['R2'],
            cph(Surv(stroke_5y.t,stroke_5y) ~ age_quartile,data=stroke_set)$stats['R2']-1.96*sd(boot_age),
            cph(Surv(stroke_5y.t,stroke_5y) ~ age_quartile,data=stroke_set)$stats['R2']+1.96*sd(boot_age))
AIC(mod_age_stroke)

# Predict-AF
boot_predict <- boot_r2(time='stroke_5y.t',status='stroke_5y',response='predict_af_pred5_quartile',data=stroke_set,runs=200,size=nrow(stroke_set))
## Values
cstat_predict <- c(summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ predict_af_pred5_quartile,data=stroke_set))$concordance[1],
                   summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ predict_af_pred5_quartile,data=stroke_set))$concordance[1]-1.96*summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ predict_af_pred5_quartile,data=stroke_set))$concordance[2],
                   summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ predict_af_pred5_quartile,data=stroke_set))$concordance[1]+1.96*summary(coxph(Surv(stroke_5y.t,stroke_5y) ~ predict_af_pred5_quartile,data=stroke_set))$concordance[2])
r2_predict <- c(cph(Surv(stroke_5y.t,stroke_5y) ~ predict_af_pred5_quartile,data=stroke_set)$stats['R2'],
                cph(Surv(stroke_5y.t,stroke_5y) ~ predict_af_pred5_quartile,data=stroke_set)$stats['R2']-1.96*sd(boot_predict),
                cph(Surv(stroke_5y.t,stroke_5y) ~ predict_af_pred5_quartile,data=stroke_set)$stats['R2']+1.96*sd(boot_predict))
AIC(mod_predict_stroke)

############################################# STEP 13: EUR subset
# Define white_noprevaf
eur_set <- af_set[genetic_white==1]

### Re-standardization
eur_set[,':='(charge_std = (charge-mean(charge))/sd(charge),
                     predict_af_std = (predict_af-mean(predict_af))/sd(predict_af),
                     prs_std = (prs_norm-mean(prs_norm))/sd(prs_norm))]

# Associations
## AF
charge_af_std <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge_std, data=eur_set)
prs_af_std <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_std, data=eur_set)
predict_af_std <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_std, data=eur_set)

# AF
# CHARGE
boot_charge <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='charge',data=eur_set,runs=200,size=nrow(eur_set))
## Values
cstat_charge <- c(summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge,data=eur_set))$concordance[1],
                  summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge,data=eur_set))$concordance[1]-1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge,data=eur_set))$concordance[2],
                  summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge,data=eur_set))$concordance[1]+1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge,data=eur_set))$concordance[2])
r2_charge <- c(cph(Surv(incd_af_5y.t,incd_af_5y) ~ charge,data=eur_set)$stats['R2'],
               cph(Surv(incd_af_5y.t,incd_af_5y) ~ charge,data=eur_set)$stats['R2']-1.96*sd(boot_charge),
               cph(Surv(incd_af_5y.t,incd_af_5y) ~ charge,data=eur_set)$stats['R2']+1.96*sd(boot_charge))
AIC(charge_af_std)

# PRS
boot_prs <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='prs_norm',data=eur_set,runs=200,size=nrow(eur_set))
## Values
mod_prs <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_std,data=eur_set)
cstat_prs <- c(summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_norm,data=eur_set))$concordance[1],
               summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_norm,data=eur_set))$concordance[1]-1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_norm,data=eur_set))$concordance[2],
               summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_norm,data=eur_set))$concordance[1]+1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_norm,data=eur_set))$concordance[2])
r2_prs <- c(cph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_norm,data=eur_set)$stats['R2'],
            cph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_norm,data=eur_set)$stats['R2']-1.96*sd(boot_prs),
            cph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_norm,data=eur_set)$stats['R2']+1.96*sd(boot_prs))
AIC(prs_af_std)

# Predict-AF
boot_predict <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='predict_af',data=eur_set,runs=200,size=nrow(eur_set))
## Values
mod_predict <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_std,data=eur_set)
cstat_predict <- c(summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=eur_set))$concordance[1],
                   summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=eur_set))$concordance[1]-1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=eur_set))$concordance[2],
                   summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=eur_set))$concordance[1]+1.96*summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=eur_set))$concordance[2])
r2_predict <- c(cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=eur_set)$stats['R2'],
                cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=eur_set)$stats['R2']-1.96*sd(boot_predict),
                cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=eur_set)$stats['R2']+1.96*sd(boot_predict))
AIC(predict_af_std)

############################################# STEP 14: Merits OAC at AF analysis
# Define subset
needs_oac <- af_set[needs_oac_visit0==1]

# Cuminc
age <- cuminc(data=needs_oac[age65_visit0==1],time='incd_af_5y.t',status='incd_af_5y')
charge <- cuminc(data=needs_oac[charge_pred5_above5==1],time='incd_af_5y.t',status='incd_af_5y')
predict <- cuminc(data=needs_oac[predict_af_pred5_above5==1],time='incd_af_5y.t',status='incd_af_5y')

# Binary models
mod_charge <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge_pred5_above5,data=needs_oac)
mod_predict <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_pred5_above5,data=needs_oac)
mod_age <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ age65_visit0,data=needs_oac)

############################################# STEP 15: ROCs and best threshold
# Plot ROCs
##### ALL
## Models
charge <- timeROC(T=af_set$incd_af_5y.t, delta=af_set$incd_af_5y,
                  marker=af_set$charge,cause=1,times=c(0,1,2,3,4,4.999))

predict <- timeROC(T=af_set$incd_af_5y.t, delta=af_set$incd_af_5y,
                   marker=af_set$predict_af,cause=1,times=c(0,1,2,3,4,4.999))

age <- timeROC(T=af_set$incd_af_5y.t, delta=af_set$incd_af_5y,
               marker=af_set$agevisit0,cause=1,times=c(0,1,2,3,4,4.999))

## Plotting
pdf(file='rocs.pdf',height=3,width=3,
    pointsize=3)
par(oma=c(3,3,1,1))
par(mar=c(4,4,1,1))
plot.new()

# Plot 1 by 1
plot(charge,4.999,add=T,col='#2b83ba',lwd=0.8)
par(new=TRUE)
plot(predict,4.999,add=T,col='#d7191c',lwd=0.8)
par(new=TRUE)
plot(age,4.999,add=T,col='#abdda4',lwd=0.8)

## Axes
axis(1,at=seq(0,1,0.1),cex.axis=1.4)
axis(2,at=seq(0,1,0.1),cex.axis=1.4)

## Labels
title(xlab='1-specificity',line=2.5,cex.lab=1.5)
title(ylab='sensitivity',line=2.5,cex.lab=1.5)

## Legend
legend(0.5,0.25,legend=c('Predict-AF','CHARGE-AF','Age'),col=c('#d7191c','#2b83ba','#abdda4'),
       lty=1,lwd=0.8,pch=1,bty='n',cex=1.2)
## Stop
dev.off()

##### AGE < 65
## Models
roc_charge_less65 <- timeROC(T=af_set$incd_af_before65_5y.t, delta=af_set$incd_af_before65_5y,
                             marker=af_set$charge,cause=1,times=c(0,1,2,3,4,4.999))

roc_predict_less65 <- timeROC(T=af_set$incd_af_before65_5y.t, delta=af_set$incd_af_before65_5y,
                              marker=af_set$predict_af,cause=1,times=c(0,1,2,3,4,4.999))

## Plotting
pdf(file='rocs_less65.pdf',height=3,width=3,
    pointsize=3)
par(oma=c(3,3,1,1))
par(mar=c(4,4,1,1))
plot.new()

# Plot 1 by 1
plot(roc_charge_less65,4.999,add=T,col='#2b83ba',lwd=0.8)
par(new=TRUE)
plot(roc_predict_less65,4.999,add=T,col='#d7191c',lwd=0.8)

## Axes
axis(1,at=seq(0,1,0.1),cex.axis=1.4)
axis(2,at=seq(0,1,0.1),cex.axis=1.4)

## Labels
title(xlab='1-specificity',line=2.5,cex.lab=1.5)
title(ylab='sensitivity',line=2.5,cex.lab=1.5)

## Legend
legend(0.5,0.25,legend=c('Predict-AF','CHARGE-AF'),col=c('#d7191c','#2b83ba'),
       lty=1,lwd=0.8,pch=1,bty='n',cex=1.2)
## Stop
dev.off()

# Find best thresholds (this assumes complete follow-up, so should be done on subset with complete follow-up at 5 years)
complete <- af_set[incd_af_5y==1 | (incd_af_5y==0 & incd_af_5y.t==5)]

predict_roc <- roc(response=complete$incd_af_5y,predictor=complete$predict_af,auc=T,ci=T)
predict_thres <- coords(roc=predict_roc,input="threshold",x='best',ret=c("threshold","sensitivity","specificity")) # best threshold is 21.67

charge_roc <- roc(response=complete$incd_af_5y,predictor=complete$charge,auc=T,ci=T)
charge_thres <- coords(roc=charge_roc,input="threshold",x='best',ret=c("threshold","sensitivity","specificity")) # best threshold is 12.14

age_roc <- roc(response=complete$incd_af_5y,predictor=complete$agevisit0,auc=T,ci=T)
age_thres <- coords(roc=age_roc,input="threshold",x='best',ret=c("threshold","sensitivity","specificity")) # best threshold is 60.69

# "Best Threshold"
complete[,':='(charge_pred5_best = ifelse(charge >= as.numeric(charge_thres[1]),1,0),
               predict_af_best = ifelse(predict_af >= as.numeric(predict_thres[1]),1,0))]

# Cuminc
age <- cuminc(data=complete[age65_visit0==1],time='incd_af_5y.t',status='incd_af_5y')
charge <- cuminc(data=complete[charge_pred5_best==1],time='incd_af_5y.t',status='incd_af_5y')
predict <- cuminc(data=complete[predict_af_best==1],time='incd_af_5y.t',status='incd_af_5y')

# Binary models
mod_charge <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge_pred5_best,data=complete)
mod_predict <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af_best,data=complete)
mod_age <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ age65_visit0,data=complete)

############################################# STEP 16: Calibration plots
####### Predict-AF
survmod<-with(af_set,Surv(incd_af_5y.t*365.25,incd_af_5y))
fit<-cph(survmod~predict_af,data=af_set,surv=TRUE,time.inc=5*365.25,u=5*365.25,x=T,y=T)
fit2 <- coxph(survmod~predict_af,data=af_set)
cal_predict <- c(fit2$coefficients[1],confint(fit2,"predict_af")[1],confint(fit2,"predict_af")[2]) # beta and 95%CI corresponds to calibration slope

# plot calibrations
cal.predict<-calibrate(fit,u=5*365.25,cmethod='hare',B=200) # ultimately want B=200 but takes a while, B=10 good enough to see how it looks

pdf('predict_cal_rms.pdf',height=3,width=3,pointsize=3)
par(oma=c(1,1,1,1))
par(oma=c(1,1,1,1))

col=paste(rgb(252,146,114,maxColorValue=255),sep="")

plot(cal.predict,scat1d.opts=list(frac=0.1,side=1),xlim=c(1,0.5),ylim=c(1,0.5),
     xaxt="n",yaxt="n", xlab="Predicted 5-year risk of AF", ylab="Proportion with AF at 5 years", subtitles=F, 
     par.corrected=list(col=col, lty=1, lwd=1.2), bty='n', cex.lab=1.25)
axis(1,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))
axis(2,at=seq(0,1,0.1),labels=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0), las=1)
mtext(text="Predict-AF", side=3, line=-0.5, adj=0.5, cex=1.2)
mtext(text="B", side=3, line=-1, at=1.4, cex=1.3, padj=-2)
mtext(text="Calibration slope", side=3, line=-2.5, at=0.9, cex=0.9)
mtext(text=paste0(format(round(cal_predict[1],2),nsmall=2),' (95%CI ',
                  format(round(cal_predict[2],2),nsmall=2),'-',
                  format(round(cal_predict[3],2),nsmall=2),')'), 
      side=3, line=-3.5, at=0.9, cex=0.9)

legend(0.7,0.9,c('optimal','observed','optimism-corrected'),col=c('darkgray','black',col),lty=1,lwd=1.5,bty='n',cex=0.75)
segments(1,1,0,0,col='darkgray',lty=1)

dev.off()

############################################# STEP 17: Proportionality assumption
# Make AF models
predict <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=af_set)
charge <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge,data=af_set)
prs <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs,data=af_set)

# Make stroke models
predict_stroke <- coxph(Surv(stroke_5y.t,stroke_5y) ~ predict_af,data=stroke_set)
charge_stroke <- coxph(Surv(stroke_5y.t,stroke_5y) ~ charge,data=stroke_set)
prs_stroke <- coxph(Surv(stroke_5y.t,stroke_5y) ~ prs,data=stroke_set)

# Standard
ph_predict <- cox.zph(predict)
ph_charge <- cox.zph(charge)
ph_prs <- cox.zph(prs)

# Stroke
ph_predict_stroke <- cox.zph(predict_stroke)
ph_charge_stroke <- cox.zph(charge_stroke)
ph_prs_stroke <- cox.zph(prs_stroke)

# Plots
plot(cox.zph(predict)[1])
abline(h=0, lty=2)
abline(h=coef(predict)[1],lty=3)

plot(cox.zph(charge)[1])
abline(h=0, lty=2)
abline(h=coef(charge)[1],lty=3)

plot(cox.zph(prs)[1])
abline(h=0, lty=2)
abline(h=coef(prs)[1],lty=3)

############################################# STEP 18: Decision curve using STDCA
# Needs DF
setDF(af_set)

# Needs predictions in 0-1 range
af_set$predict_af_pred5_fraction <- af_set$predict_af_pred5/100
af_set$charge_pred5_fraction <- af_set$charge_pred5/100
af_set$age_pred5_fraction <- af_set$age_pred5/100

# Generate data
stdca_std <- stdca(data=af_set,outcome='incd_af_5y',ttoutcome='incd_af_5y.t',
                   timepoint=4.999,predictors=c('predict_af_pred5_fraction','charge_pred5_fraction','age_pred5_fraction'),
                   xstart=0,xstop=0.15,xby=0.005,probability=c(TRUE,TRUE,TRUE),graph=TRUE)

# Plot standard decision curve
pdf('decision_stdca.pdf',height=4,width=4.5,pointsize=5)
par(mar=c(4.5,3,1,1),oma=c(4.5,3,1,1))

y1 <- stdca_std$net.benefit[,'all'][!is.na(stdca_std$net.benefit[,'all'])]/max(stdca_std$net.benefit[,'all'])
y2 <- stdca_std$net.benefit[,'none'][!is.na(stdca_std$net.benefit[,'none'])]/max(stdca_std$net.benefit[,'all'])
y3 <- stdca_std$net.benefit[,'predict_af_pred5_fraction'][!is.na(stdca_std$net.benefit[,'predict_af_pred5_fraction'])]/max(stdca_std$net.benefit[,'all'])
y4 <- stdca_std$net.benefit[,'charge_pred5_fraction'][!is.na(stdca_std$net.benefit[,'charge_pred5_fraction'])]/max(stdca_std$net.benefit[,'all'])
y5 <- stdca_std$net.benefit[,'age_pred5_fraction'][!is.na(stdca_std$net.benefit[,'age_pred5_fraction'])]/max(stdca_std$net.benefit[,'all'])

x1 <- stdca_std$net.benefit[,'threshold'][!is.na(stdca_std$net.benefit[,'all'])]
x2 <- stdca_std$net.benefit[,'threshold'][!is.na(stdca_std$net.benefit[,'none'])]
x3 <- stdca_std$net.benefit[,'threshold'][!is.na(stdca_std$net.benefit[,'predict_af_pred5_fraction'])]
x4 <- stdca_std$net.benefit[,'threshold'][!is.na(stdca_std$net.benefit[,'charge_pred5_fraction'])]
x5 <- stdca_std$net.benefit[,'threshold'][!is.na(stdca_std$net.benefit[,'age_pred5_fraction'])]

plot(x1,y1,type='l',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',
     ylim=c(-0.02,1),xlim=c(0,0.15),col='darkgray',lwd=1.2)
axis(2,seq(0,1,0.1),las=2,cex.axis=1.3)
axis(1,seq(0,0.15,0.025),pos=-0.07,cex.axis=1.3)
axis(1,at=c(1/101,1/26,1/11,1/7),labels=c('1:100','1:25','1:10','1:6'),pos=-0.24,
     cex.axis=1.3)

mtext('Threshold',1,line=2.8,cex=1.5)
mtext('Cost:Benefit ratio',1,line=7.5,cex=1.5)
mtext('Standardized Net Benefit',2,line=3.3,cex=1.5)

legend(x=0.09,y=1,c('Predict-AF','CHARGE-AF','Age','All','None'),
       col=c('#9e0142','#66c2a5','#5e4fa2','darkgray','black'),
       bty='n',lty=1,cex=1.4,lwd=1.5)

par(new=TRUE)
plot(x3,y3,type='l',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',ylim=c(-0.02,1),
     xlim=c(0,0.15),col='#9e0142',lwd=1.2)

par(new=TRUE)
plot(x4,y4,type='l',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',ylim=c(-0.02,1),
     xlim=c(0,0.15),col='#66c2a5',lwd=1.2)

par(new=TRUE)
plot(x5,y5,type='l',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',ylim=c(-0.02,1),
     xlim=c(0,0.15),col='#5e4fa2',lwd=1.2)

par(new=TRUE)
plot(x2,y2,type='l',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',ylim=c(-0.02,1),
     xlim=c(0,0.15),col='black',lwd=1.2)

dev.off()

############################################# STEP 19: Interventions avoided 
# Generate data
stdca <- stdca(data=af_set,outcome='incd_af_5y',ttoutcome='incd_af_5y.t',
               timepoint=4.999,predictors=c('predict_af_pred5_fraction','charge_pred5_fraction','age_pred5_fraction'),
               xstart=0,xstop=0.10,xby=0.005,probability=c(TRUE,TRUE,TRUE),intervention=TRUE,graph=FALSE)
stdca2 <- stdca(data=af_set[af_set$age65_visit0==1,],outcome='incd_af_5y',ttoutcome='incd_af_5y.t',
                timepoint=4.999,predictors=c('predict_af_pred5_fraction','charge_pred5_fraction'),
                xstart=0,xstop=0.15,xby=0.005,probability=c(TRUE,TRUE),intervention=TRUE,graph=FALSE)

# Plot all strategies in all population
pdf('intervention_avoided.pdf',height=3,width=3,pointsize=3)
par(mar=c(2,3,1,1),oma=c(2,3,1,1))

y1 <- stdca$interventions.avoided[,'predict_af_pred5_fraction'][!is.na(stdca$interventions.avoided[,'predict_af_pred5_fraction'])]
y2 <- stdca$interventions.avoided[,'charge_pred5_fraction'][!is.na(stdca$interventions.avoided[,'charge_pred5_fraction'])]
y3 <- stdca$interventions.avoided[,'age_pred5_fraction'][!is.na(stdca$interventions.avoided[,'age_pred5_fraction'])]

x1 <- stdca$interventions.avoided[,'threshold'][!is.na(stdca$interventions.avoided[,'predict_af_pred5_fraction'])]
x2 <- stdca$interventions.avoided[,'threshold'][!is.na(stdca$interventions.avoided[,'charge_pred5_fraction'])]
x3 <- stdca$interventions.avoided[,'threshold'][!is.na(stdca$interventions.avoided[,'age_pred5_fraction'])]

plot(x1,y1,type='l',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',
     ylim=c(0,90),xlim=c(0,0.08),col='#9e0142',lwd=1.5)
axis(2,seq(0,90,10),las=2,cex.axis=1.2)
axis(1,seq(0,0.08,0.01),pos=0,cex.axis=1.2)

mtext('Threshold',1,line=2,cex=1.5)
mtext('Interventions avoided per 100 persons',2,line=3,cex=1.5)

legend(x=0.045,y=22,c('Predict-AF','CHARGE-AF','Age'),col=c('#9e0142','#abdda4','#5e4fa2'),
       bty='n',lty=1,cex=1.2)

par(new=TRUE)
plot(x2,y2,type='l',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',ylim=c(0,90),
     xlim=c(0,0.08),col='#abdda4',lwd=1.5)

par(new=TRUE)
plot(x3,y3,type='l',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',ylim=c(0,90),
     xlim=c(0,0.08),col='#5e4fa2',lwd=1.5)

dev.off()

# Plot risk in age >= 65 population
pdf('intervention_avoided_65.pdf',height=3,width=3,pointsize=3)
par(mar=c(2,3,1,1),oma=c(2,3,1,1))

y1 <- stdca2$interventions.avoided[,'predict_af_pred5_fraction'][c(!is.na(stdca$interventions.avoided[,'predict_af_pred5_fraction'])
                                                                   & (stdca$interventions.avoided[,'predict_af_pred5_fraction'] > 0))]
y2 <- stdca2$interventions.avoided[,'charge_pred5_fraction'][c(!is.na(stdca$interventions.avoided[,'charge_pred5_fraction'])
                                                               & (stdca$interventions.avoided[,'charge_pred5_fraction'] > 0))]

x1 <- stdca2$interventions.avoided[,'threshold'][c(!is.na(stdca$interventions.avoided[,'predict_af_pred5_fraction'])
                                                   & (stdca$interventions.avoided[,'predict_af_pred5_fraction'] > 0))]
x2 <- stdca2$interventions.avoided[,'threshold'][c(!is.na(stdca$interventions.avoided[,'charge_pred5_fraction'])
                                                   & (stdca$interventions.avoided[,'charge_pred5_fraction'] > 0))]

plot(x1,y1,type='l',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',
     ylim=c(0,70),xlim=c(0,0.10),col='#9e0142',lwd=1.5)
axis(2,seq(0,70,10),las=2,cex.axis=1.2)
axis(1,seq(0,0.10,0.01),pos=0,cex.axis=1.2)

mtext('Threshold',1,line=2,cex=1.5)
mtext('Interventions avoided per 100 persons',2,line=3,cex=1.5)

legend(x=0.055,y=15,c('Predict-AF','CHARGE-AF'),col=c('#9e0142','#abdda4'),
       bty='n',lty=1)

par(new=TRUE)
plot(x2,y2,type='l',xaxt='n',yaxt='n',xlab='',ylab='',bty='n',ylim=c(0,70),
     xlim=c(0,0.10),col='#abdda4',lwd=1.5)

dev.off()

setDT(af_set)

############################################# STEP 20: NNS
# 1-year AF for all pop
age <- cuminc(data=af_set[age65_visit0==1],time='incd_af_1y.t',status='incd_af_1y')
charge <- cuminc(data=af_set[charge_pred5_above5==1],time='incd_af_1y.t',status='incd_af_1y')
predict <- cuminc(data=af_set[predict_af_pred5_above5==1],time='incd_af_1y.t',status='incd_af_1y')

nns_age <- c(ceiling(1/(age[,c(3,5,4)]/100)/0.9/0.3),ceiling(1/(age[,c(3,5,4)]/100)/0.9/0.35),ceiling(1/(age[,c(3,5,4)]/100)/0.9/0.40))
nns_charge <- c(ceiling(1/(charge[,c(3,5,4)]/100)/0.9/0.3),ceiling(1/(charge[,c(3,5,4)]/100)/0.9/0.35),ceiling(1/(charge[,c(3,5,4)]/100)/0.9/0.40))
nns_predict <- c(ceiling(1/(predict[,c(3,5,4)]/100)/0.9/0.3),ceiling(1/(predict[,c(3,5,4)]/100)/0.9/0.35),ceiling(1/(predict[,c(3,5,4)]/100)/0.9/0.40))

# 1-year AF for age < 65
charge_before65 <- cuminc(data=af_set[charge_pred5_above5==1 & !is.na(incd_af_before65_1y.t)],time='incd_af_before65_1y.t',status='incd_af_before65_1y')
predict_before65 <- cuminc(data=af_set[predict_af_pred5_above5==1 & !is.na(incd_af_before65_1y.t)],time='incd_af_before65_1y.t',status='incd_af_before65_1y')

nns_charge_before65 <- c(ceiling(1/(charge_before65[,c(3,5,4)]/100)/0.9/0.3),ceiling(1/(charge_before65[,c(3,5,4)]/100)/0.9/0.35),ceiling(1/(charge_before65[,c(3,5,4)]/100)/0.9/0.40))
nns_predict_before65 <- c(ceiling(1/(predict_before65[,c(3,5,4)]/100)/0.9/0.3),ceiling(1/(predict_before65[,c(3,5,4)]/100)/0.9/0.35),ceiling(1/(predict_before65[,c(3,5,4)]/100)/0.9/0.40))

# 1-year AF for age >= 65
charge_after65 <- cuminc(data=af_set[charge_pred5_above5==1 & !is.na(incd_af_after65_1y.t)],time='incd_af_after65_1y.t',status='incd_af_after65_1y')
predict_after65 <- cuminc(data=af_set[predict_af_pred5_above5==1 & !is.na(incd_af_after65_1y.t)],time='incd_af_after65_1y.t',status='incd_af_after65_1y')

nns_charge_after65 <- c(ceiling(1/(charge_after65[,c(3,5,4)]/100)/0.9/0.3),ceiling(1/(charge_after65[,c(3,5,4)]/100)/0.9/0.35),ceiling(1/(charge_after65[,c(3,5,4)]/100)/0.9/0.40))
nns_predict_after65 <- c(ceiling(1/(predict_after65[,c(3,5,4)]/100)/0.9/0.3),ceiling(1/(predict_after65[,c(3,5,4)]/100)/0.9/0.35),ceiling(1/(predict_after65[,c(3,5,4)]/100)/0.9/0.40))

##### PLOTS
### CHARGE
af1_charge <- prodlim(Hist(incd_af_1y.t,incd_af_1y)~1,data=af_set[charge_pred5_above5==1])
af1_charge_sub1 <- af1_charge; af1_charge_sub1$surv <- 1-((1-af1_charge_sub1$surv)*0.30)
af1_charge_sub2 <- af1_charge; af1_charge_sub2$surv <- 1-((1-af1_charge_sub2$surv)*0.35)
af1_charge_sub3 <- af1_charge; af1_charge_sub3$surv <- 1-((1-af1_charge_sub3$surv)*0.40)

CairoPDF(file='af_charge_nns.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(5,6,3,1),mar=c(5,6,3,1))
plot(af1_charge,"cuminc",ylim=c(0,0.012),
     axis2.at=seq(0,0.012,0.001),axis2.las=2,lwd=1.5,background=F,
     atrisk=F,col=c("#253494"),atrisk.col='black',confint=FALSE,axis1.at=seq(0,1,0.25),
     legend.x=0,legend.y=0.012,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     xlab='',ylab='')
mtext("Cumulative incidence of AF (%)",side=2,line=9,at=0.006,cex=2.5)
mtext("Years",side=1, line=5,cex=2.5)
mtext("CHARGE-AF 5-yr risk \u2265 5%",side=3, line=1.5,cex=2.5,at=0.42)
legend(0,0.012,c('Observed AF','Subclinical AF (rate 30%)','Subclinical AF (rate 35%)','Subclinical AF (rate 40%)'),
       col=c('#253494','#66c2a5','#1a9850','#f46d43'),bty='n',cex=1.5,lty=c(1,2,4,6),lwd=1.5,seg.len=3)
plot(af1_charge_sub1,"cuminc",add=TRUE,confint=FALSE,atrisk=F,lwd=1.5,lty=3,col='#f46d43')
plot(af1_charge_sub2,"cuminc",add=TRUE,confint=FALSE,atrisk=F,lwd=1.5,lty=2,col='#1a9850')
plot(af1_charge_sub3,"cuminc",add=TRUE,confint=FALSE,atrisk=F,lwd=1.5,lty=2,col='#d53e4f')
dev.off()

### PREDICT
af1_predict <- prodlim(Hist(incd_af_1y.t,incd_af_1y)~1,data=af_set[predict_af_pred5_above5==1])
af1_predict_sub1 <- af1_predict; af1_predict_sub1$surv <- 1-((1-af1_predict_sub1$surv)*0.30)
af1_predict_sub2 <- af1_predict; af1_predict_sub2$surv <- 1-((1-af1_predict_sub2$surv)*0.35)
af1_predict_sub3 <- af1_predict; af1_predict_sub3$surv <- 1-((1-af1_predict_sub3$surv)*0.40)

CairoPDF(file='af_predict_nns.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(5,6,3,1),mar=c(5,6,3,1))
plot(af1_predict,"cuminc",ylim=c(0,0.012),
     axis2.at=seq(0,0.012,0.001),axis2.las=2,lwd=1.5,background=F,
     atrisk=F,col=c("#253494"),atrisk.col='black',confint=FALSE,axis1.at=seq(0,1,0.25),
     legend.x=0,legend.y=0.012,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     xlab='',ylab='')
mtext("Cumulative incidence of AF (%)",side=2,line=9,at=0.006,cex=2.5)
mtext("Years",side=1, line=5,cex=2.5)
mtext("Predict-AF 5-yr risk \u2265 5%",side=3, line=1.5,cex=2.5,at=0.42)
legend(0,0.012,c('Observed AF','Subclinical AF (rate 30%)','Subclinical AF (rate 35%)','Subclinical AF (rate 40%)'),
       col=c('#253494','#66c2a5','#1a9850','#f46d43'),bty='n',cex=1.5,lty=c(1,2,4,6),lwd=1.5,seg.len=3)
plot(af1_predict_sub1,"cuminc",add=TRUE,confint=FALSE,atrisk=F,lwd=1.5,lty=3,col='#f46d43')
plot(af1_predict_sub2,"cuminc",add=TRUE,confint=FALSE,atrisk=F,lwd=1.5,lty=3,col='#1a9850')
plot(af1_predict_sub3,"cuminc",add=TRUE,confint=FALSE,atrisk=F,lwd=1.5,lty=2,col='#d53e4f')
dev.off()

### AGE
af1_age <- prodlim(Hist(incd_af_1y.t,incd_af_1y)~1,data=af_set[age65_visit0==1])
af1_age_sub1 <- af1_age; af1_age_sub1$surv <- 1-((1-af1_age_sub1$surv)*0.30)
af1_age_sub2 <- af1_age; af1_age_sub2$surv <- 1-((1-af1_age_sub2$surv)*0.35)
af1_age_sub3 <- af1_age; af1_age_sub3$surv <- 1-((1-af1_age_sub3$surv)*0.40)

CairoPDF(file='af_age_nns.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(5,6,3,1),mar=c(5,6,3,1))
plot(af1_age,"cuminc",ylim=c(0,0.012),
     axis2.at=seq(0,0.012,0.001),axis2.las=2,lwd=1.5,background=F,
     atrisk=F,col=c("#253494"),atrisk.col='black',confint=FALSE,axis1.at=seq(0,1,0.25),
     legend.x=0,legend.y=0.012,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     xlab='',ylab='')
mtext("Cumulative incidence of AF (%)",side=2,line=9,at=0.006,cex=2.5)
mtext("Years",side=1, line=5,cex=2.5)
mtext("Age \u2265 65",side=3, line=1.5,cex=2.5,at=0.42)
legend(0,0.012,c('Observed AF','Subclinical AF (rate 30%)','Subclinical AF (rate 35%)','Subclinical AF (rate 40%)'),
       col=c('#253494','#d53e4f','#1a9850','#f46d43'),bty='n',cex=1.5,lty=c(1,2,4,6),lwd=1.5,seg.len=3)
plot(af1_age_sub1,"cuminc",add=TRUE,confint=FALSE,atrisk=F,lwd=1.5,lty=3,col='#f46d43')
plot(af1_age_sub2,"cuminc",add=TRUE,confint=FALSE,atrisk=F,lwd=1.5,lty=2,col='#1a9850')
plot(af1_age_sub3,"cuminc",add=TRUE,confint=FALSE,atrisk=F,lwd=1.5,lty=2,col='#d53e4f')
dev.off()

############################################# STEP 21: Bar plots
########### AF
## Overall
ukbb <- c(3.61,8.15,7.63); ukbb_lb <- c(3.48,7.67,7.30); ukbb_ub <- c(3.73,8.62,7.97)
fg <- c(4.90,6.06,6.45); fg_lb <- c(3.85,4.77,5.20); fg_ub <- c(5.93,7.33,7.69)
geisinger <- c(8.56,7.95,7.65); geisinger_lb <- c(7.94,7.45,7.19); geisinger_ub <- c(9.17,8.45,8.12)
fhs <- c(11.06,10.71,11.09); fhs_lb <- c(9.29,9.08,9.42); fhs_ub <- c(12.80,12.31,12.72)

CairoPDF(file='bar_plot_cohort_all.pdf',height=3,width=5.5,
         pointsize=5)
par(oma=c(2.5,3,2.5,1),mar=c(2.5,3,2.5,1))
coords <- barplot(c(ukbb,NA,fg,NA,geisinger,NA,fhs),
                  col=rep(c('darkgray','#92c5de','#2166ac',NA),4),
                  yaxt='n',bty='n',border = NA,ylim=c(0,15))
par(xpd=TRUE)
legend(14,16.5,c('Age \u2265 65','CHARGE-AF \u2265 5%','Predict-AF \u2265 5%'),bty='n',lty=1,lwd=2,
       col=c('darkgray','#92c5de','#2166ac'),cex=1.3)

# Y-axis
axis(2,at=seq(0,16,2.5),las=2,cex.axis=1.3,labels=c('0','2.5','5','7.5','10','12.5','15'))
mtext('5-year cumulative incidence of AF (%)',2,line=4,at=7.5,cex=1.5)

# Error bars and counts
lb_all <- c(ukbb_lb,fg_lb,geisinger_lb,fhs_lb)
ub_all <- c(ukbb_ub,fg_ub,geisinger_ub,fhs_ub)

# Looping function that calls barplot coords
n<-1
events <- c('3148','1039','1835',
            '81','81','97',
            '681','889','953',
            '134','152','155')
at_risk <- c('88572','13131','24573',
             '1695','1381','1555',
             '7958','11183','12453',
             '1305','1502','1479')
for (i in coords[c(1:3,5:7,9:11,13:15)]){
  arrows(x0=i,y0=lb_all[n],x1=i,y1=ub_all[n],angle=90,length=0.03,lwd=0.8,code=3)
  mtext(events[n],1,at=i,line=0.8,cex=1.1)
  mtext(at_risk[n],1,at=i,line=1.7,cex=1.1)
  n <- n+1}

mtext('UK Biobank',1,at=coords[2],line=3.2,cex=1.5)
mtext('FINRISK',1,at=coords[6],line=3.2,cex=1.5)
mtext('Geisinger',1,at=coords[10],line=3.2,cex=1.5)
mtext('Framingham',1,at=coords[14],line=3.2,cex=1.5)
mtext('N high risk',1,at=-0.9,line=1.7,cex=1.1)
mtext('N AF events',1,at=-0.9,line=0.8,cex=1.1)

dev.off()

## Age < 65
ukbb <- c(0,11.17,8.41); ukbb_lb <- c(0,7.85,6.91); ukbb_ub <- c(0,14.37,9.88)
fg <- c(0,6.96,6.61); fg_lb <- c(0,2.57,3.38); fg_ub <- c(0,11.16,9.73)
geisinger <- c(0,5.54,5.65); geisinger_lb <- c(0,4.68,4.91); geisinger_ub <- c(0,6.39,6.38)
fhs <- c(0,7.63,9.09); fhs_lb <- c(0,2.87,4.30); fhs_ub <- c(0,12.16,13.64)

CairoPDF(file='bar_plot_cohort_less65.pdf',height=3,width=5.5,
         pointsize=5)
par(oma=c(2.5,3,2.5,1),mar=c(2.5,3,2.5,1))
coords <- barplot(c(ukbb,NA,fg,NA,geisinger,NA,fhs),col=rep(c('darkgray','#92c5de','#2166ac',NA),4),
                  yaxt='n',bty='n',border = NA,ylim=c(0,15))
par(xpd=TRUE)
legend(14,16.5,c('Age \u2265 65','CHARGE-AF \u2265 5%','Predict-AF \u2265 5%'),bty='n',lty=1,lwd=2,
       col=c('darkgray','#92c5de','#2166ac'),cex=1.3)

# Y-axis
axis(2,at=seq(0,16,2.5),las=2,cex.axis=1.3,labels=c('0','2.5','5','7.5','10','12.5','15'))
mtext('5-year cumulative incidence of AF (%)',2,line=4,at=7.5,cex=1.5)

# Error bars and counts
lb_all <- c(ukbb_lb,fg_lb,geisinger_lb,fhs_lb)
ub_all <- c(ukbb_ub,fg_ub,geisinger_ub,fhs_ub)

# Looping function that calls barplot coords
n<-1
events <- c('0','68','178',
            '0','11','17',
            '0','165','225',
            '0','11','15')
at_risk <- c('0','1952','5553',
             '0','281','457',
             '0','3915','5094',
             '0','271','275')
for (i in coords[c(1:3,5:7,9:11,13:15)]){
  arrows(x0=i,y0=lb_all[n],x1=i,y1=ub_all[n],angle=90,length=0.03,lwd=0.8,code=3)
  mtext(events[n],1,at=i,line=0.8,cex=1.1)
  mtext(at_risk[n],1,at=i,line=1.7,cex=1.1)
  n <- n+1}

mtext('UK Biobank',1,at=coords[2],line=3.2,cex=1.5)
mtext('FINRISK',1,at=coords[6],line=3.2,cex=1.5)
mtext('Geisinger',1,at=coords[10],line=3.2,cex=1.5)
mtext('Framingham',1,at=coords[14],line=3.2,cex=1.5)
mtext('N high risk',1,at=-0.9,line=1.7,cex=1.1)
mtext('N AF events',1,at=-0.9,line=0.8,cex=1.1)

dev.off()

## Age >= 65
ukbb <- c(3.61,7.86,7.53); ukbb_lb <- c(3.48,7.35,7.15); ukbb_ub <- c(3.73,8.36,7.91)
fg <- c(4.90,5.81,6.38); fg_lb <- c(3.85,4.40,4.90); fg_ub <- c(5.93,7.21,7.83)
geisinger <- c(8.56,9.15,9.04); geisinger_lb <- c(7.94,8.49,8.38); geisinger_ub <- c(9.17,9.81,9.69)
fhs <- c(11.06,11.25,11.44); fhs_lb <- c(9.29,9.47,9.63); fhs_ub <- c(12.80,12.99,13.22)

CairoPDF(file='bar_plot_cohort_over65.pdf',height=3,width=5.5,
         pointsize=5)
par(oma=c(2.5,3,2.5,1),mar=c(2.5,3,2.5,1))
coords <- barplot(c(ukbb,NA,fg,NA,geisinger,NA,fhs),col=rep(c('darkgray','#92c5de','#2166ac',NA),3),
                  yaxt='n',bty='n',border = NA,ylim=c(0,15))
par(xpd=TRUE)
legend(14,16.5,c('Age \u2265 65','CHARGE-AF \u2265 5%','Predict-AF \u2265 5%'),bty='n',lty=1,lwd=2,
       col=c('darkgray','#92c5de','#2166ac'),cex=1.3)

# Y-axis
axis(2,at=seq(0,16,2.5),las=2,cex.axis=1.3,labels=c('0','2.5','5','7.5','10','12.5','15'))
mtext('5-year cumulative incidence of AF (%)',2,line=4,at=7.5,cex=1.5)

# Error bars and counts
lb_all <- c(ukbb_lb,fg_lb,geisinger_lb,fhs_lb)
ub_all <- c(ukbb_ub,fg_ub,geisinger_ub,fhs_ub)

# Looping function that calls barplot coords
n<-1
events <- c('3148','854','1400',
            '81','62','68',
            '681','665','665',
            '134','141','140')
at_risk <- c('88572','11179','19020',
             '1695','1100','1098',
             '7958','7268','7359',
             '1305','1417','1386')
for (i in coords[c(1:3,5:7,9:11,13:15)]){
  arrows(x0=i,y0=lb_all[n],x1=i,y1=ub_all[n],angle=90,length=0.03,lwd=0.8,code=3)
  mtext(events[n],1,at=i,line=0.8,cex=1.1)
  mtext(at_risk[n],1,at=i,line=1.7,cex=1.1)
  n <- n+1}

mtext('UK Biobank',1,at=coords[2],line=3.2,cex=1.5)
mtext('FINRISK',1,at=coords[6],line=3.2,cex=1.5)
mtext('Geisinger',1,at=coords[10],line=3.2,cex=1.5)
mtext('Framingham',1,at=coords[14],line=3.2,cex=1.5)
mtext('N high risk',1,at=-0.9,line=1.7,cex=1.1)
mtext('N AF events',1,at=-0.9,line=0.8,cex=1.1)

dev.off()

########### Stroke
## Overall
ukbb <- c(1.13,2.04,1.86); ukbb_lb <- c(1.06,1.79,1.68); ukbb_ub <- c(1.20,2.29,2.03)
fg <- c(4.90,6.06,6.45); fg_lb <- c(3.85,4.77,5.20); fg_ub <- c(5.93,7.33,7.69)
geisinger <- c(1.99,2.00,1.89); geisinger_lb <- c(1.68,1.74,1.65); geisinger_ub <- c(2.29,2.26,2.13)
fhs <- c(3.09,3.27,3.25); fhs_lb <- c(2.11,2.33,2.30); fhs_ub <- c(4.07,4.21,4.19)

CairoPDF(file='bar_plot_cohort_stroke.pdf',height=3,width=5.5,
         pointsize=5)
par(oma=c(2.5,3,2.5,1),mar=c(2.5,3,2.5,1))
coords <- barplot(c(ukbb,NA,fg,NA,geisinger,NA,fhs),
                  col=rep(c('darkgray','#92c5de','#2166ac',NA),4),
                  yaxt='n',bty='n',border = NA,ylim=c(0,10))
par(xpd=TRUE)
legend(14,11,c('Age \u2265 65','CHARGE-AF \u2265 5%','Predict-AF \u2265 5%'),bty='n',lty=1,lwd=2,
       col=c('darkgray','#92c5de','#2166ac'),cex=1.3)

# Y-axis
axis(2,at=seq(0,10,1),las=2,cex.axis=1.3,labels=c('0','1','2','3','4','5','6','7','8','9','10'))
mtext('5-year cumulative incidence of stroke (%)',2,line=4,at=5,cex=1.5)

# Error bars and counts
lb_all <- c(ukbb_lb,fg_lb,geisinger_lb,fhs_lb)
ub_all <- c(ukbb_ub,fg_ub,geisinger_ub,fhs_ub)

# Looping function that calls barplot coords
n<-1
events <- c('957','247','427',
            '84','75','85',
            '157','222','234',
            '37','45','44')
at_risk <- c('86284','12472','23544',
             '1639','1325','1492',
             '7911','11115','12379',
             '1264','1455','1432')
for (i in coords[c(1:3,5:7,9:11,13:15)]){
  arrows(x0=i,y0=lb_all[n],x1=i,y1=ub_all[n],angle=90,length=0.03,lwd=0.8,code=3)
  mtext(events[n],1,at=i,line=0.8,cex=1.1)
  mtext(at_risk[n],1,at=i,line=1.7,cex=1.1)
  n <- n+1}

mtext('UK Biobank',1,at=coords[2],line=3.2,cex=1.5)
mtext('FINRISK',1,at=coords[6],line=3.2,cex=1.5)
mtext('Geisinger',1,at=coords[10],line=3.2,cex=1.5)
mtext('Framingham',1,at=coords[14],line=3.2,cex=1.5)
mtext('N high risk',1,at=-0.9,line=1.7,cex=1.1)
mtext('N stroke events',1,at=-0.9,line=0.8,cex=1.1)

dev.off()

## Age < 65
ukbb <- c(0,3.74,2.08); ukbb_lb <- c(0,1.82,1.26); ukbb_ub <- c(0,5.62,2.89)
fg <- c(0,5.43,4.25); fg_lb <- c(0,1.35,1.74); fg_ub <- c(0,9.35,6.70)
geisinger <- c(0,1.77,1.57); geisinger_lb <- c(0,1.29,1.18); geisinger_ub <- c(0,2.24,1.95)
fhs <- c(0,4.57,2.79); fhs_lb <- c(0,0.55,0); fhs_ub <- c(0,8.43,5.63)

CairoPDF(file='bar_plot_cohort_less65_stroke.pdf',height=3,width=5.5,
         pointsize=5)
par(oma=c(2.5,3,2.5,1),mar=c(2.5,3,2.5,1))
coords <- barplot(c(ukbb,NA,fg,NA,geisinger,NA,fhs),col=rep(c('darkgray','#92c5de','#2166ac',NA),4),
                  yaxt='n',bty='n',border = NA,ylim=c(0,10))
par(xpd=TRUE)
legend(14,11,c('Age \u2265 65','CHARGE-AF \u2265 5%','Predict-AF \u2265 5%'),bty='n',lty=1,lwd=2,
       col=c('darkgray','#92c5de','#2166ac'),cex=1.3)

# Y-axis
axis(2,at=seq(0,10,1),las=2,cex.axis=1.3,labels=c('0','1','2','3','4','5','6','7','8','9','10'))
mtext('5-year cumulative incidence of stroke (%)',2,line=4,at=5,cex=1.5)

# Error bars and counts
lb_all <- c(ukbb_lb,fg_lb,geisinger_lb,fhs_lb)
ub_all <- c(ukbb_ub,fg_ub,geisinger_ub,fhs_ub)

# Looping function that calls barplot coords
n<-1
events <- c('0','22','38',
            '0','8','12',
            '0','51','60',
            '0','4','6')
at_risk <- c('0','1839','5318',
             '0','264','432',
             '0','3894','5067',
             '0','265','268')
for (i in coords[c(1:3,5:7,9:11,13:15)]){
  arrows(x0=i,y0=lb_all[n],x1=i,y1=ub_all[n],angle=90,length=0.03,lwd=0.8,code=3)
  mtext(events[n],1,at=i,line=0.8,cex=1.1)
  mtext(at_risk[n],1,at=i,line=1.7,cex=1.1)
  n <- n+1}

mtext('UK Biobank',1,at=coords[2],line=3.2,cex=1.5)
mtext('FINRISK',1,at=coords[6],line=3.2,cex=1.5)
mtext('Geisinger',1,at=coords[10],line=3.2,cex=1.5)
mtext('Framingham',1,at=coords[14],line=3.2,cex=1.5)
mtext('N high risk',1,at=-0.9,line=1.7,cex=1.1)
mtext('N stroke events',1,at=-0.9,line=0.8,cex=1.1)

dev.off()

## Age >= 65
ukbb <- c(1.13,1.94,1.90); ukbb_lb <- c(1.05,1.68,1.70); ukbb_ub <- c(1.20,2.21,2.11)
fg <- c(5.28,6.15,6.25); fg_lb <- c(4.18,4.67,4.75); fg_ub <- c(6.38,7.61,7.72)
geisinger <- c(1.99,2.12,2.09); geisinger_lb <- c(1.68,1.79,1.76); geisinger_ub <- c(2.29,2.45,2.42)
fhs <- c(3.09,3.30,3.38); fhs_lb <- c(2.11,2.25,2.34); fhs_ub <- c(4.07,4.34,4.44)

CairoPDF(file='bar_plot_cohort_over65_stroke.pdf',height=3,width=5.5,
         pointsize=5)
par(oma=c(2.5,3,2.5,1),mar=c(2.5,3,2.5,1))
coords <- barplot(c(ukbb,NA,fg,NA,geisinger,NA,fhs),col=rep(c('darkgray','#92c5de','#2166ac',NA),3),
                  yaxt='n',bty='n',border = NA,ylim=c(0,10))
par(xpd=TRUE)
legend(14,11,c('Age \u2265 65','CHARGE-AF \u2265 5%','Predict-AF \u2265 5%'),bty='n',lty=1,lwd=2,
       col=c('darkgray','#92c5de','#2166ac'),cex=1.3)

# Y-axis
axis(2,at=seq(0,10,1),las=2,cex.axis=1.3,labels=c('0','1','2','3','4','5','6','7','8','9','10'))
mtext('5-year cumulative incidence of stroke (%)',2,line=4,at=5,cex=1.5)

# Error bars and counts
lb_all <- c(ukbb_lb,fg_lb,geisinger_lb,fhs_lb)
ub_all <- c(ukbb_ub,fg_ub,geisinger_ub,fhs_ub)

# Looping function that calls barplot coords
n<-1
events <- c('957','201','339',
            '84','63','64',
            '157','156','156',
            '37','37','37')
at_risk <- c('86284','10633','18226',
             '1639','1062','1060',
             '7911','7221','7312',
             '1264','1190','1164')
for (i in coords[c(1:3,5:7,9:11,13:15)]){
  arrows(x0=i,y0=lb_all[n],x1=i,y1=ub_all[n],angle=90,length=0.03,lwd=0.8,code=3)
  mtext(events[n],1,at=i,line=0.8,cex=1.1)
  mtext(at_risk[n],1,at=i,line=1.7,cex=1.1)
  n <- n+1}

mtext('UK Biobank',1,at=coords[2],line=3.2,cex=1.5)
mtext('FINRISK',1,at=coords[6],line=3.2,cex=1.5)
mtext('Geisinger',1,at=coords[10],line=3.2,cex=1.5)
mtext('Framingham',1,at=coords[14],line=3.2,cex=1.5)
mtext('N high risk',1,at=-0.9,line=1.7,cex=1.1)
mtext('N stroke events',1,at=-0.9,line=0.8,cex=1.1)

dev.off()

################################################ Step 22: Targeted tallies
### Over 65 but low risk
over65_predict_low <- af_set[age65_visit0==1 & predict_af_pred5_above5==0]
nrow(over65_predict_low); nrow(over65_predict_low[incd_af_after65_5y==1])
count(over65_predict_low[incd_af_after65_5y==1]$needs_oac_at_af)

### Less 65 but high risk
less65_predict_high <- af_set[age65_visit0==0 & predict_af_pred5_above5==1]
nrow(less65_predict_high); nrow(less65_predict_high[incd_af_5y==1])
count(less65_predict_high[incd_af_before65_5y==1]$needs_oac_at_af)

################################################ Step 23: Density plots in age subgroup
### Density Plots
over65 <- af_set[age65_visit0==1]

# Create separate data.tables for cases/controls
gp_incident_af <- over65[incd_af_5y==1]
gp_noaf <- over65[incd_af_5y==0]

# CHARGE VS Predict-AF (controls)
x <- list(v1=gp_noaf$charge_pred5,v2=gp_noaf$predict_af_pred5)
data <- melt(x)

# Plot distribution
plot1 <- ggplot(data) + 
  geom_area(aes(x=value,y=..count..,fill=L1),stat='bin',position='identity',bins=2*nrow(gp_noaf)^(1/3)) +
  scale_y_continuous(breaks=seq(0,5500,500),expand=c(0,0),limits=c(0,5500)) +
  scale_x_continuous(breaks=seq(0,15,3),expand=c(0,0),limits=c(0,15),position='top') +
  scale_y_reverse(limits=c(5500,0),expand=c(0,0),breaks=seq(5500,0,-500)) +
  scale_fill_manual(values=alpha(c("#addd8e","#00441b"),0.65),name='',labels=c('CHARGE-AF','Predict-AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.8,0.20),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.text=element_text(size=20)) +
  geom_vline(xintercept=5,linetype='dashed')+ 
  labs(x='Predicted 5-year AF risk (%)',y='Count')
ggsave(filename='over65_density_ctrl.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

# CHARGE VS Predict-AF (cases)
x <- list(v1=gp_incident_af$charge_pred5,v2=gp_incident_af$predict_af_pred5)
data <- melt(x)

# Plot distribution
ggplot(data) +
  geom_area(aes(x=value,y=..count..,fill=L1),stat='bin',position='identity',bins=2*nrow(gp_incident_af)^(1/3)) +
  scale_y_continuous(breaks=seq(0,5500,500),expand=c(0,0),limits=c(0,5500)) +
  scale_x_continuous(breaks=seq(0,15,3),expand=c(0,0),limits=c(0,15)) +
  scale_fill_manual(values=alpha(c("#f16913","#b2182b"),0.65),name='',labels=c('CHARGE-AF','Predict-AF')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.8,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),
        legend.text=element_text(size=20)) +
  geom_vline(xintercept=5,linetype='dashed')+
  labs(x='Predicted 5-year AF risk (%)',y='Count')
ggsave(filename='over65_density_case.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

############################################# Step 24: Comparing LDpred-PRS to PT_PRS
# Load LD
ld_prs <- fread(file='UKBB.LDpred.score')

# Add LD
setkey(af_set,ID); setkey(ld_prs,IID)
af_set[ld_prs,ld_norm := i.ALLELIC_SCORE]

# PRS only
mod_ld <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ ld_norm,data=af_set)

# Predict-AF
mod_predict_ld <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ charge + ld_norm,data=af_set)

## Create scores
af_set[,':='(predict_ld = charge*mod_predict_ld$coefficients[1] + ld_norm*mod_predict_ld$coefficients[2],
                           ld_score = ld_norm*mod_ld$coefficients[1])]

############################################# Calculate predicted risk
### LD Score
ld_avgbeta <- mean(af_set$ld_score)
res <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ ld_score, data=af_set)
km <- survfit(res, data=data.frame(x1=mean(ld_score)),type="kaplan-meier")
ld_s0 <- summary(km, times=5, extend=TRUE)$surv
af_set[,ld_pred5 := (1-(ld_s0)^exp(ld_score - (ld_avgbeta)))*100]

### Predict-LD
predict_ld_avgbeta <- mean(af_set$predict_ld)
res <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_ld, data=af_set)
km <- survfit(res, data=data.frame(x1=mean(predict_ld)),type="kaplan-meier")
predict_ld_s0 <- summary(km, times=5, extend=TRUE)$surv
af_set[,predict_ld_pred5 := (1-(predict_ld_s0)^exp(predict_ld - (predict_ld_avgbeta)))*100]

############################################# Compare model fits and concordance
# Create standardized scores to facilitate comparison
af_set[,':='(ld_score_std = (ld_score - mean(ld_score))/sd(ld_score),
                           predict_ld_std = (predict_ld - mean(predict_ld))/sd(predict_ld))]

################################# AF
# LD PRS
boot_ld <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='ld_score_std',data=af_set,runs=200)
ld_prs_cox <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ ld_score_std,data=af_set)
ld_prs_cph <- cph(Surv(incd_af_5y.t,incd_af_5y) ~ ld_score_std,data=af_set)

## Values
cstat_ld <- c(summary(ld_prs_cox)$concordance[1],
              summary(ld_prs_cox)$concordance[1]-1.96*summary(ld_prs_cox)$concordance[2],
              summary(ld_prs_cox)$concordance[1]+1.96*summary(ld_prs_cox)$concordance[2])
r2_ld <- c(ld_prs_cph$stats['R2'],ld_prs_cph$stats['R2']-1.96*sd(boot_ld),ld_prs_cph$stats['R2']+1.96*sd(boot_ld))
hr_ld <- c(exp(ld_prs_cox$coefficients[1]),exp(confint(ld_prs_cox)[1]),exp(confint(ld_prs_cox)[2]))
aic_ld <- AIC(ld_prs_cox)

# Predict-LD
boot_predict_ld <- boot_r2(time='incd_af_5y.t',status='incd_af_5y',response='predict_ld_std',data=af_set,runs=200)
predict_ld_cox <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_ld_std,data=af_set)
predict_ld_cph <- cph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_ld_std,data=af_set)

## Values
cstat_predict_ld <- c(summary(predict_ld_cox)$concordance[1],
                      summary(predict_ld_cox)$concordance[1]-1.96*summary(predict_ld_cox)$concordance[2],
                      summary(predict_ld_cox)$concordance[1]+1.96*summary(predict_ld_cox)$concordance[2])
r2_predict_ld <- c(predict_ld_cph$stats['R2'],predict_ld_cph$stats['R2']-1.96*sd(boot_predict_ld),predict_ld_cph$stats['R2']+1.96*sd(boot_predict_ld))
hr_predict_ld <- c(exp(predict_ld_cox$coefficients[1]),exp(confint(predict_ld_cox)[1]),exp(confint(predict_ld_cox)[2]))
aic_predict_ld <- AIC(predict_ld_cox)

# C-stat comparisons
## Estimate standard errors using bootstrapping
ld_pt <- boot_compare(time="incd_af_5y.t",status='incd_af_5y',data=af_set,
                      response1='ld_score',response2='prs_score',runs=200)
predict_predict <- boot_compare(time="incd_af_5y.t",status='incd_af_5y',data=af_set,
                                response1='predict_ld',response2='predict_af',runs=200)

# Can now calculate 95% CI around difference in c-stat
c_predict <- summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_af,data=af_set))$concordance[1]
c_predict_ld <- summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_ld,data=af_set))$concordance[1]
c_pt <- summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ prs_score,data=af_set))$concordance[1]
c_ld <- summary(coxph(Surv(incd_af_5y.t,incd_af_5y) ~ ld_score,data=af_set))$concordance[1]

# Differences (with 95% CI)
c_diff_ld_pt <- c(c_ld-c_pt,(c_ld-c_pt)-1.96*sd(ld_pt),(c_ld-c_pt)+1.96*sd(ld_pt))
c_diff_predict_predict <- c(c_predict-c_predict,(c_predict-c_predict)-1.96*sd(predict_predict),(c_predict-c_predict)+1.96*sd(predict_predict))

# P-values
p_ld_pt <- 2*(1-pnorm(abs(c_diff_ld_pt[1]/sd(ld_pt))))
p_predict_predict <- 2*(1-pnorm(abs(c_diff_predict_predict[1]/sd(predict_predict))))

############################################# Calibration
#### Slopes
# Fit models
mod_ld <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ ld_score,data=af_set)
mod_predict_ld <- coxph(Surv(incd_af_5y.t,incd_af_5y) ~ predict_ld,data=af_set)

# Coefficients are slopes
slope_ld <- c(mod_ld$coefficients[1],confint(mod_ld)[1],confint(mod_ld)[2])
slope_predict_ld <- c(mod_predict_ld$coefficients[1],confint(mod_predict_ld)[1],confint(mod_predict_ld)[2])

# Correct slopes for internal estimates
## Estimate distributions
boot_ld <- boot_cal(time='incd_af_5y.t',status='incd_af_5y',response='ld_score',data=af_set,runs=200)
boot_predict_ld <- boot_cal(time='incd_af_5y.t',status='incd_af_5y',response='predict_ld',data=af_set,runs=200)

## Corrected estimates
corrected_ld <- slope_ld[1]*2 - mean(boot_ld)
corrected_predict_ld <- slope_predict_ld[1]*2 - mean(boot_predict_ld)

############################################# Step 25: Standardized NB example
# Determine 5-year AF rate in sample
af_rate <- cuminc(data=af_set,time='incd_af_5y.t',status='incd_af_5y')

# Determine sens/spec of Predict-AF at 5% risk threshold
chars <- SeSpPPVNPV(cutpoint=5, T=af_set$incd_af_5y.t, delta=af_set$incd_af_5y, 
                    marker=af_set$predict_af_pred5, cause=1,
                    weighting = "marginal", times=4.999)

# Calculate SNB
nb <- (af_rate[,3]/100)*chars$TP[2]-0.05/0.95*(1-(af_rate[,3]/100))*chars$FP[2]
max_nb <- af_rate[,3]/100
snb <- nb/max_nb