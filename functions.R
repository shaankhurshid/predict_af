# F1: R2 bootstrap
boot_r2 <- function(time,status,response,data,runs,size=nrow(data)){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE)]
    out[i] <- cph(Surv(sample[,get(time)],sample[,get(status)]) ~ sample[,get(response)],data=sample)$stats['R2']
    if (i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(out)
}

# F2: Interval validation of calibration slope
boot_cal <- function(time,status,response,data,runs,size=nrow(data)){
  out <- rep(NA,times=runs)
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE)]
    model <- coxph(Surv(sample[,get(time)],sample[,get(status)]) ~ sample[,get(response)],data=sample)
    out[i] <- model$coefficients[1]
    if (i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(out)
}

# F3: Diffrence in c-index bootstrap
boot_compare <- function(time,status,data,response1,response2,runs,size=nrow(data)){
  diff <- list()
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE),]
    c1 <- summary(coxph(Surv(sample[,get(time)],sample[,get(status)]) ~ sample[,get(response1)]))$concordance[1]
    c2 <- summary(coxph(Surv(sample[,get(time)],sample[,get(status)]) ~ sample[,get(response2)]))$concordance[1]
    diff[[i]] <- abs(c1-c2)
    if (i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(do.call('rbind',diff))
}

# F4: Cumulative incidence of event 
cuminc <- function(data,time,status,response){
  obj <- survfit(Surv(data[,get(time)],data[,get(status)]) ~ 1)
  ci <- c((1-obj$surv[length(obj$surv)])*100,
          (1-obj$upper[length(obj$upper)])*100,
          (1-obj$lower[length(obj$lower)])*100)
  n_event <- sum(data[,get(status)])
  n_total <- nrow(data)
  out <- data.table(n_total=n_total,n_event=n_event,ci=ci[1],ci_lower=ci[2],ci_upper=ci[3])
  return(out)
}

# F5: Bootstrap function to obtain distribution of differences in event rates
boot_cuminc <- function(time,status,data,response1,response2,runs,size=nrow(data)){
  diff <- list()
  for (i in 1:runs){
    sample <- data[sample(1:nrow(data),size=size,replace=TRUE)]
    subset1 <- sample[get(response1)==1]
    subset2 <- sample[get(response2)==1]
    est1 <- survfit(Surv(subset1[,get(time)],subset1[,get(status)]) ~ 1)$surv
    est2 <- survfit(Surv(subset2[,get(time)],subset2[,get(status)]) ~ 1)$surv
    diff[[i]] <- c((1-est1[length(est1)]), (1-est2[length(est2)]), ((1-est1[length(est1)]) - (1-est2[length(est2)])))
    if (i %% 50 == 0){print(paste0('run ',i,' complete'))}
  }
  return(do.call('rbind',diff))
}