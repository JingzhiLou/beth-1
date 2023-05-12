
## Here, we introduce the pipline to estimate transition time tau and take data from Hong Kong for example.

library(plyr)

setwd("C:/Users/jingzhi/Desktop/Site-basedPredictionModel/data/transition time estimation/HK")

###################### calculate g-measure and transition time under different theta and h #######################
## The function 'Getgmeature' allows to input prevalence table and output g-measure and average transition time

Getgmeature <- function(input_prev,theta_range,h_range){
  myh <- input_prev
  rowname <- myh[,1]
  myh <- myh[,-1]
  gmeasure <- matrix(ncol=0,nrow=length(myh[,1]))
  colname <- c()
  transition_time <- c()
  for (theta in theta_range){
    theta <- theta/10
    for (h in h_range){
      mut <- matrix(0,nrow=length(myh[,1]),ncol=length(myh))
      trans_time <- c()
      for (i in 1:length(myh[1,])){
        start <- 1
        r<-1
        while(r <= length(myh[,1])){
          # detect effective mutation
          if (myh[r,i]>=theta && any(myh[start:r,i]==0)){
            c <- r # record fisrt effective mutation site#
            start <- r+1 # prepare to detect next effective mutation#
            for (s in 1:r){
              if (myh[s,i]==0){
                a <- s           }
            }# record last zero#
            mut[(a+1):c,i] <- 1
            trans_time <- c(trans_time,c-a)
            # extend h years
            if (h!=0){
              if (c+h > length(myh[,1])){
                fakeh <- length(myh[,1])-c 
                if (fakeh!=0){
                  for (j in 1:fakeh){
                    if (myh[(c+j),i] >= theta){
                      mut[(c+j),i] <- 1       }else
                      {break}    
                  }
                }  
              }else{
                for (j in 1:h){
                  if (myh[(c+j),i] >= theta){
                    mut[(c+j),i] <-1}else
                    {break}}}
            }}
          r=r+1               }
      }
      transition_time <- c(transition_time,mean(trans_time))
      product <- c()
      for (y in 1:length(myh)){
        for(x in 1:length(myh[,1])){
          product <- c(product, as.numeric(myh[x,y])*mut[x,y])
        }
      }
      myy <- matrix(product, ncol=length(myh),nrow=length(myh[,1]),byrow=F)
      gsum <- rowSums(myy)
      gmeasure <- cbind(gmeasure,gsum)
      colname <- c(colname,paste0("theta=",theta,",h=",h))
    }
  }
  rownames(gmeasure)<-rowname
  colnames(gmeasure)<-colname
  #gmeasure <- rbind(gmeasure,transition_time)
  return(gmeasure)
}

ha_gmeasure <- Getgmeature(input_prev = read.csv("HKH3_HA_prevalence.csv"), theta_range = c(1:9),h_range = c(0:9))
na_gmeasure <- Getgmeature(input_prev = read.csv("HKH3_NA_prevalence.csv"), theta_range = c(1:9),h_range = c(0:9))

## please note that for coding convenience, the theta_range here is the actual theta*10. For example, theta_range = c(5:10) means theta from 0.5 to 1.

## output: ha_gmeasure and na_gmeasure 

write.csv(ha_gmeasure,file="HKH3_HA_gmeasure.csv")
write.csv(na_gmeasure,file="HKH3_NA_gmeasure.csv")


###################### regression fit #######################
## fit regression between epidemic trend and g-measure to find the optimal (theta,h)
## The function 'FitRegression' outputs optimal (theta,h) and corresponding R-squared, and the result is put in 'summary of regression fit.xlsx'
## The epidemic and meterological data are retrieved and saved in file 'XXXX_regression fit.csv'. The first column is the response variable, 
# and here we use sero-positivity rate of pH1N1 or H3N2 subtype. The following columns include year, average temperature, average humidity, 
# average absolute humidity, HA g-measure and NA g-measure under different theta and h.

FitRegression <- function(input_epi_data){
  myv <- input_epi_data
  ha_lm_coeff <- matrix(0,nrow=4,ncol=90)
  for (i in 1:90){
    res <- NULL
    res <- summary(m1 <- lm(log(myv[,1])~year+temp+ab_humid+myv[,i+5],data=myv))
    ha_lm_coeff[1,i] <- res$coefficients[5,1]
    ha_lm_coeff[2,i] <- res$coefficients[5,4]
    ha_lm_coeff[3,i] <- res$fstatistic[1]
    ha_lm_coeff[4,i] <- res$r.squared
  }
  print(paste0("HA protein  Optimal theta and h: ",colnames(myv)[which.max(ha_lm_coeff[4,])+5], 
               "  R-squared=",ha_lm_coeff[4,which.max(ha_lm_coeff[4,])]))
  
  na_lm_coeff <- matrix(0,nrow=4,ncol=90)
  for (i in 1:90){
    res <- NULL
    res <- summary(m1 <- lm(positive_rate~year+temp+ab_humid+myv[,i+95],data=myv))
    na_lm_coeff[1,i] <- res$coefficients[5,1]
    na_lm_coeff[2,i] <- res$coefficients[5,4]
    na_lm_coeff[3,i] <- res$fstatistic[1]
    na_lm_coeff[4,i] <- res$r.squared
  }
  print(paste0("NA protein  Optimal theta and h: ",colnames(myv)[which.max(na_lm_coeff[4,])+5], 
               "  R-squared=",na_lm_coeff[4,which.max(na_lm_coeff[4,])]))
}

regression_result <- FitRegression(input_epi_data = read.csv("HKH3_regression fit.csv"))


###################### EM and average transition time under optimal (theta,h) #######################
## The function EffMutation obtains EMs and tau under sepcific theta and h, and the result is put in 'summary of EM prevalence table.xlsx'

EffMutation <- function(input_prev_table, theta, year){
  myh <- input_prev_table
  myh <- myh[which(myh$season %in% year),]
  rowname <- myh[,1]
  myh <- myh[,-1]
  sel_col <- c()
  trans_time <- c()
  for(i in 1:length(myh)){
    start <- 1
    r <- 1
    while(r <= length(myh[,1])){
      if (myh[r,i]>=theta && any(myh[start:r,i]==0)){
        c <- r # record fisrt effective mutation time#
        start <- r+1 # prepare to detect next effective mutation in the same substitution
        for (s in 1:r){
          if (myh[s,i]==0){
            a <- s           }
        }#record last zero#
        trans_time <- c(trans_time,(c-a))
        sel_col <- cbind(sel_col,i)                }
      r=r+1          }
  }
  effm <- myh[,sel_col] 
  effm <- rbind(effm,trans_time)
  row.names(effm) <- c(rowname,"tau")
  return(effm)
}

HA_EMsite <- EffMutation(input_prev_table = read.csv("HKH3_HA_prevalence.csv"), theta = 0.8, year=c(2005:2019))
NA_EMsite <- EffMutation(input_prev_table = read.csv("HKH3_NA_prevalence.csv"), theta = 0.9, year=c(2006:2019))

## output
write.csv(HA_EMsite,file="HKH3_HA_EM sites.csv")
write.csv(NA_EMsite,file="HKH3_NA_EM sites.csv")


