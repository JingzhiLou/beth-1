# Re-estimate theta and h recursively at time T using toydata as example

library(lubridate)

# Introduce function 'SitePrev', 'Getgmeasure' and 'FitRegression'

SitePrev <- function(input_csv_seq){
  seq <- input_csv_seq
  year <- sort(unique(seq$season))
  year <- as.vector(year)
  codata <- matrix(nrow=length(year),ncol=0)
  for (i in 1:(length(seq[1,])-4)){
    ppdata <- NULL
    for (j in year){
      data <- seq[which(seq$season == j),]
      data <- data[which(month(seq$date) < 3 | month(seq$date) > 9),]
      datat <- table(data[,i+4])
      dataf <- prop.table(datat)
      dataf <- as.data.frame(dataf)
      dataf <- t(dataf)
      colnames <- NULL
      for (k in 1:length(dataf[1,])){
        colnames <- c(colnames,paste(i,dataf[1,k],seq=""))
      }
      dimnames(dataf)=list(NULL,colnames)
      dataf <- dataf[-1,]
      dataf <- as.data.frame(t(dataf))
      ppdata <- rbind.fill(ppdata,dataf)
    }
    codata <- cbind(codata,ppdata)
  }
  codata[is.na(codata)] <- 0
  codata <- as.data.frame((lapply(codata,as.numeric)))
  rownames(codata)<-year
  return(codata)
}

Getgmeature <- function(input_prev,theta_range,h_range){
  myh <- input_prev
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
          product <- c(product, myh[x,y]*mut[x,y])
        }
      }
      myy <- matrix(product, ncol=length(myh),nrow=length(myh[,1]),byrow=F)
      gsum <- rowSums(myy)
      gmeasure <- cbind(gmeasure,gsum)
      colname <- c(colname,paste0("theta=",theta,",h=",h))
    }
  }
  rownames(gmeasure)<-rownames(myh)
  colnames(gmeasure)<-colname
  gmeasure <- rbind(gmeasure,transition_time)
  return(gmeasure)
}

FitRegression <- function(input_gmeasure,input_epi_data){
  myv_gmeasure <- input_gmeasure
  myv_epi <- input_epi_data
  myv_seq <- myv_gmeasure[-(length(myv_gmeasure[,1])),]
  myv <- cbind(myv_epi,myv_seq)
  lm_coeff <- matrix(0,nrow=4,ncol=length(myv_seq[1,]))
  for (i in 1:length(myv_seq[1,])){
    res <- NULL
    res <- summary(m1 <- lm(h3rate~season+temperature+ab_humidity+myv[,i+4],data=myv))
    lm_coeff[1,i] <- res$coefficients[5,1]
    lm_coeff[2,i] <- res$coefficients[5,4]
    lm_coeff[3,i] <- res$fstatistic[1]
    lm_coeff[4,i] <- res$r.squared
  }
  transtime <- round(myv_gmeasure[length(myv_gmeasure[,1]),which.min(lm_coeff[2,])],digits = 2)
  print(paste0("Optimal theta and h: ",colnames(myv)[which.max(lm_coeff[4,])+5], 
               "  R-squared=",round(lm_coeff[4,which.max(lm_coeff[4,])],digit = 2),
               "  tau=",transtime," year"))
  return(round(transtime))
}


# Re-estimate theta and h recursively at time T and output optimal (theta,h), corresponding R-squared and average transition time
# Run time: ~10 minutes (depends on the length of test_year)

ha_sequence <- read.csv("./toydata/H3N2_HA_HK_sequence.csv")
na_sequence <- read.csv("./toydata/H3N2_NA_HK_sequence.csv")
Epidata <- read.csv("./toydata/H3N2_HK_EpidemicData.csv")
test_year <- c(2017:2019)

for (year in test_year){
  print(paste0("test year: ",year))
  ha_seq <- ha_sequence[which(ha_sequence$date < make_date(year-1,3,1)),]
  na_seq <- na_sequence[which(na_sequence$date < make_date(year-1,3,1)),]
  ha_prev <- SitePrev(input_csv_seq = ha_seq)
  na_prev <- SitePrev(input_csv_seq = na_seq)
  ha_gmeasure <- Getgmeature(input_prev = ha_prev, theta_range = c(5:10),h_range = c(0:9))
  na_gmeasure <- Getgmeature(input_prev = na_prev, theta_range = c(5:10),h_range = c(0:9))
  test_epi_data <- Epidata[which(Epidata$season < year),]
  print("HA protein fitting result:")
  ha_tau <- FitRegression(input_gmeasure = ha_gmeasure, input_epi_data = test_epi_data)
  print("NA protein fitting result:")
  na_tau <- FitRegression(input_gmeasure = na_gmeasure, input_epi_data = test_epi_data)
}

# output: average transition time for HA and NA protein by year
# end.
