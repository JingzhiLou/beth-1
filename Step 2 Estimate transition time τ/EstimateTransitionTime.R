
# estimate transition time
# reference: Wang MH, Lou J, Cao L, et al. Characterization of key amino acid substitutions and dynamics of the influenza virus H3N2 hemagglutinin. J Infect. 2021;83(6):671-677. doi:10.1016/j.jinf.2021.09.026


# step 1: calculate site-based amino acid prevalence
# input csv format sequence data and out put site wise prevalence through time

library(plyr)
inputdata_name <- c("H3N2_HA_HK_sequence.csv")
outputdata_name <- c("H3N2_HA_HK_prevalence.csv")

seq <- read.csv(inputdata_name,sep=",")
year <- sort(unique(seq$season))
year <- as.vector(year)
codata <- matrix(nrow=length(year),ncol=0)
for (i in 1:(length(seq[1,])-4)){
  ppdata <- NULL
  for (j in year){
    data <- seq[which(seq$season == j),]
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
rownames(codata)<-year
write.csv(codata,file=outputdata_name)



# step 2: calculate g-measure and transition time under different theta and h
# input site wise prevalence and output g-measure and average transition time

inputdata_name <- c("H3N2_HA_HK_prevalence.csv")
outputdata_name <- c("H3N2_HA_HK_gmeasure.csv")

myh <- read.csv(inputdata_name)
myh <- myh[,-1]
gmeasure <- matrix(ncol=0,nrow=length(myh[,1]))
colname <- c()
transition_time <- c()
for (theta in c(5:10)){
  theta <- theta/10
  for (h in c(0:9)){
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
} # here we test theta ranged from 0.5 to 1, and h ranged from 0 to 9 years
rownames(gmeasure)<-year
colnames(gmeasure)<-colname
gmeasure <- rbind(gmeasure,transition_time)
write.csv(gmeasure,file=outputdata_name)   



# step 3: fit regression between epidemic trend and g-measure to find the optimal (theta,h) and corresponding transition time
# output optimal (theta,h), corresponding p-value and average transition time

inputdata_gmeasure_name <- c("H3N2_HA_HK_gmeasure.csv")
inputdata_epidemicdata_name <- c("H3N2_HK_EpidemicData.csv")

myv_seq <- read.csv(inputdata_gmeasure_name)
myv_epi <- read.csv(inputdata_epidemicdata_name)
myv_seq <- myv_seq[-(length(myv_seq[,1])),]
myv <- cbind(myv_epi,myv_seq)
lm_coeff <- matrix(0,nrow=2,ncol=60)
for (i in 1:60){
  res <- NULL
  res <- summary(m1 <- lm(h3rate~season+temperature+ab_humidity+myv[,i+5],data=myv))
  lm_coeff[1,i] <- res$coefficients[5,1]
  lm_coeff[2,i] <- res$coefficients[5,4]
}
colnames(lm_coeff) <- colname
print(c(colname[which.min(lm_coeff[2,])],paste0("p-value=",lm_coeff[2,which.min(lm_coeff[2,])]),
      paste0("transition time=",transition_time[which.min(lm_coeff[2,])])))



