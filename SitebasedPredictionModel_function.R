#install packages:

#biostring 
#typical installation time: 5-6 minutes 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

#lubridate
install.packages("tidyverse")

#plyr
install.packages('plyr')


###########################################################
# Step 1: prepare protein sequence
# The function 'Prepareseq' allows to input an aligned protein sequence data and output a csv format data with flu season

library(Biostrings)
library(lubridate)

PrepareSeq <- function(input_fasta_seq){
  fasta_seq <- readAAStringSet(input_fasta_seq)
  name <- strsplit(fasta_seq@ranges@NAMES,split=" | ",fixed=T)
  seq <- strsplit(as.character(tolower(fasta_seq)),split="")
  aa_name <- c()
  complete_data <- matrix(0,ncol=(length(fasta_seq[[1]])+4),nrow=length(fasta_seq))
  for (i in 1:length(fasta_seq)){
    complete_data[i,1]<-name[[i]][1] #accession number#
    complete_data[i,2]<-name[[i]][3] #isolate name#
    complete_data[i,3]<-name[[i]][2] #date
    if(month(name[[i]][2])<10){
      complete_data[i,4]<-year(name[[i]][2])}else{
        complete_data[i,4]<-year(name[[i]][2])+1
      }
    for (j in 1:length(fasta_seq[[1]])){
      complete_data[i,j+4]<-seq[[i]][j]
    }
  }
  for (k in 1:length(fasta_seq[[1]])){
    aa_name <- c(aa_name,paste("x",k,sep="",collapse=""))
  }
  colnames(complete_data) <- c("accession number","name","date","season",aa_name)
  complete_data <- as.data.frame(complete_data)
  complete_data[complete_data=="-"] <- NA
  return(complete_data)
}

ha_seq <- PrepareSeq(input_fasta_seq = "H3N2_HA_HK_sequence.fasta")
na_seq <- PrepareSeq(input_fasta_seq = "H3N2_NA_HK_sequence.fasta")

# output: ha_seq and na_seq


###########################################################
# step 2: calculate site wise amino acid prevalence
# The function 'SitePrev' allows to input csv format sequence data and output site wise prevalence through time
# Run time: 3-5 minutes

library(plyr)

SitePrev <- function(input_csv_seq){
  seq <- input_csv_seq
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
  codata <- as.data.frame((lapply(codata,as.numeric)))
  rownames(codata)<-year
  return(codata)
}

ha_prev <- SitePrev(input_csv_seq = ha_seq)
na_prev <- SitePrev(input_csv_seq = na_seq)

# output: ha_prev and na_prev


###########################################################
# step 3: calculate g-measure and transition time under different theta and h
# The function 'Getgmeature' allows to input site wise prevalence and output g-measure and average transition time
# run time: 2-3 minutes
# reference: Wang MH, Lou J, Cao L, et al. Characterization of key amino acid substitutions and dynamics of the influenza virus H3N2 hemagglutinin. J Infect. 2021;83(6):671-677. doi:10.1016/j.jinf.2021.09.026

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

ha_gmeasure <- Getgmeature(input_prev = ha_prev, theta_range = c(5:10),h_range = c(0:9))
na_gmeasure <- Getgmeature(input_prev = na_prev, theta_range = c(5:10),h_range = c(0:9))

# please note that for coding convenience, the theta_range here is the actual theta*10. For example, theta_range = c(5:10) means theta from 0.5 to 1.
# output: ha_gmeasure and na_gmeasure 


####################################################################################
# step 4: fit regression between epidemic trend and g-measure to find the optimal (theta,h) and corresponding transition time
# The function 'FitRegression' outputs optimal (theta,h), corresponding p-value and average transition time

FitRegression <- function(input_gmeasure,input_epi_data){
  myv_gmeasure <- input_gmeasure
  myv_epi <- input_epi_data
  myv_seq <- myv_gmeasure[-(length(myv_gmeasure[,1])),]
  myv <- cbind(myv_epi,myv_seq)
  lm_coeff <- matrix(0,nrow=2,ncol=length(myv_seq[1,]))
  for (i in 1:length(myv_seq[1,])){
    res <- NULL
    res <- summary(m1 <- lm(h3rate~season+temperature+ab_humidity+myv[,i+4],data=myv))
    lm_coeff[1,i] <- res$coefficients[5,1]
    lm_coeff[2,i] <- res$coefficients[5,4]
  }
  print(c(colnames(myv_gmeasure)[which.min(lm_coeff[2,])],paste0("p-value=",lm_coeff[2,which.min(lm_coeff[2,])]),
          paste0("transition time=",myv_gmeasure[length(myv_gmeasure[,1]),which.min(lm_coeff[2,])])))
  transtime <- round(myv_gmeasure[length(myv_gmeasure[,1]),which.min(lm_coeff[2,])])
  return(transtime)
}

ha_tau <- FitRegression(input_gmeasure = ha_gmeasure, input_epi_data = read.csv("H3N2_HK_EpidemicData.csv"))
na_tau <- FitRegression(input_gmeasure = na_gmeasure, input_epi_data = read.csv("H3N2_HK_EpidemicData.csv"))

# output: ha_tau and na_tau


####################################################################################
# step 5: prediction of future mutation prevalence and consensus sequence
# The function 'PredictConseq' allows to input site wise prevalence table and output yearly predicted consensus sequence

PredictConseq <- function(input_prev,transition_time,predict_year){
  Prev_Table <- input_prev
  transtime <- transition_time
  Pred_Prevalence <- matrix(0,ncol=length(Prev_Table))
  colnames(Pred_Prevalence) <- colnames(Prev_Table)
  for (year in predict_year){
    pred_prev <- c()
    Test_PrevTable <- Prev_Table[which(rownames(Prev_Table)<year),]
    t <- length(Test_PrevTable[,1]) 
    #Trend Prediction#
    for(i in 1:length(Test_PrevTable[1,])){
      tzero_list <- which(Test_PrevTable[1:(t-1),i]==0)
      if (length(tzero_list)!=0){
        tzero <- max(which(Test_PrevTable[1:(t-1),i]==0)) 
        r <- t-max(tzero, t-transtime)
        delta <- (Test_PrevTable[t,i]-Test_PrevTable[t-r,i])/r
        if (delta <=0){
          pred_prev <- c(pred_prev,max(0,(Test_PrevTable[t,i]+delta)))}
        else{pred_prev <- c(pred_prev,min(1,(Test_PrevTable[t,i]+delta)))}
      }
      if (length(tzero_list)==0){
        r <- transtime
        delta <- (Test_PrevTable[t,i]-Test_PrevTable[t-r,i])/r	
        if (delta <=0){
          pred_prev <- c(pred_prev,max(0,(Test_PrevTable[t,i]+delta)))}
        else{pred_prev <- c(pred_prev,min(1,(Test_PrevTable[t,i]+delta)))}
      }
    }
    Pred_Prevalence <- rbind(Pred_Prevalence,pred_prev)
  }
  Pred_Prevalence <- Pred_Prevalence[-1,]
  colname <- colnames(Prev_Table)
  test <- strsplit(colname,split=".",fixed=T)
  xname <- c() 
  yname <- c() 
  for (j in 1:length(test)){
    xname <- c(xname,test[[j]][1])
    yname <- c(yname,test[[j]][2])}
  xuniq <- unique(xname)
  Pred_Sequence <- matrix(0,ncol=length(xuniq))
  for (k in 1:length(predict_year)){
    preseq <- rbind(xname,yname,Pred_Prevalence[k,])
    pred_seq <- matrix(0,nrow=3)
    for (l in xuniq){
      xpos <- which(xname==l)
      posmax <- max(Pred_Prevalence[k,min(xpos):max(xpos)])
      for (m in xpos){
        if (preseq[3,m]==posmax){
          pred_seq <- cbind(pred_seq,preseq[,m])
          break 
        }}}
    pred_seq <- pred_seq[,-1]
    Pred_Sequence <- rbind(Pred_Sequence,pred_seq)}
  Pred_Sequence <- Pred_Sequence[seq(3,3*length(predict_year),3),]
  Pred_Sequence <- as.data.frame(Pred_Sequence)
  rownames(Pred_Sequence) <- predict_year
  return(Pred_Sequence)
}

ha_predicted_conseq <- PredictConseq(input_prev = ha_prev, transition_time = ha_tau, predict_year = c(2012:2021))
na_predicted_conseq <- PredictConseq(input_prev = na_prev, transition_time = na_tau, predict_year = c(2012:2021))

# output: ha_predicted_conseq and na_predicted_conseq


##########################################################################
# step 6: select wild type strain based on dynamic predictor sites set W(t)
# The function allpes to input yearly predicted consensus strains and output wildtype sequences
# Note: the EM sites and predictor sites include the signal piptide and start with M(Methionine)

WTselection <- function(HA_sequence,NA_sequence,HA_conseq,NA_conseq,EMset,predict_year){
  ha_seq <- HA_sequence
  na_seq <- NA_sequence
  ha_pre_seq <- cbind(predict_year,HA_conseq)
  na_pre_seq <- cbind(predict_year,NA_conseq)
  All_emset <- EMset
  colnames(ha_pre_seq)[1] <- "season"
  colnames(na_pre_seq)[1] <- "season"
  h3ha_as <- c(122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198,44,45,46,47,48,50,51,53,54,273,275,276,278,279,280,294,297,299,300,304,305,307,308,309,310,311,312,96,102,103,117,121,167,170,171,172,173,174,175,176,177,179,182,201,203,205,207,208,209,212,213,214,215,216,217,218,219,226,227,228,229,230,240,242,244,246,247,248,57,59,62,63,67,75,78,80,81,82,83,86,87,88,91,92,94,109,260,261,262,265)
  h3na_as <- c(150,198,199,220,221,253,329,334,344,368,370,403)
  ha_beta <- 4.17 # please refer reference 1 & 2
  na_beta <- 12.28
  alpha <- ha_beta/(ha_beta+na_beta)
  ha_allem <- All_emset[,c(1:2)] 
  na_allem <- All_emset[,c(3:4)]
  ha_as <- h3ha_as 
  na_as <- h3na_as
  pre_year <- predict_year
  wildtype_ha <- matrix(ncol=length(ha_seq[1,]),nrow=0)
  wildtype_na <- matrix(ncol=length(na_seq[1,]),nrow=0)
  for (s in pre_year){
    ha_test <- ha_seq[which(ha_seq$season==s-1),]
    na_test <- na_seq[which(na_seq$season==s-1),]
    pre_hatest <- ha_pre_seq[which(ha_pre_seq$season==s),]
    pre_natest <- na_pre_seq[which(na_pre_seq$season==s),]
    ha_em <- ha_allem[which(ha_allem[,2]==s-1),1]
    na_em <- na_allem[which(na_allem[,2]==s-1),1]
    total_emd <- c()
    EMD_ha <- c()
    EMD_na <- c()
    for (k in 1:length(ha_test[,1])){
      name_test <- as.vector(ha_test[k,2])
      na_corresp <- na_test[which(na_test$name==name_test),]
      ha_emd <- 0
      na_emd <- 0
      for (m in union(ha_em,ha_as)){# genetic distance in dynamic predictor set
        if (is.na(ha_test[k,m+4])){ha_emd <- ha_emd +1}else{
          if (ha_test[k,m+4]==pre_hatest[1,m+1]){}else{
            ha_emd <- ha_emd+1
          }
        }
      } 
      for (n in union(na_em,na_as)){
        if (is.na(na_corresp[1,n+4])){}else{
          if (na_corresp[1,n+4]==pre_natest[1,n+1]){}else{
            na_emd <- na_emd+1
          }
        }
      } 
      emd <- alpha*ha_emd+(1-alpha)*na_emd
      total_emd <- c(total_emd,emd)
      EMD_ha <- c(EMD_ha,ha_emd)
      EMD_na <- c(EMD_na,na_emd)
    }
    position <- which(total_emd == min(total_emd))
    ha_test_second <- ha_test[position,]
    na_test_second <- na_test[match(as.vector(ha_test_second$name),as.vector(na_test$name)),]
    total_Hamming_Dis <- c()
    ha_site <- setdiff(c(1:566),union(ha_allem,ha_as))
    na_site <- setdiff(c(1:469),union(na_allem,na_as))
    for (k in 1:length(ha_test_second[,1])){
      ha_hamming_d <- 0
      na_hamming_d <- 0
      for (m in ha_site){# genetic distance in other sites
        if (is.na(ha_test_second[k,m+4])){ha_hamming_d <- ha_hamming_d+1}else{
          if (ha_test_second[k,m+4]==pre_hatest[1,m+1]){}else{
            ha_hamming_d <- ha_hamming_d+1
          }
        }
      }
      for (n in na_site){
        if (is.na(na_test_second[k,n+4])){na_hamming_d <- na_hamming_d +1}else{
          if (na_test_second[k,n+4]==pre_natest[1,n+1]){}else{
            na_hamming_d <- na_hamming_d+1
          }
        }
      }
      hamming_dis <- ha_hamming_d+na_hamming_d
      total_Hamming_Dis <- c(total_Hamming_Dis,hamming_dis)
    } 
    position <- which.min(total_Hamming_Dis)
    wildtype_ha <- rbind(wildtype_ha,ha_test_second[position,])
    wildtype_ha_name <- as.vector(ha_test_second[position,2])
    pre_naseq <- na_test_second[which(na_test_second$name==wildtype_ha_name),]
    wildtype_na <- rbind(wildtype_na,pre_naseq[1,])
  }
  colnames(wildtype_ha)[4] <- "predicted season"
  colnames(wildtype_na)[4] <- "predicted season"
  wildtype_ha[,4] <- predict_year
  wildtype_na[,4] <- predict_year
  list(wildtype_ha=wildtype_ha,wildtype_na=wildtype_na)
}

wildtype_seq <- WTselection(HA_sequence = ha_seq, NA_sequence = na_seq, HA_conseq = ha_predicted_conseq, NA_conseq = na_predicted_conseq, 
                        EMset = read.csv("H3N2_HK_EMset.csv"),predict_year = c(2012:2021))

wildtype_ha <- wildtype_seq$wildtype_ha
wildtype_na <- wildtype_seq$wildtype_na

#output: wildtype_ha and wildtype_na


##################################################################################3
# Step 7: calculate epitope distance between predicted vaccine strains and circulating strains
# input vaccine strain and circulating strains, and output genetic mismatch

h3ha_epitope <- c(122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198,44,45,46,47,48,50,51,53,54,273,275,276,278,279,280,294,297,299,300,304,305,307,308,309,310,311,312,96,102,103,117,121,167,170,171,172,173,174,175,176,177,179,182,201,203,205,207,208,209,212,213,214,215,216,217,218,219,226,227,228,229,230,240,242,244,246,247,248,57,59,62,63,67,75,78,80,81,82,83,86,87,88,91,92,94,109,260,261,262,265)+16
h3na_epitope <- c(150,198,199,220,221,253,329,334,344,368,370,403)

GetMismatch <- function(input_pre_seq, input_cir_seq, epitope_site, cir_year){
  cir_seq <- input_cir_seq
  pre_seq <- input_pre_seq
  pre_seq <- as.matrix(pre_seq)
  start_season <- cir_year[1]
  emsite <- epitope_site
  ave_mismatch <- c()
  sd <- c()
  for (z in 1:length(cir_year)){
    obs_seq <- subset(cir_seq,season==z+start_season-1)
    obs_seq <- as.matrix(obs_seq)
    diff <- c()
    for (i in 1:length(obs_seq[,1])){
      d <- 0
      site <- c()
      for (j in emsite){
        if (is.na(obs_seq[i,j+4])){}else{
          if (pre_seq[z,j+4]==obs_seq[i,j+4]){}else{
            d <- d+1
            site <- c(site,j)}}}
      diff <- c(diff,d)}
    ave_mismatch <- c(ave_mismatch,mean(diff))
    sd <- c(sd,sd(diff))}
  Mismatch <- cbind(cir_year,ave_mismatch,sd)
}

ha_mismatch <- GetMismatch(input_pre_seq = wildtype_ha, input_cir_seq = ha_seq, epitope_site = h3ha_epitope, cir_year = c(2012:2020))
na_mismatch <- GetMismatch(input_pre_seq = wildtype_na, input_cir_seq = na_seq, epitope_site = h3na_epitope, cir_year = c(2012:2020))

## output: epitope mismatch tables

## end.
  
