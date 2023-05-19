# This R script is conducted using toydata in Hong Kong as example.

#install packages:
#typical installation time: 5-6 minutes 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
install.packages("lubridate")
install.packages('plyr')
install.packages('stringr')

# set working dictionary and library packages
setwd("./Site-basedPredictionModel/toydata")
library(Biostrings)
library(lubridate)
library(plyr)
library(stringr)

###########################################################
# Step 1: prepare protein sequence
# The function 'Prepareseq' allows to input an aligned protein sequence data and output a csv format data with flu season

PrepareSeq <- function(input_fasta_seq){
  fasta_seq <- readAAStringSet(input_fasta_seq)
  name <- strsplit(fasta_seq@ranges@NAMES,split="_|_",fixed=T)
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

ha_seq <- PrepareSeq(input_fasta_seq = "H3N2_HA_HK_sequence.fas")
na_seq <- PrepareSeq(input_fasta_seq = "H3N2_NA_HK_sequence.fas")

# write.csv(ha_seq,file="HKH3_HA_sequence.csv",row.names = F)
# write.csv(na_seq,file="HKH3_NA_sequence.csv",row.names = F)

# Output: HA and NA sequences with sequence name, accession number, collection date and flu season


###########################################################
# step 2: calculate site wise amino acid prevalence
# The function 'SitePrev' allows to input csv format sequence data and output site wise prevalence through time
# Run time: 3-5 minutes

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

ha_prev <- SitePrev(input_csv_seq = ha_seq)
na_prev <- SitePrev(input_csv_seq = na_seq)

# write.csv(ha_prev,file="./toydata_result/H3N2_HA_HK_prevalence.csv",row.names = F)
# write.csv(na_prev,file="./toydata_result/H3N2_HA_HK_prevalence.csv",row.names = F)

# Output: HA and NA yearly site prevalence table


###########################################################
# step 3: calculate g-measure and transition time under different theta and h
# The function 'Getgmeature' allows to input site wise prevalence and output g-measure and average transition time
# Run time: 2-3 minutes
# Reference: Wang MH, Lou J, Cao L, et al. Characterization of key amino acid substitutions and dynamics of the influenza virus H3N2 hemagglutinin. J Infect. 2021;83(6):671-677. doi:10.1016/j.jinf.2021.09.026

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

# write.csv(ha_gmeasure,file="./toydata_result/H3N2_HA_HK_gmeasure.csv",row.names = F)
# write.csv(na_gmeasure,file="./toydata_result/H3N2_NA_HK_gmeasure.csv",row.names = F)

# Output: HA and NA yearly g-measure


####################################################################################
# step 4: fit regression between epidemic trend and g-measure to find the optimal (theta,h) and corresponding transition time
# The function 'FitRegression' outputs optimal (theta,h), corresponding R-squared and average transition time

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

ha_tau <- FitRegression(input_gmeasure = ha_gmeasure, input_epi_data = read.csv("H3N2_HK_EpidemicData.csv"))
na_tau <- FitRegression(input_gmeasure = na_gmeasure, input_epi_data = read.csv("H3N2_HK_EpidemicData.csv"))

# output: average transition time for HA and NA protein


####################################################################################
# step 5: Obtain history EMs under optimal (theta,h)
# The function 'EffMutation' obtains history EMs under specific theta and h

EffMutation <- function(input_prev_table, theta){
  myh <- input_prev_table
  myh$season <- row.names(input_prev_table)
  rowname <- myh$season
  colname <- colnames(myh)
  em_site <- c()
  em_season <- c()
  em_tau <- c()
  for(i in 1:length(myh)){
    start <- 1
    r <- 1
    while(r <= length(myh[,1])){
      if (myh[r,i]>=theta && any(myh[start:r,i]==0)){
        c <- r # record first effective mutation time#
        start <- r+1 # prepare to detect next effective mutation in the same substitution
        for (s in 1:r){
          if (myh[s,i]==0){
            a <- s           }
        }
        site_split <- strsplit(colname[i],split=".",fixed=T)[[1]][1]
        site_position <- as.numeric(strsplit(site_split,split="X",fixed=T)[[1]][2])
        em_site <- c(em_site, site_position)
        em_season <- c(em_season,as.numeric(rowname[c]))
        em_tau <- c(em_tau,(c-a))
                     }
      r=r+1          }
  }
  effm <- cbind(em_site,em_season,em_tau)
  return(effm)
}

HA_EMsite <- EffMutation(input_prev_table = ha_prev, theta = 0.7)
NA_EMsite <- EffMutation(input_prev_table = na_prev, theta = 0.6)

# write.csv(HA_EMsite,file="./toydata_result/H3N2_HA_HK_EMsites.csv",row.names = F)
# write.csv(NA_EMsite,file="./toydata_result/H3N2_NA_HK_EMsites.csv",row.names = F)

# output: history EMs with their transition time


####################################################################################
# Step 6: beth-1 prediction
# The beth-1 consists of future mutation prevalence prediction and vaccine strain selection


# Step 6.1: The function PreYearlyStrain allows to calculate one-step prediction using beth-1 and history tau, and output predicted consensus strain and predicted EM set
# Run time: 5-7 minutes

PreYearlyStrain <- function(input_csv_seq, predict_season, transition_time, his_EMset){
  Pred_Conseq <- matrix(0,ncol=length(input_csv_seq[1,])-4)
  Pred_EM <- matrix (0,ncol=2)
  for (pre_season in predict_season){
    print(pre_season)
    seq <- input_csv_seq[which(input_csv_seq$date >= make_date(1998,10,1) & input_csv_seq$date < make_date(pre_season-1,3,1)),] 
    year <- sort(unique(seq$season))
    year <- as.vector(year)
    codata <- matrix(nrow=length(year),ncol=0) # get site-wise prevalence table
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
    
    #beth-1 Prediction#
    Test_his_emset <- his_EMset[which(his_EMset[,2]<pre_season),]
    Pred_Prevalence <- matrix(0,ncol=length(codata))
    pred_prev <- c()
    Test_PrevTable <- codata
    t <- length(Test_PrevTable[,1]) 
    for(i in 1:length(Test_PrevTable[1,])){
      if (i %in% Test_his_emset[,1]){
        position <- max(which(Test_his_emset[,1]==i))
        transtime <- Test_his_emset[position,3]
      }else{transtime <- transition_time}
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
    colnames(Pred_Prevalence) <- colnames(Test_PrevTable)
    
    colname <- colnames(Test_PrevTable)
    test <- strsplit(colname,split=".",fixed=T)
    xname <- c() 
    yname <- c() 
    for (j in 1:length(test)){
      xname <- c(xname,test[[j]][1])
      yname <- c(yname,test[[j]][2])}
    xuniq <- unique(xname)
    Pred_Sequence <- matrix(0,ncol=length(xuniq))
    preseq <- rbind(xname,yname,Pred_Prevalence[2,])
    pred_seq <- matrix(0,nrow=3)
    for (l in xuniq){
      xpos <- which(xname==l)
      posmax <- max(Pred_Prevalence[2,min(xpos):max(xpos)])
      for (m in xpos){
        if (preseq[3,m]==posmax){
          pred_seq <- cbind(pred_seq,preseq[,m])
          break 
        }}}
    pred_seq <- pred_seq[,-1]
    Pred_Conseq <- rbind(Pred_Conseq,pred_seq[2,])
    
    #detect predicted EMs#
    Pre_codata <- rbind(codata, Pred_Prevalence[2,])
    for( col_k in 1:length(Pre_codata[1,])){
      start_year <- 1
      theta_year <- 1
      while(theta_year <= length(Pre_codata[,1])){
        if (Pre_codata[theta_year,col_k]>=0.5 & any(Pre_codata[start_year:theta_year,col_k]==0)){
          start_year <- theta_year+1
          pred_em_site <- c(colnames(Pre_codata)[col_k],pre_season)
          Pred_EM <- rbind(Pred_EM,pred_em_site)}
        theta_year <- theta_year+1
      }
    }
  }
  Pred_Conseq <- Pred_Conseq[-1,]
  rownames(Pred_Conseq) <- predict_season
  Pred_EM <- Pred_EM[-1,]
  return(list(beth_conseq = Pred_Conseq, Predicted_EM = Pred_EM, PrevTable = codata))
}

ha_predicted_conseq <- PreYearlyStrain(input_csv_seq = ha_seq, transition_time = ha_tau, predict_season = c(2012:2020), his_EMset = HA_EMsite)
na_predicted_conseq <- PreYearlyStrain(input_csv_seq = na_seq, transition_time = na_tau, predict_season = c(2012:2020), his_EMset = NA_EMsite)

# write.csv(ha_predicted_conseq$beth_conseq,file="./toydata_result/BethH3_HA_predicted consensus seq.csv")
# write.csv(na_predicted_conseq$beth_conseq,file="./toydata_result/BethH3_NA_predicted consensus seq.csv")
# write.csv(ha_predicted_conseq$Predicted_EM,file="./toydata_result/HKH3_predicted_EMset_HA.csv")
# write.csv(na_predicted_conseq$Predicted_EM,file="./toydata_result/HKH3_predicted_EMset_NA.csv")

# output: predicted consensus sequence and predicted EM set
# The output predicted consensus sequence are "BethH3_HA_predicted consensus seq.csv" and "BethH3_NA_predicted consensus seq
# The output predicted EM set is summarized in "HKH3_predicted_EMset.csv"


# Step 6.2: The function WTselection allows to select wildtype strain integrating mutation information on both HA and NA

## H3N2 epitopes
h3ha_as_ab <- c(122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198)+16
h3ha_as <- c(122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198,44,45,46,47,48,50,51,53,54,273,275,276,278,279,280,294,297,299,300,304,305,307,308,309,310,311,312,96,102,103,117,121,167,170,171,172,173,174,175,176,177,179,182,201,203,205,207,208,209,212,213,214,215,216,217,218,219,226,227,228,229,230,240,242,244,246,247,248,57,59,62,63,67,75,78,80,81,82,83,86,87,88,91,92,94,109,260,261,262,265)+16
h3na_as <- c(150,198,199,220,221,253,329,334,344,368,370,403)  

WTselection <- function(input_ha_seq, input_na_seq, pred_emset, ha_pre_seq, na_pre_seq, ha_beta, na_beta, ha_as, na_as, pre_year){
  seq1 <- input_ha_seq
  seq2 <- input_na_seq
  name1 <- str_trim(as.vector(seq1$name))
  name2 <- str_trim(as.vector(seq2$name))
  name <- intersect(name1,name2)
  ha_seq <- seq1[which(seq1$name %in% name),]
  na_seq <- seq2[which(seq2$name %in% name),]
  ha_predem <- pred_emset[,c(1:2)] 
  na_predem <- pred_emset[,c(3:4)]
  colnames(ha_pre_seq)[1] <- "season"
  colnames(na_pre_seq)[1] <- "season"
  wildtype_ha <- matrix(ncol=length(ha_seq[1,]),nrow=0)
  wildtype_na <- matrix(ncol=length(na_seq[1,]),nrow=0)
  for (s in pre_year){
    print(s)
    ha_test <- ha_seq[which(ha_seq$date<make_date(s-1,3,1)),] #s-1
    na_test <- na_seq[which(na_seq$date<make_date(s-1,3,1)),] #s-1
    pre_hatest <- ha_pre_seq[which(ha_pre_seq$season==s),]
    pre_natest <- na_pre_seq[which(na_pre_seq$season==s),]
    ha_em <- unique(ha_predem[which(ha_predem[,2]<=s),1])
    na_em <- unique(na_predem[which(na_predem[,2]<=s),1])
    total_VE_EMD <- c()
    EMD_ha <- c()
    EMD_na <- c()
    ha_ve_site <- intersect(ha_em,h3ha_as_ab) # H1N1: EMs on antigenic sites, H3N2: EMs on antigenic sites A & B  ¡Ì¡Ì¡Ì
    na_ve_site <- na_as # H1N1: EMs on antigenic sites, H3N2: antigenic sites  ¡Ì¡Ì¡Ì
    
    for (k in 1:length(ha_test[,1])){ # genetic distance in dynamic predictor set
      name_test <- as.vector(ha_test[k,2])
      na_corresp <- na_test[which(na_test$name==name_test),]
      ha_emd <- 0
      na_emd <- 0
      for (m in ha_ve_site){ 
        if (is.na(ha_test[k,m+4])){ha_emd <- ha_emd +1}else{
          if (ha_test[k,m+4]==pre_hatest[1,m+1]){}else{
            ha_emd <- ha_emd+1
          }
        }
      } 
      
      for (n in na_ve_site){
        if (is.na(na_corresp[1,n+4])){na_emd <- na_emd +1}else{
          if (na_corresp[1,n+4]==pre_natest[1,n+1]){}else{
            na_emd <- na_emd+1
          }
        }
      } 
      
      ve_emd <- ha_beta*ha_emd+na_beta*na_emd
      total_VE_EMD <- c(total_VE_EMD,ve_emd)
      EMD_ha <- c(EMD_ha,ha_emd)
      EMD_na <- c(EMD_na,na_emd)
    } #select the wildtype sequence by maximizing predicted VE
    
    position <- which(total_VE_EMD == max(total_VE_EMD))
    ha_test_second <- ha_test[position,]
    na_test_second <- na_test[match(as.vector(ha_test_second$name),as.vector(na_test$name)),]
    
    # when multiple wild-type strains give tied result, they are prioritized by 
    # (1) minimizing EM or antigenic sites mismatch; 
    # (2) minimizing mismatch on the remaining sites of full-length strain 
    
    total_antigenic_dis <- c()
    ha_site <- setdiff(union(ha_em,ha_as),ha_ve_site) 
    na_site <- setdiff(union(na_em,na_as),na_ve_site) 
    
    for (k in 1:length(ha_test_second[,1])){# genetic distance in EM and antigenic sites
      ha_antigenic_d <- 0
      na_antigenic_d <- 0
      for (m in ha_site){# distance on other site
        if (is.na(ha_test_second[k,m+4])){ha_antigenic_d <- ha_antigenic_d +1}else{
          if (ha_test_second[k,m+4]==pre_hatest[1,m+1]){}else{
            ha_antigenic_d <- ha_antigenic_d+1
          }
        }
      }
      for (n in na_site){# distance on other site
        if (is.na(na_test_second[k,n+4])){na_antigenic_d <- na_antigenic_d +1}else{
          if (na_test_second[k,n+4]==pre_natest[1,n+1]){}else{
            na_antigenic_d <- na_antigenic_d+1
          }
        }
      }
      antigenic_dis <- ha_antigenic_d+na_antigenic_d
      total_antigenic_dis <- c(total_antigenic_dis,antigenic_dis)
    } 
    
    position <- which(total_antigenic_dis == min(total_antigenic_dis))
    ha_test_third <- ha_test_second[position,]
    na_test_third <- na_test_second[match(as.vector(ha_test_third$name),as.vector(na_test_second$name)),]
    total_hamming_dis <- c()
    ha_rest_site <- setdiff(c(1:566),union(ha_site,ha_ve_site)) 
    na_rest_site <- setdiff(c(1:469),union(na_site,na_ve_site))
    
    for (k in 1:length(ha_test_third[,1])){# genetic distance in the rest sites
      ha_hamming_d <- 0
      na_hamming_d <- 0
      for (m in ha_rest_site){# distance on other site
        if (is.na(ha_test_third[k,m+4])){ha_hamming_d <- ha_hamming_d +1}else{
          if (ha_test_third[k,m+4]==pre_hatest[1,m+1]){}else{
            ha_hamming_d <- ha_hamming_d+1
          }
        }
      }
      for (n in na_rest_site){# distance on other site
        if (is.na(na_test_third[k,n+4])){na_hamming_d <- na_hamming_d +1}else{
          if (na_test_third[k,n+4]==pre_natest[1,n+1]){}else{
            na_hamming_d <- na_hamming_d+1
          }
        }
      }
      hamming_dis <- ha_hamming_d+na_hamming_d
      total_hamming_dis <- c(total_hamming_dis,hamming_dis)
    } 
    position <- which(total_hamming_dis == min(total_hamming_dis))
    final_position <- position[length(position)]
    wildtype_ha <- rbind(wildtype_ha,ha_test_third[final_position,])
    wildtype_ha_name <- as.vector(ha_test_third[final_position,2])
    pre_naseq <- na_test_third[which(na_test_third$name==wildtype_ha_name),]
    wildtype_na <- rbind(wildtype_na,pre_naseq[1,])
  }
  wildtype_ha$season <- pre_year
  wildtype_na$season <- pre_year
  return(list(beth_ha = wildtype_ha, beth_na = wildtype_na))
}

h3n2_wildtype_result <- WTselection(input_ha_seq = ha_seq, input_na_seq = na_seq, pred_emset = read.csv("./toydata_result/HKH3_predicted_EMset.csv"),
                                    ha_pre_seq = read.csv("./toydata_result/BethH3_HA_predicted consensus seq.csv"),
                                    na_pre_seq = read.csv("./toydata_result/BethH3_NA_predicted consensus seq.csv"),
                                    ha_beta = -11, na_beta = -19.8, ha_as <- h3ha_as, na_as <- h3na_as, pre_year <- c(2012:2020))

# write.csv(h3n2_wildtype_result$beth_ha,file="./toydata_result/BethH3_HA_Wildtype_history tau.csv",row.names = F)
# write.csv(h3n2_wildtype_result$beth_na,file="./toydata_result/BethH3_NA_Wildtype_history tau.csv",row.names = F)

# output: wildtype vaccine strains


##################################################################################3
# Step 7: calculate average genetic mismatch between predicted vaccine strains and circulating strains
# input vaccine strain and circulating strains, and output genetic mismatch

GetMismatch <- function(input_pre_seq, input_cir_seq, test_site, test_year){
  cir_seq <- input_cir_seq
  pre_seq <- input_pre_seq
  pre_seq <- as.matrix(pre_seq)
  emsite <- test_site
  ave_mismatch <- c()
  sd <- c()
  for (z in 1:length(test_year)){
    obs_seq <- subset(cir_seq,season==test_year[z])
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
  Mismatch <- cbind(ave_mismatch,sd)
  rownames(Mismatch)<-test_year
  return(Mismatch)
}

ha_mismatch <- GetMismatch(input_pre_seq = h3n2_wildtype_result$beth_ha, input_cir_seq = ha_seq, test_site = c(1:566), test_year = c(2012:2020))
na_mismatch <- GetMismatch(input_pre_seq = h3n2_wildtype_result$beth_na, input_cir_seq = na_seq, test_site = c(1:469), test_year = c(2012:2020))

## output: average genetic mismatch tables

## end.
  
