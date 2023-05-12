
setwd("C:/Users/jingzhi/Desktop/Site-basedPredictionModel/data/training set/SEA/H3N2")

###################### Prediction of future mutation prevalence and consensus strain #######################
## The function PreYearlyStrain allows to calculate one-step prediction using beth-1 and history tau, and output predicted consensus strain and predicted EM set

library(plyr)
library(lubridate)

PreYearlyStrain <- function(input_csv_seq, predict_season,transition_time,his_EMset){
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

## for H1N1, predict_year=2012-2019 ha_tau=2 na_tau=2
## the history EM set summarizing history transition events is obtained by optimal (theta,h) 

ha_seq <- read.csv("SEAH1_HA_sequence.csv")
na_seq <- read.csv("SEAH1_NA_sequence.csv")
his_emset <- read.csv("SEAH1_history_EMset.csv")

ha_predicted_conseq <- PreYearlyStrain(input_csv_seq = ha_seq, transition_time = 2, predict_season = c(2012:2019), his_EMset = his_emset[,c(1:3)])
na_predicted_conseq <- PreYearlyStrain(input_csv_seq = na_seq, transition_time = 2, predict_season = c(2012:2019), his_EMset = his_emset[,c(4:6)])

## output: predicted consensus sequence and predicted EM set

write.csv(ha_predicted_conseq$beth_conseq,file="BethH1_HA_predicted consensus seq_history tau.csv")
write.csv(na_predicted_conseq$beth_conseq,file="BethH1_NA_predicted consensus seq_history tau.csv")
write.csv(ha_predicted_conseq$Predicted_EM,file="SEAH1_predicted_EMset_HA.csv")
write.csv(na_predicted_conseq$Predicted_EM,file="SEAH1_predicted_EMset_NA.csv")


## for H3N2, predict_year=2002-2019 ha_tau=2 na_tau=2

ha_seq <- read.csv("SEAH3_HA_sequence.csv")
na_seq <- read.csv("SEAH3_NA_sequence.csv")
his_emset <- read.csv("SEAH3_history_EMset.csv")

ha_predicted_conseq <- PreYearlyStrain(input_csv_seq = ha_seq, transition_time = 2, predict_season = c(2002:2019), his_EMset = his_emset[,c(1:3)])
na_predicted_conseq <- PreYearlyStrain(input_csv_seq = na_seq, transition_time = 2, predict_season = c(2002:2019), his_EMset = his_emset[,c(4:6)])

## output: predicted consensus sequence and predicted EM set

write.csv(ha_predicted_conseq$beth_conseq,file="BethH3_HA_predicted consensus seq_history tau.csv")
write.csv(na_predicted_conseq$beth_conseq,file="BethH3_NA_predicted consensus seq_history tau.csv")
write.csv(ha_predicted_conseq$Predicted_EM,file="SEAH3_predicted_EMset_HA.csv")
write.csv(na_predicted_conseq$Predicted_EM,file="SEAH3_predicted_EMset_NA.csv")



###################### Wildtype vaccine strain selection based on VE-GD model #######################
## The function WTselection allows to select wildtype strain integrating mutation information on both HA and NA

library(stringr)

## H1N1 epitopes
h1ha_as <- c(141,142,c(170:174),c(176:181),c(201:212),c(183:187),220,221,222,252,253,254,c(154:159),238,239,87,88,90,91,92,132)
h1na_as <- c(93,94,95,216,217,219,220,221,250,251,252,254,c(262:268),270,355,358,375,377,378,388,389,449,450,451)
## H3N2 epitopes
h3ha_as_ab <- c(122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198)+16
h3ha_as <- c(122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198,44,45,46,47,48,50,51,53,54,273,275,276,278,279,280,294,297,299,300,304,305,307,308,309,310,311,312,96,102,103,117,121,167,170,171,172,173,174,175,176,177,179,182,201,203,205,207,208,209,212,213,214,215,216,217,218,219,226,227,228,229,230,240,242,244,246,247,248,57,59,62,63,67,75,78,80,81,82,83,86,87,88,91,92,94,109,260,261,262,265)+16
h3na_as <- c(150,198,199,220,221,253,329,334,344,368,370,403)  


## for H1N1 predict_year 2012-2019 ha_beta = -8.43, na_beta = -10.6
H1N1_WTselection <- function(input_ha_seq, input_na_seq, his_emset, pred_emset, ha_pre_seq, na_pre_seq, ha_beta, na_beta, ha_as, na_as, pre_year){
  seq1 <- input_ha_seq
  seq2 <- input_na_seq
  name1 <- str_trim(as.vector(seq1$name))
  name2 <- str_trim(as.vector(seq2$name))
  name <- intersect(name1,name2)
  ha_seq <- seq1[which(seq1$name %in% name),]
  na_seq <- seq2[which(seq2$name %in% name),]
  ha_predem <- pred_emset[,c(1:2)] 
  na_predem <- pred_emset[,c(3:4)]
  wildtype_ha <- matrix(ncol=length(ha_seq[1,]),nrow=0)
  wildtype_na <- matrix(ncol=length(na_seq[1,]),nrow=0)
  for (s in pre_year){
    print(s)
    ha_test <- ha_seq[which(ha_seq$date<make_date(s-1,3,1)),] 
    na_test <- na_seq[which(na_seq$date<make_date(s-1,3,1)),] 
    pre_hatest <- ha_pre_seq[which(ha_pre_seq$season==s),]
    pre_natest <- na_pre_seq[which(na_pre_seq$season==s),]
    ha_em <- unique(ha_predem[which(ha_predem[,2]<=s),1])
    na_em <- unique(na_predem[which(na_predem[,2]<=s),1])
    total_VE_EMD <- c()
    EMD_ha <- c()
    EMD_na <- c()
    ha_ve_site <- intersect(ha_em,ha_as) # H1N1: EMs on antigenic sites, H3N2: EMs on antigenic sites A & B  ¡Ì¡Ì¡Ì
    na_ve_site <- intersect(na_em,na_as) # H1N1: EMs on antigenic sites, H3N2: antigenic sites  ¡Ì¡Ì¡Ì
    
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

h1n1_wildtype_result <- H1N1_WTselection(input_ha_seq = read.csv("SEAH1_HA_sequence.csv"), input_na_seq = read.csv("SEAH1_NA_sequence.csv"),
                               his_emset = read.csv("SEAH1_history_EMset.csv"), pred_emset = read.csv("SEAH1_predicted_EMset.csv"),
                               ha_pre_seq = read.csv("BethH1_HA_predicted consensus seq_history tau.csv"),
                               na_pre_seq = read.csv("BethH1_NA_predicted consensus seq_history tau.csv"),
                               ha_beta = -8.43, na_beta = -10.6, ha_as <- h1ha_as, na_as <- h1na_as, pre_year <- c(2012:2019))

## output: wildtype strain

write.csv(h1n1_wildtype_result$beth_ha,file="BethH1N1_Wildtype_HA_history tau.csv",row.names = F)
write.csv(h1n1_wildtype_result$beth_na,file="BethH1N1_Wildtype_NA_history tau.csv",row.names = F)


## for H3N2 predict_year 2002-2019 ha_beta = -11, na_beta = -19.8
H3N2_WTselection <- function(input_ha_seq, input_na_seq, his_emset, pred_emset, ha_pre_seq, na_pre_seq, ha_beta, na_beta, ha_as, na_as, pre_year){
  seq1 <- input_ha_seq
  seq2 <- input_na_seq
  name1 <- str_trim(as.vector(seq1$name))
  name2 <- str_trim(as.vector(seq2$name))
  name <- intersect(name1,name2)
  ha_seq <- seq1[which(seq1$name %in% name),]
  na_seq <- seq2[which(seq2$name %in% name),]
  ha_predem <- pred_emset[,c(1:2)] 
  na_predem <- pred_emset[,c(3:4)]
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

h3n2_wildtype_result <- H3N2_WTselection(input_ha_seq = read.csv("SEAH3_HA_sequence.csv"), input_na_seq = read.csv("SEAH3_NA_sequence.csv"),
                                    his_emset = read.csv("SEAH3_history_EMset.csv"), pred_emset = read.csv("SEAH3_predicted_EMset.csv"),
                                    ha_pre_seq = read.csv("BethH3_HA_predicted consensus seq_history tau.csv"),
                                    na_pre_seq = read.csv("BethH3_NA_predicted consensus seq_history tau.csv"),
                                    ha_beta = -11, na_beta = -19.8, ha_as <- h3ha_as, na_as <- h3na_as, pre_year <- c(2002:2019))

## output: wildtype strain

write.csv(h3n2_wildtype_result$beth_ha,file="BethH3N2_Wildtype_HA_history tau.csv",row.names = F)
write.csv(h3n2_wildtype_result$beth_na,file="BethH3N2_Wildtype_NA_history tau.csv",row.names = F)




