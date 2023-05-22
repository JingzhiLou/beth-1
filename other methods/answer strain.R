
###########################################################
# The function GetAnswerCon allows to calculate consensus sequence in T+1 flu season

GetAnswerCon <- function(input_csv_seq, predict_season){
  Pred_Conseq <- matrix(0,ncol=length(input_csv_seq[1,])-4)
  for (pre_season in predict_season){
    print(pre_season)
    seq <- input_csv_seq[which(input_csv_seq$season == pre_season),]
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
    rownames(codata)<-year # get site-wise prevalence table
    
    colname <- colnames(codata)
    test <- strsplit(colname,split=".",fixed=T)
    xname <- c() 
    yname <- c() 
    for (j in 1:length(test)){
      xname <- c(xname,test[[j]][1])
      yname <- c(yname,test[[j]][2])}
    xuniq <- unique(xname)
    preseq <- rbind(xname,yname,codata[1,])
    pred_seq <- matrix(0,nrow=3)
    for (l in xuniq){
      xpos <- which(xname==l)
      posmax <- max(codata[1,min(xpos):max(xpos)])
      for (m in xpos){
        if (preseq[3,m]==posmax){
          pred_seq <- cbind(pred_seq,preseq[,m])
          break 
        }}}
    pred_seq <- pred_seq[,-1]
    Pred_Conseq <- rbind(Pred_Conseq,pred_seq[2,])
  }
  Pred_Conseq <- Pred_Conseq[-1,]
  rownames(Pred_Conseq) <- predict_season
  return(Pred_Conseq)
}

ha_answer_conseq <- GetAnswerCon(input_csv_seq = read.csv("H3N2_HA_HK_sequence.csv"), predict_season = c(2012:2020))
na_answer_conseq <- GetAnswerCon(input_csv_seq = read.csv("H3N2_NA_HK_sequence.csv"), predict_season = c(2012:2020))

# Output: consensus answer strain in T+1 flu season

write.csv(ha_answer_conseq,file="H3N2_HA_HK_answer consensus.csv")
write.csv(na_answer_conseq,file="H3N2_NA_HK_answer consensus.csv")


# The function GetAnswerWT allows to select wildtype answer strain by minimizing Hamming distance of both HA and NA to the consensus sequence in time T+1

GetAnswerWT<- function(input_ha_seq, input_na_seq, ha_pre_seq, na_pre_seq, pre_year){
  seq1 <- input_ha_seq
  seq2 <- input_na_seq
  name1 <- str_trim(as.vector(seq1$name))
  name2 <- str_trim(as.vector(seq2$name))
  name <- intersect(name1,name2)
  ha_seq <- seq1[which(seq1$name %in% name),]
  na_seq <- seq2[which(seq2$name %in% name),]
  colnames(ha_pre_seq)[1] <- "season"
  colnames(na_pre_seq)[1] <- "season"
  alpha <- 0.5 # equal wight for HA and NA distance
  wildtype_ha <- matrix(ncol=length(ha_seq[1,]),nrow=0)
  wildtype_na <- matrix(ncol=length(na_seq[1,]),nrow=0)
  EMD_total <- matrix(ncol=3,nrow=0)
  for (s in pre_year){
    print(s)
    ha_test <- ha_seq[which(ha_seq$season==s),]
    na_test <- na_seq[which(na_seq$season==s),]
    pre_hatest <- ha_pre_seq[which(ha_pre_seq$season==s),]
    pre_natest <- na_pre_seq[which(na_pre_seq$season==s),]
    total_emd <- c()
    EMD_ha <- c()
    EMD_na <- c()
    site_ha <- c()
    site_na <- c()
    for (k in 1:length(ha_test[,1])){
      name_test <- as.vector(ha_test[k,2])
      na_corresp <- na_test[which(na_test$name==name_test),]
      ha_emd <- 0
      na_emd <- 0
      for (m in c(1:566)){
        if (is.na(ha_test[k,m+4])){ha_emd <- ha_emd +1}else{
          if (ha_test[k,m+4]==pre_hatest[1,m+1]){}else{
            ha_emd <- ha_emd+1
            site_ha <- c(site_ha,m)
          }
        }
      }
      for (n in c(1:469)){
        if (is.na(na_corresp[1,n+4])){}else{
          if (na_corresp[1,n+4]==pre_natest[1,n+1]){}else{
            na_emd <- na_emd+1
            site_ha <- c(site_ha,n)
          }
        }
      }
      emd <- alpha*ha_emd+(1-alpha)*na_emd
      total_emd <- c(total_emd,emd)
      EMD_ha <- c(EMD_ha,ha_emd)
      EMD_na <- c(EMD_na,na_emd)
    }
    position <- which.min(total_emd)
    wildtype_ha <- rbind(wildtype_ha,ha_test[position,])
    wildtype_ha_name <- as.vector(ha_test[position,2])
    pre_naseq <- na_test[which(na_test$name==wildtype_ha_name),]
    wildtype_na <- rbind(wildtype_na,pre_naseq[1,])
    EMD <- c(EMD_ha[position],EMD_na[position],total_emd[position])
    EMD_total <- rbind(EMD_total,EMD)
  }
  return(list(answer_ha = wildtype_ha, answer_na = wildtype_na))
}


answer_wildtype_result <- GetAnswerWT(input_ha_seq = read.csv("H3N2_HA_HK_sequence.csv"), 
                                      input_na_seq = read.csv("H3N2_NA_HK_sequence.csv"), 
                                      ha_pre_seq <- read.csv("H3N2_HA_HK_answer consensus.csv"),
                                      na_pre_seq <- read.csv("H3N2_NA_HK_answer consensus.csv"),
                                      pre_year = c(2012:2020))

# Output: wildtype answer strain in T+1 flu season

write.csv(answer_wildtype_result$answer_ha,file="H3N2_HA_HK_answer wildtype.csv",row.names = F)
write.csv(answer_wildtype_result$answer_na,file="H3N2_NA_HK_answer wildtype.csv",row.names = F)
