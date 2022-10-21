
# site-based prevalence prediction and wild type selection
# reference: 1. Cao L, Lou J, Zhao S, et al. In silico prediction of influenza vaccine effectiveness by sequence analysis. Vaccine. 2021;39(7):1030-1034. doi:10.1016/j.vaccine.2021.01.006
#            2. Cao L, Zhao S, Lou J, et al. Differential Influence of Age on the Relationship between Genetic Mismatch and A(H1N1)pdm09 Vaccine Effectiveness. Viruses. 2021;13(4):619. Published 2021 Apr 4. doi:10.3390/v13040619


# step 1: prediction of future mutation prevalence and consensus strain
# input site wise prevalence table and output yearly predicted consensus strains

inputdata_name <- c("H3N2_HA_HK_prevalence.csv")
outputdata_name <- c("HA_predicted_consensus_sequence.csv")
transtime <- 2 # transition time is an integer and is adopted from the last step EstimateTransitionTime
pred_year <- 10 # define # of predicted years, and it requires (transition time + 1) years data to burn in

Prev_Table <- read.csv(inputdata_name)
last_year <- Prev_Table[length(Prev_Table[,1]),1]
Prev_Table <- Prev_Table[,-1]
Pred_Prevalence <- matrix(0,ncol=length(Prev_Table))
colnames(Pred_Prevalence) <- colnames(Prev_Table)
for (year in pred_year:1){
  pred_prev <- c()
  Test_PrevTable <- Prev_Table[1:(length(Prev_Table[,1])-year+1),]
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
for (k in 1:pred_year){
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
Pred_Sequence <- Pred_Sequence[seq(3,3*pred_year,3),]
rownames(Pred_Sequence) <- sort(c((last_year+1):(last_year-pred_year+2)))
write.csv(Pred_Sequence,file=outputdata_name)

# Note: if you want to output predicted site wise prevalence, please output matrix Pred_Prevalence



# step 2: select wild type strain based on dynamic predictor sites set W(t)
# input yearly predicted consensus strains and output wildtype sequences
# Note: the EM sites and predictor sites include the signal piptide and start with M(Methionine)

ha_seq <- read.csv("H3N2_HA_HK_sequence.csv")
na_seq <- read.csv("H3N2_NA_HK_sequence.csv")
ha_pre_seq <- read.csv("HA_predicted_consensus_sequence.csv")
na_pre_seq <- read.csv("NA_predicted_consensus_sequence.csv")
All_emset <- read.csv("H3N2_HK_EMset.csv")
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
pre_year <- c(2012:2021)
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

    for (n in union(na_em,na_as)){# current EM set of NA
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
wildtype_ha[,4] <- wildtype_ha[,4]+1
wildtype_na[,4] <- wildtype_na[,4]+1
write.csv(wildtype_ha,file="CombinedH3N2_HA_Wildtype.csv",row.names = F)
write.csv(wildtype_na,file="CombinedH3N2_NA_Wildtype.csv",row.names = F)



