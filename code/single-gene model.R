
setwd("C:/Users/jingzhi/Desktop/Site-basedPredictionModel/data/training set/SEA/H3N2")

###################### single-protein model #######################
## The function HAonlyWT and NAonlyWT allow to select wildtype by Hamming distance of single gene

library(stringr)

## H1N1 epitopes
h1ha_as <- c(141,142,c(170:174),c(176:181),c(201:212),c(183:187),220,221,222,252,253,254,c(154:159),238,239,87,88,90,91,92,132)
h1na_as <- c(93,94,95,216,217,219,220,221,250,251,252,254,c(262:268),270,355,358,375,377,378,388,389,449,450,451)
## H3N2 epitopes
h3ha_as_ab <- c(122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198)+16
h3ha_as <- c(122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198,44,45,46,47,48,50,51,53,54,273,275,276,278,279,280,294,297,299,300,304,305,307,308,309,310,311,312,96,102,103,117,121,167,170,171,172,173,174,175,176,177,179,182,201,203,205,207,208,209,212,213,214,215,216,217,218,219,226,227,228,229,230,240,242,244,246,247,248,57,59,62,63,67,75,78,80,81,82,83,86,87,88,91,92,94,109,260,261,262,265)+16
h3na_as <- c(150,198,199,220,221,253,329,334,344,368,370,403)  

HAonlyWT <- function(input_ha_seq, input_na_seq, his_emset, pred_emset, ha_pre_seq, na_pre_seq, ha_as, na_as, pre_year,subtype){
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
    pre_hatest <- ha_pre_seq[which(ha_pre_seq$season==s),]
    ha_em <- unique(ha_predem[which(ha_predem[,2]<=s),1])
    total_EMD <- c()
    if(subtype=="H1N1"){
      ha_ve_site <- intersect(ha_em,ha_as)
    }else{
      ha_ve_site <- intersect(ha_em,h3ha_as_ab)
    }
    
    for (k in 1:length(ha_test[,1])){ # genetic distance in dynamic predictor set
      ha_emd <- 0
      for (m in ha_ve_site){ 
        if (is.na(ha_test[k,m+4])){ha_emd <- ha_emd +1}else{
          if (ha_test[k,m+4]==pre_hatest[1,m+1]){}else{
            ha_emd <- ha_emd+1
          }
        }
      } 
      emd <- ha_emd
      total_EMD <- c(total_EMD,emd)
    } #select the wildtype sequence based on genetic distance on EMs
    
    position <- which(total_EMD == min(total_EMD))
    ha_test_second <- ha_test[position,]
    total_antigenic_dis <- c()
    ha_site <- setdiff(union(ha_em,ha_as),ha_ve_site) 
    
    for (k in 1:length(ha_test_second[,1])){# genetic distance in EM and antigenic sites
      ha_antigenic_d <- 0
      for (m in ha_site){# distance on other site
        if (is.na(ha_test_second[k,m+4])){ha_antigenic_d <- ha_antigenic_d +1}else{
          if (ha_test_second[k,m+4]==pre_hatest[1,m+1]){}else{
            ha_antigenic_d <- ha_antigenic_d+1
          }
        }
      }
      
      antigenic_dis <- ha_antigenic_d
      total_antigenic_dis <- c(total_antigenic_dis,antigenic_dis)
    } 
    
    position <- which(total_antigenic_dis == min(total_antigenic_dis))
    ha_test_third <- ha_test_second[position,]
    total_hamming_dis <- c()
    ha_rest_site <- setdiff(c(1:566),union(ha_site,ha_ve_site)) 
    
    for (k in 1:length(ha_test_third[,1])){# genetic distance in the rest sites
      ha_hamming_d <- 0
      for (m in ha_rest_site){# distance on other site
        if (is.na(ha_test_third[k,m+4])){ha_hamming_d <- ha_hamming_d +1}else{
          if (ha_test_third[k,m+4]==pre_hatest[1,m+1]){}else{
            ha_hamming_d <- ha_hamming_d+1
          }
        }
      }
      
      hamming_dis <- ha_hamming_d
      total_hamming_dis <- c(total_hamming_dis,hamming_dis)
    } 
    position <- which(total_hamming_dis == min(total_hamming_dis))
    final_position <- position[length(position)]
    wildtype_ha <- rbind(wildtype_ha,ha_test_third[final_position,])
    wildtype_ha_name <- as.vector(ha_test_third[final_position,2])
    pre_naseq <- na_seq[which(na_seq$name==wildtype_ha_name),]
    wildtype_na <- rbind(wildtype_na,pre_naseq[1,])
  }
  wildtype_ha$season <- pre_year
  wildtype_na$season <- pre_year
  return(list(beth_ha = wildtype_ha, beth_na = wildtype_na))
}

h1n1_haonly_result <- HAonlyWT(input_ha_seq = read.csv("SEAH1_HA_sequence.csv"), input_na_seq = read.csv("SEAH1_NA_sequence.csv"),
                                         his_emset = read.csv("SEAH1_history_EMset.csv"), pred_emset = read.csv("SEAH1_predicted_EMset.csv"),
                                         ha_pre_seq = read.csv("BethH1_HA_predicted consensus seq_history tau.csv"),
                                         na_pre_seq = read.csv("BethH1_NA_predicted consensus seq_history tau.csv"),
                                         ha_as <- h1ha_as, na_as <- h1na_as, pre_year <- c(2012:2019), subtype="H1N1")

h3n2_haonly_result <- HAonlyWT(input_ha_seq = read.csv("SEAH3_HA_sequence.csv"), input_na_seq = read.csv("SEAH3_NA_sequence.csv"),
                               his_emset = read.csv("SEAH3_history_EMset.csv"), pred_emset = read.csv("SEAH3_predicted_EMset.csv"),
                               ha_pre_seq = read.csv("BethH3_HA_predicted consensus seq_history tau.csv"),
                               na_pre_seq = read.csv("BethH3_NA_predicted consensus seq_history tau.csv"),
                               ha_as <- h3ha_as, na_as <- h3na_as, pre_year <- c(2002:2019), subtype="H3N2")

## output

write.csv(h1n1_haonly_result$beth_ha,file="BethH1N1_HAonly model_HA_wildtype.csv",row.names = F)
write.csv(h1n1_haonly_result$beth_na,file="BethH1N1_HAonly model_NA_wildtype.csv",row.names = F)
write.csv(h3n2_haonly_result$beth_ha,file="BethH3N2_HAonly model_HA_wildtype.csv",row.names = F)
write.csv(h3n2_haonly_result$beth_na,file="BethH3N2_HAonly model_NA_wildtype.csv",row.names = F)



NAonlyWT <- function(input_ha_seq, input_na_seq, his_emset, pred_emset, ha_pre_seq, na_pre_seq, ha_as, na_as, pre_year,subtype){
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
    na_test <- na_seq[which(na_seq$date<make_date(s-1,3,1)),] 
    pre_natest <- na_pre_seq[which(na_pre_seq$season==s),]
    na_em <- unique(na_predem[which(na_predem[,2]<=s),1])
    total_EMD <- c()
    if(subtype=="H1N1"){
      na_ve_site <- intersect(na_em,na_as)
    }else{
      na_ve_site <- na_as
    }
    
    for (k in 1:length(na_test[,1])){ # genetic distance in dynamic predictor set
      na_emd <- 0
      na_emd <- 0
      for (m in na_ve_site){ 
        if (is.na(na_test[k,m+4])){na_emd <- na_emd +1}else{
          if (na_test[k,m+4]==pre_natest[1,m+1]){}else{
            na_emd <- na_emd+1
          }
        }
      } 
      emd <- na_emd
      total_EMD <- c(total_EMD,emd)
    } #select the wildtype sequence based on genetic distance on EMs
    
    position <- which(total_EMD == min(total_EMD))
    na_test_second <- na_test[position,]
    total_antigenic_dis <- c()
    na_site <- setdiff(union(na_em,na_as),na_ve_site) 
    
    for (k in 1:length(na_test_second[,1])){# genetic distance in EM and antigenic sites
      na_antigenic_d <- 0
      for (m in na_site){# distance on other site
        if (is.na(na_test_second[k,m+4])){na_antigenic_d <- na_antigenic_d +1}else{
          if (na_test_second[k,m+4]==pre_natest[1,m+1]){}else{
            na_antigenic_d <- na_antigenic_d+1
          }
        }
      }
      
      antigenic_dis <- na_antigenic_d
      total_antigenic_dis <- c(total_antigenic_dis,antigenic_dis)
    } 
    
    position <- which(total_antigenic_dis == min(total_antigenic_dis))
    na_test_third <- na_test_second[position,]
    total_hamming_dis <- c()
    na_rest_site <- setdiff(c(1:469),union(na_site,na_ve_site)) 
    
    for (k in 1:length(na_test_third[,1])){# genetic distance in the rest sites
      na_hamming_d <- 0
      for (m in na_rest_site){# distance on other site
        if (is.na(na_test_third[k,m+4])){na_hamming_d <- na_hamming_d +1}else{
          if (na_test_third[k,m+4]==pre_natest[1,m+1]){}else{
            na_hamming_d <- na_hamming_d+1
          }
        }
      }
      
      hamming_dis <- na_hamming_d
      total_hamming_dis <- c(total_hamming_dis,hamming_dis)
    } 
    position <- which(total_hamming_dis == min(total_hamming_dis))
    final_position <- position[length(position)]
    wildtype_na <- rbind(wildtype_na,na_test_third[final_position,])
    wildtype_na_name <- as.vector(na_test_third[final_position,2])
    pre_haseq <- ha_seq[which(ha_seq$name==wildtype_na_name),]
    wildtype_ha <- rbind(wildtype_ha,pre_haseq[1,])
  }
  wildtype_ha$season <- pre_year
  wildtype_na$season <- pre_year
  return(list(beth_ha = wildtype_ha, beth_na = wildtype_na))
}

h1n1_naonly_result <- NAonlyWT(input_ha_seq = read.csv("SEAH1_HA_sequence.csv"), input_na_seq = read.csv("SEAH1_NA_sequence.csv"),
                               his_emset = read.csv("SEAH1_history_EMset.csv"), pred_emset = read.csv("SEAH1_predicted_EMset.csv"),
                               ha_pre_seq = read.csv("BethH1_HA_predicted consensus seq_history tau.csv"),
                               na_pre_seq = read.csv("BethH1_NA_predicted consensus seq_history tau.csv"),
                               ha_as <- h1ha_as, na_as <- h1na_as, pre_year <- c(2012:2019), subtype="H1N1")

h3n2_naonly_result <- NAonlyWT(input_ha_seq = read.csv("SEAH3_HA_sequence.csv"), input_na_seq = read.csv("SEAH3_NA_sequence.csv"),
                               his_emset = read.csv("SEAH3_history_EMset.csv"), pred_emset = read.csv("SEAH3_predicted_EMset.csv"),
                               ha_pre_seq = read.csv("BethH3_HA_predicted consensus seq_history tau.csv"),
                               na_pre_seq = read.csv("BethH3_NA_predicted consensus seq_history tau.csv"),
                               ha_as <- h3ha_as, na_as <- h3na_as, pre_year <- c(2002:2019), subtype="H3N2")

## output

write.csv(h1n1_naonly_result$beth_ha,file="BethH1N1_NAonly model_HA_wildtype.csv",row.names = F)
write.csv(h1n1_naonly_result$beth_na,file="BethH1N1_NAonly model_NA_wildtype.csv",row.names = F)
write.csv(h3n2_naonly_result$beth_ha,file="BethH3N2_NAonly model_HA_wildtype.csv",row.names = F)
write.csv(h3n2_naonly_result$beth_na,file="BethH3N2_NAonly model_NA_wildtype.csv",row.names = F)

