
setwd("C:/Users/jingzhi/Desktop/Site-basedPredictionModel/data")

###################### Genetic mismatch by strain #######################
## The function AAMismatch allows to calculate mismatch between vaccine strain and circulating strain in ten regions

AAMismatch <- function(vaccine_type, region_test, subtype, segment, start_season, pre_year, fullsite, epitope){
  samplesize <- 0
  for (re in region_test){
    sequence <- read.csv(paste0("./validation set/",subtype,"/",re,segment,"_sequence.csv"))
    sequence <- sequence[which(sequence$season>start_season & sequence$season<=start_season+pre_year),]
    samplesize <- samplesize + length(sequence[,1])
  }
  MismatchTable <- matrix(nrow=samplesize)
  for(vaccine_seq in vaccine_type){
    mismatch_table <- matrix(ncol=4)
    pre_seq <- read.csv(paste0("./training set/",vaccine_seq,".csv"))
    for (z in 1:pre_year){
      for (re in region_test){
        sequence <- read.csv(paste0("./validation set/",subtype,"/",re,segment,"_sequence.csv"))
        obs_seq <- subset(sequence,season==z+start_season)
        obs_seq <- as.matrix(obs_seq)
        predict_seq <- subset(pre_seq,season==z+start_season)
        predict_seq <- as.matrix(predict_seq)
        mismatch <- c()
        if (length(obs_seq[,1])==0){
        }else{
          for (i in 1:length(obs_seq[,1])){
            distance <- 0
            site <- c()
            for (j in fullsite){
              if (is.na(obs_seq[i,j+4])){}else{
                if (predict_seq[1,j+4]==obs_seq[i,j+4]){}else{
                  distance <- distance+1
                  site <- c(site,j)}}
            }
            mismatch <- c(mismatch,distance)}
          year <- rep(start_season+z,length(mismatch))
          region <- rep(re,length(mismatch))
          if(re %in% c("NY","CA","CAN","NH")){
            continent <- rep("North America",length(mismatch))
          }
          if(re %in% c("UK","GER","FR")){
            continent <- rep("Europe",length(mismatch))
          }
          if(re %in% c("HK","SCP","SG","JP")){
            continent <- rep("Asia",length(mismatch))
          }
          mismatch <- cbind(year, region, continent,mismatch)
          colnames(mismatch)<-c("season","region","continent","mismatch")
          mismatch_table <- rbind(mismatch_table,mismatch)
        }
        
      }
    }
    mismatch_table <- mismatch_table[-1,]
    MismatchTable <- cbind(MismatchTable, mismatch_table)
  }  # mismatch in full length HA or NA
  for(vaccine_seq in vaccine_type){
    mismatch_table <- matrix(ncol=4)
    pre_seq <- read.csv(paste0("./training set/",vaccine_seq,".csv"))
    for (z in 1:pre_year){
      for (re in region_test){
        sequence <- read.csv(paste0("./validation set/",subtype,"/",re,segment,"_sequence.csv"))
        obs_seq <- subset(sequence,season==z+start_season)
        obs_seq <- as.matrix(obs_seq)
        predict_seq <- subset(pre_seq,season==z+start_season)
        predict_seq <- as.matrix(predict_seq)
        mismatch <- c()
        if (length(obs_seq[,1])==0){
        }else{
          for (i in 1:length(obs_seq[,1])){
            distance <- 0
            site <- c()
            for (j in epitope){
              if (is.na(obs_seq[i,j+4])){}else{
                if (predict_seq[1,j+4]==obs_seq[i,j+4]){}else{
                  distance <- distance+1
                  site <- c(site,j)}}
            }
            mismatch <- c(mismatch,distance)}
          year <- rep(start_season+z,length(mismatch))
          region <- rep(re,length(mismatch))
          if(re %in% c("NY","CA","CAN","NH")){
            continent <- rep("North America",length(mismatch))
          }
          if(re %in% c("UK","GER","FR")){
            continent <- rep("Europe",length(mismatch))
          }
          if(re %in% c("HK","SCP","SG","JP")){
            continent <- rep("Asia",length(mismatch))
          }
          mismatch <- cbind(year, region, continent,mismatch)
          colnames(mismatch)<-c("season","region","continent","mismatch")
          mismatch_table <- rbind(mismatch_table,mismatch)
        }
        
      }
    }
    mismatch_table <- mismatch_table[-1,]
    MismatchTable <- cbind(MismatchTable, mismatch_table)
  }  # mismatch in epitopes of HA or NA
  MismatchTable <- MismatchTable[,c(2:5,seq(9,(8*length(vaccine_type)+1),4))]
  MismatchTable <- as.data.frame(MismatchTable)
  colnames(MismatchTable)[4:11] <- c("who_full_mismatch","beth_full_mismatch","answer_full_mismatch", "lbi_full_mismatch",
                                       "who_epitope_mismatch","beth_epitope_mismatch","answer_epitope_mismatch", "lbi_epitope_mismatch")
  return(MismatchTable)
}

## epitopes

h1ha_epitope <- c(141,142,c(170:174),c(176:181),c(201:212),c(183:187),220,221,222,253,254,252,c(154:159),238,239,87,88,90,91,92,132)
h1na_epitope <- c(93,94,95,216,217,219,220,221,250,251,252,254,c(262:268),270,355,358,375,377,378,388,389,449,450,451)
h3ha_epitope <- c(122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198,44,45,46,47,48,50,51,53,54,273,275,276,278,279,280,294,297,299,300,304,305,307,308,309,310,311,312,96,102,103,117,121,167,170,171,172,173,174,175,176,177,179,182,201,203,205,207,208,209,212,213,214,215,216,217,218,219,226,227,228,229,230,240,242,244,246,247,248,57,59,62,63,67,75,78,80,81,82,83,86,87,88,91,92,94,109,260,261,262,265)+16
h3na_epitope <- c(150,198,199,220,221,253,329,334,344,368,370,403)
ha_full <- c(1:566)
na_full <- c(1:469)

## get mismatch table

h1ha_mismatch_table <- AAMismatch(vaccine_type = c("WHOh1n1_HA","BethH1N1_Wildtype_HA_history tau","NHH1_HA_answer strain", "Neher_H1HA_sequence"), 
                                region_test = c("NY","CA","CAN","UK","GER","FR","HK","SCP","JP","SG"),
                                subtype = "H1N1", segment = "H1_HA", start_season = 2012, pre_year = 8, fullsite = ha_full, epitope= h1ha_epitope)

h1na_mismatch_table <- AAMismatch(vaccine_type = c("WHOh1n1_NA","BethH1N1_Wildtype_NA_history tau","NHH1_NA_answer strain", "Neher_H1NA_sequence"), 
                                  region_test = c("NY","CA","CAN","UK","GER","FR","HK","SCP","JP","SG"),
                                  subtype = "H1N1", segment = "H1_NA", start_season = 2012, pre_year = 8, fullsite = na_full, epitope= h1na_epitope)

h3ha_mismatch_table <- AAMismatch(vaccine_type = c("WHOh3n2_HA","BethH3N2_Wildtype_HA_history tau","NHH3_HA_answer strain", "Neher_H3HA_sequence"), 
                                  region_test = c("NY","CA","CAN","UK","GER","FR","HK","SCP","JP","SG"),
                                  subtype = "H3N2", segment = "H3_HA", start_season = 2002, pre_year = 18, fullsite = ha_full, epitope= h3ha_epitope)

h3na_mismatch_table <- AAMismatch(vaccine_type = c("WHOh3n2_NA","BethH3N2_Wildtype_NA_history tau","NHH3_NA_answer strain", "Neher_H3NA_sequence"), 
                                  region_test = c("NY","CA","CAN","UK","GER","FR","HK","SCP","JP","SG"),
                                  subtype = "H3N2", segment = "H3_NA", start_season = 2002, pre_year = 18, fullsite = na_full, epitope= h3na_epitope)

## output

write.csv(h1ha_mismatch_table,file="h1ha mismatch table.csv",row.names = F)
write.csv(h1na_mismatch_table,file="h1na mismatch table.csv",row.names = F)
write.csv(h3ha_mismatch_table,file="h3ha mismatch table.csv",row.names = F)
write.csv(h3na_mismatch_table,file="h3na mismatch table.csv",row.names = F)



###################### make average mismatch table by year and by region #######################
## The function SummaryMismatch allows to summary mismatch by year and by region

SummaryMismatch <- function(mismatch_table){
  mismatch_type <- c("who_full_mismatch","beth_full_mismatch","answer_full_mismatch", "lbi_full_mismatch",
    "who_epitope_mismatch","beth_epitope_mismatch","answer_epitope_mismatch", "lbi_epitope_mismatch")
  ave_mismatch_table <- matrix(ncol=10)
  col_number <- c(1:length(mismatch_type))
  for (col in col_number){
    ave_mismatch_table <- rbind(ave_mismatch_table,c("NY","CA","CAN","UK","GER","FR","HK","SCP","JP","SG"))
    for (year in unique(mismatch_table$season)){
      yearly_mismatch <- c()
      for(re in c("NY","CA","CAN","UK","GER","FR","HK","SCP","JP","SG")){
        test_mismatch <- mismatch_table[which(mismatch_table$season==year & mismatch_table$region==re),]
        if(length(test_mismatch[,1])==0){
          yearly_mismatch <- c(yearly_mismatch,NA)
        }else{
          yearly_mismatch <- c(yearly_mismatch,mean(test_mismatch[,col+3]))
        }
      }
      ave_mismatch_table <- rbind(ave_mismatch_table,yearly_mismatch)
    }
  }
  ave_mismatch_table <- ave_mismatch_table[-1,]
  rownames(ave_mismatch_table) <- c("who_full_mismatch",unique(mismatch_table$season),"beth_full_mismatch",unique(mismatch_table$season),
                                    "answer_full_mismatch",unique(mismatch_table$season),"lbi_full_mismatch",unique(mismatch_table$season),
                                    "who_epitope_mismatch",unique(mismatch_table$season),"beth_epitope_mismatch",unique(mismatch_table$season),
                                    "answer_epitope_mismatch",unique(mismatch_table$season),"lbi_epitope_mismatch",unique(mismatch_table$season))
  return(ave_mismatch_table)
}

## get average mismatch by year and by region

h1ha_summary_mismatch <- SummaryMismatch(mismatch_table = read.csv("h1ha mismatch table.csv"))
h1na_summary_mismatch <- SummaryMismatch(mismatch_table = read.csv("h1na mismatch table.csv"))
h3ha_summary_mismatch <- SummaryMismatch(mismatch_table = read.csv("h3ha mismatch table.csv"))
h3na_summary_mismatch <- SummaryMismatch(mismatch_table = read.csv("h3na mismatch table.csv"))

## output

write.csv(h1ha_summary_mismatch,file="h1ha_ave_mismatch_table.csv")
write.csv(h1na_summary_mismatch,file="h1na_ave_mismatch_table.csv")
write.csv(h3ha_summary_mismatch,file="h3ha_ave_mismatch_table.csv")
write.csv(h3na_summary_mismatch,file="h3na_ave_mismatch_table.csv")
