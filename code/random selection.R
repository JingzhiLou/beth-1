
setwd("C:/Users/jingzhi/Desktop/Site-basedPredictionModel/data")

###################### random selection #######################
## The function RandomPickStrain allows to randomly select a sequence in time T as a prediction and calculate genetic mismatch (repeat 10 times)

## for h1n1, random seed=38

load("./training set/NH/H1N1/NHH1_HA_sequence.Rdata")
load("./training set/NH/H1N1/NHH1_NA_sequence.Rdata")

RandomPickStrain <- function(input_csv_seq, predict_season, repeat_time, emsite, region, subtype){
  mat_yearly_mismatch <- matrix(ncol = repeat_time)
  for (pre_season in predict_season){
    print(pre_season)
    seq <- input_csv_seq[which(input_csv_seq$season == (pre_season-1) & input_csv_seq$date < make_date(pre_season-1,3,1)),]
    set.seed(38) #set a seed for random pick
    seed <- sample(1:length(seq[,1]),repeat_time,replace=F)
    ran_seq <- seq[seed,]
    ave_yearly_mismatch <- c()
    for(ran_i in 1:length(ran_seq[,1])){
      ave_difference <- c()
      for (re in region){
        sequence <- read.csv(paste0("./validation set/H1N1/",re,subtype))
        obs_seq <- subset(sequence,season==pre_season)
        obs_seq <- as.matrix(obs_seq)
        diff <- 0
        if (length(obs_seq[,1])==0){
          ave_difference <- c(ave_difference,NA)}else{
            for (i in 1:length(obs_seq[,1])){
              d <- 0
              for (j in emsite){
                if (is.na(obs_seq[i,j+4])){}else{
                  if (ran_seq[ran_i,j+4]==obs_seq[i,j+4]){}else{
                    d <- d+1}}}
              diff <- c(diff,d)}
            ave_difference <- c(ave_difference,mean(diff))
          }
      }
      ave_yearly_mismatch <- c(ave_yearly_mismatch,mean(ave_difference, na.rm=T))
    }
    mat_yearly_mismatch <- rbind(mat_yearly_mismatch,ave_yearly_mismatch)
  }
  mat_yearly_mismatch <- mat_yearly_mismatch[-1,]
  rownames(mat_yearly_mismatch) <- predict_season
  colnames(mat_yearly_mismatch) <- paste0("repeat",c(1:10))
  return(mat_yearly_mismatch)
}

ha_ran_mismatch <- RandomPickStrain(input_csv_seq = ha_allseq, predict_season = c(2012:2020), 
                                    repeat_time = 10, emsite = c(1:566), 
                                    region = c("NY","CA","CAN","UK","GER","FR","HK","SCP","JP","SG"), subtype = c("H1_HA_sequence.csv"))
na_ran_mismatch <- RandomPickStrain(input_csv_seq = na_allseq, predict_season = c(2012:2020), 
                                    repeat_time = 10, emsite = c(1:469), 
                                    region = c("NY","CA","CAN","UK","GER","FR","HK","SCP","JP","SG"), subtype = c("H1_NA_sequence.csv"))

## output: average mismatch by year

write.csv(ha_ran_mismatch,file="h1ha_random_mismatch.csv")
write.csv(na_ran_mismatch,file="h1na_random_mismatch.csv")


## for h3n2, random seed=76
load("./training set/NH/H3N2/NHH3_HA_sequence.Rdata")
load("./training set/NH/H3N2/NHH3_NA_sequence.Rdata")

RandomPickStrain <- function(input_csv_seq, predict_season, repeat_time, emsite, region, subtype){
  mat_yearly_mismatch <- matrix(ncol = repeat_time)
  for (pre_season in predict_season){
    print(pre_season)
    seq <- input_csv_seq[which(input_csv_seq$season == (pre_season-1) & input_csv_seq$date < make_date(pre_season-1,3,1)),]
    set.seed(76) #set a seed for random pick
    seed <- sample(1:length(seq[,1]),repeat_time,replace=F)
    ran_seq <- seq[seed,]
    ave_yearly_mismatch <- c()
    for(ran_i in 1:length(ran_seq[,1])){
      ave_difference <- c()
      for (re in region){
        sequence <- read.csv(paste0("./validation set/H3N2/",re,subtype))
        obs_seq <- subset(sequence,season==pre_season)
        obs_seq <- as.matrix(obs_seq)
        diff <- 0
        if (length(obs_seq[,1])==0){
          ave_difference <- c(ave_difference,NA)}else{
            for (i in 1:length(obs_seq[,1])){
              d <- 0
              for (j in emsite){
                if (is.na(obs_seq[i,j+4])){}else{
                  if (ran_seq[ran_i,j+4]==obs_seq[i,j+4]){}else{
                    d <- d+1}}}
              diff <- c(diff,d)}
            ave_difference <- c(ave_difference,mean(diff))
          }
      }
      ave_yearly_mismatch <- c(ave_yearly_mismatch,mean(ave_difference, na.rm=T))
    }
    mat_yearly_mismatch <- rbind(mat_yearly_mismatch,ave_yearly_mismatch)
  }
  mat_yearly_mismatch <- mat_yearly_mismatch[-1,]
  rownames(mat_yearly_mismatch) <- predict_season
  colnames(mat_yearly_mismatch) <- paste0("repeat",c(1:10))
  return(mat_yearly_mismatch)
}

ha_ran_mismatch <- RandomPickStrain(input_csv_seq = ha_allseq, predict_season = c(2003:2020), 
                                    repeat_time = 10, emsite = c(1:566), 
                                    region = c("NY","CA","CAN","UK","GER","FR","HK","SCP","JP","SG"), subtype = c("H3_HA_sequence.csv"))
na_ran_mismatch <- RandomPickStrain(input_csv_seq = na_allseq, predict_season = c(2003:2020), 
                                    repeat_time = 10, emsite = c(1:469), 
                                    region = c("NY","CA","CAN","UK","GER","FR","HK","SCP","JP","SG"), subtype = c("H3_NA_sequence.csv"))

## output: average mismatch by year

write.csv(ha_ran_mismatch,file="h3ha_random_mismatch.csv")
write.csv(na_ran_mismatch,file="h3na_random_mismatch.csv")


