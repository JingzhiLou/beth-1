
###########################################################
# The function RandomPickStrain allows to randomly select a sequence in time T as a prediction strain and calculate average genetic mismatch against circulating virus (repeat 10 times)

RandomPickStrain <- function(input_csv_seq, predict_season, repeat_time, emsite, val_set){
  mat_yearly_mismatch <- matrix(ncol = repeat_time)
  for (pre_season in predict_season){
    print(pre_season)
    seq <- input_csv_seq[which(input_csv_seq$season == (pre_season-1)),]
    seq <- seq[which(month(seq$date) < 3 | month(seq$date) > 9),]
    set.seed(60) #set a seed for random pick
    seed <- sample(1:length(seq[,1]),repeat_time,replace=T)
    ran_seq <- seq[seed,]
    ave_yearly_mismatch <- c()
    for(ran_i in 1:length(ran_seq[,1])){
      ave_difference <- c()
      sequence <- val_set
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
      ave_yearly_mismatch <- c(ave_yearly_mismatch,mean(ave_difference, na.rm=T))
    }
    mat_yearly_mismatch <- rbind(mat_yearly_mismatch,ave_yearly_mismatch)
  }
  mat_yearly_mismatch <- mat_yearly_mismatch[-1,]
  mat_yearly_mismatch <- cbind(mat_yearly_mismatch,rowMeans(mat_yearly_mismatch))
  rownames(mat_yearly_mismatch) <- predict_season
  colnames(mat_yearly_mismatch) <- c(paste0("repeat",c(1:10)),"mean")
  return(mat_yearly_mismatch)
}

# Random seed = 60
ha_ran_mismatch <- RandomPickStrain(input_csv_seq = read.csv("H3N2_HA_HK_sequence.csv"), predict_season = c(2012:2020), 
                                    repeat_time = 10, emsite = c(1:566), 
                                    val_set = read.csv("H3N2_HA_HK_sequence.csv"))
na_ran_mismatch <- RandomPickStrain(input_csv_seq = read.csv("H3N2_NA_HK_sequence.csv"), predict_season = c(2012:2020), 
                                    repeat_time = 10, emsite = c(1:469), 
                                    val_set = read.csv("H3N2_NA_HK_sequence.csv"))

# Output: genetic mismatch table by year

#write.csv(ha_ran_mismatch,file="h1ha_random_mismatch.csv")
#write.csv(na_ran_mismatch,file="h1na_random_mismatch.csv")
