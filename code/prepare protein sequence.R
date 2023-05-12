
setwd("C:/Users/jingzhi/Desktop/Site-basedPredictionModel/data/training set/SEA/H1N1")

###################### Prepare protein sequence data #######################
# The function 'Prepareseq' allows to input an aligned protein sequence data with fasta format and output a csv format data with flu season

library(Biostrings)
library(lubridate)

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

ha_seq <- PrepareSeq(input_fasta_seq = "SEAH1_HA_sequence.fas")
na_seq <- PrepareSeq(input_fasta_seq = "SEAH1_NA_sequence.fas")

# output: HA sequence and NA sequence data
write.csv(ha_seq,file="SEAH1_HA_sequence.csv",row.names = F)
write.csv(na_seq,file="SEAH1_NA_sequence.csv",row.names = F)