
# prepare protein sequence
# input an aligned protein sequence data and output a csv format data with flu season

library(Biostrings)
library(lubridate)

inputdata_name <- c("H3N2_HA_HK_sequence.fasta") # "H3N2_NA_HK_sequence.fasta" for NA sequence
outputdata_name <- c("H3N2_HA_HK_sequence.csv") # "H3N2_NA_HK_sequence.csv" for NA sequence

# please specify the input and output data name and run the following code for both HA and NA sequence

fasta_seq <- readAAStringSet(inputdata_name)
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
write.csv(complete_data,file=outputdata_name,row.names = F)
