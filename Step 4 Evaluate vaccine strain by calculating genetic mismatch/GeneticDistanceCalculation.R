
# calculate epitope distance between predicted vaccine strains and circulating strains
# input vaccine strain and circulating strains, and output genetic mismatch

h3ha_epitope <- c(122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198,44,45,46,47,48,50,51,53,54,273,275,276,278,279,280,294,297,299,300,304,305,307,308,309,310,311,312,96,102,103,117,121,167,170,171,172,173,174,175,176,177,179,182,201,203,205,207,208,209,212,213,214,215,216,217,218,219,226,227,228,229,230,240,242,244,246,247,248,57,59,62,63,67,75,78,80,81,82,83,86,87,88,91,92,94,109,260,261,262,265)+16
h3na_epitope <- c(150,198,199,220,221,253,329,334,344,368,370,403)

cir_seq <- read.csv("H3N2_HA_HK_sequence.csv")
pre_seq <- read.csv("CombinedH3N2_HA_Wildtype.csv")
pre_seq <- as.matrix(pre_seq)
pre_year <- 9
start_season <- 2012
emsite <- h3ha_epitope # change for h3na_epitope when compute mismatch on NA protein

ave_mismatch <- c()
sd <- c()
for (z in 1:pre_year){
    obs_seq <- subset(cir_seq,season==z+start_season-1)
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
Mismatch <- cbind(c(start_season:(start_season+pre_year-1)),ave_mismatch,sd)
print(Mismatch)
