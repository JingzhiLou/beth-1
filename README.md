# Site-basedPredictionModel
R codes for the article "Predictive evolutionary modelling for influenza virus by site-based dynamics of mutations".

R (4.1.3) or newer version is required. To download R, please see the website https://www.r-project.org/ for detailed instructions. The R packages Biostring (2.64.1), lubridate (1.8.0) and plyr (1.8.7) are required in the following analysis, and you may install them according to the instructions on the following websites.

Biostring (2.64.1): https://bioconductor.org/packages/release/bioc/html/Biostrings.html 

lubridate (1.8.0): https://www.rdocumentation.org/packages/lubridate/versions/1.8.0 

plyr (1.8.7): https://www.rdocumentation.org/packages/plyr/versions/1.8.7 

This repository includes the following four folders: Step 1 Preparing data; Step 2 Estimate transition time τ; Step 3 Prediction of future mutation prevalence and select the optimal wild-type strain; Step 4 Evaluate vaccine strain by calculating genetic mismatch. Data analysis involved in this work can be performed with these codes. 

Here, we provide a demo dataset of A/H3N2 in Hong Kong SAR, China to illustrate the prediction procedures. The dataset includes HA and NA sequences of 836 influenza virus isolates collected during Oct 2008 to Sep 2020 (H3N2_HA_HK_sequence.fasta, H3N2_NA_HK_sequence.fasta), influenza epidemic data (H3N2_HK_EpidemicData.csv)and a set of dynamic predictor codon sites (H3N2_HK_EMset.csv).
