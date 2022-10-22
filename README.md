# Site-basedPredictionModel
R codes for the article "Predictive evolutionary modelling for influenza virus by site-based dynamics of mutations".

R (4.1.3) or newer version is required. To download R, please see the website https://www.r-project.org/ for detailed instructions. The R packages Biostring (2.64.1), lubridate (1.8.0) and plyr (1.8.7) are required in the following analysis, and you may install them according to the instructions on the following websites.

Biostring (2.64.1): https://bioconductor.org/packages/release/bioc/html/Biostrings.html 

lubridate (1.8.0): https://www.rdocumentation.org/packages/lubridate/versions/1.8.0 

plyr (1.8.7): https://www.rdocumentation.org/packages/plyr/versions/1.8.7  

Here, we provide a demo dataset of A/H3N2 in Hong Kong SAR, China to illustrate the prediction procedures. The dataset includes HA and NA sequences of 836 influenza virus isolates collected during Oct 2008 to Sep 2020 (H3N2_HA_HK_sequence.fasta, H3N2_NA_HK_sequence.fasta), influenza epidemic data (H3N2_HK_EpidemicData.csv) and a set of dynamic predictor codon sites (H3N2_HK_EMset.csv). For comparison of the result of genetic mismatch evaluation, we also provide WHO recommended vaccine strain sequences (H3N2_WHOvaccine_HA.csv, H3N2_WHOvaccine_NA.csv).

Step 1 Preparing data: viral sequence data can be retrieved from public database. After completing sequence alignment, the code can be conducted to sort data and convert fasta format to csv format for further analysis. You may use the demo data to run the code for both HA and NA sequences, and two csv format sequence files will be obtained. Since the flu season is usually from October to April the next year, here 2012 season represent 2011/12 flu season and so on.

Step 2 Calculate site wise amino acid prevalence: we calculate site wise amino acid prevalence through time for both HA and NA sequence. It may require 3-5 minutes to run the code. This will output two prevalence files for HA and NA protein, respectively. 

Step 3 Calculate g-measure and transition time under different θ and h: the g-measure and transition time τ under different threshold θ and extended effective mutation period (h) are calculated (see Methods section of the paper for detailed description). Here, we test threshold θ ranged from 0.5 to 1.0, and h ranged from 0 to 9 years, which leads to in total 60 g-measure vectors and corresponding average transition time. It may require 2-3 minutes to run the code. This will output two g-measure files for HA and NA protein, respectively.

Step 4 Estimate transition time τ: in this step, we introduce method to estimate average transition time τ, which describes the duration that a mutation takes from its first emergence to reaching an influential prevalence in population. The τ may diverse across different subtypes or segments, and is also a parameter in Step 5 to predict future mutation prevalence. The optimal (θ, h) are estimated by optimizing a fitness function of the g-measure against an epidemic variable [1]. Here, we used a generalized linear model with H3 positivity rate as response, g-measure as predictor, and the covariates are yearly mean temperature and absolute humidity in Hong Kong. Alternative epidemic variable and fitting function can also be employed. 

Step 5 Prediction of future mutation prevalence and consensus strain: the prediction requires 4 year-data to burn in, so here we predict future prevalence from 2011/12 to 2020/21 season. In this step, two files containing predicted consensus sequences will be obtained. Besides, if you want to know more details about predicted yearly site wise prevalence, please output matrix Pred_Prevalence.

Step 6 Select the optimal wild-type strain: we use the prediction model for vaccine effectiveness by genetic distance (VE-GD) [2, 3] to select the available present wildtype virus that maximizes the estimated VE. A dynamic predictor codon sites set is provided based on the demo dataset, which is a feature-selected set from the EM sites by maximizing the fitness between observed VE data and the circulating viruses in corresponding epidemic seasons. Here, the relative importance of HA and NA protein to VE are quantifies by their effect size of genetic mismatch contributing to VE respectively [2, 3]. In this step, two wildtype sequence files for HA and NA protein will be output.

Step 7 Evaluate vaccine strain by calculating genetic mismatch: this step aims to compute yearly epitope mismatch of predicted vaccine strains against circulating virus. We also provide sequences of WHO recommended vaccine strains for comparison. This will output yearly average epitope mismatch and standard deviation. 

The files " H1N1_gisaid_acknowledge_table.csv" and " H3N2_gisaid_acknowledge_table.csv" shows accession numbers of influenza sequences used in this study. Viral sequence data were downloaded from the global initiative on sharing all influenza data (GISAID) at http://platform.gisaid.org/. We thank the contributions of all the health care workers and scientists, the GISAID team, and the submitting and the originating laboratories.

Reference:
1.	Wang MH, Lou J, Cao L, et al. Characterization of key amino acid substitutions and dynamics of the influenza virus H3N2 hemagglutinin. J Infect. 2021;83(6):671-677. doi:10.1016/j.jinf.2021.09.026
2.	Cao L, Lou J, Zhao S, et al. In silico prediction of influenza vaccine effectiveness by sequence analysis. Vaccine. 2021;39(7):1030-1034. doi:10.1016/j.vaccine.2021.01.006
3.	Cao L, Zhao S, Lou J, et al. Differential Influence of Age on the Relationship between Genetic Mismatch and A(H1N1)pdm09 Vaccine Effectiveness. Viruses. 2021;13(4):619. Published 2021 Apr 4. doi:10.3390/v13040619


