# Site-basedPredictionModel
R codes for the article "Predictive evolutionary modelling for influenza virus by site-based dynamics of mutations".

Input datasets:
- genetic sequences of A/H3N2 virus from Oct 2008 to Sep 2020 in Hong Kong SAR (file names: H3N2_HA_HK_sequence.fas, H3N2_NA_HK_sequence.fas)
- epidemic data (H3N2_HK_EpidemicData.csv)
- As reference group, the set of WHO recommended vaccine strains (H3N2_WHOvaccine_HA.csv, H3N2_WHOvaccine_NA.csv) 

Step 0.  Install dependent packages 
R (4.1.3) or newer version is required. To download R, please see the website https://www.r-project.org/ for detailed instructions. The R packages Biostring (2.64.1), lubridate (1.8.0), plyr (1.8.7) and stringr (1.4.0) are pre-requisites, which can be installed from the following websites.

Biostring (2.64.1): https://bioconductor.org/packages/release/bioc/html/Biostrings.html

lubridate (1.8.0): https://www.rdocumentation.org/packages/lubridate/versions/1.8.0

plyr (1.8.7): https://www.rdocumentation.org/packages/plyr/versions/1.8.7

stringr (1.4.0): https://www.rdocumentation.org/packages/stringr/versions/1.4.0

Step 1.  Preparing data:
In this step, the R function ‘Prepareseq’ reads viral sequence data in fasta format and converts the data to csv for further analysis. Two data frames will be generated, one for the HA and one for the NA sequences with sequence name, accession number, collection date and flu season. (Sequence nomenclature: flu season 2011/12 is labeled as 2012). 

Step 2.  Calculate site wise amino acid prevalence:
Site-wise amino acid prevalence through time is calculated for both HA and NA, output in two separate data frames. Run time: 3-5 minutes. 

Step 3.  Calculate g-measure and transition time under different θ and h:
The g-measure and transition time τ under different threshold (θ) and extended effective mutation period (h) are calculated (see Methods section of the paper for detailed description). Here, we test threshold θ ranged from 0.5 to 1.0, and h ranged from 0 to 9 years, which leads to in total 60 g-measure vectors and corresponding average transition time. This will output two g-measure data frames for HA and NA protein, respectively. Run time: 2-3 minutes. 

Step 4.  Estimate transition time τ:
Transition time describes the duration that a mutation takes from its first emergence to reaching an influential prevalence in population. The τ may diverse across different subtypes or segments, and is also a parameter in Step 5 to predict future mutation prevalence. The optimal (θ, h) are estimated by optimizing a fitness function of the g-measure against epidemic levels [1]. The function 'FitRegression' will output optimal (θ, h), the corresponding R-squared and average transition time.

Step 5.  Obtain history EM set under the optimal (θ, h)
The function 'EffMutation' will output history EM set with their site-specific and time-specific transition time, and it will be treated as history transition event in the following prediction procedure.

Step 6.  Beth-1 prediciton
The beth-1 consists of future mutation prevalence prediction and vaccine strain selection.

Step 6.1  Prediction of future mutation prevalence and obtain consensus strain:
The prediction requires 3 year-data to burn in. Here we output year-to-year prediction from 2011/12 to 2019/20 season. The function ‘PreYearlyStrain’ allows to input genetic sequence, predicted season, average transition time and history EM set and output predicted consensus sequence and EM set by year. The predicted EM set will be used as part of predictor codon sites to select wildtype strains in the following procedure.

Step 6.2  Select the optimal wild-type strain:
The VE-GD model [2-4] for influenza virus was best fitted when the genetic distance was evaluated on a subset of EMs residing on epitope A or B for H3N2 and the EMs on antigenic sites for pH1N1, referred to as the predictor codon sites. Here, the relative importance of HA and NA protein to VE are quantified by their effect size of genetic mismatch contributing to VE respectively [2-4]. Two wildtype sequence files for HA and NA protein will be output.

Step 7.  Evaluate vaccine strain performance by genetic mismatch:
This step aims to compute yearly epitope mismatch of predicted vaccine strains against circulating virus. We also provide sequences of WHO recommended vaccine strains for comparison. This will output yearly average epitope mismatch and standard deviation.

The files " H1N1_gisaid_acknowledge_table.csv" and " H3N2_gisaid_acknowledge_table.csv" provide accession numbers of influenza sequences used in this study. Viral sequence data were downloaded from the global initiative on sharing all influenza data (GISAID) at http://platform.gisaid.org/. We thank the contributions of all the health care workers and scientists, the GISAID team, and the submitting and the originating laboratories.

References:

[1] Wang MH, Lou J, Cao L, et al. Characterization of key amino acid substitutions and dynamics of the influenza virus H3N2 hemagglutinin. J Infect. 2021;83(6):671-677. doi:10.1016/j.jinf.2021.09.026
[2] Cao L, Lou J, Zhao S, et al. In silico prediction of influenza vaccine effectiveness by sequence analysis. Vaccine. 2021;39(7):1030-1034. doi:10.1016/j.vaccine.2021.01.006
[3] Cao L, Zhao S, Lou J, et al. Differential Influence of Age on the Relationship between Genetic Mismatch and A(H1N1)pdm09 Vaccine Effectiveness. Viruses. 2021;13(4):619. Published 2021 Apr 4. doi:10.3390/v13040619
[4] Cao, L. et al. Improving the prediction of influenza vaccine effectiveness by refined genetic distance measure. medRxiv, 2023.2002. 2014.23285900 (2023).
