# Site-basedPredictionModel
R codes for the article "Predictive evolutionary modelling for influenza virus by site-based dynamics of mutations".

Input datasets:
-  genetic sequences of A/H3N2 virus from Oct 2008 to Sep 2020 in Hong Kong SAR (file names: H3N2_HA_HK_sequence.fasta, H3N2_NA_HK_sequence.fasta)
- epidemic data (H3N2_HK_EpidemicData.csv)
- the set of dynamic predictor codon sites for HK (H3N2_HK_EMset.csv)
- As reference group, the set of WHO recommended vaccine strains (H3N2_WHOvaccine_HA.csv, H3N2_WHOvaccine_NA.csv) 

Step 0.  Install dependent packages 
R (4.1.3) or newer version is required. To download R, please see the website https://www.r-project.org/ for detailed instructions. The R packages Biostring (2.64.1), lubridate (1.8.0) and plyr (1.8.7) are pre-requisites, which can be installed from the following websites.

Biostring (2.64.1): https://bioconductor.org/packages/release/bioc/html/Biostrings.html

lubridate (1.8.0): https://www.rdocumentation.org/packages/lubridate/versions/1.8.0

plyr (1.8.7): https://www.rdocumentation.org/packages/plyr/versions/1.8.7

Step 1. Preparing data
In this step, the program reads in viral sequence data in fasta format and converts the data to csv for further analysis. Two csv files will be converted, one for the HA and one for the NA sequences. (sequence nomenclature: epidemic season 2011/12 is labeled as 2012). 

Step 2. Calculate site wise amino acid prevalence
Site-wise amino acid prevalence through time is calculated for both HA and NA, output in two separate files. Run time: 3-5 minutes. 

Step 3. Calculate g-measure and transition time under different θ and h
The g-measure and transition time τ under different threshold θ and extended effective mutation period (h) are calculated (see Methods section of the paper for detailed description). Here, we test threshold θ ranged from 0.5 to 1.0, and h ranged from 0 to 9 years, which leads to in total 60 g-measure vectors and corresponding average transition time. This will output two g-measure files for HA and NA protein, respectively. Run time: 2-3 minutes. 

Step 4. Estimate transition time τ
Transition time describes the duration that a mutation takes from its first emergence to reaching an influential prevalence in population. The τ may diverse across different subtypes or segments, and is also a parameter in Step 5 to predict future mutation prevalence. The optimal (θ, h) are estimated by optimizing a fitness function of the g-measure against epidemic levels [1]. 

Step 5. Prediction of future mutation prevalence and obtain consensus strain
The prediction requires 4 year-data to burn in. Here we output year-to-year prediction from 2011/12 to 2020/21 season. In this step, two files containing predicted consensus sequences will be obtained. 

Step 6. Select the optimal wild-type strain 
The VE-GD model [2, 3] is used to select the available present wildtype virus that maximizes the expected VE. A dynamic predictor codon sites set is provided, which is a feature-selected set from the EM sites by maximizing the fitness between observed VE data and the circulating viruses in corresponding epidemic seasons. Here, the relative importance of HA and NA protein to VE are quantified by their effect size of genetic mismatch contributing to VE respectively [2, 3]. Two wildtype sequence files for HA and NA protein are output.

Step 7. Evaluate vaccine strain performance by genetic mismatch
This step aims to compute yearly epitope mismatch of predicted vaccine strains against circulating virus. We also provide sequences of WHO recommended vaccine strains for comparison. This will output yearly average epitope mismatch and standard deviation.

The files " H1N1_gisaid_acknowledge_table.csv" and " H3N2_gisaid_acknowledge_table.csv" provide accession numbers of influenza sequences used in this study. Viral sequence data were downloaded from the global initiative on sharing all influenza data (GISAID) at http://platform.gisaid.org/. We thank the contributions of all the health care workers and scientists, the GISAID team, and the submitting and the originating laboratories.

References:

[1] Wang MH, Lou J, Cao L, et al. Characterization of key amino acid substitutions and dynamics of the influenza virus H3N2 hemagglutinin. J Infect. 2021;83(6):671-677. doi:10.1016/j.jinf.2021.09.026
[2] Cao L, Lou J, Zhao S, et al. In silico prediction of influenza vaccine effectiveness by sequence analysis. Vaccine. 2021;39(7):1030-1034. doi:10.1016/j.vaccine.2021.01.006
[3] Cao L, Zhao S, Lou J, et al. Differential Influence of Age on the Relationship between Genetic Mismatch and A(H1N1)pdm09 Vaccine Effectiveness. Viruses. 2021;13(4):619. Published 2021 Apr 4. doi:10.3390/v13040619



