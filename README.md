# Testing the existence of an unadmixed ancestor from a specific population `t` generations ago

### Abstract

  Next generation sequencing allows to determine the complete genomes of individuals. With this information together with reference panels of ancestral source populations, the ancestry of each position of the genome can be estimated (local ancestry). The information on those ancestry-specific genomic segments are commonly used to understand migration waves and admixture events. In short time scales, it is often of interest to determine the existence of the most recent unadmixed ancestor from a specific population `t` generations ago.
    We build a hypothesis test to determine if an individual has an ancestor belonging to a target ancestral population `t` generations ago based on these ancestry-specific segments at an individual level. We apply this test on a data set that includes Uruguayan admixed individuals to estimate for each one how many generations ago the most recent indigenous ancestor lived. 

### Author summary

Little is known about the Amerindian groups that populated Latin America in pre-colonial times, especially in South America. In Uruguay, the Amerindian heritage has been neglected and underestimated for decades, and it is not until the last 20 years that Charrúan (indigenous group that inhabited Uruguay) heritage has become an important part of the Uruguayan identity. The Charrúas were extinguished as a result of a genocide conducted by the government of that time in 1831. The information regarding this ethnic groups was mostly lost. A recent study of our group has used genomic data of the descendant of the Charrúas to shed some light on this chapter of Uruguayan history. 

Motivated by this new information (genomic data of the descendants) and the will of contributing with more information to the history of the country, we develop a statistical test to determine if it is likely that an individual has an ancestor belonging to a target ancestral population `t` generations ago. By applying this test to the Uruguayan genomic data we were able to assess, for each individual, how many generations ago lived the most recent indigenous ancestor. 

### Usage of the code

1) rcode.R takes the information from http://urugenomes.org/en/databases/ and returns files for the desired chromosome statistics and ancestral population. In our case, it returns est_max.csv (maximum length -in morgans- of native-american tract for each non-sexual chromosome and individual) and est_sum.csv (sum of lengths -in morgans- of native-american tracts for each non-sexual chromosome and individual). Just in case the Urugenomes data base is not available, we upload both output files, which will be used by jcode.jl.
2) largos_chr.csv contains the length -in centimorgans- of every non-sexual chromosome; this will be used by jcode.jl.
3) jcode.jl runs all simulations, tests and plots presented in the article. Some comments, examples and basic usage are shown. This code is easy to edit to use for any data base (once the statistcs are extracted from the raw data).
