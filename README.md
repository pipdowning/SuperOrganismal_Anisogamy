# SuperOrganismal_Anisogamy
Data, code, and supplementary information for SUPERORGANISMAL ANISOGAMY: A COMPARATIVE TEST OF AN EXTENDED THEORY

Number of supplementary items: four
1. Table_S1.csv
2. Table_S2.csv
3. SOA_R_code.R
4. SOA_Matlab_code.m


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: Table_S1.csv

This csv document contains:
	+ raw data for the phylogenetic comparative analysis
	+ column descriptions:
		A. subfamily = the ant subfamily to which each species belongs
		B. tip_label = latin binomial of each species, matches the phylogeny
		C. male_MW = maximum mesosoma width for male of each species
		D. male_ML = maximum mesosoma length for male of each species
		E. queen_MW = maximum mesosoma width for queen of each species
		F. queen_ML = maximum mesosoma width for queen of each species
		G. colony_size = number of workers in mature colonies
		H. queen_number_continuous = average number of queens observed per colony
		I. caste_number = binary number of worker castes (one or more than one caste)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: Table_S2.csv

This csv document contains:
	+ results from the statistical analyses using an inverse Gamma prior
	+ column descriptions:
		A. Fixed terms = the predictor variable(s) used in each model
		B. parameter = name of the parameter estimate
		C. estimate = value of the parameter estimate
		D. Lwr CI = lower 95% credible interval of the parameter estimate
		E. Upr CI = upper 95% credible interval of the parameter estimate
		F. Phylogenetic H2 = the proportion of variance explained by phylogeny


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: SOA_R_code.R

This R script contains all the code needed to replicate the comparative analyses.
Note that these analyses can take multiple days to run.

- Packages and Data (lines 14 to 108)
- Part 1: queen-male coevolution (lines 128 to 201)
- Part 2: queen size in relation to colony size, caste number, and queen number (lines 207 to 441)
- Part 3: dimorphism in relation to colony size, caste number, and queen number (lines 447 to 675)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: SOA_Matlab_code.m

This matlab script contains the code needed to replicate the theoretical analyses.

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
