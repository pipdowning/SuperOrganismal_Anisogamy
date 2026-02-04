# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #
# 																										 #
# R Code for:																							 #
# 																										 #
#				SUPERORGANISMAL ANISOGAMY: A COMPARATIVE TEST OF AN EXTENDED THEORY	 				 	 #
# 																										 #
# 																										 #
# Code by 																								 #
# 																										 #
# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

## PACKAGES AND DATA ##

library(ape)
library(phytools)
library(MCMCglmm)


## load trait data and trees ##
ant_data <- read.csv(".../Table_S1.csv")
head(ant_data)
str(ant_data)	# 732 rows, 9 variables
length(unique(ant_data$tip_label))	# 732 species
ant_trees <- read.tree(".../15K_FBD_crown_posterior.tre")
class(ant_trees)						# multiphylo
length(ant_trees)						# 100
length(ant_trees[[1]]$tip.label)		# 14594 species

## sample 50 trees at random for use in analyses ##
sample(1:100, 50)
tree_sample <- c(46,76,7,14,30,6,12,18,90,38,5,1,74,78,97,24,53,47,93,73,35,94,69,23,52,4,43,91,98,61,71,42,55,63,65,11,95,16,57,33,88,67,60,40,72,89,15,41,87,58)		# the trees we used
ant_trees <- ant_trees[tree_sample]
class(ant_trees)		# multiphylo
length(ant_trees)		# 50


## trim the trees ##
drop_tip <- ant_trees[[1]]$tip.label[which(is.na(match(ant_trees[[1]]$tip.label, ant_data$tip_label)))]
length(drop_tip)		# 13862 species to drop
ant_trees <- lapply(ant_trees, drop.tip, drop_tip, trim.internal=T)
ant_trees <- lapply(ant_trees, makeNodeLabel, method = "number")
# check match between species in data and species in tree
ant_data$tip_label[which((ant_data$tip_label %in% ant_trees[[1]]$tip.label) == FALSE)]
ant_trees[[1]]$tip.label[which((ant_trees[[1]]$tip.label %in% ant_data$tip_label) == FALSE)]
length(ant_trees[[1]]$tip.label)		# 732 species

## variable summaries and transformations ##
# note - no measurements for queens and males are missing
# colony size
length(na.omit(ant_data$colony_size))					# 154 species
range(ant_data$colony_size, na.rm=TRUE)					# 7 to 5 500 000
hist(ant_data$colony_size)								# strongly right skewed
ant_data$lg_colony_size <- log(ant_data$colony_size)		# log transform
hist(ant_data$lg_colony_size)							# normal
# queen number
length(na.omit(ant_data$queen_number_continuous))					# 51 species
range(ant_data$queen_number_continuous, na.rm=TRUE)					# 1 to 70.5
hist(ant_data$queen_number_continuous)								# strongly right skewed
ant_data$lg_queen_number <- log(ant_data$queen_number_continuous)	# log transform
hist(ant_data$lg_queen_number)										# still right skewed
# caste number
length(na.omit(ant_data$caste_number))		# 187 species
range(ant_data$caste_number, na.rm=TRUE)		# 1 to 4
table(ant_data$caste_number)				# one = 153, two = 28, three = 4, four = 2

## create queen and male thorax volume estimates ##
# queen
ant_data$Qarea <- ant_data$queen_MW * ant_data$queen_ML
ant_data$Qvol <- ant_data$Qarea^(3/2)
hist(ant_data$Qvol)							# strongly right skewed
ant_data$lg_Qvol <- log(ant_data$Qvol)
hist(ant_data$lg_Qvol)						# normal
plot(lg_Qvol ~ Qarea, ant_data)
# male
ant_data$Marea <- ant_data$male_MW * ant_data$male_ML
ant_data$Mvol <- ant_data$Marea^(3/2)
hist(ant_data$Mvol)							# strongly right skewed
ant_data$lg_Mvol <- log(ant_data$Mvol)
hist(ant_data$lg_Mvol)						# normal
plot(lg_Mvol ~ Marea, ant_data)

## queen - male thorax volume difference ##
ant_data$delta_QM_vol <- ant_data$lg_Qvol - ant_data$lg_Mvol	# 0 means =, > 0 means queen bigger
hist(ant_data$delta_QM_vol)										# normalish
length(na.omit(ant_data$delta_QM_vol))							# 732 estimates

## make caste number binary since few species have > 2 ##
ant_data$bin_caste <- as.factor(ifelse(ant_data$caste_number == 1, "one", "two"))
table(ant_data$bin_caste)		# one = 153, two = 34

## priors ##
# inverse-Gamma distribution
priorA <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
priorB <- list(R = list(V=diag(2), nu=2), G = list(G1 = list(V=diag(2), nu=2)))
# parameter expanded priors
priorC <- list(R = list(V=1, nu=1), G = list(G1 = list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))
priorD <- list(R = list(V=diag(2), nu=2), G = list(G1 = list(V=diag(2), nu=2, alpha.mu=rep(0,2), alpha.V = diag(2)*1000)))

# nitt=1100000, thin=1000, burnin=100000 from 'Model_Timings.R'
# for 6/9 parameters gives effective sample of 1000 (for 1 tree)
# for 5/11 paramters < 1000, min = 743 (for 1 tree)

# use these to create mcmc objects combining 20 values from the posterior mode per tree
# gives 1000 in total since 50 trees (i.e. 20 * 50)
store <- seq(50,1000, by = 50)		# chooses every 50th observation from each chain = 20 in total per tree
ranges <- data.frame(start = seq(1,981, by = 20), stop = seq(20,1000, by = 20))		# start and stop values for each tree


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

## sample sizes ##

# main effects
dim(na.omit(ant_data[,c("lg_Qvol", "lg_colony_size")]))							# 154 species
dim(na.omit(ant_data[,c("lg_Qvol", "bin_caste")]))								# 187 species
dim(na.omit(ant_data[,c("lg_Qvol", "lg_queen_number")]))						# 51 species
# interactions
dim(na.omit(ant_data[,c("lg_Qvol", "lg_colony_size", "lg_queen_number")]))		# 46 species
dim(na.omit(ant_data[,c("lg_Qvol", "bin_caste", "lg_queen_number")]))			# 46 species


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 1: QUEEN-MALE CO-EVOLUTION ###

## Vital Statistics ##

table(ant_data$delta_QM_vol == 0)		# 0 males = queens exactly
table(ant_data$delta_QM_vol < 0)		# 524 queens larger than males
table(ant_data$delta_QM_vol > 0)		# 208 males larger than queens


## Multi-response model to estimate phylogenetic and phenotypic correlations ##

plot(lg_Qvol ~ lg_Mvol, ant_data)
abline(lm(lg_Qvol ~ lg_Mvol, ant_data))
abline(0,1)

phylo_r_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  phylo_r_mods[[i]] <- MCMCglmm(cbind(lg_Qvol, lg_Mvol) ~ trait-1, random = ~us(trait):tip_label, rcov = ~ us(trait):units, family=c("gaussian", "gaussian"), data=ant_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorD, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(phylo_r_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
phylo_r_VCV <- mcmc(data.frame(QQphylo=1:1000, QMphylo=1:1000, MQphylo=1:1000, MMphylo=1:1000, QQunits=1:1000, QMunits=1:1000, MQunits=1:1000, MMunits=1:1000))
phylo_r_Sol <- mcmc(data.frame(Qint = 1:1000, Mint = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){phylo_r_VCV[ranges$start[i]:ranges$stop[i],] <- phylo_r_mods[[i]]$VCV[store,]}
for(i in 1:50){phylo_r_Sol[ranges$start[i]:ranges$stop[i],] <- phylo_r_mods[[i]]$Sol[store,]}

# chain convergence
plot(phylo_r_VCV)    		# all parameters well mixed
plot(phylo_r_Sol)   		# all parameters well mixed
autocorr(phylo_r_VCV)   		# correlation between successive samples < 0.1
autocorr(phylo_r_Sol)		# correlation between successive samples < 0.1
heidel.diag(phylo_r_VCV)		# all paramaters passed halfwidth
heidel.diag(phylo_r_Sol)		# intercepts passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(phylo_r_Sol),2)		# Qint = 0.77, Mint = 0.35
round(HPDinterval(phylo_r_Sol),2)		# Qint = -0.35 to 1.68, Mint = -0.52 to 1.28
round(posterior.mode(exp(phylo_r_Sol)),2)	# Qint = 1.7, Mint = 1.4
round(HPDinterval(exp(phylo_r_Sol)),2)		# Qint = 0.39 to 4.81, Mint = 0.49 to 3.27

## heritabilities ##
# queen phylo heritability
H2Qvol <- phylo_r_VCV[,"QQphylo"] / (phylo_r_VCV[,"QQphylo"] + phylo_r_VCV[,"QQunits"])
posterior.mode(H2Qvol)		# 0.82
HPDinterval(H2Qvol)			# 0.74 to 0.89
# male phylo heritability
H2Mvol <- phylo_r_VCV[,"MMphylo"] / (phylo_r_VCV[,"MMphylo"] + phylo_r_VCV[,"MMunits"])
posterior.mode(H2Mvol)		# 0.83
HPDinterval(H2Mvol)			# 0.74 to 0.89

## correlations ##
# queen:male phylo r
prQM <- phylo_r_VCV[,"QMphylo"] / sqrt(phylo_r_VCV[,"QQphylo"] * phylo_r_VCV[,"MMphylo"])
posterior.mode(prQM)		# 0.92
HPDinterval(prQM)			# 0.87 to 0.94
# queen:male pheno r
rQM <- phylo_r_VCV[,"QMunits"] / sqrt(phylo_r_VCV[,"QQunits"] * phylo_r_VCV[,"MMunits"])
posterior.mode(rQM)			# 0.68
HPDinterval(rQM)			# 0.57 to 0.74

cor.test(ant_data$lg_Qvol, ant_data$lg_Mvol)		# 0.80 < 0.82 < 0.85

## major axis phylogenetic slope estimate (Table 4, Warton et al. 2006) ##
# this treats males thorax volume as the X variable
# (1/ (2*cov(x,y)) ) * (varY - varX + sqrt(((varY - varX)^2) + (4*cov(x,y) ) ) )
posterior.mode((1/ (2*phylo_r_VCV[,"QMphylo"])) * (phylo_r_VCV[,"QQphylo"] - phylo_r_VCV[,"MMphylo"] + sqrt(((phylo_r_VCV[,"QQphylo"] - phylo_r_VCV[,"MMphylo"])^2) + (4*phylo_r_VCV[,"QMphylo"]))))
HPDinterval((1/ (2*phylo_r_VCV[,"QMphylo"])) * (phylo_r_VCV[,"QQphylo"] - phylo_r_VCV[,"MMphylo"] + sqrt(((phylo_r_VCV[,"QQphylo"] - phylo_r_VCV[,"MMphylo"])^2) + (4*phylo_r_VCV[,"QMphylo"]))))
# slope = 0.75, CI = 0.62 to 0.93


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 2: QUEEN SIZE IN RELATION TO COLONY SIZE, CASTE NUMBER, AND QUEEN NUMBER ###


## Colony Size Main Effect ##

# only include species with colony size estimates
queen_colony_data <- ant_data[-which(is.na(ant_data$lg_colony_size)),]
length(unique(queen_colony_data$tip_label))			# 154 species
plot(lg_Qvol ~ lg_colony_size, queen_colony_data)
abline(lm(lg_Qvol ~ lg_colony_size, queen_colony_data))
summary(lm(lg_Qvol ~ lg_colony_size, queen_colony_data))

qn_col_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  qn_col_mods[[i]] <- MCMCglmm(lg_Qvol ~ lg_colony_size, random = ~tip_label, family="gaussian", data=queen_colony_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorC, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(qn_col_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
qn_col_VCV <- mcmc(data.frame(phylo = 1:1000, units = 1:1000))
qn_col_Sol <- mcmc(data.frame(int = 1:1000, slope = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){qn_col_VCV[ranges$start[i]:ranges$stop[i],] <- qn_col_mods[[i]]$VCV[store,]}
for(i in 1:50){qn_col_Sol[ranges$start[i]:ranges$stop[i],] <- qn_col_mods[[i]]$Sol[store,]}

# chain convergence
plot(qn_col_VCV)    		# all parameters well mixed
plot(qn_col_Sol)   			# all parameters well mixed
autocorr(qn_col_VCV)   		# correlation between successive samples < 0.1
autocorr(qn_col_Sol)		# correlation between successive samples < 0.1
heidel.diag(qn_col_VCV)		# phylo and units passed halfwidth
heidel.diag(qn_col_Sol)		# intercept and slope passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(qn_col_Sol),2)			# intercept = -0.76, slope = 0.29
round(HPDinterval(qn_col_Sol),2)			# intercept = -2.00 to 0.34, slope = 0.21 to 0.38
round(posterior.mode(exp(qn_col_Sol)),2)	# intercept = 0.46, slope = 1.34
round(HPDinterval(exp(qn_col_Sol)),2)		# intercept = 0.11 to 1.33, slope = 1.23 to 1.46
# I^2 values
mean(qn_col_VCV[,1] / (qn_col_VCV[,1] + qn_col_VCV[,2]))		# phylogeny = 0.87
mean(qn_col_VCV[,2] / (qn_col_VCV[,1] + qn_col_VCV[,2]))		# residual = 0.13



## Caste Number Main Effect ##

# only include species with caste number estimates
queen_caste_data <- ant_data[-which(is.na(ant_data$bin_caste)),]
length(unique(queen_caste_data$tip_label))			# 187 species
plot(lg_Qvol ~ bin_caste, queen_caste_data)
abline(lm(lg_Qvol ~ bin_caste, queen_caste_data))
summary(lm(lg_Qvol ~ bin_caste, queen_caste_data))

qn_caste_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  qn_caste_mods[[i]] <- MCMCglmm(lg_Qvol ~ bin_caste, random = ~tip_label, family="gaussian", data=queen_caste_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorC, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(qn_caste_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
qn_caste_VCV <- mcmc(data.frame(phylo = 1:1000, units = 1:1000))
qn_caste_Sol <- mcmc(data.frame(int = 1:1000, slope = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){qn_caste_VCV[ranges$start[i]:ranges$stop[i],] <- qn_caste_mods[[i]]$VCV[store,]}
for(i in 1:50){qn_caste_Sol[ranges$start[i]:ranges$stop[i],] <- qn_caste_mods[[i]]$Sol[store,]}

# chain convergence
plot(qn_caste_VCV)    		# all parameters well mixed
plot(qn_caste_Sol)   		# all parameters well mixed
autocorr(qn_caste_VCV)   	# correlation between successive samples < 0.1
autocorr(qn_caste_Sol)		# correlation between successive samples < 0.1
heidel.diag(qn_caste_VCV)	# phylo and units passed halfwidth
heidel.diag(qn_caste_Sol)	# intercept and slope passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(qn_caste_Sol),2)			# intercept = 0.76, slope = 1.27
round(HPDinterval(qn_caste_Sol),2)				# intercept = -0.30 to 1.84, slope = 0.77 to 1.96
round(posterior.mode(exp(qn_caste_Sol)),2)		# intercept = 1.87, slope = 3.56
round(HPDinterval(exp(qn_caste_Sol)),2)			# intercept = 0.52 to 5.75, slope = 1.77 to 6.40
# I^2 values
mean(qn_caste_VCV[,1] / (qn_caste_VCV[,1] + qn_caste_VCV[,2]))		# phylogeny = 0.85
mean(qn_caste_VCV[,2] / (qn_caste_VCV[,1] + qn_caste_VCV[,2]))		# residual = 0.15



## Queen Number Main Effect ##

# only include species with colony size estimates
queen_number_data <- ant_data[-which(is.na(ant_data$lg_queen_number)),]
length(unique(queen_number_data$tip_label))			# 51 species
plot(lg_Qvol ~ lg_queen_number, queen_number_data)
abline(lm(lg_Qvol ~ lg_queen_number, queen_number_data))
summary(lm(lg_Qvol ~ lg_queen_number, queen_number_data))

qn_num_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  qn_num_mods[[i]] <- MCMCglmm(lg_Qvol ~ lg_queen_number, random = ~tip_label, family="gaussian", data=queen_number_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorC, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(qn_num_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
qn_num_VCV <- mcmc(data.frame(phylo = 1:1000, units = 1:1000))
qn_num_Sol <- mcmc(data.frame(int = 1:1000, slope = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){qn_num_VCV[ranges$start[i]:ranges$stop[i],] <- qn_num_mods[[i]]$VCV[store,]}
for(i in 1:50){qn_num_Sol[ranges$start[i]:ranges$stop[i],] <- qn_num_mods[[i]]$Sol[store,]}

# chain convergence
plot(qn_num_VCV)   	 		# all parameters well mixed
plot(qn_num_Sol)   			# all parameters well mixed
autocorr(qn_num_VCV)   		# correlation between successive samples < 0.1
autocorr(qn_num_Sol)		# correlation between successive samples < 0.1
heidel.diag(qn_num_VCV)		# phylo passed halfwidth
heidel.diag(qn_num_Sol)		# intercept passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(qn_num_Sol),2)				# intercept = 1.58, slope = -0.03
round(HPDinterval(qn_num_Sol),2)				# intercept = 0.38 to 2.70, slope = -0.27 to 0.33
round(posterior.mode(exp(qn_num_Sol)),2)		# intercept = 2.57, slope = 0.97
round(HPDinterval(exp(qn_num_Sol)),2)			# intercept = 0.68 to 12.29, slope = 0.74 to 1.36
# I^2 values
mean(qn_num_VCV[,1] / (qn_num_VCV[,1] + qn_num_VCV[,2]))		# phylogeny = 0.90
mean(qn_num_VCV[,2] / (qn_num_VCV[,1] + qn_num_VCV[,2]))		# residual = 0.10



## Colony Size - Queen Number Interaction ##

# only include species with colony size and queen number estimates
dim(na.omit(ant_data[,c("colony_size", "queen_number_continuous")]))		# 46 species
queen_int_data <- ant_data[-c(which(is.na(ant_data$lg_colony_size)), which(is.na(ant_data$lg_queen_number))),]
length(unique(queen_int_data$tip_label))			# 46 species
summary(lm(lg_Qvol ~ lg_colony_size * lg_queen_number, queen_int_data))

qn_int_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  qn_int_mods[[i]] <- MCMCglmm(lg_Qvol ~ lg_colony_size * lg_queen_number, random = ~tip_label, family="gaussian", data=queen_int_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorC, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(qn_int_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
qn_int_VCV <- mcmc(data.frame(phylo = 1:1000, units = 1:1000))
qn_int_Sol <- mcmc(data.frame(int = 1:1000, slope1 = 1:1000, slope2 = 1:1000, cross = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){qn_int_VCV[ranges$start[i]:ranges$stop[i],] <- qn_int_mods[[i]]$VCV[store,]}
for(i in 1:50){qn_int_Sol[ranges$start[i]:ranges$stop[i],] <- qn_int_mods[[i]]$Sol[store,]}

# chain convergence
plot(qn_int_VCV)    		# all parameters well mixed
plot(qn_int_Sol)   			# all parameters well mixed
autocorr(qn_int_VCV)   		# correlation between successive samples < 0.1
autocorr(qn_int_Sol)		# correlation between successive samples < 0.1
heidel.diag(qn_int_VCV)		# phylo passed halfwidth
heidel.diag(qn_int_Sol)		# slope1 passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(qn_int_Sol),2)
# intercept = -0.09
# slope 1 = 0.25
# slope 2 = 0.16
# interaction = 0.00
round(HPDinterval(qn_int_Sol),2)
# intercept = -1.77 to 1.35
# slope 1 = 0.06 to 0.46
# slope 2 = -0.77 to 0.80
# interaction = -0.09 to 0.11
# I^2 values
mean(qn_int_VCV[,1] / (qn_int_VCV[,1] + qn_int_VCV[,2]))		# phylogeny = 0.87
mean(qn_int_VCV[,2] / (qn_int_VCV[,1] + qn_int_VCV[,2]))		# residual = 0.13



## Caste Number - Queen Number Interaction ##

# only include species with caste number and queen number estimates
dim(na.omit(ant_data[,c("bin_caste", "queen_number_continuous")]))		# 46 species
queen_cqn_data <- ant_data[-c(which(is.na(ant_data$bin_caste)), which(is.na(ant_data$lg_queen_number))),]
length(unique(queen_cqn_data$tip_label))			# 47 species
summary(lm(lg_Qvol ~ bin_caste * lg_queen_number, queen_cqn_data))

qn_cqn_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  qn_cqn_mods[[i]] <- MCMCglmm(lg_Qvol ~ bin_caste * lg_queen_number, random = ~tip_label, family="gaussian", data=queen_cqn_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorC, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(qn_cqn_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
qn_cqn_VCV <- mcmc(data.frame(phylo = 1:1000, units = 1:1000))
qn_cqn_Sol <- mcmc(data.frame(int = 1:1000, slope1 = 1:1000, slope2 = 1:1000, cross = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){qn_cqn_VCV[ranges$start[i]:ranges$stop[i],] <- qn_cqn_mods[[i]]$VCV[store,]}
for(i in 1:50){qn_cqn_Sol[ranges$start[i]:ranges$stop[i],] <- qn_cqn_mods[[i]]$Sol[store,]}

# chain convergence
plot(qn_cqn_VCV)    		# all parameters well mixed
plot(qn_cqn_Sol)   			# all parameters well mixed
autocorr(qn_cqn_VCV)   		# correlation between successive samples < 0.1
autocorr(qn_cqn_Sol)		# correlation between successive samples < 0.1
heidel.diag(qn_cqn_VCV)		# phylo passed halfwidth
heidel.diag(qn_cqn_Sol)		# slope1 passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(qn_cqn_Sol),2)
# intercept = 1.36
# slope 1 = 0.68
# slope 2 = -0.19
# interaction = 0.31
round(HPDinterval(qn_cqn_Sol),2)	
# intercept = 0.40 to 2.28
# slope 1 = -0.06 to 2.15
# slope 2 = -0.51 to 0.27
# interaction = -0.36 to 0.84
# I^2 values
mean(qn_cqn_VCV[,1] / (qn_cqn_VCV[,1] + qn_cqn_VCV[,2]))		# phylogeny = 0.85
mean(qn_cqn_VCV[,2] / (qn_cqn_VCV[,1] + qn_cqn_VCV[,2]))		# residual = 0.15


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 3: DIMORPHISM IN RELATION TO COLONY SIZE, CASTE NUMBER, AND QUEEN NUMBER ###


## Colony Size Main Effect ##

# only include species with colony size estimates
di_colony_data <- ant_data[-which(is.na(ant_data$lg_colony_size)),]
length(unique(di_colony_data$tip_label))			# 154 species
plot(delta_QM_vol ~ lg_colony_size, di_colony_data)
abline(lm(delta_QM_vol ~ lg_colony_size, di_colony_data))
summary(lm(delta_QM_vol ~ lg_colony_size, di_colony_data))

di_col_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  di_col_mods[[i]] <- MCMCglmm(delta_QM_vol ~ lg_colony_size, random = ~tip_label, family="gaussian", data=di_colony_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorC, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(di_col_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
di_col_VCV <- mcmc(data.frame(phylo = 1:1000, units = 1:1000))
di_col_Sol <- mcmc(data.frame(int = 1:1000, slope = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){di_col_VCV[ranges$start[i]:ranges$stop[i],] <- di_col_mods[[i]]$VCV[store,]}
for(i in 1:50){di_col_Sol[ranges$start[i]:ranges$stop[i],] <- di_col_mods[[i]]$Sol[store,]}

# chain convergence
plot(di_col_VCV)    			# all parameters well mixed
plot(di_col_Sol)   				# all parameters well mixed
autocorr(di_col_VCV)   			# correlation between successive samples < 0.1
autocorr(di_col_Sol)			# correlation between successive samples < 0.1
heidel.diag(di_col_VCV)			# phylo and units passed halfwidth
heidel.diag(di_col_Sol)			# slope passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(di_col_Sol),2)			# intercept = 0.13, slope = 0.04
round(HPDinterval(di_col_Sol),2)			# intercept = -0.68 to 0.73, slope = -0.02 to 0.09
# I^2 values
mean(di_col_VCV[,1] / (di_col_VCV[,1] + di_col_VCV[,2]))		# phylogeny = 0.76
mean(di_col_VCV[,2] / (di_col_VCV[,1] + di_col_VCV[,2]))		# residual = 0.24



## Caste Number Main Effect ##

# only include species with caste size estimates
di_caste_data <- ant_data[-which(is.na(ant_data$bin_caste)),]
length(unique(di_caste_data$tip_label))			# 187 species
plot(delta_QM_vol ~ bin_caste, di_caste_data)
abline(lm(delta_QM_vol ~ bin_caste, di_caste_data))
summary(lm(delta_QM_vol ~ bin_caste, di_caste_data))

di_caste_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  di_caste_mods[[i]] <- MCMCglmm(delta_QM_vol ~ bin_caste, random = ~tip_label, family="gaussian", data=di_caste_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorC, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(di_caste_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
di_caste_VCV <- mcmc(data.frame(phylo = 1:1000, units = 1:1000))
di_caste_Sol <- mcmc(data.frame(int = 1:1000, slope = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){di_caste_VCV[ranges$start[i]:ranges$stop[i],] <- di_caste_mods[[i]]$VCV[store,]}
for(i in 1:50){di_caste_Sol[ranges$start[i]:ranges$stop[i],] <- di_caste_mods[[i]]$Sol[store,]}

# chain convergence
plot(di_caste_VCV)    			# all parameters well mixed
plot(di_caste_Sol)   			# all parameters well mixed
autocorr(di_caste_VCV)   		# correlation between successive samples < 0.1
autocorr(di_caste_Sol)			# correlation between successive samples < 0.1
heidel.diag(di_caste_VCV)		# phylo and units passed halfwidth
heidel.diag(di_caste_Sol)		# slope and intercept passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(di_caste_Sol),2)		# intercept = 0.18, slope = 0.33
round(HPDinterval(di_caste_Sol),2)			# intercept = -0.29 to 0.83, slope = 0.04 to 0.75
# I^2 values
mean(di_caste_VCV[,1] / (di_caste_VCV[,1] + di_caste_VCV[,2]))		# phylogeny = 0.75
mean(di_caste_VCV[,2] / (di_caste_VCV[,1] + di_caste_VCV[,2]))		# residual = 0.25



## Queen Number Main Effect ##

# only include species with colony size estimates
di_number_data <- ant_data[-which(is.na(ant_data$lg_queen_number)),]
length(unique(di_number_data$tip_label))			# 51 species
plot(delta_QM_vol ~ lg_queen_number, di_number_data)
abline(lm(delta_QM_vol ~ lg_queen_number, di_number_data))
summary(lm(delta_QM_vol ~ lg_queen_number, di_number_data))

di_num_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  di_num_mods[[i]] <- MCMCglmm(delta_QM_vol ~ lg_queen_number, random = ~tip_label, family="gaussian", data=di_number_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorC, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(di_num_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
di_num_VCV <- mcmc(data.frame(phylo = 1:1000, units = 1:1000))
di_num_Sol <- mcmc(data.frame(int = 1:1000, slope = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){di_num_VCV[ranges$start[i]:ranges$stop[i],] <- di_num_mods[[i]]$VCV[store,]}
for(i in 1:50){di_num_Sol[ranges$start[i]:ranges$stop[i],] <- di_num_mods[[i]]$Sol[store,]}

# chain convergence
plot(di_num_VCV)   	 		# all parameters well mixed
plot(di_num_Sol)   			# all parameters well mixed
autocorr(di_num_VCV)   		# correlation between successive samples < 0.1
autocorr(di_num_Sol	)		# correlation between successive samples < 0.1
heidel.diag(di_num_VCV)		# phylo and units passed halfwidth
heidel.diag(di_num_Sol)		# intercept and slope passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(di_num_Sol),2)			# intercept = 0.29, slope = 0.07
round(HPDinterval(di_num_Sol),2)			# intercept = -0.38 to 0.92, slope = -0.14 to 0.26
# I^2 values
mean(di_num_VCV[,1] / (di_num_VCV[,1] + di_num_VCV[,2]))		# phylogeny = 0.74
mean(di_num_VCV[,2] / (di_num_VCV[,1] + di_num_VCV[,2]))		# residual = 0.26



## Colony Size - Queen Number Interaction ##

# only include species with colony size and queen number estimates
dim(na.omit(ant_data[,c("colony_size", "queen_number_continuous")]))		# 46 species
di_int_data <- ant_data[-c(which(is.na(ant_data$lg_colony_size)), which(is.na(ant_data$lg_queen_number))),]
length(unique(di_int_data$tip_label))			# 46 species
summary(lm(delta_QM_vol ~ lg_colony_size * lg_queen_number, di_int_data))

di_int_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  di_int_mods[[i]] <- MCMCglmm(delta_QM_vol ~ lg_colony_size * lg_queen_number, random = ~tip_label, family="gaussian", data=di_int_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorC, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(di_int_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
di_int_VCV <- mcmc(data.frame(phylo = 1:1000, units = 1:1000))
di_int_Sol <- mcmc(data.frame(int = 1:1000, slope1 = 1:1000, slope2 = 1:1000, cross = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){di_int_VCV[ranges$start[i]:ranges$stop[i],] <- di_int_mods[[i]]$VCV[store,]}
for(i in 1:50){di_int_Sol[ranges$start[i]:ranges$stop[i],] <- di_int_mods[[i]]$Sol[store,]}

# chain convergence
plot(di_int_VCV)    		# all parameters well mixed
plot(di_int_Sol)   			# all parameters well mixed
autocorr(di_int_VCV)   		# correlation between successive samples < 0.1
autocorr(di_int_Sol)		# correlation between successive samples < 0.1
heidel.diag(di_int_VCV)		# phylo and units passed halfwidth
heidel.diag(di_int_Sol)		# intercept passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(di_int_Sol),2)
# intercept = 0.33
# slope 1 = -0.03
# slope 2 = 0.35
# interaction = -0.03
round(HPDinterval(di_int_Sol),2)	
# intercept = -0.67 to 1.28
# slope 1 = -0.15 to 0.13
# slope 2 = -0.14 to 0.94
# interaction = -0.11 to 0.04
# I^2 values
mean(di_int_VCV[,1] / (di_int_VCV[,1] + di_int_VCV[,2]))		# phylogeny = 0.70
mean(di_int_VCV[,2] / (di_int_VCV[,1] + di_int_VCV[,2]))		# residual = 0.30



## Caste Number - Queen Number Interaction ##

# only include species with colony size and queen number estimates
dim(na.omit(ant_data[,c("bin_caste", "queen_number_continuous")]))		# 46 species
di_cqn_data <- ant_data[-c(which(is.na(ant_data$bin_caste)), which(is.na(ant_data$lg_queen_number))),]
length(unique(di_cqn_data$tip_label))			# 46 species
summary(lm(delta_QM_vol ~ bin_caste * lg_queen_number, di_cqn_data))

di_cqn_mods <- list()          
for(i in 1:50){
  INtree <- inverseA(ant_trees[[i]], nodes="ALL")
  di_cqn_mods[[i]] <- MCMCglmm(delta_QM_vol ~ bin_caste * lg_queen_number, random = ~tip_label, family="gaussian", data=di_cqn_data, pl=TRUE, nitt=1100000, thin=1000, burnin=100000, prior=priorC, ginverse=list(tip_label=INtree$Ainv), verbose=FALSE)
}

summary(di_cqn_mods[[5]])

# create empty mcmc objects to combine posterior chains for random and fixed effect estimates
di_cqn_VCV <- mcmc(data.frame(phylo = 1:1000, units = 1:1000))
di_cqn_Sol <- mcmc(data.frame(int = 1:1000, slope1 = 1:1000, slope2 = 1:1000, cross = 1:1000))

# run loop to populate mcmc objects
for(i in 1:50){di_cqn_VCV[ranges$start[i]:ranges$stop[i],] <- di_cqn_mods[[i]]$VCV[store,]}
for(i in 1:50){di_cqn_Sol[ranges$start[i]:ranges$stop[i],] <- di_cqn_mods[[i]]$Sol[store,]}

# chain convergence
plot(di_cqn_VCV)    		# all parameters well mixed
plot(di_cqn_Sol)   			# all parameters well mixed
autocorr(di_cqn_VCV)   		# correlation between successive samples < 0.1
autocorr(di_cqn_Sol)		# correlation between successive samples < 0.1
heidel.diag(di_cqn_VCV)		# phylo and units passed halfwidth
heidel.diag(di_cqn_Sol)		# intercept passed halfwidth

# parameter estimates
# fixed effects
round(posterior.mode(di_cqn_Sol),2)
# intercept = 0.38
# slope 1 = -0.18
# slope 2 = -0.20
# interaction = 0.54
round(HPDinterval(di_cqn_Sol),2)	
# intercept = -0.05 to 1.09
# slope 1 = -0.98 to 0.52
# slope 2 = -0.48 to 0.11
# interaction = 0.20 to 1.09
# I^2 values
mean(di_cqn_VCV[,1] / (di_cqn_VCV[,1] + di_cqn_VCV[,2]))		# phylogeny = 0.67
mean(di_cqn_VCV[,2] / (di_cqn_VCV[,1] + di_cqn_VCV[,2]))		# residual = 0.33


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #
												### THE END ###
# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #