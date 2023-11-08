########################################################################
#This script runs statistical tests for Figure 2 in the following manuscript:

#Bower, S.J., Shobe, C.M., Maxwell, A.E., and Campforts, B. (2023) The 
#uncertain future of mountaintop-removal-mined landscapes 2: Modeling the
#influence of topography and vegetation. Geomorphology.

#Please cite the paper if you use this code in any way.

#Brief description: this script uses the software and methodology of 
#
#       van Doorn, J.B., Ly, A., Marsman, M., & Wagenmakers, E.-J.  (2020). 
#           Bayesian Rank-Based Hypothesis Testing for the Rank Sum Test, the Signed Rank Test, 
#           and Spearman's rho. Journal of Applied Statistics, 47, 2984-3006.
#
#as downloaded from the article's online supplement: https://osf.io/gny35/
#
#to conduct Bayesian signed-rank tests to test whether pre- and post-mining
#distributions of elevation, slope, and area--slope product are different.

########################################################################

library(HDInterval)
library(logspline)

#source needed functions
source('rankBasedCommonFunctions.R')
source('signRankSampler.R') # Wilcoxon signed-rank function

bc_pre_elev <- scan('../flat_inputs/bc_pre_elev.txt')
bc_post_elev <- scan('../flat_inputs/bc_post_elev.txt')

pre_name <- 'bc_pre_elev'
post_name <- 'bc_post_elev'

bc_elev <- data.frame(bc_pre_elev, bc_post_elev)
names(bc_elev) <- c(pre_name,post_name)

#proportion of measurements to use 
#(because the DEMs have so many pixels we subsample for speed)
pro_measurements = 0.05

set.seed(123456)
random_subsample <- bc_elev[sample(nrow(bc_elev), round(length(bc_pre_elev)*pro_measurements,0)), ]

x <- random_subsample$bc_pre_elev 
y <- random_subsample$bc_post_elev

signedRankSamples <- signRankGibbsSampler(xVals = x,
                                          yVals = y,
                                          nSamples = 5e2, 
                                          cauchyPriorParameter = 1/sqrt(2), 
                                          testValue = 0, 
                                          progBar = FALSE, 
                                          nBurnin = 1, 
                                          nGibbsIterations = 10, 
                                          nChains = 10)

# Posterior distribution
hist(signedRankSamples$deltaSamples, freq = FALSE)

# Give the posterior samples for delta to the function below to compute BF01
bayes_factor_bc_elev <- computeBayesFactorOneZero(signedRankSamples$deltaSamples, 
                          whichTest = "Wilcoxon",
                          priorParameter = 1 / sqrt(2))

#calculate 95% and 99% HPDI (highest posterior density interval) 
cred_interval_99_bc_elev <- hdi(signedRankSamples$deltaSamples, 0.99)
cred_interval_95_bc_elev <- hdi(signedRankSamples$deltaSamples, 0.95)

#posterior median
posterior_median <- median(signedRankSamples$deltaSamples)

#stuff to save!
labels_column <- c("posterior_median", 
                   "bayes_factor",
                   "cred_interval_99_min",
                   "cred_interval_99_max",
                   "cred_interval_95_min",
                   "cred_interval_95_max")
values_column <- c(posterior_median,
                   bayes_factor_bc_elev,
                   cred_interval_99_bc_elev[1],
                   cred_interval_99_bc_elev[2],
                   cred_interval_95_bc_elev[1],
                   cred_interval_95_bc_elev[2])

#build dataframe
export_df <- data.frame(labels_column, values_column)
write.csv(export_df, '../R_outputs/sam_bc_elev.csv', row.names = FALSE)
write.csv(signedRankSamples$deltaSamples, '../R_outputs/sam_bc_elev_samples.csv', row.names = FALSE)
