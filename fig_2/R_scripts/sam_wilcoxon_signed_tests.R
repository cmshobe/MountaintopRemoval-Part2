library(HDInterval)
library(logspline)
#source needed functions
path = '/home/wvu_esd/Documents/Sammy B/MTR/Scripts/bayesian_wilcoxon/'
#setwd('/home/wvu_esd/Documents/Sammy B/MTR/Scripts/bayesian_wilcoxon')
source('/home/wvu_esd/Documents/Sammy B/MTR/Scripts/bayesian_wilcoxon/scripts/rankBasedCommonFunctions.R')
source('/home/wvu_esd/Documents/Sammy B/MTR/Scripts/bayesian_wilcoxon/scripts/signRankSampler.R') # Wilcoxon signed-rank function

#import_list <- list('sam_signed_rank_data/bc_pre_elev.txt', 'sam_signed_rank_data/bc_post_elev.txt')
#bc_elev <- read.table(import_list, sep = '\t') 

bc_pre_elev <- scan('/home/wvu_esd/Documents/Sammy B/MTR/Scripts/bayesian_wilcoxon/inputs/bc_pre_elev.txt')
bc_post_elev <- scan('/home/wvu_esd/Documents/Sammy B/MTR/Scripts/bayesian_wilcoxon/inputs/bc_post_elev.txt')

pre_name <- 'bc_pre_elev'
post_name <- 'bc_post_elev'

bc_elev <- data.frame(bc_pre_elev, bc_post_elev)
names(bc_elev) <- c(pre_name,post_name)

set.seed(123456)
random_subsample <- bc_elev[sample(nrow(bc_elev), 6000), ]

# Example 2: Wilcoxon Signed Rank
# Does progabide reduce the frequency of epileptic seizures? 
x <- random_subsample$bc_pre_elev # Alcohol intake of students who failed
y <- random_subsample$bc_post_elev # Alcohol intake of students who passed

# Default Cauchy prior width is set to 1/sqrt(2) for the sampler 
signedRankSamples <- signRankGibbsSampler(xVals = x,
                                          yVals = y,
                                          nSamples = 1e1, 
                                          cauchyPriorParameter = 1/sqrt(2), 
                                          testValue = 0, 
                                          progBar = TRUE, 
                                          nBurnin = 1, 
                                          nGibbsIterations = 10, 
                                          nChains = 10)
# OR
#signedRankSamples <- signRankGibbsSampler(xVals = x - y,
#                                          testValue = 0,
#                                          nSamples = 1e1,
#                                          progBar = TRUE)


# Posterior distribution
hist(signedRankSamples$deltaSamples, freq = FALSE)

# Give the posterior samples for delta to the function below to compute BF01
bayes_factor_bc_elev <- computeBayesFactorOneZero(signedRankSamples$deltaSamples, 
                          whichTest = "Wilcoxon",
                          priorParameter = 1 / sqrt(2))
#credible intervals
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
write.csv(export_df, '/home/wvu_esd/Documents/Sammy B/MTR/Scripts/bayesian_wilcoxon/outputs/sam_bc_elev.csv', row.names = FALSE)
write.csv(signedRankSamples$deltaSamples, '/home/wvu_esd/Documents/Sammy B/MTR/Scripts/bayesian_wilcoxon/outputs/sam_bc_elev_samples.csv', row.names = FALSE)
