# load the required R packages
library(bnlearn)
library(parallel)
library(doParallel)

# source the required scripts
source("bootstrap.likelihood.fit.R")

# read the datasets
dataset_discrete = read.table(file="test_discrete_dataset_sample_size_10_noise_level_0.2.txt",sep=",",check.names=FALSE,stringsAsFactors=FALSE)
dataset_continuous = read.table(file="test_continuous_dataset_sample_size_10_noise_level_0.1.txt",sep=",",check.names=FALSE,stringsAsFactors=FALSE)

# perform the likelihood fit on the discrete dataset with all the available regularizators
results_discrete = list()
for (r in c("loglik","aic","bic")) {
    results_discrete[[r]] = likelihood.fit(dataset_discrete,regularization=r,command="hc") 
}

# perform the likelihood fit on the continuous dataset with all the available regularizators
results_continuous = list()
for (r in c("loglik-g","aic-g","bic-g")) {
    results_continuous[[r]] = likelihood.fit(dataset_continuous,regularization=r,command="hc") 
}

# perform the estimation by bootstrap
bootstrap_results_discrete = bootstrap.estimation(dataset_discrete,regularization="loglik")
bootstrap_results_continuous = bootstrap.estimation(dataset_continuous,regularization="loglik-g")
