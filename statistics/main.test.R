# load the required R packages
library(bnlearn)
library(parallel)
library(doParallel)
library(igraph)
library(infotheo)

# source the required scripts
source("bootstrap.likelihood.fit.R")
source("consensus.likelihood.fit.R")

file_true_adj_matrix = "adj_matrix.txt"
file_dataset_discrete = "dataset_sample_size_5_noise_level_0.txt"

# set the seed
set.seed(12345)

# read the true adjacency matrix
true_adj_matrix = read.table(file=file_true_adj_matrix,sep=",",check.names=FALSE,stringsAsFactors=FALSE)

# read the datasets
dataset_discrete = read.table(file=file_dataset_discrete,sep=",",check.names=FALSE,stringsAsFactors=FALSE)

# perform the estimation by bootstrap
bootstrap_results_discrete = bootstrap.estimation(dataset_discrete,regularization="bic")

# estimate the structure by consensus
results_discrete_consensus = perform.consensus.likelihood.fit(bootstrap_results_discrete,dataset_discrete,regularization="bic")

# compute the mutual information
mutual_information_discrete = mutinformation(dataset_discrete)

# RESULTS
results = list(ground.true=true_adj_matrix,dataset=dataset_discrete,bootstrap.scores=Reduce("+",bootstrap_results_discrete),mutual.information=mutual_information_discrete,reconstruction=results_discrete_consensus)

bootstrap.scores=Reduce("+",bootstrap_results_discrete)
source("guido.plot.R")
guido.plot(true_adj_matrix, true_adj_matrix, true_adj_matrix, bootstrap.scores, mutual_information_discrete, loop = NULL)

