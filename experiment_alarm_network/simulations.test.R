# set the working directoy
my.wd = "~/Desktop/experiment_alarm_network"
setwd(my.wd)

# Giulio
# setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiment_alarm_network/giulio/')

# load the required R packages
library(bnlearn)
library(parallel)
library(doParallel)
library(igraph)

# source the required scripts
source(paste0(getwd(),"/perform.bootstrap.inference.R"))

# set the seed
set.seed(12345)

# set some settings to be used in the test
regularization = "bic"
boot.first.pass = 100
boot.second.pass = 100
test.pvalue = 0.01

as.categorical.dataset <- function(dataset) {

    # Create a categorical data frame from the dataset
    data = array("missing", c(nrow(dataset), ncol(dataset)))

    for (i in 1:nrow(dataset)) {
        for (j in 1:ncol(dataset)) {
            if (dataset[i,j] == 1) {
                data[i,j] = "observed"
            }
        }
    }

    data = data.frame(data, stringsAsFactors = TRUE)
    for (n in names(data)) {
        levels(data[[n]]) = c('missing', 'observed')
    }

    # Renaming
    colnames(data) = colnames(dataset)
    rownames(data) = rownames(dataset)
    return(data)
    
}

# set the dataset and the true adjacency matrix for the test
dataset1 = read.table(file=paste0(getwd(),"/datasets_simulations_test/dataset_sample_size_5_noise_level_0.txt"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)
dataset1 = as.categorical.dataset(dataset1)
dataset2 = read.table(file=paste0(getwd(),"/datasets_simulations_test/dataset_sample_size_10_noise_level_0.txt"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)
dataset2 = as.categorical.dataset(dataset2)
dataset3 = read.table(file=paste0(getwd(),"/datasets_simulations_test/dataset_sample_size_50_noise_level_0.txt"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)
dataset3 = as.categorical.dataset(dataset3)
adj.matrix = read.table(file=paste0(getwd(),"/datasets_simulations_test/adj_matrix.txt"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)
colnames(adj.matrix) = as.character(1:ncol(adj.matrix))
rownames(adj.matrix) = as.character(1:nrow(adj.matrix))

# perform the test
results1 = perform.bootstrap.inference(dataset1,regularization,boot.first.pass,boot.second.pass,agony_files=paste0(getwd(),"/agony_files_sim1"))
results2 = perform.bootstrap.inference(dataset2,regularization,boot.first.pass,boot.second.pass,agony_files=paste0(getwd(),"/agony_files_sim2"))
results3 = perform.bootstrap.inference(dataset3,regularization,boot.first.pass,boot.second.pass,agony_files=paste0(getwd(),"/agony_files_sim3"))
