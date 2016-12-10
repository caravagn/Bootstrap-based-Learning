# set the working directoy
my.wd = "~/Desktop/experiment_simulation_networks"
setwd(my.wd)

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
agony_files = paste0(getwd(),"/agony_files")
regularization = "bic"
boot.first.pass = 100
boot.second.pass = 100
test.pvalue = 0.01

as.categorical.dataset <- function(dataset){

    ## Create a categorical data frame from the dataset
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

    ## Renaming
    colnames(data) = colnames(dataset)
    rownames(data) = rownames(dataset)
    return(data)
}

# set the dataset and the true adjacency matrix for the test
dataset1 = read.table(file=paste0(getwd(),"/datasets/dataset_sample_size_5_noise_level_0.txt"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)
dataset1 = as.categorical.dataset(dataset1)
dataset2 = read.table(file=paste0(getwd(),"/datasets/dataset_sample_size_10_noise_level_0.txt"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)
dataset2 = as.categorical.dataset(dataset2)
dataset3 = read.table(file=paste0(getwd(),"/datasets/dataset_sample_size_50_noise_level_0.txt"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)
dataset3 = as.categorical.dataset(dataset3)
adj.matrix = read.table(file=paste0(getwd(),"/datasets/adj_matrix.txt"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)
colnames(adj.matrix) = as.character(1:ncol(adj.matrix))
rownames(adj.matrix) = as.character(1:nrow(adj.matrix))

# data(alarm)
# dataset = alarm
# adj.matrix = array(0,c(ncol(dataset),ncol(dataset)))
# colnames(adj.matrix) = colnames(dataset)
# rownames(adj.matrix) = colnames(dataset)
# res = empty.graph(colnames(dataset))
# modelstring(res) = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
    # "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]", "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]", 
    # "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]", "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]", 
    # "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]", "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
# my.arcs = res$arcs
# # create the adjacency matrix of the reconstructed topology
# for (i in 1:nrow(my.arcs)) {
    # # [i,j] refers to the edge i --> j
    # adj.matrix[which(colnames(dataset)%in%as.character(my.arcs[i,1])),which(colnames(dataset)%in%as.character(my.arcs[i,2]))] = 1
# }

# perform the test
unlink(agony_files,recursive=TRUE,force=TRUE)
dir.create(agony_files,showWarnings=FALSE)
results1 = perform.bootstrap.inference(dataset1,regularization,boot.first.pass,agony_files,boot.second.pass)
unlink(agony_files,recursive=TRUE,force=TRUE)
unlink(agony_files,recursive=TRUE,force=TRUE)
dir.create(agony_files,showWarnings=FALSE)
results2 = perform.bootstrap.inference(dataset2,regularization,boot.first.pass,agony_files,boot.second.pass)
unlink(agony_files,recursive=TRUE,force=TRUE)
unlink(agony_files,recursive=TRUE,force=TRUE)
dir.create(agony_files,showWarnings=FALSE)
results3 = perform.bootstrap.inference(dataset3,regularization,boot.first.pass,agony_files,boot.second.pass)
unlink(agony_files,recursive=TRUE,force=TRUE)
