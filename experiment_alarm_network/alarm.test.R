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

# set the dataset and the true adjacency matrix for the test
data(alarm)
dataset = alarm

res = empty.graph(colnames(dataset))
modelstring(res) = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
    "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]", "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]", 
    "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]", "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]", 
    "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]", "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")

adj.matrix = amat(res)

# perform the test
results = perform.bootstrap.inference(dataset,regularization,boot.first.pass,boot.second.pass,test.pvalue)

# save the results
results_alarm = list()
results_alarm[["true_adj_matrix"]] = adj.matrix
results_alarm[["inference"]] = results
save(results_alarm,file="results_alarm.RData")
