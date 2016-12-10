# set the working directoy
# my.wd = "~/Desktop/experiment_alarm_network"
# setwd(my.wd)

# Giulio 
setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiment_alarm_network/giulio/')

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
boot.first.pass = 3
boot.second.pass = 3
test.pvalue = 0.01

# set the dataset and the true adjacency matrix for the test
data(alarm)
dataset = alarm

res = empty.graph(colnames(dataset))
modelstring(res) = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
    "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]", "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]", 
    "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]", "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]", 
    "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]", "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")

adj.matrix = amat(res) # bnlearn fun...
# library(pheatmap)
# pheatmap(adj.matrix)


# Moved inside perform.bootstrap.inference
# perform the test
# unlink(agony_files,recursive=TRUE,force=TRUE)
# dir.create(agony_files,showWarnings=FALSE)
results = perform.bootstrap.inference(dataset,regularization,boot.first.pass,boot.second.pass)
# unlink(agony_files,recursive=TRUE,force=TRUE)



# STEP 5: perform the final inference on both confidence based and agony based posets
# BINARY
confidence_based_inference_binary_example = perform.inference(confidence_bootstrap_inference_binary_example,confidence_based_poset_binary_example,test.pvalue)
agony_based_inference_binary_example = perform.inference(agony_bootstrap_inference_binary_example,agony_based_poset_binary_example,test.pvalue)
