# set the working directoy
# my.wd = "~/Desktop/experiment_alarm_network"
# setwd(my.wd)

# Giulio
setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiment_alarm_network//giulio/')

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

# generate random data from that model 
bnalarm = bn.fit(res, alarm)
dataset = rbn(bnalarm, n = 200)

#### CRASH - perform the test
results = perform.bootstrap.inference(dataset,regularization,boot.first.pass,boot.second.pass,test.pvalue)

# save the results
results_alarm = list()
results_alarm[["true_adj_matrix"]] = adj.matrix
results_alarm[["inference"]] = results
save(results_alarm,file="results_alarm.RData")

source('plot/plotter.R')

plt(results$agony.inference, results$agony.inference.pvalues, "pvalues", adj.matrix, 100, 0.01, 'Agony without MHC')
dev.copy2pdf(file='100_agony_pvalues.pdf')
plt(results$agony.inference, 'qvalues.fdr', adj.matrix, 100, 0.01, 'Agony with FDR')
dev.copy2pdf(file='100_agony_fdr.pdf')
plt(results$agony.inference, 'qvalues.holm', adj.matrix, 100, 0.01, 'Agony with Bonferroni')
dev.copy2pdf(file='100_agony_bonferroni.pdf')

plt(results$confidence.inference, "pvalues", adj.matrix, 100, 0.01, "Confidence without MHC")
dev.copy2pdf(file = "100_confidence_pvalues.pdf")
plt(results$confidence.inference, "qvalues.fdr", adj.matrix, 100, 0.01, "Confidence with FDR")
dev.copy2pdf(file = "100_confidence_fdr.pdf")
plt(results$confidence.inference, "qvalues.holm", adj.matrix, 100, 0.01, "Confidence with Bonferroni")
dev.copy2pdf(file = "100_confidence_bonferroni.pdf")


hc0 = hc(dataset, restart = 0, score = regularization)
hc100 = hc(dataset, restart = boot.first.pass + boot.second.pass, score = regularization)

h = NULL
h$hc0 = amat(hc0)
h$hc100 = amat(hc100)

plt(h, 'hc0', adj.matrix, 0, 'none', 'Hill Climbing')
dev.copy2pdf(file='0_HC.pdf')
plt(h, 'hc100', adj.matrix, 200,' none',  'Hill Climbing')
dev.copy2pdf(file='100_HC.pdf')
# plt(h, 'hc1000', 1100, 51, 'Hill Climbing')
# dev.copy2pdf(file='1000_HC.pdf')

names(results$agony.inference.pvalues$pvalues)

pvalues = results$agony.inference.pvalues$pvalues
# pvalues[pvalues == 1] = NA
ggd.qqplot(pvalues, "p-values' distribution")