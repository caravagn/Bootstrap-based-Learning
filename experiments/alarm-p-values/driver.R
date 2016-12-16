GIULIO = TRUE
DANIELE = !GIULIO

if(GIULIO)
{
	setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/alarm-p-values/')
	git = '/Volumes/Data/Github/Bootstrap-based-Learning/'
	agony.binaries = '/Volumes/Data/Github/Bootstrap-based-Learning/agony/giulio/agony'
} 
if(DANIELE)
{
	setwd('~/Desktop/')
	git = '....'
	agony.binaries = '/Volumes/Data/Github/Bootstrap-based-Learning/agony/daniele/agony'	
}

# load the required R packages
library(bnlearn)
library(parallel)
library(doParallel)
library(igraph)

# set the dataset and the true adjacency matrix for the test
data(alarm)

bnalarm = empty.graph(colnames(alarm))
modelstring(bnalarm) = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
    "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]", "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]", 
    "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]", "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]", 
    "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]", "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
bnalarm = bn.fit(bnalarm, alarm)

adj.matrix = amat(bnalarm)

# source the required scripts
source(paste0(git, "src/perform.bootstrap.inference.R"))
source(paste0(git, "utils/plotter.R"))
source(paste0(git, "utils/stat.R"))

# set the seed
set.seed(12345)

# set some settings to be used in the test
regularization = "bic"
boot.first.pass = 100
boot.second.pass = 100
numsamples = 100

res05 = list()
res01 = list()
res001 = list()

sample_size = 29
for(i in 1:sample_size)
{
	cat(i, '\n')
# generate random data from that model, ensure no NAs 
repeat {
	dataset = rbn(bnalarm, n = numsamples)
	if(!any(is.na(dataset))) break
}	

#### perform the test
results05 = perform.bootstrap.inference(
	dataset,
	regularization,
	boot.first.pass,
	boot.second.pass,
	test.pvalue = 0.05,
	agony.binaries = agony.binaries,
	do.hc = FALSE)

results01 = perform.bootstrap.inference(
	dataset,
	regularization,
	boot.first.pass,
	boot.second.pass,
	test.pvalue = 0.01,
	agony.binaries = agony.binaries,
	do.hc = FALSE)

results001 = perform.bootstrap.inference(
	dataset,
	regularization,
	boot.first.pass,
	boot.second.pass,
	test.pvalue = 0.001,
	agony.binaries = agony.binaries,
	do.hc = FALSE)


### Save results
wrapper = NULL
wrapper$results05 = results05
wrapper$results01 = results01
wrapper$results001 = results001
wrapper$params = c(regularization, boot.first.pass, boot.second.pass, 0.5, 0.01, 0.001)
wrapper$randomnet = bnalarm
wrapper$numsamples = numsamples
wrapper$dataset = dataset
file = paste0(stringi::stri_rand_strings(1, 8), '.Rdata')
print(paste('Results saved to:', file))
save(wrapper, file = file)


source('plot.R')

res05 = append(res05, list(x1))
res01 = append(res01, list(x2))
res001 = append(res001, list(x3))
}


sta05 = Reduce(cbind, res05)
sta01 = Reduce(cbind, res01)
sta001 = Reduce(cbind, res001)

library(vioplot)

par(mfrow = c(1,3))
vioplot(sta05['PPV', ], sta01['PPV', ], sta001['PPV', ], col = 'cornflowerblue', lty = 1, rectCol="gray",
  colMed = 'black', names = c('0.05', '0.01', '0.001'), pchMed = 15, horizontal = T)
title(main = 'PPV')

vioplot(sta05['TPR', ], sta01['TPR', ], sta001['TPR', ], col = 'cornflowerblue', lty = 1, rectCol="gray",
  colMed = 'black', names = c('0.05', '0.01', '0.001'), pchMed = 15, horizontal = T)
title(main = 'TPR')

dev.copy2pdf(file = 'p-values-violinplot.pdf')