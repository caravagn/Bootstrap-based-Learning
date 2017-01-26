GIULIO = TRUE
DANIELE = !GIULIO

if(GIULIO)
{
	setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/alarm-test-statistics/')
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

# set the seed
set.seed(12345)

# set some settings to be used in the test
regularization = "bic"
boot.first.pass = 100
boot.second.pass = 100
test.pvalue = 1E-2
numsamples = 100


repeat {
	dataset = rbn(bnalarm, n = numsamples)
	if(!any(is.na(dataset))) break
}


# Standard test
source(paste0(git, "src/perform.bootstrap.inference.R"))
results = perform.bootstrap.inference(
	dataset,
	regularization,
	boot.first.pass,
	boot.second.pass,
	test.pvalue,
	agony.binaries = agony.binaries)

# We create a dummy prior, and use Suppes' code
dummy = matrix(rep(1, ncol(adj.matrix) * nrow(adj.matrix)), nrow = nrow(adj.matrix), ncol = ncol(adj.matrix))
colnames(dummy) = colnames(adj.matrix)
rownames(dummy) = rownames(adj.matrix)
diag(dummy) = 0

source(paste0(git, "src/suppes.bootstrap.inference.R"))
results.dummy = suppes.bootstrap.inference(
	dataset,
	regularization,
	suppes.poset = dummy,
	nboot.second = boot.second.pass,
	test.pvalue,
	agony.binaries = agony.binaries,
	nboot.first = boot.first.pass
	)

source(paste0(git, "utils/plotter.R"))
source(paste0(git, "utils/stat.R"))


ABNF = paste0('Agony with Bonferroni [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')
NOT = paste0('No test [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')


DEV.OFF = T
STATS = TRUE

x1 = plt(results$agony.inference, results$agony.inference.pvalues, 'qvalues.holm', adj.matrix, ABNF, stats = STATS, dev.off = DEV.OFF)

x1 = plt(results$confidence.inference, results$confidence.inference.pvalues, 'qvalues.holm', adj.matrix, ABNF, stats = STATS, dev.off = DEV.OFF)


x2 = plt(results.dummy$confidence.inference, results.dummy$confidence.inference.pvalues, 'qvalues.holm', adj.matrix, NOT, stats = STATS, dev.off = DEV.OFF)

# x4 = plt(results$confidence.inference, results$confidence.inference.pvalues, "pvalues", adj.matrix, CMHC, stats = STATS, dev.off = DEV.OFF)
# x5 = plt(results$confidence.inference, results$confidence.inference.pvalues, 'qvalues.fdr', adj.matrix, CFDR, stats = STATS, dev.off = DEV.OFF)
# x6 = plt(results$confidence.inference, results$confidence.inference.pvalues, 'qvalues.holm', adj.matrix, CBNF, stats = STATS, dev.off = DEV.OFF)

# hc1 = plt(results, NULL, 'hill.climing.no.restarts.inference', adj.matrix, HC0, stats = STATS, dev.off = DEV.OFF)
# hc2 = plt(results, NULL, 'hill.climing.with.restarts.inference', adj.matrix, HC, stats = STATS, dev.off = DEV.OFF)

