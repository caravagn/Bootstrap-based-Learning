GIULIO = TRUE
DANIELE = !GIULIO

if(GIULIO)
{
	# Giulio
	setwd('~/bootstrap/12Dic/experiments/catnet')
	git = '~/bootstrap/12Dic/'
	agony.binaries = '~/bootstrap/12Dic/agony/giulio/agony'
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
library(catnet)

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

# set the seed
set.seed(12345)

# set some settings to be used in the test
regularization = "bic"
boot.first.pass = 100
boot.second.pass = 100
test.pvalue = 1E-2

each = 100
values = c(100, 1000, 10000, 20000)

res1E3 = res1E4 = res1E5 = res2E5 = NULL

TOT = each
dataset = NULL
do = function(x,y){
	cat("********** ", x, ":", y, "/", TOT, "\n")

	# generate random data from that model, ensure no NAs 
	repeat {
		dataset = rbn(bnalarm, n = x)
		if(!any(is.na(dataset))) break
	}

	return(perform.bootstrap.inference(
		dataset,
		regularization,
		boot.first.pass,
		boot.second.pass,
		test.pvalue,
		agony.binaries = agony.binaries,
		doP = ifelse(j==1,TRUE, FALSE)))
}

for(j in 1:each) res1E3 = append(res1E3, list(do(100,j)))
save(res1E3, file = 'res1E3.Rdata')

for(j in 1:each) res1E4 = append(res1E4, list(do(1000,j)))
save(res1E4, file = 'res1E4.Rdata')

for(j in 1:each) res1E5 = append(res1E5, list(do(10000,j)))
save(res1E5, file = 'res1E5.Rdata')

for(j in 1:each) res2E5 = append(res1E3, list(do(20000,j)))
save(res2E5, file = 'res2E5.Rdata')
	


#### perform 

