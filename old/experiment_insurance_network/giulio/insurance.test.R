# set the working directoy
# my.wd = "~/Desktop/experiment_alarm_network"
# setwd(my.wd)

# Giulio
setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiment_insurance_network//giulio/')

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
data(insurance)
dataset = insurance
nrow(dataset)

res = empty.graph(names(insurance))
modelstring(res) = paste("[Age][Mileage][SocioEcon|Age][GoodStudent|Age:SocioEcon]",
  "[RiskAversion|Age:SocioEcon][OtherCar|SocioEcon][VehicleYear|SocioEcon:RiskAversion]",
  "[MakeModel|SocioEcon:RiskAversion][SeniorTrain|Age:RiskAversion]",
  "[HomeBase|SocioEcon:RiskAversion][AntiTheft|SocioEcon:RiskAversion]",
  "[RuggedAuto|VehicleYear:MakeModel][Antilock|VehicleYear:MakeModel]",
  "[DrivingSkill|Age:SeniorTrain][CarValue|VehicleYear:MakeModel:Mileage]",
  "[Airbag|VehicleYear:MakeModel][DrivQuality|RiskAversion:DrivingSkill]",
  "[Theft|CarValue:HomeBase:AntiTheft][Cushioning|RuggedAuto:Airbag]",
  "[DrivHist|RiskAversion:DrivingSkill][Accident|DrivQuality:Mileage:Antilock]",
  "[ThisCarDam|RuggedAuto:Accident][OtherCarCost|RuggedAuto:Accident]",
  "[MedCost|Age:Accident:Cushioning][ILiCost|Accident]",
  "[ThisCarCost|ThisCarDam:Theft:CarValue][PropCost|ThisCarCost:OtherCarCost]",
  sep = "")

adj.matrix = amat(res)

# generate random data from that model 
# bnalarm = bn.fit(res, insurance)
# dataset = rbn(bnalarm, n = 200)

#### CRASH - perform the test
results = perform.bootstrap.inference(dataset,regularization,boot.first.pass,boot.second.pass,test.pvalue)

# save the results
results_alarm = list()
results_alarm[["true_adj_matrix"]] = adj.matrix
results_alarm[["inference"]] = results
save(results_alarm,file="results_insurance.RData")

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

pvalues = results$agony.inference.pvalues$pvalues
# pvalues[pvalues == 1] = NA
ggd.qqplot(pvalues, "p-values' distribution")