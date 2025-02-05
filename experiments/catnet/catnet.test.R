GIULIO = TRUE
DANIELE = !GIULIO

if(GIULIO)
{
	# Giulio
	setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/catnet/')
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
library(catnet)


# source the required scripts
source(paste0(git, "src/perform.bootstrap.inference.R"))

# set the seed
set.seed(12345)

# set some settings to be used in the test
regularization = "bic"
boot.first.pass = 100
boot.second.pass = 100
test.pvalue = 1E-2
numsamples = 100

# Random catnet network, with adjacency matrix and data
randomnet = cnRandomCatnet(
	numnodes = 20,
	maxParents = 4,
	numCategories = 3)

adj.matrix = t(cnMatParents(randomnet)) 
dataset = cnSamples(randomnet, numsamples = numsamples)

#### perform the test
results = perform.bootstrap.inference(
	dataset,
	regularization,
	boot.first.pass,
	boot.second.pass,
	test.pvalue,
	agony.binaries = agony.binaries)

### Save results
wrapper = NULL
wrapper$results = results
wrapper$params = c(regularization, boot.first.pass, boot.second.pass, test.pvalue)
wrapper$randomnet = randomnet
wrapper$numsamples = numsamples
wrapper$dataset = dataset
file = paste0(stringi::stri_rand_strings(1, 8), '.Rdata')
print(paste('Results saved to:', file))
save(wrapper, file = file)

### Some analysis
source(paste0(git, "utils/plotter.R"))
source(paste0(git, "utils/stat.R"))

AMHC = paste0('Agony without MHC [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')
AFDR = paste0('Agony with FDR [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')
ABNF = paste0('Agony with Bonferroni [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')

CMHC = paste0('Confidence without MHC [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')
CFDR = paste0('Confidence with FDR [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')
CBNF = paste0('Confidence with Bonferroni [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')

HC0 = paste0('Hill Climbing [', 0, ' restarts]')
HC = paste0('Hill Climbing [', boot.first.pass + boot.second.pass, ' restarts]')


DEV.OFF = FALSE

source(paste0(git, "utils/plotter.R"))
x1 = plt(results$agony.inference, results$agony.inference.pvalues, "pvalues", adj.matrix, AMHC, stats = FALSE, dev.off = DEV.OFF)
x2 = plt(results$agony.inference, results$agony.inference.pvalues, 'qvalues.fdr', adj.matrix, AFDR, dev.off = DEV.OFF)
x3 = plt(results$agony.inference, results$agony.inference.pvalues, 'qvalues.holm', adj.matrix, ABNF, dev.off = DEV.OFF)

x4 = plt(results$confidence.inference, results$confidence.inference.pvalues, "pvalues", adj.matrix, CMHC, dev.off = DEV.OFF)
x5 = plt(results$confidence.inference, results$confidence.inference.pvalues, 'qvalues.fdr', adj.matrix, CFDR, dev.off = DEV.OFF)
x6 = plt(results$confidence.inference, results$confidence.inference.pvalues, 'qvalues.holm', adj.matrix, CBNF, dev.off = DEV.OFF)

hc1 = plt(results, NULL, 'hill.climing.no.restarts.inference', adj.matrix, HC0, dev.off = DEV.OFF)
hc2 = plt(results, NULL, 'hill.climing.with.restarts', adj.matrix, HC, dev.off = DEV.OFF)


plt.stats(cbind(x1,x2,x3,x4,x5,x6, hc1, hc2), legend.cex = .6)
dev.copy2pdf(file = 'All comparison.pdf')
