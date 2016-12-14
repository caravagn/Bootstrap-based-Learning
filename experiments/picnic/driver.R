
setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/picnic/')
git = '/Volumes/Data/Github/Bootstrap-based-Learning/'
agony.binaries = '/Volumes/Data/Github/Bootstrap-based-Learning/agony/giulio/agony'

# PiCnIc processing of TCGA-COADREAD data (cancer)
source('main.R', echo = TRUE)

DOPLOTS = FALSE

# Snapshot some  output
if(DOPLOTS) {
	
tronco.plot(MSI.models, 
	 pathways = pathway.list, 
	 models = "capri_bic",
	 fontsize = 15,
	 edge.cex = 1.5,
 	 legend.cex = .5,
	 scale.nodes = .6,
	 confidence = c('tp', 'pr', 'hg'), # Display p-values 
	 pathways.color = pathways.color,
	 label.edge.size = 9,
	disconnected = F, 
	height.logic = .3)
dev.copy2pdf(file='MSI-CAPRI.pdf')

tronco.plot(MSS.models, 
	 pathways = pathway.list, 
	 models = "capri_bic",
	 fontsize = 15,
	 edge.cex = 1.5,
 	 legend.cex = .5,
	 scale.nodes = .6,
	 confidence = c('tp', 'pr', 'hg'), # Display p-values 
	 pathways.color = pathways.color,
	 label.edge.size = 9,
	disconnected = F, 
	height.logic = .3)
dev.copy2pdf(file='MSS-CAPRI.pdf')
dev.off()

dev.new(width=20, height = 20)
tronco.plot(MSS.models, 
	 pathways = pathway.list, 
	 title = "",
	 models = "capri_bic",
	 # fontsize = 15,
	 pf = TRUE,
	 # edge.cex = 1.5,
	 legend = F,
	 scale.nodes = .6,
	 pathways.color = pathways.color,
	disconnected = F, 
	expand = F)
dev.copy2pdf(file='MSS-CAPRI.pdf')
}

# Prepare data -- genotypes
dataset = keysToNames(MSS.models, as.genotypes(MSS.models))
dataset = as.genotypes(MSS.models)
str(dataset)

dataset <- as.data.frame(dataset, stringsAsFactors = TRUE) # the "[]" keeps the dataframe structure
dataset[] <- lapply( dataset, factor) # the "[]" keeps the dataframe structure
str(dataset)

# load the required R packages
library(bnlearn)
library(parallel)
library(doParallel)
library(igraph)
library(catnet)

# source the required scripts
source(paste0(git, "utils/plotter.R"))
source(paste0(git, "utils/stat.R"))

# set the seed
set.seed(12345)

# set some settings to be used in the test
boot.first.pass = 100
boot.second.pass = 100
test.pvalue = 1E-2


export.primafacie = function(TCGA, model)
{
	if(model == 'capri_bic') regularization = "bic"
	if(model == 'capri_aic') regularization = "aic"

	suppes = as.adj.matrix(TCGA, models = model, type = "pf")
	suppes = suppes[[model]]

	fit.adj = as.adj.matrix(TCGA, models = model)
	fit.adj = keysToNames(TCGA, fit.adj[[model]])


	#### Suppes
	source(paste0(git, "src/suppes.bootstrap.inference.R"))

	results = suppes.bootstrap.inference(
		dataset,
		regularization,
		suppes,
		boot.second.pass,
		test.pvalue,
		agony.binaries = agony.binaries,
		nboot.first = boot.first.pass)

	save(results, file = paste0(model, '.suppes.Rdata'))

	#### Ago, Conf
	source(paste0(git, "src/perform.bootstrap.inference.R"))
	regularization = "aic"

	# results.others = perform.bootstrap.inference(
		# dataset,
		# regularization,
		# boot.first.pass,
		# boot.second.pass,
		# test.pvalue,
		# agony.binaries = agony.binaries
	# )

	# save(results.others, file = paste0(model, '.others.Rdata'))
}

export.primafacie(MSS.models, 'capri_bic')
export.primafacie(MSS.models, 'capri_aic')

# source('plot.R', echo = TRUE)
