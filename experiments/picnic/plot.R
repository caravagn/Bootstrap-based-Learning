setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/picnic/')
git = '/Volumes/Data/Github/Bootstrap-based-Learning/'
agony.binaries = '/Volumes/Data/Github/Bootstrap-based-Learning/agony/giulio/agony'
library(bnlearn)
library(parallel)
library(doParallel)
library(igraph)
library(catnet)
library(TRONCO)
source(paste0(git, "utils/plotter.R"))
source(paste0(git, "utils/stat.R"))
set.seed(12345)
Wnt = c("APC", "CTNNB1", "DKK1", "DKK2", "DKK3", "DKK4", "LRP5", "FZD10", "FAM123B", 
        "AXIN2", "TCF7L2", "FBXW7", "ARID1A", "SOX9")
RAS = c("ERBB2", "ERBB3", "NRAS", "KRAS", "BRAF")
PI3K = c("IGF2", "IRS2", "PIK3CA", "PIK3R1", "PTEN")
TGFb = c("TGFBR1", "TGFBR2", "ACVR1B", "ACVR2A", "SMAD2", "SMAD3", "SMAD4")
P53 = c("TP53", "ATM")
pathway.genes = c(Wnt, RAS, PI3K, TGFb, P53)
pathway.names = c('Wnt', 'RAS', 'PI3K', 'TGFb', 'P53')
pathway.list = list(Wnt = Wnt, RAS = RAS, PI3K = PI3K, TGFb = TGFb, P53 = P53)
pathways.color = c('firebrick1', 'darkblue', 'darkgreen', 'darkmagenta', 'darkorange')


WORKWITH = 'MSI'
if(WORKWITH == 'MSI')
{
	load('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/picnic/MSI.models.Rdata')
	MSS.models = MSI.models
}
if(WORKWITH == 'MSS')
{
	load('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/picnic/MSS.models.Rdata')
}

# Load results from Suppes with both AIC and BIC
load('capri_bic.suppes.Rdata')
bic.suppes = results

load('capri_aic.suppes.Rdata')
aic.suppes = results

# Consider the MSS.models BIC adjacency matrix as the true model
adj.matrix = as.adj.matrix(MSS.models, models = 'capri_bic')[['capri_bic']]

DEV.OFF = TRUE
STATS = TRUE

x = plt(bic.suppes$confidence.inference, 
		bic.suppes $confidence.inference.pvalues, 
		"qvalues.holm", 
		adj.matrix, 
		"SuppesBonferroniBIC", 
		stats = STATS, 
		dev.off = DEV.OFF)

# Consider the MSS.models AIC adjacency matrix as the true model
adj.matrix = as.adj.matrix(MSS.models, models = 'capri_aic')[['capri_aic']]

x = plt(aic.suppes$confidence.inference, 
		aic.suppes$confidence.inference.pvalues, 
		"qvalues.holm", 
		adj.matrix, 
		"SuppesBonferroniAIC", 
		stats = STATS, 
		dev.off = DEV.OFF)

aux = function(x, t) {
	SB.MSS = MSS.models
	SB.MSS$model$capri_bic$adj.matrix$adj.matrix.fit = bic.suppes$confidence.inference[[x]]
	SB.MSS$model$capri_aic$adj.matrix$adj.matrix.fit = aic.suppes$confidence.inference[[x]]

	conf = aic.suppes$confidence.inference.pvalues[[x]]
	colnames(conf) = colnames(SB.MSS$model$capri_aic$adj.matrix$adj.matrix.fit)
	rownames(conf) = rownames(SB.MSS$model$capri_aic$adj.matrix$adj.matrix.fit)
	
	SUPPBONF = SB.MSS 
	save(SUPPBONF, file = 'SuppBONF.Rdata')

	dev.new(width = 4.5,height = 3)
	ggd.qqplot(conf, paste(t, "p-values"))
	dev.copy2pdf(file = paste0(t, '-pvalues.pdf'))
	dev.off()

	SB.MSS$confidence[[1]] = conf

	dev.new(width = 20, height = 20)
	tronco.plot(SB.MSS, 
		 pathways = pathway.list, 
		 title = paste("Suppes with ", x, " - MSS COADREAD (TCGA)"),
		 fontsize = 15,
		 label.edge.size = 13,
		 edge.cex = 3.5,
 		 legend.cex = .5,
 		 confidence = 'tp',
	 	scale.nodes = .6,
		 pathways.color = pathways.color,
		disconnected = T, 
		expand = T,
		height.logic = .3)		
	dev.copy2pdf(file = paste0(t, '-model.pdf'))

	tronco.plot(SB.MSS, 
		 pathways = pathway.list, 
		 title = paste("Suppes with ", x, " - MSS COADREAD (TCGA)"),
		 fontsize = 15,
		 label.edge.size = 13,
		 edge.cex = 3.5,
 		 legend.cex = .5,
 		 confidence = 'tp',
	 	scale.nodes = .6,
		 pathways.color = pathways.color,
		disconnected = T, 
		expand = FALSE,
		height.logic = .3)		
	dev.copy2pdf(file = paste0(t, '-model-compressed.pdf'))
	dev.off()

}

aux("pvalues", 'SuppNoMHC')
aux("qvalues.fdr", 'SuppFDR')
aux("qvalues.holm", 'SuppBONF')

# Load MSS obj (TRONCO)
load('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/picnic/SUPPBONF.Rdata')

par(mfrow = c(2,1))
tronco.plot(MSS.models, 
	 pathways = pathway.list, 
	 title = paste("CAPRI - MSS COADREAD (TCGA)"),
	 models = 'capri_bic',
	 fontsize = 15,
	 label.edge.size = 13,
	 edge.cex = 3.5,
	 legend.cex = .5,
	 confidence = 'tp',
	 scale.nodes = .6,
	 pathways.color = pathways.color,
	 disconnected = F, 
	 expand = T,
	 height.logic = .3)		


tronco.plot(SUPPBONF, 
	 pathways = pathway.list, 
	 title = paste("Suppes with Bonferroni BIC - MSS COADREAD (TCGA)"),
	 models = 'capri_bic',
	 fontsize = 15,
	 label.edge.size = 13,
	 edge.cex = 3.5,
	 legend.cex = .5,
	 confidence = 'tp',
	 scale.nodes = .6,
	 pathways.color = pathways.color,
	 disconnected = T, 
	 expand = T,
	 height.logic = .3)		


tronco.plot(SUPPBONF, models = 'capri_bic')
tronco.plot(SUPPBONF, models = 'capri_aic')

SUPPBONF = tronco.kfold.eloss(SUPPBONF, runs = 100)
SUPPBONF = tronco.kfold.prederr(SUPPBONF, runs = 100)


as.kfold.eloss(SUPPBONF)
as.kfold.eloss(MSS.models)



as.selective.advantage.relations(SUPPBONF)

as.kfold.prederr(MSS.models, models = 'capri_bic', events = as.events(MSS.models, genes = 'KRAS'))
24                     Mutation KRAS  0.467105263 0.000000000
23                     Mutation NRAS  0.085526316 0.000000000
as.kfold.prederr(SUPPBONF, models = 'capri_bic', events = as.events(SUPPBONF, genes = 'KRAS'))
24                     Mutation KRAS  0.534868421 0.036570900
23                     Mutation NRAS  0.085526316 0.000000000

as.kfold.posterr(SUPPBONF, models = 'capri_bic', values = T)


tabular = function(obj, M){
  tab = Reduce(
    function(...) merge(..., all = T), 
      list(
      as.selective.advantage.relations(obj)[[M]],
      # as.bootstrap.scores(obj)[[M]],
      as.kfold.prederr(obj)[[M]],
      as.kfold.posterr(obj)[[M]]
    )
  )
  # merge reverses first with second column
  tab = tab[, c(2,1,3:ncol(tab))]
  tab = tab[order(tab$MEAN.PREDERR, na.last = TRUE, decreasing = TRUE), ]
  
  return(tab)
}

MSS.bic = tabular(MSS.models, 'capri_bic')
SUPPBONF.bic = tabular(SUPPBONF, 'capri_bic')

M = merge(MSS.bic, SUPPBONF.bic, by = 'SELECTED')
BETTER = WORST = EQ = NULL
for(i in 1:nrow(M)) {
	if(M[i, 'MEAN.PREDERR.x'] > M[i, 'MEAN.PREDERR.y']) 
		BETTER = rbind(BETTER, M[i, c('SELECTED', 'SELECTS.x', 'SELECTS.y', 'MEAN.PREDERR.x', 'MEAN.PREDERR.y', 'SD.PREDERR.x', 'SD.PREDERR.y')])
	if(M[i, 'MEAN.PREDERR.x'] == M[i, 'MEAN.PREDERR.y']) 
		EQ = rbind(EQ, M[i, c('SELECTED', 'SELECTS.x', 'SELECTS.y', 'MEAN.PREDERR.x', 'MEAN.PREDERR.y', 'SD.PREDERR.x', 'SD.PREDERR.y')])
	else
		WORST = rbind(WORST, M[i, c('SELECTED', 'SELECTS.x', 'SELECTS.y', 'MEAN.PREDERR.x', 'MEAN.PREDERR.y', 'SD.PREDERR.x', 'SD.PREDERR.y')])
}
BETTER
WORST
EQ


kfoldlogLik = function(g = 10, x = 'capri_bic')
{
	data = as.genotypes(SUPPBONF)
	data <- as.data.frame(as.genotypes(SUPPBONF), stringsAsFactors = TRUE) # the "[]" keeps the dataframe structure
	data[] <- lapply(data, factor) # the "[]" keeps the dataframe structure
	# str(data)

	test = data[sample(1:nrow(data), size = nrow(data)/g), ]
	training = data[!(rownames(data) %in% rownames(test)),]

	bn.SUPPBONF = empty.graph(colnames(data))
	amat(bn.SUPPBONF) = as.adj.matrix(SUPPBONF, models = x)[[x]]

	bn.trte.SUPPBONF = bn.fit(bn.SUPPBONF, data)
	bn.tr.SUPPBONF = bn.fit(bn.SUPPBONF, training)
	
	bn.CAPRI = empty.graph(colnames(data))
	amat(bn.CAPRI) = as.adj.matrix(MSS.models, models = x)[[x]]

	bn.trte.CAPRI = bn.fit(bn.CAPRI, data)
	bn.tr.CAPRI = bn.fit(bn.CAPRI, training)

	# BIC(bn.trte.SUPPBONF, data = data)
	# BIC(bn.trte.CAPRI, data = data)
	
	if(x == 'capri_bic'){
		sspb = - BIC(bn.tr.SUPPBONF, data = test)
		scp = - BIC(bn.tr.CAPRI, data = test)	
	}

	if(x == 'capri_aic'){
		sspb = - AIC(bn.tr.SUPPBONF, data = test)
		scp = - AIC(bn.tr.CAPRI, data = test)	
	}

	sspb = - logLik(bn.tr.SUPPBONF, data = test)
	scp = - logLik(bn.tr.CAPRI, data = test)	
	
	cat('SBONF', sspb, 'CAPRI', scp, ' [', sspb - scp, ']\n')
	
	return(sspb - scp)	
}

points = sapply(1:100, kfoldlogLik)
boxplot((points), ylim = c(-10, 10))

points = sapply(1:100, kfoldlogLik, x = 'capri_aic')
points = points[!is.na(points)]
points = points[points !=  -Inf]

boxplot((points))

tabular(MSS.models, 'capri_aic')
tabular(MSI.models, 'capri_bic')
tabular(MSI.models, 'capri_bic')

