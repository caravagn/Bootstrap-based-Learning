load('capri_bic.suppes.Rdata')
bic.suppes = results

load('capri_bic.others.Rdata')
bic.others = results.others

load('capri_aic.suppes.Rdata')
aic.suppes = results

load('capri_aic.others.Rdata')
aic.others = results.others

adj.matrix = MSS.models$model$capri_bic$adj.matrix$adj.matrix.fit

SMHC = paste0('Suppes without MHC [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')
SFDR = paste0('Suppes with FDR [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')
SBNF = paste0('Suppes with Bonferroni [', boot.first.pass, '+', boot.second.pass, ' npb, p ', test.pvalue, ']')

DEV.OFF = TRUE
STATS = TRUE

x4 = plt(results$confidence.inference, results$confidence.inference.pvalues, "pvalues", adj.matrix, SMHC, stats = STATS, dev.off = DEV.OFF)
x5 = plt(results$confidence.inference, results$confidence.inference.pvalues, 'qvalues.fdr', adj.matrix, SFDR, stats = STATS, dev.off = DEV.OFF)
x6 = plt(results$confidence.inference, results$confidence.inference.pvalues, 'qvalues.holm', adj.matrix, SBNF, stats = STATS, dev.off = DEV.OFF)


SB.MSS = MSS.models
SB.MSS$model$capri_bic$adj.matrix$adj.matrix.fit = bic.others $confidence.inference$qvalues.holm
SB.MSS$model$capri_aic$adj.matrix$adj.matrix.fit = aic.suppes$confidence.inference$qvalues.holm

conf = aic.suppes$confidence.inference.pvalues$qvalues.holm
colnames(conf) = colnames(SB.MSS$model$capri_aic$adj.matrix$adj.matrix.fit)
rownames(conf) = rownames(SB.MSS$model$capri_aic$adj.matrix$adj.matrix.fit)

dev.new(width = 4.5,height = 3)
ggd.qqplot(conf, "P-values")
dev.copy2pdf(file = 'SB-pval.pdf')
dev.off()

SB.MSS$confidence[[1]] = conf

dev.new(width = 20, height = 20)
tronco.plot(SB.MSS, 
	 pathways = pathway.list, 
	 title = "Suppes with Bonferroni - MSS COADREAD (TCGA)",
	 # models = "capri_bic",
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
dev.copy2pdf(file = 'SB.pdf')

AB.MSS = MSS.models
AB.MSS $model$capri_bic$adj.matrix$adj.matrix.fit = bic.others$agony.inference$qvalues.holm
AB.MSS $model$capri_aic$adj.matrix$adj.matrix.fit = aic.others$agony.inference $qvalues.holm

conf = aic.others $confidence.inference.pvalues$qvalues.holm
colnames(conf) = colnames(SB.MSS$model$capri_aic$adj.matrix$adj.matrix.fit)
rownames(conf) = rownames(SB.MSS$model$capri_aic$adj.matrix$adj.matrix.fit)

dev.new(width = 4.5,height = 3)
ggd.qqplot(conf, "P-values")
dev.copy2pdf(file = 'AB-pval.pdf')
dev.off()

AB.MSS$confidence[[1]] = conf

tronco.plot(AB.MSS, 
	 pathways = pathway.list, 
	 title = "Agony with Bonferroni - MSS COADREAD (TCGA)",
	 # models = "capri_bic",
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
dev.copy2pdf(file = 'AB.pdf')

