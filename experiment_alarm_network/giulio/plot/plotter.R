library(Rgraphviz)
library(bnlearn)
library(RColorBrewer)


plt = function(results, pvalues, algorithm, adj.true, nboot, p, title, stats = TRUE, save.images = TRUE) {
		
	# ratio = c(5,5,1)	
	# L1 = matrix(1, nrow=ratio[1], ncol = ratio[2])
	# L2 = matrix(2, nrow=ratio[1], ncol = ratio[2])
	# L3 = matrix(3, nrow=ratio[3], ncol = ratio[2])
	# L4 = matrix(4, nrow=ratio[3], ncol = ratio[2])

	# L = cbind(L1,L2)
	# L = rbind(L,cbind(L3,L4))
	# layout(L)
	par(mfrow = c(1,2))
	
	X = adj.true
	Xgraph = new("graphAM", adjMat = X, edgemode = "directed")

	Y = results[[algorithm]]
	colnames(Y) = colnames(adj.true)
	rownames(Y) = colnames(adj.true)
	Ygraph = new("graphAM", adjMat = Y, edgemode = "directed")

	colors = c('firebrick2', 'orange', 'dodgerblue3', 'lightblue')
	add.alpha <- function(col, alpha=1){
	apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
	}
	colors[4] = add.alpha(colors[4], .3)

    source('plot/stat.R')
    classes = classifyedges(X, Y)
	print(classes)

	edgesX = edgesY = NULL
	colorX = c(
		rep(colors[1], classes$nreversed), 
		rep(colors[2], classes$nmadeup), 
		rep(colors[3], classes$nmissing))
	names(colorX) = c(
		unlist(classes$reversedX), 
		unlist(classes$madeup), 
		unlist(classes$missing))	
		
	colorY = c(
		rep(colors[1], classes$nreversed), 
		rep(colors[2], classes$nmadeup), 
		rep(colors[3], classes$nmissing))
	names(colorY) = c(
		unlist(classes$reversedY), 
		unlist(classes$madeup), 
		unlist(classes$missing))	


	edgesX$color = colorX	
	edgesY$color = colorY	

	plot(Xgraph, attrs = list(node = list(fillcolor = colors[4])), edgeAttrs = edgesX)

	stat = stats(X, Y)
	print(stat)

	legend("topleft", title = "True model",
		c(
			paste('', classes$nreversed), 
			paste('', classes$nmissing), 
			paste('', classes$nmadeup), 
			paste('', stat$tp)), 
		cex=0.8, 
		col=c(colors[1],colors[3],colors[2], 'black'), 
		lty=1, 
		lwd = 2,
		bty='n')

	plot(Ygraph, attrs = list(node = list(fillcolor = colors[4])), edgeAttrs = edgesY)

	legend("topleft", title = title, 
		c(
		paste('k', nboot), 
		paste('p <', p), 
		paste(''),
		paste('pr.', round(stat$PPV, 3)), 
		paste('re.', round(stat$TPR, 3))), 
		cex=0.8, 
		bty='n')

	if(save.images) dev.copy2pdf(file=paste(nboot, p, title, 'networks.pdf', sep ='_'))

	if(stats)
	{
		dev.new(width = 8, height = 3.5)
		par(mfrow = c(1,2))

		# library(gridExtra)
		# library(gridBase)

		# vps <- baseViewports()
		# pushViewport(vps$figure)
		# vp1 <-plotViewport()
	
		# stat = lapply(stat, round, digits = 3)
		# stat = t(stat)
		stat = t(as.data.frame(stat))
	
		# mytheme <- ttheme_default(
    		# core = list(fg_params=list(cex = 0.5)),
    		# colhead = list(fg_params=list(cex = 0.5)),
    		# rowhead = list(fg_params=list(cex = 0.5)))

		# g <- tableGrob(df, theme = mytheme)
		# grid.draw(g)
		# popViewport()
		
		dfbp = stat[5:length(stat)]
		# print(dfbp)
		colors = colorRampPalette(brewer.pal(6, 'Dark2'))(length(dfbp))
		
		# print(colors)
		barplot(dfbp, 
			main = "Scores",
			# col = "cornflowerblue",
			col =colors,
			names.arg= rownames(stat)[5:length(stat)],
			# cex.names = .7,
			las = 1,
			horiz = TRUE) 

		# print(is.null(pvalues))
		if(!is.null(pvalues)) ggd.qqplot(pvalues[[algorithm]], "P-values")
		if(save.images) dev.copy2pdf(file=paste(nboot, p, title, 'stats.pdf', sep ='_'))

		colnames(stat) = title
		return(stat)
	}
	
	
}



ggd.qqplot =function(pvector, main=NULL, ...) {
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        xlim=c(0,max(e)), ylim=c(0,max(o)))
    lines(e,e,col="red")
}
		

# }



# load("100-alarm-agony.Rdata")
# plt(agony_model, "pvalues", adj.true, 100, 0.01, 'Agony without MHC')
# dev.copy2pdf(file='../100_agony_pvalues.pdf')
# plt(agony_model, 'qvalues.fdr', adj.true, 100, 0.01, 'Agony with FDR')
# dev.copy2pdf(file='../100_agony_fdr.pdf')
# plt(agony_model, 'qvalues.holm', adj.true, 100, 0.01, 'Agony with Bonferroni')
# dev.copy2pdf(file='../100_agony_bonferroni.pdf')

# load("100-alarm-confidence.Rdata")
# plt(confidence_model, "pvalues", adj.true, 100, 0.01, "Confidence without MHC")
# dev.copy2pdf(file = "../100_confidence_pvalues.pdf")
# plt(confidence_model, "qvalues.fdr", adj.true, 100, 0.01, "Confidence with FDR")
# dev.copy2pdf(file = "../100_confidence_fdr.pdf")
# plt(confidence_model, "qvalues.holm", adj.true, 100, 0.01, "Confidence with Bonferroni")
# dev.copy2pdf(file = "../100_confidence_bonferroni.pdf")


# hc0 = hc(data, restart = 0, score = reg)
# hc100 = hc(data, restart = 200, , score = reg)

# h = NULL
# h$hc0 = amat(hc0)
# h$hc100 = amat(hc100)

# plt(h, 'hc0', adj.true, 0, 51, 'Hill Climbing')
# dev.copy2pdf(file='../0_HC.pdf')
# plt(h, 'hc100', adj.true, 200, 51,  'Hill Climbing')
# dev.copy2pdf(file='../100_HC.pdf')
# # plt(h, 'hc1000', 1100, 51, 'Hill Climbing')
# # dev.copy2pdf(file='1000_HC.pdf')

