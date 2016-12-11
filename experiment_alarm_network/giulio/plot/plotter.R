
library(Rgraphviz)
library(bnlearn)
library(RColorBrewer)


plt = function(allresults, algorithm, adj.true, nboot, p, title) {
	par(mfrow = c(1, 2))
	

	X = adj.true
	Xgraph = new("graphAM", adjMat = X, edgemode = "directed")

	Y = allresults[[algorithm]]
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
    print(t(as.data.frame(stat)))

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
		paste('pr.', round(stat$precision, 3)), 
		paste('re.', round(stat$sensitivity, 3))), 
		cex=0.8, 
		bty='n')

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

