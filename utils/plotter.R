library(Rgraphviz)
library(bnlearn)
library(RColorBrewer)


# plt = function(results, pvalues, algorithm, adj.true, nboot, p, title, stats = TRUE, save.images = TRUE) {
plt = function(results, pvalues, algorithm, adj.true, title, stats = TRUE, save.images = TRUE, dev.off = FALSE) {
				
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

    classes = classifyedges(X, Y)
	# print(classes)

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
		legend = "",
		cex=0.8, 
		bty='n')

	legend("bottom", title = paste('PPV', round(stat$PPV, 3), 'TPR', round(stat$TPR, 3)), 
		legend = "",
		cex=0.8, 
		bty='n')
	
	if(save.images && !stats) dev.copy2pdf(file=paste(gsub(' ', '_',title), '.pdf', sep =''))
	if(save.images && stats) dev.copy2pdf(file=paste('tmp1.pdf', sep =''))
	if(dev.off) dev.off()
	
	
	if(stats)
	{
		dev.new(width = 8, height = 3.5)
		par(mfrow = c(1,2))
		stat = t(as.data.frame(stat))
		
		dfbp = stat[5:length(stat)]
		colors = colorRampPalette(brewer.pal(6, 'Dark2'))(length(dfbp))
		
		barplot(dfbp, 
			main = "Scores",
			# col = "cornflowerblue",
			col = colors,
			names.arg= rownames(stat)[5:length(stat)],
			# cex.names = .7,
			las = 1,
			horiz = TRUE) 

		if(!is.null(pvalues)) ggd.qqplot(pvalues[[algorithm]], "P-values")
		# if(save.images) dev.copy2pdf(file=paste('[Stats] ', title, '.pdf', sep =''))
	
		if(save.images) dev.copy2pdf(file=paste('tmp2.pdf', sep =''))

		# print(gsub(' ', '_',title))
		mergePDF(
			'tmp1.pdf',
			'tmp2.pdf',
			file = paste(gsub(' ', '_',title), '.pdf', sep ='')
			)
		file.remove('tmp1.pdf', 'tmp2.pdf')
   	    if(dev.off) dev.off()

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
		


# FROM https://github.com/trinker/plotflow/blob/master/R/mergePDF.R
mergePDF = function(..., file, gsversion = NULL, in.file = NULL) {
    if (is.null(in.file)) in.file <- substitute(...())

    infiles <- paste(unlist(lapply(in.file, function(y) as.character(y))), collapse = " ")

    if (is.null(gsversion)) {

      gsversion <- names(which(Sys.which(c("gs", "gswin32c", "gswin64c")) != ""))

      if (length(gsversion) == 0) 

        stop("Please install Ghostscript and ensure it is in your PATH")

      if (length(gsversion) > 1)

        stop("More than one Ghostscript executable was found:", 

             paste(gsversion, collapse = " "), 

             ". Please specify which version should be used with the gsversion argument")

    }   

    pre = " -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="

    system(paste(paste(gsversion, pre, file, sep = ""), infiles, collapse = " "))

}



plt.stats = function(all, palette = 'Paired', cols = 8, legend.cex = .5)
{
	k = nrow(all)
	par(mfrow = c(sqrt(k)+1, sqrt(k)+1),
	    oma = c(5,4,2,2) + 0.4,
        mar = c(2,2,2,2) + 0.4)
        
   # m = matrix(c(1:4,5:8,9:12, 13:14, 0, 0, rep(15,4)), ncol = 4, byrow = TRUE)
	# layout(mat = m, heights = c(rep(0.4, nrow(m) - 1),0.2))

	
	colors = colorRampPalette(brewer.pal(cols, palette))(ncol(all))
	names(colors) = colnames(all)
	
	
	# print(colnames(all))
	# a.colors = colorRampPalette(brewer.pal(cols, palette))(ncol(all))
	
	scores = rownames(all)
	algos = colnames(all)
	for(i in 1:nrow(all))
	{
		# print(i)
		data = sort(all[i, , drop = T])
				
		# print(colors)
		# print(colors)

		names = rep("", length(algos))
		# if(i == 1) names = colnames(data)
		
		barplot(data, 
			main = scores[i],
			col = colors[names(data)],
			names.arg = names,
			# xlim = c(0,1),
			cex.names = .7,
			las = 1,
			horiz = TRUE) 
	}
	
	plot(1, type="n", axes=FALSE, xlab="", ylab="")

	legend("center", 
		title = "Algorithms",
		algos,
		cex = legend.cex, 
		col = colors, 
		pch = 15,
		# lty=1, 
		# lwd = 2,
		bty='n')

}


