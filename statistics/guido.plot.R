
guido.plot = function(true.model, consensus.model, standard.model, X.scores, Y.scores, labels = c('NPB','MI'), loop = NULL)
{
	if(!is.matrix(true.model)) true.model = as.matrix(true_adj_matrix)
	if(!is.matrix(X.scores)) X.scores = as.matrix(X.scores)
	if(!is.matrix(Y.scores)) Y.scores = as.matrix(Y.scores)
	if(!is.matrix(consensus.model)) consensus.model = as.matrix(consensus.model)
	if(!is.matrix(standard.model)) standard.model = as.matrix(standard.model)
	
	ordering = order(X.scores)

	#plot.new()
	# par(mfrow = c(2,1))
	
	# X.scores and Y.scores
	data = NULL	
	for(i in 1:length(ordering))
		data = rbind(data, c(i, X.scores[ordering[i]], Y.scores[ordering[i]]))
	
	colnames(data) = c('X', 'X.score', 'Y.score')
	
	# ggplot2 data format	
	dframe = data.frame(row.names = 1:length(ordering))
	dframe$X = as.numeric(data[, 1])
	dframe$X.score = as.numeric(data[, 2])
	dframe$Y.score = as.numeric(data[, 3])

	dframe$true = rep('spurious', length(ordering))
	dframe$model = rep('', length(ordering))

	# edges as points with special mark
	for(i in 1:length(ordering))
	{
		if(true.model[ordering[i]] == 1)
			dframe[i, 'true'] = 'true'
		
		if(consensus.model[ordering[i]] == 1 && standard.model[ordering[i]] == 1)
			dframe[i, 'model'] = '*'
			
		if(consensus.model[ordering[i]] == 1 && standard.model[ordering[i]] == 0)
			dframe[i, 'model'] = 'C'

		if(consensus.model[ordering[i]] == 0 && standard.model[ordering[i]] == 1)
			dframe[i, 'model'] = 'S'
		
	}
				
	library(ggplot2)
	library(gridExtra)

    mycolours = c("true" = "brown3", "spurious" = "black")
    myshape = c("yes" = "a", "no" = "b")

	top.frame = dframe[, c('X', 'X.score', 'true', 'model')]
	mid.frame = dframe[, c('X', 'Y.score', 'true', 'model')]
	
	p.top = ggplot(top.frame, aes(x = X, y = X.score, fill = labels[1])) +
	ggtitle(paste(labels[1],'score'))+
	xlab(paste('X-axis ordered according to the ', labels[1],'score (*: both, C: consensus, S: standard)')) +guides(fill=FALSE) +
    geom_point(size = 2,  aes(colour = true)) +
    scale_color_manual("Edge", values = mycolours) +
    geom_text(data = top.frame, size = 3, aes(x = X + 0.05, y = X.score + 5, label = top.frame$model)) 

    # scale_color_manual("Consensus", values = myshape) 
     
     # p.top = ggplot() +
	# # scale_x_discrete(name="") +
	# # scale_y_continuous(limits=c(0,1), breaks=NA, name="") +
	# # scale_shape_discrete(solid=T, legend=F) +
	# geom_point(data= top.frame, mapping=aes(x=X, y=X.score, shape= myshape), size=10)



	p.mid = ggplot(mid.frame, aes(x = X, y = Y.score, fill = labels[2])) +
	xlab(paste('X-axis ordered as in the top panel (', labels[1], 'annotated)')) +
	guides(fill=FALSE) +
	ggtitle(paste(labels[2],'score'))+
    geom_point(size = 2, aes(colour = true)) +
    scale_color_manual("Edge", values = mycolours) +
    geom_text(data = mid.frame, size = 3, aes(x = X + 0.05, y = Y.score + 0.05, label = top.frame$X.score)) 

	grid.arrange(p.top, p.mid)
	
}


guido.wblist.plot = function(true.model, consensus.model, consensus.w.model, consensus.b.model,  standard.model, X.scores, Y.scores, labels = c('NPB','MI'), loop = NULL)
{
	if(!is.matrix(true.model)) true.model = as.matrix(true_adj_matrix)
	if(!is.matrix(X.scores)) X.scores = as.matrix(X.scores)
	if(!is.matrix(Y.scores)) Y.scores = as.matrix(Y.scores)
	if(!is.matrix(consensus.model)) consensus.model = as.matrix(consensus.model)
	if(!is.matrix(standard.model)) standard.model = as.matrix(standard.model)
	if(!is.matrix(consensus.w.model)) consensus.w.model = as.matrix(consensus.w.model)
	if(!is.matrix(consensus.b.model)) consensus.b.model = as.matrix(consensus.b.model)
	
	ordering = order(X.scores)

	#plot.new()
	# par(mfrow = c(2,1))
	
	# X.scores and Y.scores
	data = NULL	
	for(i in 1:length(ordering))
		data = rbind(data, c(i, X.scores[ordering[i]], Y.scores[ordering[i]]))
	
	colnames(data) = c('X', 'X.score', 'Y.score')
	
	# ggplot2 data format	
	dframe = data.frame(row.names = 1:length(ordering))
	dframe$X = as.numeric(data[, 1])
	dframe$X.score = as.numeric(data[, 2])
	dframe$Y.score = as.numeric(data[, 3])

	dframe$true = rep('spurious', length(ordering))
	dframe$model = rep('', length(ordering))

	# edges as points with special mark
	for(i in 1:length(ordering))
	{
		if(true.model[ordering[i]] == 1)
			dframe[i, 'true'] = 'true'
		
		if(consensus.model[ordering[i]] == 1 && standard.model[ordering[i]] == 1)
			dframe[i, 'model'] = '*'
			
		if(consensus.model[ordering[i]] == 1 && standard.model[ordering[i]] == 0)
			dframe[i, 'model'] = 'C'

		if(consensus.model[ordering[i]] == 0 && standard.model[ordering[i]] == 1)
			dframe[i, 'model'] = 'S'
			
		if( consensus.w.model[ordering[i]] == 1 )	
			dframe[i, 'model'] = paste(dframe[i, 'model'] , ' w', sep='')
	
		if( consensus.b.model[ordering[i]] == 1 )	
			dframe[i, 'model'] = paste(dframe[i, 'model'] , ' b', sep='')
	}
				
	library(ggplot2)
	library(gridExtra)

    mycolours = c("true" = "brown3", "spurious" = "black")
    myshape = c("yes" = "a", "no" = "b")

	top.frame = dframe[, c('X', 'X.score', 'true', 'model')]
	mid.frame = dframe[, c('X', 'Y.score', 'true', 'model')]
	
	p.top = ggplot(top.frame, aes(x = X, y = X.score, fill = labels[1])) +
	ggtitle(paste(labels[1],'score'))+
	xlab(paste('X-axis ordered according to the ', labels[1],'score (*: both, C: consensus, S: standard)')) +guides(fill=FALSE) +
    geom_point(size = 2,  aes(colour = true)) +
    scale_color_manual("Edge", values = mycolours) +
    geom_text(data = top.frame, size = 2, aes(x = X + 0.05, y = X.score + 5, label = top.frame$model)) 

    # scale_color_manual("Consensus", values = myshape) 
     
     # p.top = ggplot() +
	# # scale_x_discrete(name="") +
	# # scale_y_continuous(limits=c(0,1), breaks=NA, name="") +
	# # scale_shape_discrete(solid=T, legend=F) +
	# geom_point(data= top.frame, mapping=aes(x=X, y=X.score, shape= myshape), size=10)



	p.mid = ggplot(mid.frame, aes(x = X, y = Y.score, fill = labels[2])) +
	xlab(paste('X-axis ordered as in the top panel (', labels[1], 'annotated)')) +
	guides(fill=FALSE) +
	ggtitle(paste(labels[2],'score'))+
    geom_point(size = 2, aes(colour = true)) +
    scale_color_manual("Edge", values = mycolours) +
    geom_text(data = mid.frame, size = 3, aes(x = X + 0.05, y = Y.score + 0.05, label = top.frame$X.score)) 

	grid.arrange(p.top, p.mid)
	
}



guido.graph.plot = function(true.model, consensus.model, standard.model, X.scores, Y.scores, labels = c('NPB','MI'), loop = NULL)
{
	if(!is.matrix(true.model)) true.model = as.matrix(true_adj_matrix)
	if(!is.matrix(X.scores)) X.scores = as.matrix(X.scores)
	if(!is.matrix(Y.scores)) Y.scores = as.matrix(Y.scores)
	if(!is.matrix(consensus.model)) consensus.model = as.matrix(consensus.model)
	if(!is.matrix(standard.model)) standard.model = as.matrix(standard.model)
	
	colnames(true.model) = 1:ncol(true.model)
	rownames(true.model) = 1:ncol(true.model)
	colnames(consensus.model) = 1:ncol(true.model)
	rownames(consensus.model) = 1:ncol(true.model)
	# colnames(true.model) = 1:ncol(true.model)
	# colnames(true.model) = 1:ncol(true.model)
	# colnames(true.model) = 1:ncol(true.model)
	# colnames(true.model) = 1:ncol(true.model)
	# colnames(true.model) = 1:ncol(true.model)

	library(Rgraphviz)
	D.true.model <-new("graphAM", adjMat= true.model, edgemode="directed")
	D.consensus.model <-new("graphAM", adjMat= consensus.model, edgemode="directed")

	# par(mfrow = c(2,1))
		

	npb.as.label = function(m)
	{
		r = NULL
		for(i in 1:nrow(m))
			for(j in 1:ncol(m))
				{
					if(consensus.model[i, j] == 1 && standard.model[i,j] == 1)
									r = c(r, paste(m[i,j], '(*)'))
					if(consensus.model[i, j] == 1 && standard.model[i,j] == 0)
									r = c(r, paste(m[i,j], '(C)'))
					if(consensus.model[i, j] == 0 && standard.model[i,j] == 1)
									r = c(r, paste(m[i,j], '(S)'))
					if(consensus.model[i, j] == 0 && standard.model[i,j] == 0)
									r = c(r, paste(m[i,j]))
									
					names(r)[length(r)] = paste(rownames(m)[i], '~' ,  colnames(m)[j], sep = '')
				}
				
		return(r)
	}
	
	cns.as.color = function(m)
	{
		cols = NULL
		 # matrix('black', nrow = nrow(true.model), ncol = ncol(true.model))
		# colnames(cols) = colnames(true.model)
		# rownames(cols) = rownames(true.model) 

		for(i in 1:nrow(m))
			for(j in 1:ncol(m))
				{
					if(consensus.model[i,j] == 1)
					{
															cols = c(cols, paste('red'))
									
						names(cols)[length(cols)] = paste(rownames(m)[i], '~' ,  colnames(m)[j], sep = '')
					}
								}
				
		return(cols)
	}
	
	
	npb.as.lwd = function(m)
	{
		print(m)
		r = NULL
		for(i in 1:nrow(m))
			for(j in 1:ncol(m))
				{
					r = c(r, as.numeric(paste(m[i,j])))
					names(r)[length(r)] = paste(rownames(m)[i], '~' ,  colnames(m)[j], sep = '')
				}

				
		return(r)
	}


	attrs <- list(node=list(shape="ellipse", fixedsize=FALSE))


	eAttrs  = list()
	eAttrs$label <- npb.as.label(X.scores)
	eAttrs$color <- cns.as.color(X.scores)
	# eAttrs$lwd <-  npb.as.lwd(X.scores)

print(eAttrs)
	plot(D.true.model, attrs = list(node = list(fillcolor = "lightblue"),
                                edge = list(arrowsize=0.5, fontsize = 4)), edgeAttrs=eAttrs)


	# plot(D.consensus.model, attrs = list(node = list(fillcolor = "lightblue"),
                                # edge = list(arrowsize=0.5)))

	
}

guido.graphdense.plot = function(true.model, consensus.model, standard.model, X.scores, Y.scores, labels = c('NPB','MI'), loop = NULL)
{
	if(!is.matrix(true.model)) true.model = as.matrix(true_adj_matrix)
	if(!is.matrix(X.scores)) X.scores = as.matrix(X.scores)
	if(!is.matrix(Y.scores)) Y.scores = as.matrix(Y.scores)
	if(!is.matrix(consensus.model)) consensus.model = as.matrix(consensus.model)
	if(!is.matrix(standard.model)) standard.model = as.matrix(standard.model)
	
	colnames(true.model) = 1:ncol(true.model)
	rownames(true.model) = 1:ncol(true.model)
	colnames(consensus.model) = 1:ncol(true.model)
	rownames(consensus.model) = 1:ncol(true.model)
	# colnames(true.model) = 1:ncol(true.model)
	# colnames(true.model) = 1:ncol(true.model)
	# colnames(true.model) = 1:ncol(true.model)
	# colnames(true.model) = 1:ncol(true.model)
	# colnames(true.model) = 1:ncol(true.model)

	dense.model = diag(0, nrow = nrow(true.model), ncol = ncol(true.model))
	colnames(dense.model) = 1:ncol(true.model)
	rownames(dense.model) = 1:ncol(true.model)
	for(i in 1:nrow(dense.model))
		for(j in 1:ncol(dense.model))
			if(i!=j && X.scores[i,j] > 0) dense.model[i,j] = 1


	library(Rgraphviz)
	D.true.model <-new("graphAM", adjMat= dense.model, edgemode="directed")

	# par(mfrow = c(2,1))
		

	npb.as.label = function(m)
	{
		r = NULL
		for(i in 1:nrow(m))
			for(j in 1:ncol(m))
				{
					if(consensus.model[i, j] == 1 && standard.model[i,j] == 1)
									r = c(r, paste(m[i,j], '(*)'))
					if(consensus.model[i, j] == 1 && standard.model[i,j] == 0)
									r = c(r, paste(m[i,j], '(C)'))
					if(consensus.model[i, j] == 0 && standard.model[i,j] == 1)
									r = c(r, paste(m[i,j], '(S)'))
					if(consensus.model[i, j] == 0 && standard.model[i,j] == 0)
									r = c(r, paste(m[i,j]))
									
					names(r)[length(r)] = paste(rownames(m)[i], '~' ,  colnames(m)[j], sep = '')
				}
				
		return(r)
	}
	
	cns.as.color = function(m)
	{
		cols = NULL
		 # matrix('black', nrow = nrow(true.model), ncol = ncol(true.model))
		# colnames(cols) = colnames(true.model)
		# rownames(cols) = rownames(true.model) 

		for(i in 1:nrow(m))
			for(j in 1:ncol(m))
				{
					if(consensus.model[i,j] == 1)
					{
															cols = c(cols, paste('red'))
									
						names(cols)[length(cols)] = paste(rownames(m)[i], '~' ,  colnames(m)[j], sep = '')
					}
					if(true.model[i,j] == 0)
					{
															cols = c(cols, paste('gray'))
									
						names(cols)[length(cols)] = paste(rownames(m)[i], '~' ,  colnames(m)[j], sep = '')
					}

								}
				
		return(cols)
	}
	
	
	npb.as.lwd = function(m)
	{
		print(m)
		r = NULL
		for(i in 1:nrow(m))
			for(j in 1:ncol(m))
				{
					r = c(r, as.numeric(paste(m[i,j])))
					names(r)[length(r)] = paste(rownames(m)[i], '~' ,  colnames(m)[j], sep = '')
				}

				
		return(r)
	}


	attrs <- list(node=list(shape="ellipse", fixedsize=FALSE))


	eAttrs  = list()
	eAttrs$label <- npb.as.label(X.scores)
	eAttrs$color <- cns.as.color(X.scores)
	# eAttrs$lwd <-  npb.as.lwd(X.scores)

print(eAttrs)
	plot(D.true.model, attrs = list(node = list(fillcolor = "lightblue"),
                                edge = list(arrowsize=0.5, fontsize = 12)), edgeAttrs=eAttrs)


	# plot(D.consensus.model, attrs = list(node = list(fillcolor = "lightblue"),
                                # edge = list(arrowsize=0.5)))

	
}