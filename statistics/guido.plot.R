
guido.plot = function(true.model, consensus.model, X.scores, Y.scores, labels = c('NPB','MI'), loop = NULL)
{
	if(!is.matrix(true.model)) true.model = as.matrix(true_adj_matrix)
	if(!is.matrix(X.scores)) X.scores = as.matrix(X.scores)
	if(!is.matrix(Y.scores)) Y.scores = as.matrix(Y.scores)
	if(!is.matrix(consensus.model)) consensus.model = as.matrix(consensus.model)
	# if(!is.matrix(aic.model)) aic.model = as.matrix(aic.model)
	
	ordering = order(X.scores)

	#plot.new()
	par(mfrow = c(2,1))
	
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
	dframe$consensus = rep('', length(ordering))

	# edges as points with special mark
	for(i in 1:length(ordering))
	{
		if(true.model[ordering[i]] == 1)
			dframe[i, 'true'] = 'true'
		
		if(consensus.model[ordering[i]] == 1)
			dframe[i, 'consensus'] = 'C'
	}
				
	library(ggplot2)
	library(gridExtra)

    mycolours = c("true" = "brown3", "spurious" = "black")
    myshape = c("yes" = "a", "no" = "b")

	top.frame = dframe[, c('X', 'X.score', 'true', 'consensus')]
	mid.frame = dframe[, c('X', 'Y.score', 'true', 'consensus')]
	
	p.top = ggplot(top.frame, aes(x = X, y = X.score, fill = labels[1])) +
	ggtitle(paste(labels[1],'score'))+
	xlab(paste('X-axis ordered according to the ', labels[1],'score')) +guides(fill=FALSE) +
    geom_point(size = 2,  aes(colour = true)) +
    scale_color_manual("Edge", values = mycolours) +
    geom_text(data = top.frame, size = 3, aes(x = X + 0.05, y = X.score + 5, label = top.frame$consensus)) 

    # scale_color_manual("Consensus", values = myshape) 
     
     # p.top = ggplot() +
	# # scale_x_discrete(name="") +
	# # scale_y_continuous(limits=c(0,1), breaks=NA, name="") +
	# # scale_shape_discrete(solid=T, legend=F) +
	# geom_point(data= top.frame, mapping=aes(x=X, y=X.score, shape= myshape), size=10)



	p.mid = ggplot(mid.frame, aes(x = X, y = Y.score, fill = labels[2])) +
	xlab(paste('X-axis ordered as in the top panel')) +
	guides(fill=FALSE) +
	ggtitle(paste(labels[2],'score'))+
    geom_point(size = 2, aes(colour = true)) +
    scale_color_manual("Edge", values = mycolours) +
    geom_text(data = mid.frame, size = 3, aes(x = X + 0.05, y = Y.score + 0.05, label = top.frame$X.score)) 

	grid.arrange(p.top, p.mid)
	
}
