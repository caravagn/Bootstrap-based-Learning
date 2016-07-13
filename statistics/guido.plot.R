
guido.plot = function(true.model, X.scores, Y.scores, labels = c('NPB','MI'), loop = NULL)
{
	if(!is.matrix(true.model)) true.model = as.matrix(true_adj_matrix)
	if(!is.matrix(X.scores)) X.scores = as.matrix(X.scores)
	if(!is.matrix(Y.scores)) Y.scores = as.matrix(Y.scores)
	
	ordering = order(X.scores)

	plot.new()
	par(mfrow = c(2,1))
	
	# X.scores and Y.scores
	data = NULL	
	for(i in 1:length(ordering))
		data = rbind(data, c(i, X.scores[ordering[i]], Y.scores[ordering[i]]))
	
	colnames(data) = c('X', 'X.score', 'Y.score')
	
	# # true edges as points
	# true.data = NULL	
	# for(i in 1:length(ordering))
		# if(true.model[ordering[i]] == 1)
			# true.data = rbind(true.data, c(i, X.scores[ordering[i]], Y.scores[ordering[i]]))
	# colnames(true.data) = c('X', 'X.score', 'Y.score')
	
	# # data frames
	# tdframe = data.frame(row.names = 1:nrow(true.data))
	# tdframe$X = as.numeric(true.data[, 1])
	# tdframe$X.score = as.numeric(true.data[, 2])
	# tdframe$Y.score = as.numeric(true.data[, 3])
	
	# print(head(tdframe))
	
	dframe = data.frame(row.names = 1:length(ordering))
	dframe$X = as.numeric(data[, 1])
	dframe$X.score = as.numeric(data[, 2])
	dframe$Y.score = as.numeric(data[, 3])
	dframe$true = rep('spurious', length(ordering))

	# true edges as points with special mark
	for(i in 1:length(ordering))
		if(true.model[ordering[i]] == 1)
			dframe[i, 'true'] = 'true'
				
	library(ggplot2)
	library(gridExtra)

    mycolours = c("true" = "brown3", "spurious" = "black")

	top.frame = dframe[, c('X', 'X.score', 'true')]
	mid.frame = dframe[, c('X', 'Y.score', 'true')]
	
	p.top = ggplot(top.frame, aes(x = X, y = X.score, fill = labels[1])) +
	ggtitle(paste(labels[1],'score'))+
	xlab(paste('X-axis ordered according to the ', labels[1],'score')) +guides(fill=FALSE) +
    geom_point(size = 2, aes(colour = true)) +
    scale_color_manual("Edge", values = mycolours) 

	p.mid = ggplot(mid.frame, aes(x = X, y = Y.score, fill = labels[2])) +
	xlab(paste('X-axis ordered as in the top panel')) +
	guides(fill=FALSE) +
	ggtitle(paste(labels[2],'score'))+
    geom_point(size = 2, aes(colour = true)) +
    scale_color_manual("Edge", values = mycolours) +
    geom_text(data = mid.frame, size = 3, aes(x = X + 0.05, y = Y.score + 0.05, label = top.frame$X.score)) 



	
	# p.top = ggplot(top.frame, aes(x = X, y = X.score, fill = labels[1])) +
	# geom_jitter(pch = 21, position = position_jitterdodge(jitter.height = 0.2, jitter.width = 2))+
	# ggtitle(paste(labels[1],'score'))+
	# xlab(paste('X-axis ordered according to the ', labels[1],'score')) +guides(fill=FALSE) 

	
	# p.mid = ggplot(mid.frame, aes(x = X, y = Y.score, fill = labels[2])) +
	# geom_jitter(pch = 21, position = position_jitterdodge(jitter.height = 0.2, jitter.width = 2))+
	# xlab(paste('X-axis ordered as in the top panel')) +
	# guides(fill=FALSE) +
	# ggtitle(paste(labels[2],'score'))

	

	
	grid.arrange(p.top, p.mid)
	
}
