library(bnlearn)
library(Rgraphviz)
library(igraph)
library(RColorBrewer)
library(plyr)
library(relations)

## Empty adjacency Matrix
emptyadjmat = function() {
	mat = matrix(0, nrow = N, ncol = N)
	colnames(mat) = NAMES
	rownames(mat) = NAMES
	return(mat)
}

## Matrix to binary string
mat2code = function(M) {
	return(paste(as.vector(M), collapse = ""))
}

## Binary string to matrix
code2mat = function(C) {
	M = matrix(as.integer(strsplit(C, split = "")[[1]]), nrow = N, ncol = N)
	colnames(M) = NAMES
	rownames(M) = NAMES
	return(M)
}

## In the search, check if M is scheduled to visit or as visited 
iscached = function(M) {
	code = mat2code(M)
	if (code %in% CACHE) 
		return(TRUE)
	if (code %in% lapply(Q, mat2code)) 
		return(TRUE)
	return(FALSE)
}

## Add alpha to a color
add.alpha <- function(col, alpha = 1) {
	apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}

## Generate a random true model for a test
create.true.model = function() {
	TRUEMODEL = NULL

	# we sample a random structure and accept it if it has <= 2 parents per node
	# -- this makes CPTs' generation easier	
	TRUEMODEL$bn = empty.graph(NAMES)
	while (narcs(TRUEMODEL$bn) == 0) {
		# model = random.graph(NAMES)
		model = random.graph(NAMES, num = 1, method = "melancon")

		# count the size of each parent set
		if (all(lapply(sapply(nodes(model), incoming.arcs, x = model), nrow) <= 2)) 
			TRUEMODEL$bn = model
	}
	TRUEMODEL$bn

	TRUEMODEL$code = mat2code(amat(TRUEMODEL$bn))
	if(DOPLOTS) plot(TRUEMODEL$bn)

	# random params -- create CPTs
	CPTs = list()

	parents.set = lapply(sapply(nodes(TRUEMODEL$bn), incoming.arcs, x = TRUEMODEL$bn), function(z) return(z[, 
		"from"]))

	CPTs = sapply(names(parents.set), function(z) {
		create.random.CPT(parents.set[[z]], z)
	})

	TRUEMODEL$bn = custom.fit(TRUEMODEL$bn, CPTs)

	# generate big dataset
	TRUEMODEL$data = rbn(TRUEMODEL$bn, n = NOBS)

	return(TRUEMODEL)
}

# Create a random CPT for edge j. Parent set i can be empty, or have 1 or 2 nodes at most
create.random.CPT = function(i, j) {
	Y <- c("yes", "no")

	# root node -- easy
	if (length(i) == 0) {
		p = runif(1)
		table = matrix(c(p, 1 - p), ncol = 2, dimnames = list(NULL, Y))

		return(table)
	}

	## i --> j edge
	if (length(i) == 1) {
		aux = list(Y, Y)
		names(aux) = c(i, j)
		p = runif(2)
		table <- matrix(c(p[1], 1 - p[1], p[2], 1 - p[2]), ncol = 2, dimnames = aux)
		return(table)
	}

	## i1/i2 --> j edge	if(length(i) == 2)  
	{
		X1 = i[1]
		X2 = i[2]

		p = runif(4)
		table = c(p[1], 1 - p[1], p[2], 1 - p[2], p[3], 1 - p[3], p[4], 1 - p[4])
		dim(table) = c(2, 2, 2)

		aux = list(Y, Y, Y)
		names(aux) = c(j, X1, X2)
		dimnames(table) = aux

		return(table)
	}
}

# Compute a score for a network
score.network = function(M, score = SCORE) {
	BN = empty.graph(NAMES)
	amat(BN) = M

	BN = bn.fit(BN, TRUEMODEL$data)

	ret = NULL
	if (score == "BIC") 
		ret = BIC(BN, TRUEMODEL$data)
	if (score == "AIC") 
		ret = AIC(BN, TRUEMODEL$data)
	if (score == "logLik") 
		ret = logLik(BN, TRUEMODEL$data)

	return(ret)
}


# Plot all DAGs
plot.DAGS = function(X) {
	nAttrs <- list()
	nAttrs$color <- brewer.pal(N, "Pastel1")
	nAttrs$fillcolor <- brewer.pal(N, "Pastel1")
	nAttrs$label <- rep("", N)
	nAttrs$height = rep(0.2, N)
	nAttrs$width = rep(0.2, N)


	names(nAttrs$color) = NAMES
	names(nAttrs$fillcolor) = NAMES
	names(nAttrs$label) = NAMES
	names(nAttrs$height) = NAMES
	names(nAttrs$width) = NAMES

	size = 10
	if (length(X) > 100) 
		size = floor((length(Q))) * 0.25

	dev.new(width = size, height = size)
	par(mfrow = c(floor(sqrt(length(X))) + 1, floor(sqrt(length(X))) + 1))

	lapply(X, FUN = function(Y) {
		# par(bg = add.alpha('blue', runif(1)))
		# if(runif(1) > 0.5) par(bg ='red')
# else par(bg ='blue')


		DAG = new("graphAM", adjMat = Y, edgemode = "directed")
		plot(DAG, nodeAttrs = nAttrs)
	})
}

# Plot state space
plot.SSPACE = function(show.scores = FALSE, layout = "twopi", which = "SSPACE") {
	# transform weights to log-space
	# SSPACE$weight = log(SSPACE$weight)

	# create graph from data.frame, and get its adjacency matrix
	if (which == "SSPACE") 
		SS = graph.data.frame(SSPACE, directed = TRUE)
	if (which == "OSSPACE") 
		SS = graph.data.frame(OSSPACE, directed = TRUE)


	SS = as.matrix(get.adjacency(SS))

	# compute scores for all graphs
	states = colnames(SS)
	scores = NULL
	for (s in 1:length(states)) scores[s] = score.network(code2mat(states[s]))

	# Color gradient for the score (100 values)
	DeltaColor = max(scores) - min(scores)
	intColor = DeltaColor/100

	colorRamp = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

	nAttrs <- list()
	nAttrs$color = rep("white", length(states))
	nAttrs$fillcolor = rep("white", length(states))
	nAttrs$label = rep("", length(states))
	nAttrs$height = rep(0.2, length(states))
	nAttrs$width = rep(0.2, length(states))
	nAttrs$shape = rep("circle", length(states))

	names(nAttrs$color) = states
	names(nAttrs$fillcolor) = states
	names(nAttrs$label) = states
	names(nAttrs$height) = states
	names(nAttrs$width) = states
	names(nAttrs$shape) = states


	# Scores (usually useless)
	if (show.scores) {
		for (s in 1:length(states)) nAttrs$label[s] = scores[s]
	}

	# The true model is boxed
	mat2code(amat(TRUEMODEL$bn)) %in% states
	nAttrs$shape[mat2code(amat(TRUEMODEL$bn))] = "box"

	# Set the color according to the score 
	for (s in 1:length(states)) {
		colorIndx = floor((scores[s] - min(scores))/intColor)
		if (colorIndx == 0) 
			colorIndx = 1

		nAttrs$color[s] = colorRamp[colorIndx]
		nAttrs$fillcolor[s] = colorRamp[colorIndx]

	}

	DAG = new("graphAM", adjMat = SS, edgemode = "directed")
	dev.new(width = 10, height = 10)
	plot(DAG, nodeAttrs = nAttrs, layout)
}

# plot all (ID = 0) local optima, or one specific (ID > 0)
plot.local.optima = function(G, TMODEL, show.scores = FALSE, layout = "twopi", scale = 0.2, ID = 0) {
	SS = graph.data.frame(G, directed = TRUE)
	L = decompose.graph(SS)

	if(ID == 0)
	{
		k = floor(sqrt(length(L)))

		dev.new(width = 10, height = 10)
		par(mfrow = c(k + 1, k))
	}
	else L = L[ID]
	
	lapply(L, function(X) {
		
		X = as.matrix(get.adjacency(X))

		# compute scores for all graphs
		states = colnames(X)
		scores = NULL
		for (s in 1:length(states)) 
			scores[s] = score.network(code2mat(states[s]))

		# Color gradient for the score (100 values)
		DeltaColor = max(scores) - min(scores)
		intColor = DeltaColor/100

		colorRamp = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

		nAttrs <- list()
		nAttrs$color = rep("white", length(states))
		nAttrs$fillcolor = rep("white", length(states))
		nAttrs$label = rep("", length(states))
		nAttrs$height = rep(scale, length(states))
		nAttrs$width = rep(scale, length(states))
		nAttrs$shape = rep("circle", length(states))

		names(nAttrs$color) = states
		names(nAttrs$fillcolor) = states
		names(nAttrs$label) = states
		names(nAttrs$height) = states
		names(nAttrs$width) = states
		names(nAttrs$shape) = states


		# Scores (usually useless)
		if (show.scores) {
			for (s in 1:length(states)) nAttrs$label[s] = scores[s]
		}
		
		# The true model is boxed
		mat2code(amat(TMODEL)) %in% states
		nAttrs$shape[mat2code(amat(TMODEL))] = "box"

		# print(TMODEL)

		# Set the color according to the score 
		for (s in 1:length(states)) {
			colorIndx = floor((scores[s] - min(scores))/intColor)
			if (colorIndx == 0) 
				colorIndx = 1

			nAttrs$color[s] = colorRamp[colorIndx]
			nAttrs$fillcolor[s] = colorRamp[colorIndx]

		}

		DAG = new("graphAM", adjMat = X, edgemode = "directed")
		plot(DAG, nodeAttrs = nAttrs, "twopi") #circo
	})

}

# test transition M to N for PORDER compliance 
is.order.compliant = function(M, N) {
	if (any(PORDER - M < 0)) 
		return(FALSE)
	if (any(PORDER - N < 0)) 
		return(FALSE)
	return(TRUE)
}

# save random data file
random.save <- function() {
	r <- c(1:1)
	for (i in 1:1) {
		r[i] <- paste(sample(c(0:9, letters, LETTERS), 12, replace = TRUE), collapse = "")
	}
	save(SSPACE, OSSPACE, SCORE, N, NOBS, TRUEMODEL, PORDER, file = paste(r, ".Rdata", sep = ""))
}

# num of local optima
nloptima <- function(X) {
	X = graph.data.frame(X, directed = TRUE)
	return(clusters(X)$no)
}

plot.model = function(M, ORDER, node.size = 0.2)
{
	nAttrs <- list()
	nAttrs$color <- brewer.pal(N, "Pastel1")
	nAttrs$fillcolor <- brewer.pal(N, "Pastel1")
	nAttrs$label <- rep("", N)
	nAttrs$height = rep(node.size, N)
	nAttrs$width = rep(node.size, N)

	names(nAttrs$color) = NAMES
	names(nAttrs$fillcolor) = NAMES
	names(nAttrs$label) = NAMES
	names(nAttrs$height) = NAMES
	names(nAttrs$width) = NAMES
	
	eAttrs = list()
	eAttrs$color = list()
	which = ORDER - M
	ecolors = NULL
	for(i in 1:nrow(M)) { 
		for(j in 1:nrow(M)) {
			if(which[i,j] < 0)
			{
				ecolors = c(ecolors, 'red')
				names(ecolors)[length(ecolors)] = paste(colnames(M)[i], colnames(M)[j], sep = '~')
			}
		}
	}
	
	print(ecolors)
	eAttrs$color = ecolors
	print(eAttrs)
	
	DAG = new("graphAM", adjMat = M, edgemode = "directed")
	plot(DAG, nodeAttrs = nAttrs, edgeAttrs = eAttrs)
	# title('True model')	
	box('figure')

}

# G is either SSPACE/OSSPACE (list of edges)
subgraph.plot = function(G, TMODEL, scale = 0.2, ORDER, file = '')
{
	SS = graph.data.frame(G, directed = TRUE)
	L = decompose.graph(SS)

	sapply(1:length(L), function(X) {
		ID = X
		X = L[[X]]
		X = as.matrix(get.adjacency(X))
		
		# A dataframe of all the states in G, w
		states = data.frame(state = character(0), score = integer(0), stringsAsFactors = FALSE)
		colnames(states) = c("state", "score")
	
		for(i in 1:nrow(X))
		{
			states = rbind(states, 
				data.frame(state = rownames(X)[i], 
					score = score.network(code2mat(rownames(X)[i])), 
					stringsAsFactors = FALSE))
		}
		states = states[order(states$score, decreasing = TRUE),]
		# B = states$state[1:5]
		# states$score = sapply(states $score, round, digits = 2)
		
		# layout(matrix(c(1,2,2,3,2,2,4,5,6),nrow=3, ncol = 3), widths = c(1.5,1,0.5), heights = c(1,1.5,1.5))
		
		# plot.model(code2mat(B[1]), ORDER, node.size = scale)
		# title(paste('Optimum [', states$score[1], ']', sep = '') )
		# plot.local.optima(G = G, TMODEL = TMODEL, ID = ID, scale = scale, show.scores = T)
		# title(paste('Landscape'))
		# plot.model(code2mat(B[2]), ORDER, node.size = scale)
		# title(paste('2nd [', states$score[2], ']', sep = '') )
		# plot.model(code2mat(B[3]), ORDER, node.size = scale )
		# title(paste('3rd [', states$score[3], ']', sep = '') )
		# plot.model(code2mat(B[4]), ORDER, node.size = scale)
		# title(paste('4th [', states$score[4], ']', sep = '') )
		# plot.model(code2mat(B[5]), ORDER, node.size = scale)
		# title(paste('5th [', states$score[5], ']', sep = '') )	
		B = states$state[1:3]
		states$score = sapply(states $score, round, digits = 2)
		
		layout(matrix(c(1,1,1,1,1,1,2,3,4),nrow=3, ncol = 3))
		# , widths = c(1.5,1,0.5), heights = c(1,1.5,1.5))
		
		plot.local.optima(G = G, TMODEL = TMODEL, ID = ID, scale = scale, show.scores = T)
		title(paste('Landscape'))
		plot.model(code2mat(B[1]), ORDER, node.size = scale)
		title(paste('Optimum [', states$score[1], ']', sep = '') )
		plot.model(code2mat(B[2]), ORDER, node.size = scale)
		title(paste('2nd [', states$score[2], ']', sep = '') )
		plot.model(code2mat(B[3]), ORDER, node.size = scale )
		title(paste('3rd [', states$score[3], ']', sep = '') )
		dev.copy2pdf(file=paste(file, '/Optimum-', ID, '.pdf', sep =''))
		}) # lapply	
}


# G is either SSPACE/OSSPACE (list of edges)
compare.optima = function(G, TMODEL, ORDER, scale = 0.2)
{
	SS = graph.data.frame(G, directed = TRUE)
	L = decompose.graph(SS)

	idx = unlist(lapply(L, function(x) optimum.scores(x)))

	k = ceiling(length(L)/3)

	# dev.new(width = 10, height = 10)
	par(mfrow = c(k, 3))

	sapply(order(idx, decreasing = TRUE), function(X) {
		M = L[[X]]
		M = as.matrix(get.adjacency(M))
		
		# A dataframe of all the states in G, w
		states = data.frame(state = character(0), score = integer(0), stringsAsFactors = FALSE)
		colnames(states) = c("state", "score")
	
		for(i in 1:nrow(M))
		{
			states = rbind(states, 
				data.frame(state = rownames(M)[i], 
					score = score.network(code2mat(rownames(M)[i])), 
					stringsAsFactors = FALSE))
		}
		states = states[order(states$score, decreasing = TRUE),]
		B = states$state[1:5]
				
		plot.model(code2mat(B[1]), ORDER)
		tit = paste('Optimum [',round(states$score[1], 2), ']', sep='') 
		if(mat2code(amat(TMODEL)) == B[1]) tit = paste(tit, ' TM')
		
		title(tit)	
		 
	}) # lapply	
}

optimum.scores = function(G)
{
	M = as.matrix(get.adjacency(G))
		
	# A dataframe of all the states in G, w
	states = data.frame(state = character(0), score = integer(0), stringsAsFactors = FALSE)
	colnames(states) = c("state", "score")
	
	for(i in 1:nrow(M))
	{
		states = rbind(states, 
			data.frame(state = rownames(M)[i], 
				score = score.network(code2mat(rownames(M)[i])), 
				stringsAsFactors = FALSE))
	}
	states = states[order(states$score, decreasing = TRUE),]
	return(states$score[1])
}

