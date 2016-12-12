source('Functions.R')

# Model size, sample size and score
N = 3
NOBS = 10000
SCORE = "BIC"
DOPLOTS = FALSE

# Variables
NAMES = paste("n", 1:N, sep = "")

# Create a true model, and sample from data a big dataset (NOBS)
TRUEMODEL = create.true.model()

# The true partial ordering
PORDER = amat(TRUEMODEL$bn)
PORDER = relation_incidence(transitive_closure(as.relation(PORDER)))
PORDER = matrix(unlist(PORDER[1:(N * N)]), nrow = N, ncol = N)
colnames(PORDER) = NAMES
rownames(PORDER) = NAMES

# State space is a graph
SSPACE = data.frame(from = character(0), to = character(0), weight = integer(0), stringsAsFactors = FALSE)
colnames(SSPACE) = c("from", "to", "weight")

OSSPACE = data.frame(from = character(0), to = character(0), weight = integer(0), stringsAsFactors = FALSE)
colnames(OSSPACE) = c("from", "to", "weight")

# Data structures for recursive visits
Q = list(emptyadjmat())
R = list()
CACHE = list()
while (length(Q) > 0) {
	
	print(paste('#Q =', length(Q), '#CACHE =', length(CACHE), '#R =', length(R)))
	head = Q[[1]]

	# Signature for the head, that we add to the CACHE of visited structures
	code = mat2code(head)
	CACHE[[length(CACHE) + 1]] = code

	# Score function for the head
	head_score = score.network(head)

	# 1-edge move
	for (i in 1:N) {
		for (j in 1:N) {
			M = head
			if (i != j && head[i, j] == 0) { # new edge, non-self
				M[i, j] = 1
				GM = graph_from_adjacency_matrix(M)
				DAG_flag = is_dag(GM)

				# If this is a DAG -- I want it to be connected to head in SSPACE
				if (DAG_flag) {
					# Invariant, if the weight > 0, than moving to M increases the score function
					from = code
					to = mat2code(M)
					weight = score.network(M) - head_score

					SSPACE = rbind(SSPACE, data.frame(from, to, weight, stringsAsFactors = FALSE))

					# check it against the PRODER, if it is also valid store it
					if (is.order.compliant(head, M)) 
						OSSPACE = rbind(OSSPACE, data.frame(from, to, weight, stringsAsFactors = FALSE))
				}

				if (DAG_flag && !iscached(M)) { # is acyclic and we do not have visited/scheduled that
					Q[[length(Q) + 1]] = M
					R[[length(R) + 1]] = M
				}
			}
		}
	}

	if (length(Q) > 1) 
		Q = Q[-1]
	else Q = list()
}
# Q
# R
# CACHE
# SSPACE
# OSSPACE
# SSPACE == OSSPACE

if(DOPLOTS) plot.DAGS(R)

# The state space so far is undirected, but we want to maintain the only edges 
# that reflect the optimization gradient. We compute now the correct direction 
# for each edge  
for (i in 1:nrow(SSPACE)) {
	## weight < 0, I backflip from with to
	if (SSPACE[i, "weight"] < 0) {
		x = SSPACE[i, "from"]
		SSPACE[i, "from"] = SSPACE[i, "to"]
		SSPACE[i, "to"] = x
		SSPACE[i, "weight"] = SSPACE[i, "weight"] * -1
	}
}
for (i in 1:nrow(OSSPACE)) {
	if (OSSPACE[i, "weight"] < 0) {
		x = OSSPACE[i, "from"]
		OSSPACE[i, "from"] = OSSPACE[i, "to"]
		OSSPACE[i, "to"] = x
		OSSPACE[i, "weight"] = OSSPACE[i, "weight"] * -1
	}

}

# Then, HC moves towards the maximum of the gradient, so for every edge move
# I will maintain in SSPACE only those that maximize the score
# SSPACE[order(SSPACE$from),]
SSPACE = ddply(SSPACE, .(from), function(x) x[which(x$weight == max(x$weight)), ])
OSSPACE = ddply(OSSPACE, .(from), function(x) x[which(x$weight == max(x$weight)), ])

# SSPACE
# OSSPACE
print(paste('Original State-space (number of states):', nrow(SSPACE)))
print(paste('Original State-space (optima):', nloptima(SSPACE)))
print(paste(
	' Reduced State-space (number of states):', nrow(OSSPACE), 
	'[deflated =', (nrow(SSPACE) - nrow(OSSPACE))/nrow(SSPACE), '%]'))
print(paste('Reduced State-space (optima):', nloptima(OSSPACE)))

random.save()

if(DOPLOTS) {
	plot.SSPACE()
	plot.SSPACE(which = "OSSPACE")
}

# plot.SSPACE(show.scores = T)
# plot.local.optima(show.scores = T)



