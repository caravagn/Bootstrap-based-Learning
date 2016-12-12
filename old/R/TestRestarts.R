bin2categoric = function(dataset)
{
	M = nrow(dataset)
	N = ncol(dataset)
    
    data = array("missing",c(M,N))
	data[dataset == 1] = "observed"

	data = as.data.frame(data, stringsAsFactors = TRUE)
     
     for (n in names(data)) {
            levels(data[[n]]) = c('missing', 'observed')
     }

    # renaming
    colnames(data) = colnames(dataset)
    rownames(data) = rownames(dataset)

	return(data)        
}


bn2adjmatrix = function(bnobj, dataset)
{
    my.arcs = bnobj$arcs
    
    # structure to save the adjacency matrix of the reconstructed topolgy
    adj.matrix = array(0, c(ncol(dataset), ncol(dataset)))
    rownames(adj.matrix) = colnames(dataset)
    colnames(adj.matrix) = colnames(dataset)

    # create the adjacency matrix of the reconstructed topology
    if(length(nrow(my.arcs))>0 && nrow(my.arcs)>0) {
        for (i in 1:nrow(my.arcs)) {
            # [i,j] refers to the edge i --> j
            adj.matrix[as.numeric(my.arcs[i,1]),as.numeric(my.arcs[i,2])] = 1
        }
    }
    
    return(adj.matrix)
}

# parameter 1: dataset
# parameter 2: maximum number of restarts (step: 100)
# parameter 3: target adjacency matrix
# .. others as usual
restartstest = function(dataset, N = 10^5, target, cores.ratio =0.5, score = 'aic', debug = FALSE)
{	
	require(bnlearn)
	require(parallel)
	require(doParallel)
		
	# tested restarts
	tests = seq(1, N, N/10)
	results = NULL
	
	if(debug) print(paste('** Testing up to ', N, 'restarts with ', score, '(', length(tests),' tests)'))
	
	# data
	catdata = bin2categoric(dataset)

    # set the seed to be used for the random samplings
    set.seed(12345)

    # set the number of cores to be used in the parallelization
    cores = as.integer(cores.ratio * (detectCores() - 1))
    if (cores < 1) cores = 1
    
    # setup the parallelization
    cl = makeCluster(cores)
    registerDoParallel(cl)
    if(debug) print(paste('** Using ', cores, ' cores.'))

    # perform tests
    ntests = length(tests)
    r = foreach(num = 1:ntests, .packages='bnlearn', .export='bn2adjmatrix') %dopar% {

		if(debug) print(paste('** Testing ', tests[num], 'restarts'))
        my.net = bn2adjmatrix(
        	hc(catdata, score = score, restart = tests[num]),
        	catdata)        
    }

    stopCluster(cl)
    results = r
        
    stats = lapply(results, function(x){
    	if(all(x ==target)) return(1) else return(0)
    })
    stats = sum(unlist(stats))
    
	if(debug) print(paste('** Tests completed, found: ', stats/ntests))
	return(stats/ntests)
 }
 
 
# random data
m = 31
n = 12
x = sample.int (2, m*n, TRUE)-1L; 
dim(x) <- c(m,n)
colnames(x) = 1:ncol(x)

# our hypothetical adjancency matrix
o = sample.int (2, n*n, TRUE)-1L; 
dim(o) <- c(n,n)

restartstest(
	dataset = x, 
	N=10^5, 
	target = o, 
	debug = TRUE)