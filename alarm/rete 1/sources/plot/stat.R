stats = function(X,Y){
	
	res = NULL
	tp = tn = fp = fn = 0
	
	for (i in 1:nrow(X)) {
		for (j in 1:ncol(X)) {
			if(X[i,j] == Y[i,j] && X[i,j] == 0) tn = tn+1
			if(X[i,j] == Y[i,j] && X[i,j] == 1) tp = tp+1
			if(X[i,j] != Y[i,j] && X[i,j] == 1) fn = fn+1
			if(X[i,j] != Y[i,j] && X[i,j] == 0) fp = fp+1			
		}
	}	

	res$tp = tp
	res$fp = fp
	res$fn = fn
	res$tn = tn

	res$sensitivity = tp/(tp+fn)
	res$specificity = tn/(tn+fp)
	res$precision = tp/(tp+fp)
	res$negativepredictivevalue = tn/(tn+fn)
	res$falsepositiverate = 1-res$sensitivity
	res$falsediscoveryrate = 1- res$precision
	res$accuracy = (tp+tn)/(tp+fp+fn+tn)
	res$f1score = 2*tp/(2*tp+fp+fn)
	res$matthewcorrelation = (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

	return(res)

}

classifyedges = function(X, Y)
{
	reversed = list()
	madeup = list()
	missing = list()

	for (i in 1:nrow(X)) {
		for (j in 1:ncol(Y)) {
			from = rownames(X)[i]
			to = rownames(Y)[j]
			e = paste(from, "~", to, sep = "")

			# Reversed
			if (X[i, j] == 1 && Y[j, i] == 1) reversed = append(reversed, e)

			# Made up
			if (X[i, j] == 0 && Y[i,j] == 1) madeup = append(madeup, e)
			
			# Missing
			if (X[i, j] == 1 && Y[i,j] == 0) missing = append(missing, e)
		}
	}
	
	res = NULL
	res$reversed = unlist(reversed)
	res$madeup = unlist(madeup)
	res$missing = unlist(missing)
	
	res$nreversed = length(res$reversed)
	res$nmadeup = length(res$madeup)
	res$nmissing = length(res$missing)
	
	return(res)
}

