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

	res$TPR = tp/(tp+fn)
	res$SPC = tn/(tn+fp)
	res$PPV = tp/(tp+fp)
	res$NPV = tn/(tn+fn)
	res$FPR = 1-res$SPC
	res$FNR = 1-res$TPR
	res$FDR = 1- res$PPV
	res$ACC = (tp+tn)/(tp+fp+fn+tn)
	res$F1 = 2*tp/(2*tp+fp+fn)
	res$MW = (tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
	
	# res = t(data.frame(res, stringsAsFactors = FALSE))
	# print(res)

	return(res)
}

classifyedges = function(X, Y)
{
	reversedX = reversedY =list()
	madeup = list()
	missing = list()

	for (i in 1:nrow(X)) {
		for (j in 1:ncol(Y)) {
			from = rownames(X)[i]
			to = rownames(Y)[j]
			e =

			# Reversed
			if (X[i, j] == 1 && Y[j, i] == 1) {
				reversedX = append(reversedX, paste(from, "~", to, sep = ""))
				reversedY = append(reversedY, paste(to, "~", from, sep = ""))
			}
			
			# Made up
			if (X[i,j] == 0 && X[j,i] == 0 && Y[i,j] == 1) madeup = append(madeup, paste(from, "~", to, sep = ""))
			
			# Missing
			if (X[i, j] == 1 && Y[i,j] == 0  && Y[j,i] == 0) missing = append(missing, paste(from, "~", to, sep = ""))
		}
	}
	
	res = NULL
	res$reversedX = unlist(reversedX)
	res$reversedY = unlist(reversedY)
	res$madeup = unlist(madeup)
	res$missing = unlist(missing)
	
	res$nreversed = length(res$reversedX)
	res$nmadeup = length(res$madeup)
	res$nmissing = length(res$missing)
	
	return(res)
}

