# function to compute the results at each step
"getStats" = function( true_matrix, inferred_matrix ) {
    
    # compute the statistics
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    hamming_distance = 0
    for (i in 1:nrow(inferred_matrix)) {
        for (j in i:ncol(inferred_matrix)) {
            if(i!=j) {
                if(true_matrix[i,j]==0 && inferred_matrix[i,j]==0) {
                    tn = tn + 1
                }
                else if(true_matrix[i,j]==0 && inferred_matrix[i,j]==1) {
                    fp = fp + 1
                    hamming_distance = hamming_distance + 1
                }
                else if(true_matrix[i,j]==1 && inferred_matrix[i,j]==0) {
                    fn = fn + 1
                    hamming_distance = hamming_distance + 1
                }
                else if(true_matrix[i,j]==1 && inferred_matrix[i,j]==1) {
                    tp = tp + 1
                }
            }
        }
    }
    
    # compute the statistics
    accuracy = (tp+tn)/(tp+tn+fp+fn)
    sensitivity = (tp)/(tp+fn)
    specificity = (tn)/(fp+tn)
    
    # return the results
    results_values = list(accuracy=accuracy,sensitivity=sensitivity,specificity=specificity,hamming_distance=hamming_distance)
    return(results_values)
    
}
