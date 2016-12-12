# set the files for agony computation
get.edge.list = function( bootstrap_results, curr_file ) {
    for(i in 1:length(bootstrap_results)) {
        curr_adj.matrix = bootstrap_results[[i]]
        for(j in 1:nrow(curr_adj.matrix)) {
            for(k in 1:ncol(curr_adj.matrix)) {
                if(curr_adj.matrix[j,k]==1) {
                    cat(as.character(j),as.character(k),"\n",file=curr_file,append=TRUE)
                }
            }
        }
    }
}
