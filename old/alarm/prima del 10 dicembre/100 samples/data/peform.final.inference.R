perform.inference = function( bootstrap_results, curr.poset, test_pvalue ) {
	pvalues_estimates = get.pvalue.estimate(bootstrap_results,curr.poset)
	for(vals in names(pvalues_estimates)) {
		curr_res = pvalues_estimates[[vals]]
		curr_res[which(curr_res<test_pvalue,arr.ind=TRUE)] = -1
		curr_res[which(curr_res>=test_pvalue,arr.ind=TRUE)] = 0
		curr_res = curr_res * -1
		colnames(curr_res) = colnames(curr.poset)
		rownames(curr_res) = rownames(curr.poset)
		pvalues_estimates[[vals]] = curr_res
	}
	return(pvalues_estimates)
}

get.pvalue.estimate = function( bootstrap_results, curr.poset ) {

	bootstrap_res = get.bootstrap.second.pass.res(bootstrap_results)
	bootstrap_model = bootstrap_res[["model"]]
	bootstrap_null = bootstrap_res[["null"]]

	curr_confidence = array(0,c(nrow(bootstrap_model),ncol(bootstrap_model)))
    for(p in 1:nrow(bootstrap_model)) {
        for(q in 1:ncol(bootstrap_model)) {
            curr_confidence[p,q] = wilcox.test(unlist(bootstrap_model[p,q]),unlist(bootstrap_null[p,q]),alternative="greater",mu=0)$p.value
        }
    }
    curr_valid_arcs = which(curr.poset==1,arr.ind=TRUE)
    # perform pvalue correction by FDR
    curr_confidence_fdr_adj = p.adjust(as.vector(curr_confidence[curr_valid_arcs]),method="fdr")
    curr_confidence_fdr = curr_confidence
    curr_confidence_fdr[curr_valid_arcs] = curr_confidence_fdr_adj
    # perform pvalue correction by Holm
    curr_confidence_holm_adj = p.adjust(as.vector(curr_confidence[curr_valid_arcs]),method="holm")
    curr_confidence_holm = curr_confidence
    curr_confidence_holm[curr_valid_arcs] = curr_confidence_holm_adj

    return(list(pvalues=curr_confidence,qvalues.fdr=curr_confidence_fdr,qvalues.holm=curr_confidence_holm))
}

get.bootstrap.second.pass.res = function( bootstrap_second_pass ) {
    
    boot_confidence_1 = array(list(),c(nrow(bootstrap_second_pass[[1]][[1]]),ncol(bootstrap_second_pass[[1]][[1]])))
    boot_null_confidence_1 = array(list(),c(nrow(bootstrap_second_pass[[1]][[1]]),ncol(bootstrap_second_pass[[1]][[1]])))
    for(m in 1:length(bootstrap_second_pass)) {
        curr_boot_confidence_1 = bootstrap_second_pass[[m]][["model"]]
        curr_boot_null_confidence_1 = bootstrap_second_pass[[m]][["null"]]
        for(n in 1:nrow(curr_boot_confidence_1)) {
            for(o in 1:ncol(curr_boot_confidence_1)) {
                boot_confidence_1[n,o] = list(c(unlist(boot_confidence_1[n,o]),curr_boot_confidence_1[n,o]))
                boot_null_confidence_1[n,o] = list(c(unlist(boot_null_confidence_1[n,o]),curr_boot_null_confidence_1[n,o]))
            }
        }
    }
    res1 = list()
    res1[["model"]] = boot_confidence_1
    res1[["null"]] = boot_null_confidence_1
    
    bootstrap_confidence = res1
    
    return(bootstrap_confidence)
}
