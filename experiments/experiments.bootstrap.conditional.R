# set the working directory
setwd("~/Desktop/results_bootstrap_conditional")

# load the needed scripts and libraries
source("build.consensus.structure.R")
source("build.conditional.structure.R")
source("statistics.R")
library(igraph)

# load the data
load(paste0(getwd(),"/results_bootstrap.RData"))

# read all the data and get the results
results_bootstrap_conditional = NULL
results_bootstrap = get("results_bootstrap",globalenv())
cont = 0
for (runs in names(results_bootstrap[["run"]])) {
    cont = cont + 1
    cat(cont/length(names(results_bootstrap[["run"]])),"\n")
    for (var in names(results_bootstrap[["run"]][[runs]])) {
        for (den in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]])) {
            curr_ground_true = results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["adj.matrix"]][["adj.matrix.true"]]
            for (den in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]])) {
                for (samples in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["reconstructions"]])) {
                    for (noise_levels in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["reconstructions"]][[samples]])) {
                        for (reg in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["reconstructions"]][[samples]][[noise_levels]])) {
                            curr_bootstrap_results = results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["reconstructions"]][[samples]][[noise_levels]][[reg]][["boot"]]
                            curr_cardinality = results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["reconstructions"]][[samples]][[noise_levels]][[reg]][["card"]]
                            curr_consensus_estimation = build.conditional.consensus(curr_bootstrap_results,curr_cardinality)
                            curr_stats = getStats(curr_ground_true,curr_consensus_estimation)
                            # save the results
                            results_bootstrap_conditional[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["adj.matrix"]][["adj.matrix.true"]] = curr_ground_true
                            results_bootstrap_conditional[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["reconstructions"]][[samples]][[noise_levels]][[reg]] = curr_consensus_estimation
                            results_bootstrap_conditional[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["stats"]][[samples]][[noise_levels]][[reg]] = curr_stats
                        }
                    }
                }
            }
        }
    }
}

# save the results
save(results_bootstrap_conditional,file="results_bootstrap_conditional.RData")
