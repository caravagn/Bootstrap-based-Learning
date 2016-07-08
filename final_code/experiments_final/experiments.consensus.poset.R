# set the working directory
setwd("~/Desktop/results_consensus_poset")

# load the needed scripts and libraries
source("consensus.likelihood.fit.R")
source("statistics.R")
library(bnlearn)
library(igraph)

# set the seed
set.seed(4444444)

# load the data
load(paste0(getwd(),"/results_bootstrap.RData"))

# enumerate all the files to be read
all.files = list.files(paste0(getwd(),"/data"),recursive=TRUE)

# read all the data and get the results
results_consensus_poset = NULL
results_bootstrap = get("results_bootstrap",globalenv())
cont = 0
for (runs in names(results_bootstrap[["run"]])) {
    cont = cont + 1
    cat(cont/length(names(results_bootstrap[["run"]])),"\n")
    for (var in names(results_bootstrap[["run"]][[runs]])) {
        for (n_s in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]])) {
            for (den in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]])) {
                curr_ground_true = results_bootstrap[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]][[den]][["adj.matrix"]][["adj.matrix.true"]]
                #for (den in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]])) {
                    for (samples in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]][[den]][["reconstructions"]])) {
                        for (noise_levels in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]][[den]][["reconstructions"]][[samples]])) {
                            for (reg in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]][[den]][["reconstructions"]][[samples]][[noise_levels]])) {
                                my.curr.reg = gsub(".hc","",gsub("adj.matrix.","",reg))
                                if(var=="continuous") {
                                    my.curr.reg = paste0(my.curr.reg,"-g")
                                }
                                curr_dataset = read.table(file=paste0(getwd(),"/data/run/",runs,"/",var,"/","node_size/",as.character(n_s),"/density/",den,"/dataset_sample_size_",samples,"_noise_level_",noise_levels,".txt"),header=FALSE,sep=",",check.names=FALSE,stringsAsFactors=FALSE)
                                curr_bootstrap_results = results_bootstrap[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]][[den]][["reconstructions"]][[samples]][[noise_levels]][[reg]][["boot"]]
                                curr_cardinality = results_bootstrap[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]][[den]][["reconstructions"]][[samples]][[noise_levels]][[reg]][["card"]]
                                curr_consensus_estimation = perform.consensus.likelihood.fit(curr_bootstrap_results,curr_dataset,my.curr.reg)
                                curr_stats = getStats(curr_ground_true,curr_consensus_estimation)
                                # save the results
                                results_consensus_poset[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]][[den]][["adj.matrix"]][["adj.matrix.true"]] = curr_ground_true
                                results_consensus_poset[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]][[den]][["reconstructions"]][[samples]][[noise_levels]][[reg]] = curr_consensus_estimation
                                results_consensus_poset[["run"]][[runs]][[var]][["node_size"]][[n_s]][["density"]][[den]][["stats"]][[samples]][[noise_levels]][[reg]] = curr_stats
                            }
                        }
                    }
                #}
            }
        }
    }
}

# save the results
save(results_consensus_poset,file="results_consensus_poset.RData")
