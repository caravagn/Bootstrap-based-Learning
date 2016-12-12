# set the working directory
setwd("~/Desktop/results_bootstrap_agony_poset")

# load the needed scripts and libraries
source("build.agony.structure.R")
source("statistics.R")
library(bnlearn)
library(igraph)

# load the data
load(paste0(getwd(),"/results_bootstrap.RData"))
load(paste0(getwd(),"/agony_rankings.RData"))

# read all the data and get the results
results_bootstrap_agony_poset = NULL
results_bootstrap = get("results_bootstrap",globalenv())
agony_rankings = get("agony_rankings",globalenv())
cont = 0
for (runs in names(results_bootstrap[["run"]])) {
    cont = cont + 1
    cat(cont/length(names(results_bootstrap[["run"]])),"\n")
    for (var in names(results_bootstrap[["run"]][[runs]])) {
        for (den in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]])) {
            curr_ground_true = results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["adj.matrix"]][["adj.matrix.true"]]
            #for (den in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]])) {
                for (samples in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["reconstructions"]])) {
                    for (noise_levels in names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["reconstructions"]][[samples]])) {
                        for (reg in c("no_reg",names(results_bootstrap[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["reconstructions"]][[samples]][[noise_levels]]))) {
                            my.curr.reg = gsub(".hc","",gsub("adj.matrix.","",reg))
                            if(var=="continuous") {
                                my.curr.reg = paste0(my.curr.reg,"-g")
                            }
                            curr_dataset = read.table(file=paste0(getwd(),"/data/run/",runs,"/",var,"/","node_size/10/density/",den,"/dataset_sample_size_",samples,"_noise_level_",noise_levels,".txt"),header=FALSE,sep=",",check.names=FALSE,stringsAsFactors=FALSE)
                            
                            # get the agony results
                            curr_rank_loose = agony_rankings[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["ranks"]][[samples]][[noise_levels]][["all"]]
                            curr_rank_strict = agony_rankings[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["ranks"]][[samples]][[noise_levels]][["fixed"]]
                            
                            if(reg=="no_reg") {
                                curr_agony_estimation = curr_rank_strict
                            }
                            else {
                                curr_agony_estimation = build.agony.structure(curr_dataset,curr_rank_loose,my.curr.reg)
                            }
                            
                            curr_stats = getStats(curr_ground_true,curr_agony_estimation)
                            # save the results
                            results_bootstrap_agony_poset[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["adj.matrix"]][["adj.matrix.true"]] = curr_ground_true
                            results_bootstrap_agony_poset[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["reconstructions"]][[samples]][[noise_levels]][[reg]] = curr_agony_estimation
                            results_bootstrap_agony_poset[["run"]][[runs]][[var]][["node_size"]][["10"]][["density"]][[den]][["stats"]][[samples]][[noise_levels]][[reg]] = curr_stats
                        }
                    }
                }
            #}
        }
    }
}

# save the results
save(results_bootstrap_agony_poset,file="results_bootstrap_agony_poset.RData")
