# set the working directory
setwd("~/Desktop/results_hill_climbing_tabu")

# load the needed scripts and libraries
source("perform.likelihood.fit.R")
source("statistics.R")
library("bnlearn")
library("parallel")
library("doParallel")

# set the seed
set.seed(2222222)

# enumerate all the files to be read
all.files = list.files(paste0(getwd(),"/data"),recursive=TRUE)

# read all the data and get the results
restuls_hc_tabu = NULL
adj.matrix = NULL
name.res = NULL
cont = 0
for (i in all.files) {
    # read the data
    curr_split = unlist(strsplit(i,"/"))
    if(curr_split[8]=="adj_matrix.txt") {
        adj.matrix = read.table(file=paste0(getwd(),"/data/",i),header=FALSE,sep=",",check.names=FALSE,stringsAsFactors=FALSE)
        name.res = "adj.matrix"
        curr_result = list(adj.matrix.true=adj.matrix)
    }
    else {
        # perform the likelihood fit
        dataset = read.table(file=paste0(getwd(),"/data/",i),header=FALSE,sep=",",check.names=FALSE,stringsAsFactors=FALSE)
        if(curr_split[3]=="continuous") {
            adj.matrix.loglik.hc = perform.my.fit(dataset,"hc","loglik-g")
            adj.matrix.loglik.tabu = perform.my.fit(dataset,"tabu","loglik-g")
            adj.matrix.aic.hc = perform.my.fit(dataset,"hc","aic-g")
            adj.matrix.aic.tabu = perform.my.fit(dataset,"tabu","aic-g")
            adj.matrix.bic.hc = perform.my.fit(dataset,"hc","bic-g")
            adj.matrix.bic.tabu = perform.my.fit(dataset,"tabu","bic-g")
        }
        else {
            adj.matrix.loglik.hc = perform.my.fit(dataset,"hc","loglik")
            adj.matrix.loglik.tabu = perform.my.fit(dataset,"tabu","loglik")
            adj.matrix.aic.hc = perform.my.fit(dataset,"hc","aic")
            adj.matrix.aic.tabu = perform.my.fit(dataset,"tabu","aic")
            adj.matrix.bic.hc = perform.my.fit(dataset,"hc","bic")
            adj.matrix.bic.tabu = perform.my.fit(dataset,"tabu","bic")
        }
        name.res = "reconstructions"
        curr_result = list(adj.matrix.loglik.hc=adj.matrix.loglik.hc,adj.matrix.loglik.tabu=adj.matrix.loglik.tabu,adj.matrix.aic.hc=adj.matrix.aic.hc,adj.matrix.aic.tabu=adj.matrix.aic.tabu,adj.matrix.bic.hc=adj.matrix.bic.hc, adj.matrix.bic.tabu=adj.matrix.bic.tabu)
    }
    restuls_hc_tabu[[curr_split[1]]][[curr_split[2]]][[curr_split[3]]][[curr_split[4]]][[curr_split[5]]][[curr_split[6]]][[curr_split[7]]][name.res] = list(curr_result)
    
    # compute the statistics
    if(name.res=="reconstructions") {
        curr_stats = NULL
        for (j in names(curr_result)) {
            curr_stats[[j]] = list(getStats(adj.matrix,curr_result[[j]]))
        }
        restuls_hc_tabu[[curr_split[1]]][[curr_split[2]]][[curr_split[3]]][[curr_split[4]]][[curr_split[5]]][[curr_split[6]]][[curr_split[7]]]["statistics"] = list(curr_stats)
    }
    
    cont = cont + 1
    print(cont/length(all.files))
}

# save the results
save(restuls_hc_tabu,file="restuls_hc_tabu.RData")
