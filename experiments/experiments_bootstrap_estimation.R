# set the working directory
setwd("~/Desktop/results_bootstrap_estimation")

# load the needed scripts and libraries
source("bootstrap.likelihood.fit.R")
library("bnlearn")
library("parallel")
library("doParallel")

# set the seed
set.seed(3333333)

# enumerate all the files to be read
all.files = list.files(paste0(getwd(),"/data"),recursive=TRUE)

# read all the data and get the results
restuls_bootstrap = NULL
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
        # perform the boostrap estimation
        dataset = read.table(file=paste0(getwd(),"/data/",i),header=FALSE,sep=",",check.names=FALSE,stringsAsFactors=FALSE)
        if(curr_split[3]=="continuous") {
            adj.matrix.loglik.hc = bootstrap.estimation(dataset,"loglik-g","hc")
            adj.matrix.loglik.hc = list(boot=adj.matrix.loglik.hc,card=cardinality.parent.set(adj.matrix.loglik.hc))
            #adj.matrix.loglik.tabu = bootstrap.estimation(dataset,"loglik-g","tabu")
            #adj.matrix.loglik.tabu = list(boot=adj.matrix.loglik.tabu,card=cardinality.parent.set(adj.matrix.loglik.tabu))
            adj.matrix.aic.hc = bootstrap.estimation(dataset,"aic-g","hc")
            adj.matrix.aic.hc = list(boot=adj.matrix.aic.hc,card=cardinality.parent.set(adj.matrix.aic.hc))
            #adj.matrix.aic.tabu = bootstrap.estimation(dataset,"aic-g","tabu")
            #adj.matrix.aic.tabu = list(boot=adj.matrix.aic.tabu,card=cardinality.parent.set(adj.matrix.aic.tabu))
            adj.matrix.bic.hc = bootstrap.estimation(dataset,"bic-g","hc")
            adj.matrix.bic.hc = list(boot=adj.matrix.bic.hc,card=cardinality.parent.set(adj.matrix.bic.hc))
            #adj.matrix.bic.tabu = bootstrap.estimation(dataset,"bic-g","tabu")
            #adj.matrix.bic.tabu = list(boot=adj.matrix.bic.tabu,card=cardinality.parent.set(adj.matrix.bic.tabu))
        }
        else {
            adj.matrix.loglik.hc = bootstrap.estimation(dataset,"loglik","hc")
            adj.matrix.loglik.hc = list(boot=adj.matrix.loglik.hc,card=cardinality.parent.set(adj.matrix.loglik.hc))
            #adj.matrix.loglik.tabu = bootstrap.estimation(dataset,"loglik","tabu")
            #adj.matrix.loglik.tabu = list(boot=adj.matrix.loglik.tabu,card=cardinality.parent.set(adj.matrix.loglik.tabu))
            adj.matrix.aic.hc = bootstrap.estimation(dataset,"aic","hc")
            adj.matrix.aic.hc = list(boot=adj.matrix.aic.hc,card=cardinality.parent.set(adj.matrix.aic.hc))
            #adj.matrix.aic.tabu = bootstrap.estimation(dataset,"aic","tabu")
            #adj.matrix.aic.tabu = list(boot=adj.matrix.aic.tabu,card=cardinality.parent.set(adj.matrix.aic.tabu))
            adj.matrix.bic.hc = bootstrap.estimation(dataset,"bic","hc")
            adj.matrix.bic.hc = list(boot=adj.matrix.bic.hc,card=cardinality.parent.set(adj.matrix.bic.hc))
            #adj.matrix.bic.tabu = bootstrap.estimation(dataset,"bic","tabu")
            #adj.matrix.bic.tabu = list(boot=adj.matrix.bic.tabu,card=cardinality.parent.set(adj.matrix.bic.tabu))
        }
        name.res = "reconstructions"
        #curr_result = list(adj.matrix.loglik.hc=adj.matrix.loglik.hc,adj.matrix.loglik.tabu=adj.matrix.loglik.tabu,adj.matrix.aic.hc=adj.matrix.aic.hc,adj.matrix.aic.tabu=adj.matrix.aic.tabu,adj.matrix.bic.hc=adj.matrix.bic.hc, adj.matrix.bic.tabu=adj.matrix.bic.tabu)
        curr_result = list(adj.matrix.loglik.hc=adj.matrix.loglik.hc,adj.matrix.aic.hc=adj.matrix.aic.hc,adj.matrix.bic.hc=adj.matrix.bic.hc)
    }
    if(name.res=="adj.matrix") {
        restuls_bootstrap[[curr_split[1]]][[curr_split[2]]][[curr_split[3]]][[curr_split[4]]][[curr_split[5]]][[curr_split[6]]][[curr_split[7]]][name.res] = list(curr_result)
    }
    else {
        restuls_bootstrap[[curr_split[1]]][[curr_split[2]]][[curr_split[3]]][[curr_split[4]]][[curr_split[5]]][[curr_split[6]]][[curr_split[7]]][[name.res]][[gsub(".txt","",unlist(strsplit(curr_split[8],"_")))[4]]][gsub(".txt","",unlist(strsplit(curr_split[8],"_")))[7]] = list(curr_result)
    }
    
    # compute the statistics
    # if(name.res=="reconstructions") {
        # curr_stats = NULL
        # for (j in names(curr_result)) {
            # curr_stats[[j]] = list(getStats(adj.matrix,curr_result[[j]]))
        # }
        # restuls_bootstrap[[curr_split[1]]][[curr_split[2]]][[curr_split[3]]][[curr_split[4]]][[curr_split[5]]][[curr_split[6]]][[curr_split[7]]][["statistics"]][[gsub(".txt","",unlist(strsplit(curr_split[8],"_")))[4]]][gsub(".txt","",unlist(strsplit(curr_split[8],"_")))[7]] = list(curr_stats)
    # }
    
    cont = cont + 1
    print(cont/length(all.files))
}

# save the results
save(restuls_bootstrap,file="restuls_bootstrap.RData")
