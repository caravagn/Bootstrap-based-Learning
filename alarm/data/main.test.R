# set the working directoy
#my.wd = "~/Desktop/full_run_real_data"
#setwd(my.wd)

# load the required R packages
library(bnlearn)
library(parallel)
library(doParallel)
library(igraph)

# source the required scripts
source(paste0(getwd(),"/bootstrap.likelihood.fit.R"))
source(paste0(getwd(),"/consensus.likelihood.fit.R"))
source(paste0(getwd(),"/set.agony.inputs.files.R"))
source(paste0(getwd(),"/build.agony.poset.R"))
source(paste0(getwd(),"/bootstrap.likelihood.fit.constrained.R"))
source(paste0(getwd(),"/peform.final.inference.R"))

# set the seed
set.seed(12345)

# set the pvalue to be used in the test
test.pvalue = 0.01


# data(alarm)
# levels = unique(unlist(lapply(alarm,function(x)levels(x))))
# numbers = 1:length(levels)
# names(numbers) = levels

# dataset_example_discrete = matrix(0, nrow = nrow(alarm), ncol = ncol(alarm))
# str(dataset_example_discrete)
# colnames(dataset_example_discrete) = colnames(alarm)

# for(i in 1: nrow(alarm))
# 	for(j in 1: ncol(alarm))
# 	{
# 		dataset_example_discrete[i,j] = numbers[alarm[i,j]]
# 	}

# head(dataset_example_discrete)

dataset_example_discrete = alarm

# STEP 1: perform bootstrap to estimate the poset
# BINARY

bootstrap_poset_binary_example = bootstrap.estimation(dataset_example_discrete,"bic", nboot = nboot)

# STEP 2: perform the estimation of the poset based on confidence
# BINARY
confidence_based_poset_binary_example = perform.consensus.likelihood.fit(bootstrap_poset_binary_example,dataset_example_discrete)

# STEP 3: perform the estimation of the poset based on agony
unlink(paste0(getwd(),"/agony_files"),recursive=TRUE,force=TRUE)
dir.create(paste0(getwd(),"/agony_files"),showWarnings=FALSE)
# BINARY
get.edge.list(bootstrap_poset_binary_example,paste0(getwd(),"/agony_files/agony_binary_inputs.txt"))
system(paste0("./agony ",getwd(),"/agony_files/agony_binary_inputs.txt ",getwd(),"/agony_files/agony_binary_outputs.txt"))
agony_based_poset_binary_example = build_agony_poset(paste0(getwd(),"/agony_files/agony_binary_outputs.txt"),dataset_example_discrete,bootstrap_poset_binary_example)

# STEP 4: perform bootstrap to estimate the DAG both with confidence and agony posets
# BINARY
confidence_bootstrap_inference_binary_example = bootstrap.estimation_constrained(dataset_example_discrete,"bic",confidence_based_poset_binary_example)
agony_bootstrap_inference_binary_example = bootstrap.estimation_constrained(dataset_example_discrete,"bic",agony_based_poset_binary_example)

# STEP 5: perform the final inference on both confidence based and agony based posets
# BINARY
confidence_based_inference_binary_example = perform.inference(confidence_bootstrap_inference_binary_example,confidence_based_poset_binary_example,test.pvalue)
agony_based_inference_binary_example = perform.inference(agony_bootstrap_inference_binary_example,agony_based_poset_binary_example,test.pvalue)


save(confidence_based_inference_binary_example, file = paste(fname, 'alarm-confidence.Rdata', sep =''))
save(agony_based_inference_binary_example, file = paste(fname, 'alarm-agony.Rdata', sep =''))


# # 
# qvalues.fdr
# pvalues
# qvalues.holm 

# A = agony_based_inference_binary_example$qvalues.holm
# colnames(A) = colnames(alarm)
# rownames(A) = colnames(alarm)
# Abn = empty.graph(colnames(alarm))
# amat(Abn) = A

# graphviz.plot(Abn)
