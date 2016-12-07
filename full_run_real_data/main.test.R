# set the working directoy
my.wd = "~/Desktop/full_run_real_data"
setwd(my.wd)

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

# STEP 0: read the example datasets
# BINARY
dataset_example_discrete = read.table(file=paste0(getwd(),"/example_datasets/dataset_binary_sample_size_10_noise_level_0.txt"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)
# CONTINUOUS
dataset_example_continuous = read.table(file=paste0(getwd(),"/example_datasets/","dataset_continuous_sample_size_10_noise_level_0.txt"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)

# STEP 1: perform bootstrap to estimate the poset
# BINARY
bootstrap_poset_binary_example = bootstrap.estimation(dataset_example_discrete,"bic")
# CONTINUOUS
bootstrap_poset_continuous_example = bootstrap.estimation(dataset_example_continuous,"bic-g")

# STEP 2: perform the estimation of the poset based on confidence
# BINARY
confidence_based_poset_binary_example = perform.consensus.likelihood.fit(bootstrap_poset_binary_example,dataset_example_discrete)
# CONTINUOUS
confidence_based_poset_continuous_example = perform.consensus.likelihood.fit(bootstrap_poset_continuous_example,dataset_example_continuous)

# STEP 3: perform the estimation of the poset based on agony
unlink(paste0(getwd(),"/agony_files"),recursive=TRUE,force=TRUE)
dir.create(paste0(getwd(),"/agony_files"),showWarnings=FALSE)
# BINARY
get.edge.list(bootstrap_poset_binary_example,paste0(getwd(),"/agony_files/agony_binary_inputs.txt"))
system(paste0("./agony ",getwd(),"/agony_files/agony_binary_inputs.txt ",getwd(),"/agony_files/agony_binary_outputs.txt"))
agony_based_poset_binary_example = build_agony_poset(paste0(getwd(),"/agony_files/agony_binary_outputs.txt"),dataset_example_discrete,bootstrap_poset_binary_example)
# CONTINUOUS
get.edge.list(bootstrap_poset_continuous_example,paste0(getwd(),"/agony_files/agony_continuous_inputs.txt"))
system(paste0("./agony ",getwd(),"/agony_files/agony_continuous_inputs.txt ",getwd(),"/agony_files/agony_continuous_outputs.txt"))
agony_based_poset_continuous_example = build_agony_poset(paste0(getwd(),"/agony_files/agony_continuous_outputs.txt"),dataset_example_continuous,bootstrap_poset_continuous_example)

# STEP 4: perform bootstrap to estimate the DAG both with confidence and agony posets
# BINARY
confidence_bootstrap_inference_binary_example = bootstrap.estimation_constrained(dataset_example_discrete,"bic",confidence_based_poset_binary_example)
agony_bootstrap_inference_binary_example = bootstrap.estimation_constrained(dataset_example_discrete,"bic",agony_based_poset_binary_example)
# CONTINUOUS
confidence_bootstrap_inference_continuous_example = bootstrap.estimation_constrained(dataset_example_continuous,"bic-g",confidence_based_poset_continuous_example)
agony_bootstrap_inference_continuous_example = bootstrap.estimation_constrained(dataset_example_continuous,"bic-g",agony_based_poset_continuous_example)

# STEP 5: perform the final inference on both confidence based and agony based posets
# BINARY
confidence_based_inference_binary_example = perform.inference(confidence_bootstrap_inference_binary_example,confidence_based_poset_binary_example,test.pvalue)
agony_based_inference_binary_example = perform.inference(agony_bootstrap_inference_binary_example,agony_based_poset_binary_example,test.pvalue)
# CONTINUOUS
confidence_based_inference_continuous_example = perform.inference(confidence_bootstrap_inference_continuous_example,confidence_based_poset_continuous_example,test.pvalue)
agony_based_inference_continuous_example = perform.inference(agony_bootstrap_inference_continuous_example,agony_based_poset_continuous_example,test.pvalue)
