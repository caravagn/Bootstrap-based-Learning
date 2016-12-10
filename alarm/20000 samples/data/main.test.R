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

AGOFOLDER = "./agony_files/"
AGOINPUT = paste(AGOFOLDER, nboot, "_", test.pvalue, "_inputs.txt", sep ='')
AGOUTPUT = paste(AGOFOLDER, nboot, "_", test.pvalue, "_outputs.txt", sep ='')

print(paste('Nboot', nboot, 'p-value', test.pvalue))

data = alarm

# \Pi_boot: perform first bootstrap to estimate the poset
Pi_boot = bootstrap.estimation(data,"bic", nboot = nboot)
head(Pi_boot)

print('Pi_boot: OK.')

# \Pi: maximize confidence in \Pi_boot
Pi_confidence = perform.consensus.likelihood.fit(Pi_boot, data)
head(Pi_confidence)

print('Pi_confidence: OK.')

# \Pi: perform the estimation using agony
# unlink(AGOFOLDER, recursive = TRUE,force = TRUE)
# dir.create(AGOFOLDER, showWarnings=FALSE)

print(paste0("Agony: ./agony ", AGOINPUT, ' ', AGOUTPUT, sep =' '))
get.edge.list(Pi_boot, AGOINPUT)
system(paste0("./agony ", AGOINPUT, ' ', AGOUTPUT, sep =' '))
Pi_agony = build_agony_poset(AGOUTPUT, data, Pi_boot)
head(Pi_agony)

print('Pi_agony: OK.')

# Bootstrap to estimate the DAGs, Gamma etc
Gamma_confidence = bootstrap.estimation_constrained(data, "bic", Pi_confidence)
head(Gamma_confidence)

Gamma_agony = bootstrap.estimation_constrained(data, "bic", Pi_agony)
head(Pi_agony)
Gamma_agony
print('Gamma: OK.')

# Final inference on both posets
confidence_model = perform.inference(Gamma_confidence, Pi_confidence, test.pvalue)
agony_model = perform.inference(Gamma_agony, Pi_agony, test.pvalue)

print('Models: OK.')

save(confidence_model, file = paste(fname, 'alarm-confidence.Rdata', sep =''))
save(agony_model, file = paste(fname, 'alarm-agony.Rdata', sep =''))
