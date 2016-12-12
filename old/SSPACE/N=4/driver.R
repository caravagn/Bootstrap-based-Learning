library(doParallel} 
registerDoParallel(cores=2) 

foreach(i=1:5) %dopar% source('SSPACE.R')


