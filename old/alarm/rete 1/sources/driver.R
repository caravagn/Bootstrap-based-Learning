setwd('/Volumes/Data/Github/Bootstrap-based-Learning/alarm/rete 1/sources')

# set the seed
set.seed(12345)

data = read.table('dataset_sample_size_5_noise_level_0.txt', sep = ',')
data = data.frame(sapply(1:ncol(data), function(x) { data.frame(as.factor(data[,x])) }))
colnames(data) = paste("X", 1:ncol(data), sep ='')
str(data)

adj.true = as.matrix(read.table('adj_matrix.txt', sep = ',', stringsAsFactors = F))
colnames(adj.true) = colnames(data) 
rownames(adj.true) = colnames(data) 

reg = 'bic'
test.pvalue = 0.01
nboot = 100
fname = '100-'
source('main.test.R')

source('plot/plotter.R')


# nboot = 1000
# fname = '1000-'
# source('main.test.R')

# nboot = 3000
# fname = '3000-'
# source('main.test.R')

# nboot = 6000
# fname = '6000-'
# source('main.test.R')

# nboot = 10000
# fname = '10000-'
# source('main.test.R')


# test.pvalue = 0.05
# nboot = 100
# fname = '100-0.05-'
# source('main.test.R')
# test.pvalue = 0.001
# nboot = 100
# fname = '100-0.001-'
# source('main.test.R')