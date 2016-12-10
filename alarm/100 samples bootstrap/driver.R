# set the seed
# set.seed(12345)

library(bnlearn)

data(alarm)
data = alarm[1:100, ]
reg = 'aic'

test.pvalue = 0.01
nboot = 100
fname = '100-'
source('main.test.R')

nboot = 1000
fname = '1000-'
source('main.test.R')

nboot = 3000
fname = '3000-'
source('main.test.R')

nboot = 6000
fname = '6000-'
source('main.test.R')

nboot = 10000
fname = '10000-'
source('main.test.R')


# test.pvalue = 0.05
# nboot = 100
# fname = '100-0.05-'
# source('main.test.R')
# test.pvalue = 0.001
# nboot = 100
# fname = '100-0.001-'
# source('main.test.R')