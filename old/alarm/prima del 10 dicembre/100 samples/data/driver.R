test.pvalue = 0.01
nboot = 100
fname = '100-'
source('main.test.R')
nboot = 1000
fname = '1000-'
source('main.test.R')
test.pvalue = 0.05
nboot = 100
fname = '100-0.05-'
source('main.test.R')
test.pvalue = 0.001
nboot = 100
fname = '100-0.001-'
source('main.test.R')