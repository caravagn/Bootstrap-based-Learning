library(bnlearn)

# Setup a true model with some custom parameters

# A, marginal node
cptA = matrix(c(0.4, 0.6), ncol = 2, dimnames = list(NULL, c("LOW", "HIGH")))
cptA

# B, marginal node
cptB = matrix(c(0.8, 0.2), ncol = 2, dimnames = list(NULL, c("GOOD", "BAD")))
cptB

# C, marginal node
cptC = matrix(c(0.65, 0.35), ncol = 2, dimnames = list(NULL, c("GOOD", "BAD")))
cptC

# D | A, B, C 
cptD = c(0.5, 0.5, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8, 0.1, 0.9, 0.2, 0.8, 0.5, 0.5, 0.4, 0.6)
dim(cptD) = c(2, 2, 2, 2)
dimnames(cptD) = list("D" = c("TRUE", "FALSE"), "A" =  c("LOW", "HIGH"), "B" = c("GOOD", "BAD"), "C" =  c("GOOD", "BAD"))
cptD

# Bayesian Network model
net = model2network("[A][B][C][D|A:B:C]")
net = custom.fit(net, dist = list(A = cptA, B = cptB, C = cptC, D = cptD))
net

# A function that taken two BNs, returns their scores and the differential between post and pre
# The scores are computed by first fitting MLE parameters to a dataset
delta = function(pre, post, data) {
  
  # MLE fits of the parameters
  pre = bn.fit(pre, data)
  post = bn.fit(post, data)
  

  # auxiliary function that computes a score "f" for BN "x"
  f.marginal = function(f,x){
    f(x, data = data)
  }
  
  # auxiliary function that computes a differential between score "f" for BN "post" and "pre" 
  f.diff = function(f){
    f.marginal(f, post) - f.marginal(f, pre)
  }
  
  # scores that we test are functions
  functions = c(logLik, AIC, BIC)

  # get differnetial scores
  values = sapply(functions, f.diff)
  
  # get marginal scores
  marg.pre = sapply(functions, f.marginal, x = pre)
  marg.post = sapply(functions, f.marginal, x = post)
  
  names(values) = 
    names(marg.pre) =
    names(marg.post) = c('logLik', 'AIC', 'BIC')

  list(differential = values, pre = marg.pre, post = marg.post)
}


# M1 vs M1,2
M1 = model2network("[A][B][C][D|A]")
M1.2 = model2network("[A][B][C][D|A:B]")

# N1 vs N1,2
N1 = model2network("[A][B][C][D|A:C]")
N1.2 = model2network("[A][B][C][D|A:B:C]")

par(mfrow = c(3,2))
graphviz.plot(M1, main = 'Smaller parent set, without the true edge', highlight = list(nodes = c("A"), fill = 'red'))
graphviz.plot(M1.2, main = 'Smaller parent set, plus the true edge', highlight = list(nodes = c("A", "B"), fill = 'red'))

graphviz.plot(N1, main = 'Larger parent set, without the true edge ', highlight = list(nodes = c("A", "C"), fill = 'red'))
graphviz.plot(N1.2, main = 'Larger parent set, plus the true edge ', highlight = list(nodes = c("A", "C", "B"), fill = 'red'))

dev.copy2pdf(file = "schema.pdf")

# Sampler function -- reapeat K times
# 1. generate N samples data from "net"
# 2. compute scores via "delta" with M1, M1.2, N1, N1.2
# 3. check the submodular property
sampler = function(N, K){
  
  # one liner sample
  one.sample = function(w)
  {
    data = rbn(net, n = N)
    
    delta.small = delta(pre = M1, post = M1.2, data = data)
    delta.large = delta(pre = N1, post = N1.2, data = data)
    
    submodularity = delta.small$differential > delta.large$differential
    
    c(submodularity, `N` = paste(N))
  }
  
  Reduce(
    rbind,
    lapply(1:K, one.sample)
  )
}

POINTS_PER_TEST = 100 # repetitions 
DATA_SIZES = c(10, 100, 1000, 10000)

# actual test results
points =
  Reduce(
    rbind,
    lapply(DATA_SIZES, sampler, K = POINTS_PER_TEST)
  )
rownames(points) = NULL


# reshape and plot
library(reshape2)
library(ggplot2)

df.points = melt(data.frame(points), "N")

ggplot(df.points, aes(x = value, fill = variable)) +
  geom_bar(position="dodge") +
  facet_wrap(~N)

