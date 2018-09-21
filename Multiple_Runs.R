library(bnlearn)
library(gRain)
library(reshape2)
library(ggplot2)

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

sampler = function(TMODEL, M1, M1.2, N1, N1.2, N, K){
  
  cat('Sampling N =', N, ' - K =', K, '\n') 
  
  # one liner sample
  one.sample = function(w)
  {
    # sample a new dataset of N observations
    data = rbn(TMODEL, n = N)
    
    delta.small = delta(pre = M1, post = M1.2, data = data)
    delta.large = delta(pre = N1, post = N1.2, data = data)
    
    submodularity = delta.small$differential >= delta.large$differential
    
    c(submodularity, `N` = N)
  }
  
  # Repeat K times the sampler
  points = 
    Reduce(
      rbind,
      lapply(1:K, one.sample)
    )
  
  df.points = data.frame(points)
  rownames(df.points) = NULL
  
  entry = points[1, , drop = FALSE]
  entry[1, 1:3] = colSums(df.points[, 1:3]) / K
  
  print(entry)
  
  
  # # Convert points
  # df.points = data.frame(points)
  # df.points$N = factor(df.points$N, levels = DATA_SIZES)
  # 
  # df.points[, 1:3] = apply(df.points[, 1:3], 2, as.logical)
  # df.points[, 1:3] = apply(df.points[, 1:3], 2, as.numeric)
  # 
  # df.points = split(df.points, f = df.points[ ,'N'])
  # df.points = lapply(df.points, function(x) colSums(x[, 1:3]))
  
  entry
}

test = function(
  TMODEL,
  M1, M1.2, N1, N1.2,
  POINTS_PER_TEST = 1000,
  DATA_SIZES = c(10, 50, 100, 200, 500, 800, 1000, 5000, 10000)
)
{
  points =
    Reduce(
      rbind,
      lapply(DATA_SIZES, sampler, 
             K = POINTS_PER_TEST,
             TMODEL = TMODEL,
             M1 = M1, 
             M1.2 = M1.2, 
             N1 = N1, 
             N1.2 = N1.2)
    )
  rownames(points) = NULL
  
  points
}


cpt = function(LABELS)
{
  # cat("CPT for ", LABELS, '\n')
  formula = as.formula(
    paste('~ ', 
          paste(LABELS, collapse = ' + ')
    )
  )
  
  cptable(
    formula,
    values = VGAM::rdiric(1, rep(1, 2^length(LABELS))), 
    levels = Y)
}

getBN = function(adj_matrix)
{
  pset.nodes = colSums(adj_matrix)
  
  # Marginal and non marginal nodes
  marginal.nodes = pset.nodes[pset.nodes == 0]
  
  nonmarginal.nodes = pset.nodes[pset.nodes > 0]
  nonmarginal.nodes = adj_matrix[, names(nonmarginal.nodes)]
  
  # CPTs for marginal nodes
  marginal.cpts = lapply(names(marginal.nodes), cpt)
  
  # CPTs for non marginal nodes
  nonmarginal.cpts = 
    lapply(
      colnames(nonmarginal.nodes), 
      function(x) 
      {
        to = x
        from = nonmarginal.nodes[, x]
        from = names(from[from > 0])
        
        cpt(c(to, from))
      }
    )
  
  # gRain CPTs, and network
  tables = compileCPT(append(marginal.cpts, nonmarginal.cpts))
  NET = grain(tables)
  
  # bnlearn converstion
  BN = as.bn.fit(NET)
  BN.noparams = as.bn(NET)
  
  ########### We select now the true edge for which we compute the differential 

  list(BN = BN, BN.noparams = BN.noparams, NET = NET, tables = tables, adj_matrix = adj_matrix)
}


createTest = function(NODES, plots = TRUE)
{
  # TRUE MODEL -- fix a certain structure (random, but with at least one edge with NODES/3 incoming edges)
  repeat {
    adj_matrix = bnlearn::random.graph(LETTERS[1:NODES])
    adj_matrix = revolver:::DataFrameToMatrix(adj_matrix$arcs)
    
    pset.nodes = colSums(adj_matrix)
    if(max(pset.nodes) > round(NODES/3)) break
  }
  
  print(adj_matrix)
  TMODEL = getBN(adj_matrix)
  
  # graphviz.plot(TMODEL$BN)
  
  # Select the true edge for which we compute the differential 
  to = names(which.max(pset.nodes)) 

  # Define its parent set
  true.parents = adj_matrix[, to]
  true.parents = true.parents[true.parents > 0]
  true.parents = names(true.parents)
  
  true.from = sample(true.parents, 1)
  true.parents = setdiff(true.parents, true.from)
  
  cat("Differential edge: ", true.from, " --> ", to, "\n")
  
  # nested parents sets from smaller to larger models
  small = sample(true.parents, 1)
  large = true.parents
  
  cat("Small parent set: ", small, "\n")
  cat("Large parent set: ", large, "\n")
  
  
  adj_matrix[, to] = 0
  
  # Small model
  small.adj_matrix.pre = adj_matrix
  small.adj_matrix.pre[small, to] = 1
  
  small.adj_matrix.post = small.adj_matrix.pre
  small.adj_matrix.post[true.from, to] = 1

  # Larger model
  large.adj_matrix.pre = small.adj_matrix.pre
  large.adj_matrix.pre[large, to] = 1
  
  large.adj_matrix.post = large.adj_matrix.pre
  large.adj_matrix.post[true.from, to] = 1
  
  # BN models to test
  small.pre = getBN(small.adj_matrix.pre)$BN.noparams
  small.post = getBN(small.adj_matrix.post)$BN.noparams
  
  large.pre = getBN(large.adj_matrix.pre)$BN.noparams
  large.post = getBN(large.adj_matrix.post)$BN.noparams
  
  if(plots)
  {
    par(mfrow = c(2,2))
    graphviz.plot(small.pre, main = 'M1 without the true edge (pre)', highlight = list(nodes = small, fill = 'red'))
    graphviz.plot(small.post, main = 'M1,2 plus the true edge (post)', highlight = list(nodes = c(small, true.from), fill = 'red'))
    
    graphviz.plot(large.pre, main = 'N1 without the true edge (pre)', highlight = list(nodes = large, fill = 'red'))
    graphviz.plot(large.post, main = 'N1,2 plus the true edge (post)', highlight = list(nodes = c(large, true.from), fill = 'red'))
  }  
  
  list(TRUTH = TMODEL$BN, 
       edge = c(true.from, to),
       small.parentSet = small,
       large.parentSet = large,
       M1 = small.pre, 
       M1.2 = small.post,
       N1 = large.pre,
       N1.2 = large.post)
}


plotTest = function(BN.test, points) {
  
  # layout(matrix(c(1,1,1,1,2:5), ncol = 2, byrow = T))  

  # layout(matrix(c(1,1,1,1,2:5), ncol = 2, byrow = T))  
  layout(matrix(c(1:6), ncol = 2, byrow = T))  
  
  graphviz.plot(BN.test$TRUTH, main = 'Ground thruth')
  
  # scores
  plot(x = points[, 'N'], y = points[, 'BIC'], bg = 'steelblue', pch = 21, ylim = c(0,1),
       xlab = 'Sample size', ylab = 'Score', cex = 2, border = NA, log = 'x')
  lines(x = points[, 'N'], y = points[, 'BIC'], col = 'steelblue', log = 'x')
  
  points(x = points[, 'N'], y = points[, 'AIC'], bg = 'red', pch = 22, xlab = 'Sample size', ylab = 'Score',  cex = 2, border = NA, log = 'x')
  lines(x = points[, 'N'], y = points[, 'AIC'], col = 'red', log = 'x')
  
  points(x = points[, 'N'], y = points[, 'logLik'], bg = 'orange', pch = 24, xlab = 'Sample size', ylab = 'Score', cex = 2, border = NA, log = 'x')
  lines(x = points[, 'N'], y = points[, 'logLik'], col = 'orange', log = 'x')
  
  legend('topright', legend = c("logLik", "AIC", "BIC"), col = c('orange', 'red', 'steelblue'), pch = 21, cex = 2, bty = 'n')
  
  
  # M1 and N1
  graphviz.plot(BN.test$M1, main = 'M1 (pre)', highlight = list(nodes = BN.test$small.parentSet, fill = 'red'))
  graphviz.plot(BN.test$M1.2, main = 'M1 (post)', highlight = list(nodes = c(BN.test$small.parentSet, BN.test$edge[1]), fill = 'red'))
  
  graphviz.plot(BN.test$N1, main = 'N1 (pre)', highlight = list(nodes = BN.test$large.parentSet, fill = 'red'))
  graphviz.plot(BN.test$N1.2, main = 'N1 (post)', highlight = list(nodes = c(BN.test$large.parentSet, BN.test$edge[1]), fill = 'red'))
 
 
}



# ....
Y = c("YES", "NO")


pdf("tests.pdf")


for(i in 1:30)
{
  pio::pioTit(i)
  BN.test = createTest(NODES = sample(c(5, 8, 12), 1), plots = FALSE)
  
  points = test(
    TMODEL = BN.test$TRUTH,
    M1 = BN.test$M1,
    M1.2 = BN.test$M1.2,
    N1 = BN.test$N1,
    N1.2 = BN.test$N1.2,
    POINTS_PER_TEST = 100,
    # DATA_SIZES = c(10, 50, 100, 200, 500, 800, 1000, 5000, 10000, 50000, 100000)
    DATA_SIZES = c(100, 1000, 5000, 10000)
    # POINTS_PER_TEST = 100,
    # DATA_SIZES =  c(10, 50, 100)
  )
    
  
  plotTest(BN.test, points)
  
  
  
  # code = paste0(sample(LETTERS, 10), collapse = '')
  # code = paste0(code, '.pdf')
  # pdf(code)
  # plotTest(BN.test, points)
  # dev.off()
}

dev.off()
