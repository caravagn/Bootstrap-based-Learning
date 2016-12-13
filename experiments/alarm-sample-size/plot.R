GIULIO = TRUE
DANIELE = !GIULIO

if(GIULIO)
{
	# Giulio
	setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/alarm-sample-size')
	git = '/Volumes/Data/Github/Bootstrap-based-Learning/'
	agony.binaries = '/Volumes/Data/Github/Bootstrap-based-Learning/agony/giulio/agony'
} 
if(DANIELE)
{
	setwd('~/Desktop/')
	git = '....'
	agony.binaries = '/Volumes/Data/Github/Bootstrap-based-Learning/agony/daniele/agony'	
}

# load the required R packages
library(bnlearn)
library(igraph)
source(paste0(git, "utils/plotter.R"))
source(paste0(git, "utils/stat.R"))

# Results from Tom & Jerry
load('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/alarm-sample-size/res1E3.Rdata')
load('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/alarm-sample-size/res1E4.Rdata')
load('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/alarm-sample-size/res1E5.Rdata')
load('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/alarm-sample-size/res2E5.Rdata')

# set the dataset and the true adjacency matrix for the test
data(alarm)

bnalarm = empty.graph(colnames(alarm))
modelstring(bnalarm) = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
    "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]", "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]", 
    "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]", "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]", 
    "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]", "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
bnalarm = bn.fit(bnalarm, alarm)

adj.matrix = amat(bnalarm)


points = data.frame(stringsAsFactors = FALSE)

aux = function(res, which, algo, metric, m, title)
{
	g = lapply(1:length(res), function(i){
		r = res[[i]]
		
		
		if(!is.null(which)) q = r[[which]]
		else q = r
		
		# q = r[[which]]
		
		s = stats(adj.matrix, q[[algo]])	
		
		# print(adj.matrix)
		# print(q[[algo]])
		
		# print(s)
	})
	
	# print(g)
	
	vvals = unlist(lapply(g, function(x){return(x[[metric]])}))
	vsampl = rep(m, length(res))
	algorithm = rep(title, length(res))
	vscore = rep(metric, length(res))
	
	return(data.frame(vvals, algorithm, vscore, vsampl, stringsAsFactors = TRUE))
}

# scores = c('tp', 'fp', 'fn', 'tn', 'TPR', 'SPC', 'PPV', 'NPV', 'FPR', 'FNR', 'FDR', 'ACC', 'F1', 'MW')

### Agony
points = rbind(points, 
	aux(res = res1E3, m = 100, which ='agony.inference', algo = 'pvalues', metric = 'PPV', title = 'Agony with no MHC'),
	aux(res = res1E3, m = 100, which ='agony.inference', algo = 'pvalues', metric = 'TPR', title = 'Agony with no MHC'),
	aux(res = res1E3, m = 100, which ='agony.inference', algo = 'qvalues.fdr', metric = 'PPV', title = 'Agony with FDR'),
	aux(res = res1E3, m = 100, which ='agony.inference', algo = 'qvalues.fdr', metric = 'TPR', title = 'Agony with FDR'),
	aux(res = res1E3, m = 100, which ='agony.inference', algo = 'qvalues.holm', metric = 'PPV', title = 'Agony with Bonferroni'),
	aux(res = res1E3, m = 100, which ='agony.inference', algo = 'qvalues.holm', metric = 'TPR', title = 'Agony with Bonferroni')
	)

points = rbind(points, 
	aux(res = res1E4, m = 1000, which ='agony.inference', algo = 'pvalues', metric = 'PPV', title = 'Agony with no MHC'),
	aux(res = res1E4, m = 1000, which ='agony.inference', algo = 'pvalues', metric = 'TPR', title = 'Agony with no MHC'),
	aux(res = res1E4, m = 1000, which ='agony.inference', algo = 'qvalues.fdr', metric = 'PPV', title = 'Agony with FDR'),
	aux(res = res1E4, m = 1000, which ='agony.inference', algo = 'qvalues.fdr', metric = 'TPR', title = 'Agony with FDR'),
	aux(res = res1E4, m = 1000, which ='agony.inference', algo = 'qvalues.holm', metric = 'PPV', title = 'Agony with Bonferroni'),
	aux(res = res1E4, m = 1000, which ='agony.inference', algo = 'qvalues.holm', metric = 'TPR', title = 'Agony with Bonferroni')
	)

points = rbind(points, 
	aux(res = res1E5, m = 10000, which ='agony.inference', algo = 'pvalues', metric = 'PPV', title = 'Agony with no MHC'),
	aux(res = res1E5, m = 10000, which ='agony.inference', algo = 'pvalues', metric = 'TPR', title = 'Agony with no MHC'),
	aux(res = res1E5, m = 10000, which ='agony.inference', algo = 'qvalues.fdr', metric = 'PPV', title = 'Agony with FDR'),
	aux(res = res1E5, m = 10000, which ='agony.inference', algo = 'qvalues.fdr', metric = 'TPR', title = 'Agony with FDR'),
	aux(res = res1E5, m = 10000, which ='agony.inference', algo = 'qvalues.holm', metric = 'PPV', title = 'Agony with Bonferroni'),
	aux(res = res1E5, m = 10000, which ='agony.inference', algo = 'qvalues.holm', metric = 'TPR', title = 'Agony with Bonferroni')
	)

points = rbind(points, 
	aux(res = res2E5, m = 20000, which ='agony.inference', algo = 'pvalues', metric = 'PPV', title = 'Agony with no MHC'),
	aux(res = res2E5, m = 20000, which ='agony.inference', algo = 'pvalues', metric = 'TPR', title = 'Agony with no MHC'),
	aux(res = res2E5, m = 20000, which ='agony.inference', algo = 'qvalues.fdr', metric = 'PPV', title = 'Agony with FDR'),
	aux(res = res2E5, m = 20000, which ='agony.inference', algo = 'qvalues.fdr', metric = 'TPR', title = 'Agony with FDR'),
	aux(res = res2E5, m = 20000, which ='agony.inference', algo = 'qvalues.holm', metric = 'PPV', title = 'Agony with Bonferroni'),
	aux(res = res2E5, m = 20000, which ='agony.inference', algo = 'qvalues.holm', metric = 'TPR', title = 'Agony with Bonferroni')
	)
	
### Confidence
points = rbind(points, 
	aux(res = res1E3, m = 100, which ='confidence.inference', algo = 'pvalues', metric = 'PPV', title = 'Confidence with no MHC'),
	aux(res = res1E3, m = 100, which ='confidence.inference', algo = 'pvalues', metric = 'TPR', title = 'Confidence with no MHC'),
	aux(res = res1E3, m = 100, which ='confidence.inference', algo = 'qvalues.fdr', metric = 'PPV', title = 'Confidence with FDR'),
	aux(res = res1E3, m = 100, which ='confidence.inference', algo = 'qvalues.fdr', metric = 'TPR', title = 'Confidence with FDR'),
	aux(res = res1E3, m = 100, which ='confidence.inference', algo = 'qvalues.holm', metric = 'PPV', title = 'Confidence with Bonferroni'),
	aux(res = res1E3, m = 100, which ='confidence.inference', algo = 'qvalues.holm', metric = 'TPR', title = 'Confidence with Bonferroni')
	)

points = rbind(points, 
	aux(res = res1E4, m = 1000, which ='confidence.inference', algo = 'pvalues', metric = 'PPV', title = 'Confidence with no MHC'),
	aux(res = res1E4, m = 1000, which ='confidence.inference', algo = 'pvalues', metric = 'TPR', title = 'Confidence with no MHC'),
	aux(res = res1E4, m = 1000, which ='confidence.inference', algo = 'qvalues.fdr', metric = 'PPV', title = 'Confidence with FDR'),
	aux(res = res1E4, m = 1000, which ='confidence.inference', algo = 'qvalues.fdr', metric = 'TPR', title = 'Confidence with FDR'),
	aux(res = res1E4, m = 1000, which ='confidence.inference', algo = 'qvalues.holm', metric = 'PPV', title = 'Confidence with Bonferroni'),
	aux(res = res1E4, m = 1000, which ='confidence.inference', algo = 'qvalues.holm', metric = 'TPR', title = 'Confidence with Bonferroni')
	)

points = rbind(points, 
	aux(res = res1E5, m = 10000, which ='confidence.inference', algo = 'pvalues', metric = 'PPV', title = 'Confidence with no MHC'),
	aux(res = res1E5, m = 10000, which ='confidence.inference', algo = 'pvalues', metric = 'TPR', title = 'Confidence with no MHC'),
	aux(res = res1E5, m = 10000, which ='confidence.inference', algo = 'qvalues.fdr', metric = 'PPV', title = 'Confidence with FDR'),
	aux(res = res1E5, m = 10000, which ='confidence.inference', algo = 'qvalues.fdr', metric = 'TPR', title = 'Confidence with FDR'),
	aux(res = res1E5, m = 10000, which ='confidence.inference', algo = 'qvalues.holm', metric = 'PPV', title = 'Confidence with Bonferroni'),
	aux(res = res1E5, m = 10000, which ='confidence.inference', algo = 'qvalues.holm', metric = 'TPR', title = 'Confidence with Bonferroni')
	)

points = rbind(points, 
	aux(res = res2E5, m = 20000, which ='confidence.inference', algo = 'pvalues', metric = 'PPV', title = 'Confidence with no MHC'),
	aux(res = res2E5, m = 20000, which ='confidence.inference', algo = 'pvalues', metric = 'TPR', title = 'Confidence with no MHC'),
	aux(res = res2E5, m = 20000, which ='confidence.inference', algo = 'qvalues.fdr', metric = 'PPV', title = 'Confidence with FDR'),
	aux(res = res2E5, m = 20000, which ='confidence.inference', algo = 'qvalues.fdr', metric = 'TPR', title = 'Confidence with FDR'),
	aux(res = res2E5, m = 20000, which ='confidence.inference', algo = 'qvalues.holm', metric = 'PPV', title = 'Confidence with Bonferroni'),
	aux(res = res2E5, m = 20000, which ='confidence.inference', algo = 'qvalues.holm', metric = 'TPR', title = 'Confidence with Bonferroni')
	)
	
### HC
points = rbind(points, 
	aux(res = res1E3, m = 100, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = 'PPV', title = 'Hill Climbing without restarts'),
	aux(res = res1E3, m = 100, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = 'TPR', title = 'Hill Climbing without restarts'),
	aux(res = res1E4, m = 1000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = 'PPV', title = 'Hill Climbing without restarts'),
	aux(res = res1E4, m = 1000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = 'TPR', title = 'Hill Climbing without restarts'),
	aux(res = res1E5, m = 10000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = 'PPV', title = 'Hill Climbing without restarts'),
	aux(res = res1E5, m = 10000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = 'TPR', title = 'Hill Climbing without restarts'),
	aux(res = res2E5, m = 20000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = 'PPV', title = 'Hill Climbing without restarts'),
	aux(res = res2E5, m = 20000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = 'TPR', title = 'Hill Climbing without restarts')
	)

# points = rbind(points, 
	# aux(res = res1E3, m = 100, which = NULL, algo = 'hill.climing.with.restarts', metric = 'PPV', title = 'Hill Climbing with restarts'),
	# aux(res = res1E3, m = 100, which = NULL, algo = 'hill.climing.with.restarts', metric = 'TPR', title = 'Hill Climbing with restarts'),
	# aux(res = res1E4, m = 1000, which = NULL, algo = 'hill.climing.with.restarts', metric = 'PPV', title = 'Hill Climbing with restarts'),
	# aux(res = res1E4, m = 1000, which = NULL, algo = 'hill.climing.with.restarts', metric = 'TPR', title = 'Hill Climbing with restarts'),
	# aux(res = res1E5, m = 10000, which = NULL, algo = 'hill.climing.with.restarts', metric = 'PPV', title = 'Hill Climbing with restarts'),
	# aux(res = res1E5, m = 10000, which = NULL, algo = 'hill.climing.with.restarts', metric = 'TPR', title = 'Hill Climbing with restarts'),
	# aux(res = res2E5, m = 20000, which = NULL, algo = 'hill.climing.with.restarts', metric = 'PPV', title = 'Hill Climbing with restarts'),
	# aux(res = res2E5, m = 20000, which = NULL, algo = 'hill.climing.with.restarts', metric = 'TPR', title = 'Hill Climbing with restarts')
	# )

points = rbind(points, 
	aux(res = res1E3, m = 100, which = NULL, algo = 'hill.climing.with.restarts', metric = 'PPV', title = 'Hill Climbing with restarts'),
	aux(res = res1E3, m = 100, which = NULL, algo = 'hill.climing.with.restarts', metric = 'TPR', title = 'Hill Climbing with restarts'),
	aux(res = res1E4, m = 1000, which = NULL, algo = 'hill.climing.with.restarts', metric = 'PPV', title = 'Hill Climbing with restarts'),
	aux(res = res1E4, m = 1000, which = NULL, algo = 'hill.climing.with.restarts', metric = 'TPR', title = 'Hill Climbing with restarts')
	)

for(i in c('tp', 'fp', 'fn', 'tn', 'TPR', 'SPC', 'PPV', 'NPV', 'FPR', 'FNR', 'FDR', 'ACC', 'F1', 'MW'))
{
	points = data.frame(stringsAsFactors = FALSE)
	
	### Agony
	points = rbind(points, 
	aux(res = res1E3, m = 100, which ='agony.inference', algo = 'pvalues', metric = i, title = 'Agony with no MHC'),
	aux(res = res1E3, m = 100, which ='agony.inference', algo = 'qvalues.fdr', metric = i, title = 'Agony with FDR'),
	aux(res = res1E3, m = 100, which ='agony.inference', algo = 'qvalues.holm', metric = i, title = 'Agony with Bonferroni'),
	aux(res = res1E4, m = 1000, which ='agony.inference', algo = 'pvalues', metric = 'PPV', title = 'Agony with no MHC'),
	aux(res = res1E4, m = 1000, which ='agony.inference', algo = 'pvalues', metric = i, title = 'Agony with no MHC'),
	aux(res = res1E4, m = 1000, which ='agony.inference', algo = 'qvalues.fdr', metric = i, title = 'Agony with FDR'),
	aux(res = res1E4, m = 1000, which ='agony.inference', algo = 'qvalues.holm', metric = i, title = 'Agony with Bonferroni'),
	aux(res = res1E5, m = 10000, which ='agony.inference', algo = 'pvalues', metric = i, title = 'Agony with no MHC'),
	aux(res = res1E5, m = 10000, which ='agony.inference', algo = 'qvalues.fdr', metric = i, title = 'Agony with FDR'),
	aux(res = res1E5, m = 10000, which ='agony.inference', algo = 'qvalues.holm', metric = i, title = 'Agony with Bonferroni'),
	# aux(res = res2E5, m = 20000, which ='agony.inference', algo = 'pvalues', metric = i, title = 'Agony with no MHC'),
	# aux(res = res2E5, m = 20000, which ='agony.inference', algo = 'qvalues.fdr', metric = i, title = 'Agony with FDR'),
	# aux(res = res2E5, m = 20000, which ='agony.inference', algo = 'qvalues.holm', metric = i, title = 'Agony with Bonferroni'),
	aux(res = res1E3, m = 100, which ='confidence.inference', algo = 'pvalues', metric = i, title = 'Confidence with no MHC'),
	aux(res = res1E3, m = 100, which ='confidence.inference', algo = 'qvalues.fdr', metric = i, title = 'Confidence with FDR'),
	aux(res = res1E3, m = 100, which ='confidence.inference', algo = 'qvalues.holm', metric = i, title = 'Confidence with Bonferroni'),
	aux(res = res1E4, m = 1000, which ='confidence.inference', algo = 'pvalues', metric = i, title = 'Confidence with no MHC'),
	aux(res = res1E4, m = 1000, which ='confidence.inference', algo = 'qvalues.fdr', metric = i, title = 'Confidence with FDR'),
	aux(res = res1E4, m = 1000, which ='confidence.inference', algo = 'qvalues.holm', metric = i, title = 'Confidence with Bonferroni'),
	aux(res = res1E5, m = 10000, which ='confidence.inference', algo = 'pvalues', metric = i, title = 'Confidence with no MHC'),
	aux(res = res1E5, m = 10000, which ='confidence.inference', algo = 'qvalues.fdr', metric = i, title = 'Confidence with FDR'),
	aux(res = res1E5, m = 10000, which ='confidence.inference', algo = 'qvalues.holm', metric = i, title = 'Confidence with Bonferroni'),
	# aux(res = res2E5, m = 20000, which ='confidence.inference', algo = 'pvalues', metric = i, title = 'Confidence with no MHC'),
	# aux(res = res2E5, m = 20000, which ='confidence.inference', algo = 'qvalues.fdr', metric = i, title = 'Confidence with FDR'),
	# aux(res = res2E5, m = 20000, which ='confidence.inference', algo = 'qvalues.holm', metric = i, title = 'Confidence with Bonferroni'),
	aux(res = res1E3, m = 100, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = i, title = 'Hill Climbing without restarts'),
	aux(res = res1E4, m = 1000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = i, title = 'Hill Climbing without restarts'),
	aux(res = res1E5, m = 10000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = i, title = 'Hill Climbing without restarts'),
	aux(res = res2E5, m = 20000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = i, title = 'Hill Climbing without restarts'),
	aux(res = res1E3, m = 100, which = NULL, algo = 'hill.climing.with.restarts', metric = i, title = 'Hill Climbing with restarts'),
	aux(res = res1E4, m = 1000, which = NULL, algo = 'hill.climing.with.restarts', metric = i, title = 'Hill Climbing with restarts'),
	aux(res = res1E5, m = 10000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = i, title = 'Hill Climbing without restarts')
	# aux(res = res2E5, m = 20000, which = NULL, algo = 'hill.climing.no.restarts.inference', metric = i, title = 'Hill Climbing without restarts')
	)
	
	print('Data loaded')

	pl = ggplot(points, aes(x = algorithm, y = vvals))  +
		geom_boxplot(aes(color = algorithm)) +
		facet_wrap(~ vsampl)  +
		scale_color_brewer(palette="Dark2") +
		theme(
			# axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8, face = "italic"), 
			legend.position = "bottom",
			legend.title = element_text(face = "italic"),
			axis.text.x=element_blank())+
		labs(title = paste(i, "- Alarm Network with 100+100 npb, p<0.01")) +
		xlab("") +
		ylab("")
	pl
	dev.copy2pdf(file = paste0(i,'.pdf'))
}



# # r = unlist(res2E5[[1]], recursive = FALSE)
# str(r)


# colnames(adj.matrix)
# plt(r, NULL, 'hill.climing.no.restarts.inference', adj.matrix, "Hill Climbing without restarts")


# for(i in 1:1)
# {
	# x = res1E3[[i]]
	# q = x[['confidence.inference']]
	# p = x[['confidence.inference.pvalues']]
	# print(q)
	# print(p)
	# # s = stats(adj.matrix, q)	

	# # x2 = plt(results$agony.inference, results$agony.inference.pvalues, 'qvalues.fdr', adj.matrix, AFDR, dev.off = DEV.OFF)

	# plt(q, P, 'qvalues.holm', adj.matrix, "Confidence with Bonferroni", dev.off = TRUE)
# }

# head(points)


library(ggplot2)

dev.new()

for(i in c('tp', 'fp', 'fn', 'tn', 'TPR', 'SPC', 'PPV', 'NPV', 'FPR', 'FNR', 'FDR', 'ACC', 'F1', 'MW'))
{
	subset = points[points$vscore == i, ]
	print(subset)

	pl = ggplot(subset, aes(x = algorithm, y = vvals))  +
		geom_boxplot(aes(color = algorithm)) +
		facet_wrap(~ vsampl)  +
		scale_color_brewer(palette="Dark2") +
		theme(
			# axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8, face = "italic"), 
			legend.position = "bottom",
			legend.title = element_text(face = "italic"),
			axis.text.x=element_blank())+
		labs(title = paste(i, "- Alarm Network with 100+100 npb, p<0.01")) +
		xlab("") +
		ylab("")
	pl
	# dev.copy2pdf(file = paste0(i,'.pdf'))
}

head(points)

pl = ggplot(points[points$vscore == 'TPR', ], aes(x = algorithm, y = vvals))  +
	geom_boxplot(aes(color = algorithm)) +
	facet_wrap(~ vsampl)  +
	scale_color_brewer(palette="Dark2") +
	theme(
	# axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8, face = "italic"), 
		legend.position = "bottom",
		legend.title = element_text(face = "italic"),
		axis.text.x=element_blank())+
	labs(title = "True Positive Rate (sensitivity/recall) - Alarm Network with 100+100 npb, p<0.01") +
	xlab("") +
	ylab("")
pl




# # 

# # set some settings to be used in the test
# regularization = "bic"
# boot.first.pass = 100
# boot.second.pass = 100
# test.pvalue = 1E-2

# each = 100
# values = c(100, 1000, 10000, 20000)

# res1E3 = res1E4 = res1E5 = res2E5 = NULL

# TOT = each
# dataset = NULL
# do = function(x,y){
	# cat("********** ", x, ":", y, "/", TOT, "\n")

	# # generate random data from that model, ensure no NAs 
	# repeat {
		# dataset = rbn(bnalarm, n = x)
		# if(!any(is.na(dataset))) break
	# }

	# return(perform.bootstrap.inference(
		# dataset,
		# regularization,
		# boot.first.pass,
		# boot.second.pass,
		# test.pvalue,
		# agony.binaries = agony.binaries,
		# doP = ifelse(j==1,TRUE, FALSE)))
# }

# for(j in 1:each) res1E3 = append(res1E3, list(do(100,j)))
# save(res1E3, file = 'res1E3.Rdata')

# for(j in 1:each) res1E4 = append(res1E4, list(do(1000,j)))
# save(res1E4, file = 'res1E4.Rdata')

# for(j in 1:each) res1E5 = append(res1E5, list(do(10000,j)))
# save(res1E5, file = 'res1E5.Rdata')

# for(j in 1:each) res2E5 = append(res1E3, list(do(20000,j)))
# save(res2E5, file = 'res2E5.Rdata')
	


# #### perform 

