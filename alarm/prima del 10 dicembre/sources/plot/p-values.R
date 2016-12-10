
library(Rgraphviz)
library(bnlearn)
library(RColorBrewer)


plt = function(model, X, nboot, p, title) {
	data(alarm)

	res = empty.graph(names(alarm))
	modelstring(res) = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]", 
		"[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]", "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]", 
		"[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]", 
		"[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]", "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]", 
		"[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")


	par(mfrow = c(1, 2))

	al.adj = amat(res)


	# plot X	
	A = model[[X]]
	colnames(A) = colnames(alarm)
	rownames(A) = colnames(alarm)

	reversed = 0
	missing_alarm = spurious = 0


	for (i in 1:nrow(A)) {
		for (j in 1:ncol(A)) {
			from = rownames(A)[i]
			to = rownames(A)[j]
			e = paste(from, "~", to, sep = "")

			# Reversed edge
			if (A[i, j] == 1 && al.adj[i, j] == 0 && al.adj[j, i] == 1) {
				reversed = reversed + 1
			}

			# Made up edge
			if (A[i, j] == 1 && al.adj[i, j] == 0 && al.adj[j, i] == 0) {
				spurious = spurious + 1
			}
		}
	}

	tn = 0
	
	for (i in 1:nrow(al.adj)) {
		for (j in 1:ncol(al.adj)) {
			from = rownames(al.adj)[i]
			to = rownames(al.adj)[j]
			e = paste(from, "~", to, sep = "")
			
			if(al.adj[i, j] == 0 && A[i, j] == 0) tn = tn+1

			# false neg
			if (al.adj[i, j] == 1 && (A[i, j] == 0 && A[j, i] == 0)) {
				missing_alarm = missing_alarm + 1
			}
		}
	}

	true_edges = sum(al.adj) - reversed - missing_alarm
	
	# fnr = m_a/(bl+m_a)
	# fdr = m/(m+bl)
	precision = true_edges/(true_edges + reversed + spurious)
	recall = true_edges/(true_edges + reversed + missing_alarm)
	
	return(c(precision, recall))
}


rows = c('Agony without MHC', 'Agony with FDR', 'Agony with Bonferroni',
	'Confidence without MHC', 'Confidence with FDR', 'Confidence with Bonferroni',
	'Hill Climbing (k=0)', 'Hill Climbing (k=200)', 'Hill Climbing (k=1100)')

pvaltprecision = data.frame(row.names = rows, stringsAsFactors = FALSE)
pvaltrecall = data.frame(row.names = rows, stringsAsFactors = FALSE)


pvaltprecision$'0.001' = 0
pvaltprecision$'0.01' = 0
pvaltprecision$'0.05' = 0

pvaltrecall $'0.001' = 0
pvaltrecall $'0.01' = 0
pvaltrecall $'0.05' = 0


load("data/100-alarm-agony.Rdata")
x = plt(agony_model, "pvalues", 100, 0.01, 'Agony without MHC')
pvaltprecision[ 'Agony without MHC', '0.01'] = x[1]
pvaltrecall[ 'Agony without MHC', '0.01'] = x[2]

x = plt(agony_model, 'qvalues.fdr', 100, 0.01, 'Agony with FDR')
pvaltprecision[ 'Agony with FDR', '0.01'] = x[1]
pvaltrecall[ 'Agony with FDR', '0.01'] = x[2]

x = plt(agony_model, 'qvalues.holm', 100, 0.01, 'Agony with Bonferroni')
pvaltprecision[ 'Agony with Bonferroni', '0.01'] = x[1]
pvaltrecall[ 'Agony with Bonferroni', '0.01'] = x[2]

load('data/100-0.001-alarm-agony.Rdata')
x = plt(agony_model, 'pvalues', 100, 0.001,'Agony without MHC')
pvaltprecision[ 'Agony without MHC', '0.001'] = x[1]
pvaltrecall[ 'Agony without MHC', '0.001'] = x[2]

x = plt(agony_model, 'qvalues.fdr', 100, 0.001,'Agony with FDR')
pvaltprecision[ 'Agony with FDR', '0.001'] = x[1]
pvaltrecall[ 'Agony with FDR', '0.001'] = x[2]

x = plt(agony_model, 'qvalues.holm', 100, 0.001,'Agony with Bonferroni')
pvaltprecision[ 'Agony with Bonferroni', '0.001'] = x[1]
pvaltrecall[ 'Agony with Bonferroni', '0.001'] = x[2]


load('data/100-0.05-alarm-agony.Rdata')
x = plt(agony_model, 'pvalues', 100, 0.05,'Agony without MHC')
pvaltprecision[ 'Agony without MHC', '0.05'] = x[1]
pvaltrecall[ 'Agony without MHC', '0.05'] = x[2]

x = plt(agony_model, 'qvalues.fdr',100, 0.05, 'Agony with FDR')
pvaltprecision[ 'Agony with FDR', '0.05'] = x[1]
pvaltrecall[ 'Agony with FDR', '0.05'] = x[2]

x = plt(agony_model, 'qvalues.holm',100, 0.05, 'Agony with Bonferroni')
pvaltprecision[ 'Agony with Bonferroni', '0.05'] = x[1]
pvaltrecall[ 'Agony with Bonferroni', '0.05'] = x[2]


load("data/100-alarm-confidence.Rdata")
x = plt(confidence_model, "pvalues", 100, 0.01, 'Confidence without MHC')
pvaltprecision[ 'Confidence without MHC', '0.01'] = x[1]
pvaltrecall[ 'Confidence without MHC', '0.01'] = x[2]

x = plt(confidence_model, 'qvalues.fdr', 100, 0.01, 'Confidence with FDR')
pvaltprecision[ 'Confidence with FDR', '0.01'] = x[1]
pvaltrecall[ 'Confidence with FDR', '0.01'] = x[2]

x = plt(confidence_model, 'qvalues.holm', 100, 0.01, 'Confidence with Bonferroni')
pvaltprecision[ 'Confidence with Bonferroni', '0.01'] = x[1]
pvaltrecall[ 'Confidence with Bonferroni', '0.01'] = x[2]

load('data/100-0.001-alarm-confidence.Rdata')
x = plt(confidence_model, 'pvalues', 100, 0.001,'Confidence without MHC')
pvaltprecision[ 'Confidence without MHC', '0.001'] = x[1]
pvaltrecall[ 'Confidence without MHC', '0.001'] = x[2]

x = plt(confidence_model, 'qvalues.fdr', 100, 0.001,'Confidence with FDR')
pvaltprecision[ 'Confidence with FDR', '0.001'] = x[1]
pvaltrecall[ 'Confidence with FDR', '0.001'] = x[2]

x = plt(confidence_model, 'qvalues.holm', 100, 0.001,'Confidence with Bonferroni')
pvaltprecision[ 'Confidence with Bonferroni', '0.001'] = x[1]
pvaltrecall[ 'Confidence with Bonferroni', '0.001'] = x[2]

load('data/100-0.05-alarm-confidence.Rdata')
x = plt(confidence_model, 'pvalues', 100, 0.05,'Confidence without MHC')
pvaltprecision[ 'Confidence without MHC', '0.05'] = x[1]
pvaltrecall[ 'Confidence without MHC', '0.05'] = x[2]

x = plt(confidence_model, 'qvalues.fdr',100, 0.05, 'Confidence with FDR')
pvaltprecision[ 'Confidence with FDR', '0.05'] = x[1]
pvaltrecall[ 'Confidence with FDR', '0.05'] = x[2]

x = plt(confidence_model, 'qvalues.holm',100, 0.05, 'Confidence with Bonferroni')
pvaltprecision[ 'Confidence with Bonferroni', '0.05'] = x[1]
pvaltrecall[ 'Confidence with Bonferroni', '0.05'] = x[2]

alarm = alarm[1:100, ]

hc0 = hc(alarm, restart = 0, score ='bic')
hc100 = hc(alarm, restart = 200, , score ='bic')
hc1000 = hc(alarm, restart = 1100, , score ='bic')

h = NULL
h$hc0 = amat(hc0)
h$hc100 = amat(hc100)
h$hc1000 = amat(hc1000)

x = plt(h, 'hc0', 0, 51, 'Hill Climbing')
pvaltprecision[ 'Hill Climbing (k=0)', ] = x[1]
pvaltrecall[ 'Hill Climbing (k=0)', ] = x[2]

x = plt(h, 'hc100', 0, 51, 'Hill Climbing')
pvaltprecision[ 'Hill Climbing (k=200)', ] = x[1]
pvaltrecall[ 'Hill Climbing (k=200)', ] = x[2]

x = plt(h, 'hc1000', 0, 51, 'Hill Climbing')
pvaltprecision[ 'Hill Climbing (k=1100)', ] = x[1]
pvaltrecall[ 'Hill Climbing (k=1100)', ] = x[2]





exploded = data.frame(stringsAsFactors = TRUE)
for(i in 1:nrow(pvaltprecision))
{	for(j in 1:ncol(pvaltprecision))
	{
		tmp = data.frame(colnames(pvaltprecision)[j], pvaltprecision[i,j], rownames(pvaltprecision)[i])
		colnames(tmp) = c('p', 'Precision', 'Algorithm')
		exploded = rbind(exploded, tmp)
	}
}
colnames(exploded) = c('p', 'Precision', 'Algorithm')
rownames(exploded) = NULL
str(exploded)

library(ggplot2)

colourCount = length(unique(exploded $Algorithm))
colors = colorRampPalette(brewer.pal(9, "Set1"))


pl = ggplot(data= exploded, aes(x=p, y=Precision, group=Algorithm)) +
  geom_line(aes(color = Algorithm))+
  geom_point()
pl
pl + scale_color_brewer(palette="Set1") + theme_classic() + xlab('p-value')
dev.copy2pdf(file='p-value-change-precision.pdf')

exploded = data.frame(stringsAsFactors = TRUE)
for(i in 1:nrow(pvaltrecall))
{	for(j in 1:ncol(pvaltrecall))
	{
		tmp = data.frame(colnames(pvaltrecall)[j], pvaltrecall[i,j], rownames(pvaltrecall)[i])
		colnames(tmp) = c('p', 'Recall', 'Algorithm')
		exploded = rbind(exploded, tmp)
	}
}
colnames(exploded) = c('p', 'Recall', 'Algorithm')
rownames(exploded) = NULL
str(exploded)

library(ggplot2)

colourCount = length(unique(exploded $Algorithm))
colors = colorRampPalette(brewer.pal(9, "Set1"))


pl = ggplot(data= exploded, aes(x=p, y=Recall, group=Algorithm)) +
  geom_line(aes(color = Algorithm))+
  geom_point()
pl
pl + scale_color_brewer(palette="Set1") + theme_classic() + xlab('p-value')
dev.copy2pdf(file='p-value-change-recall.pdf')

