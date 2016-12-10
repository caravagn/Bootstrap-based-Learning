
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
	al = new("graphAM", adjMat = al.adj, edgemode = "directed")


	# plot X	
	A = model[[X]]
	colnames(A) = colnames(alarm)
	rownames(A) = colnames(alarm)
	Abn = empty.graph(colnames(alarm))
	amat(Abn) = A

	edge = list()
	edge_a = list()

	color = list()
	lty = list()

	color_a = list()

	reversed = 0
	missing_alarm = spurious = 0

	# colors = brewer.pal(4, 'Accent')
	colors = c('firebrick2', 'orange', 'dodgerblue3', 'lightblue')

	add.alpha <- function(col, alpha=1){
	apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
	}
	
	colors[4] = add.alpha(colors[4], .3)


	for (i in 1:nrow(A)) {
		for (j in 1:ncol(A)) {
			from = rownames(A)[i]
			to = rownames(A)[j]
			e = paste(from, "~", to, sep = "")

			# Reversed edge
			if (A[i, j] == 1 && al.adj[i, j] == 0 && al.adj[j, i] == 1) {
				color = append(color, colors[1])
				names(color)[length(color)] = e

				reversed = reversed + 1
			}

			# Made up edge
			if (A[i, j] == 1 && al.adj[i, j] == 0 && al.adj[j, i] == 0) {
				color = append(color, colors[2])
				names(color)[length(color)] = e

				lty = append(lty, "dashed")
				names(lty)[length(lty)] = e

				spurious = spurious + 1

			}

			# # false neg
			# if(al.adj[i,j] == 1 && (A[i,j] == 0 && A[j,i] == 0))
# {
# color_a = append(color_a, 'blue')
# names(color_a)[length(color_a)] = e
# }
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
				color_a = append(color_a, colors[3])
				names(color_a)[length(color_a)] = e

				missing_alarm = missing_alarm + 1
			}

			# false neg reversed
			if (al.adj[i, j] == 1 && (A[i, j] == 0 && A[j, i] == 1)) {
				color_a = append(color_a, colors[1])
				names(color_a)[length(color_a)] = e
				# r_a = r_a + 1

			}

		}
	}

	edge$color = color
	edge_a$color = color_a
	# edge$lty = lty
	


	# print(edge)
	#print(edge_a$color)

	# plot first graph (alarm)	
	plot(al, attrs = list(node = list(fillcolor = colors[4])), edgeAttrs = edge_a)

	true_edges = sum(al.adj) - reversed - missing_alarm
	
	legend("topleft", title = "Alarm (cat.)",
		c(paste('', reversed), paste('', missing_alarm), 
		paste('', spurious), paste('', true_edges)), 
		cex=0.8, 
		col=c(colors[1],colors[3],colors[2], 'black'), 
		lty=1, 
		lwd = 2,
		bty='n')

	A = new("graphAM", adjMat = A, edgemode = "directed")

	
	plot(A, attrs = list(node = list(fillcolor = colors[4])), edgeAttrs = edge)
    
	# fnr = m_a/(bl+m_a)
	# fdr = m/(m+bl)
	precision = true_edges/(true_edges + reversed + spurious)
	recall = true_edges/(true_edges + reversed + missing_alarm)
	legend("topleft", title = title, 
		c(
		paste('k', nboot), 
		paste('p <', p), 
		paste(''),
		paste('pr.', round(precision, 3)), 
		paste('re.', round(recall, 3))), 
		cex=0.8, 
		bty='n')
		
	# precision[ 'Agony', paste(nboot,', ', p, sep = '')] = precision

		
	return(c(precision, recall))
}

tprecision$'100, 0.01' = 0
tprecision$'100, 0.001' = 0
tprecision$'100, 0.05' = 0
tprecision$'1000, 0.01' = 0



load("data/100-alarm-agony.Rdata")
x = plt(agony_model, "pvalues", 100, 0.01, 'Agony without MHC')
dev.copy2pdf(file='100_agony_pvalues.pdf')

plt(agony_model, 'qvalues.fdr', 100, 0.01, 'Agony with FDR')
dev.copy2pdf(file='100_agony_fdr.pdf')

plt(agony_model, 'qvalues.holm', 100, 0.01, 'Agony with Bonferroni')
dev.copy2pdf(file='100_agony_bonferroni.pdf')


hc0 = hc(alarm[1:100, ], restart = 0, score ='aic')
hc100 = hc(alarm[1:100, ], restart = 200, , score ='aic')

h = NULL
h$hc0 = amat(hc0)
h$hc100 = amat(hc100)

plt(h, 'hc0', 0, 51, 'Hill Climbing')
dev.copy2pdf(file='0_HC.pdf')
plt(h, 'hc100', 200, 51, 'Hill Climbing')
dev.copy2pdf(file='100_HC.pdf')

