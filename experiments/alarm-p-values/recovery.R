files <- list.files(pattern='*.Rdata')

bnalarm = empty.graph(colnames(alarm))
modelstring(bnalarm) = paste("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",
    "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA]", "[HRSA|ERCA:HR][ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK]", 
    "[MINV|INT:VLNG][FIO2][PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB]", "[SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS]", 
    "[VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV]", "[CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")
bnalarm = bn.fit(bnalarm, alarm)

adj.matrix = amat(bnalarm)

# source the required scripts
source(paste0(git, "src/perform.bootstrap.inference.R"))
source(paste0(git, "utils/plotter.R"))
source(paste0(git, "utils/stat.R"))

res05 = list()
res01 = list()
res001 = list()

for(i in 1:length(files))
{
	load(files[i])
		
	results05 = wrapper$results05
	results01 = wrapper$results01
	results001 = wrapper$results001

	source('plot.R')

	res05 = append(res05, list(x1))
	res01 = append(res01, list(x2))
	res001 = append(res001, list(x3))
}


sta05 = Reduce(cbind, res05)
sta01 = Reduce(cbind, res01)
sta001 = Reduce(cbind, res001)

library(vioplot)

par(mfrow = c(1,2))
vioplot(sta05['PPV', ], sta01['PPV', ], sta001['PPV', ], col = 'cornflowerblue', lty = 1, rectCol="gray",
  colMed = 'black', names = c('0.05', '0.01', '0.001'), pchMed = 15, horizontal = F)
title(main = 'PPV', xlab = 'p-value')


vioplot(sta05['TPR', ], sta01['TPR', ], sta001['TPR', ], col = 'cornflowerblue', lty = 1, rectCol="gray",
  colMed = 'black', names = c('0.05', '0.01', '0.001'), pchMed = 15, horizontal = F)
title(main = 'TPR', xlab = 'p-value')

dev.copy2pdf(file = 'p-values-violinplot.pdf')