GIULIO = TRUE
DANIELE = !GIULIO

if(GIULIO)
{
	setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/alarm-p-values/')
	git = '/Volumes/Data/Github/Bootstrap-based-Learning/'
	agony.binaries = '/Volumes/Data/Github/Bootstrap-based-Learning/agony/giulio/agony'
} 
if(DANIELE)
{
	setwd('~/Desktop/')
	git = '....'
	agony.binaries = '/Volumes/Data/Github/Bootstrap-based-Learning/agony/daniele/agony'	
}

AB05 = paste0('Agony with Bonferroni p 0.05')
AB01 = paste0('Agony with Bonferroni p 0.01')
AB001 = paste0('Agony with Bonferroni p 0.001')

# HC0 = paste0('Hill Climbing [', 0, ' restarts]')
# HC = paste0('Hill Climbing [', boot.first.pass + boot.second.pass, ' restarts]')


DEV.OFF = T
STATS = TRUE

x1 = plt(results05$agony.inference, results05$agony.inference.pvalues, "qvalues.holm", adj.matrix, AB05, stats = STATS, dev.off = DEV.OFF)
x2 = plt(results01$agony.inference, results01$agony.inference.pvalues, "qvalues.holm", adj.matrix, AB01, stats = STATS, dev.off = DEV.OFF)
x3 = plt(results001$agony.inference, results001$agony.inference.pvalues, "qvalues.holm", adj.matrix, AB001, stats = STATS, dev.off = DEV.OFF)



# # hc1 = plt(results05, NULL, 'hill.climing.no.restarts.inference', adj.matrix, HC0, stats = STATS, dev.off = DEV.OFF)
# # hc2 = plt(results05, NULL, 'hill.climing.with.restarts.inference', adj.matrix, HC, stats = STATS, dev.off = DEV.OFF)


# plt.stats(cbind(x1,x2,x3, hc1, hc2), legend.cex = .6)
# dev.copy2pdf(file = 'All comparison.pdf')
