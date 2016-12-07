source('Functions.R')

plotter = function(file) { 
	dir.create(file.path('./', file), showWarnings = FALSE)

	subgraph.plot(OSSPACE, TRUEMODEL$bn, scale = .2, ORDER = PORDER)
	dev.copy2pdf(file = paste(file, '/State_space_preorder.pdf', sep =''))

# subgraph.plot(SSPACE, TRUEMODEL$bn, ORDER = PORDER)

	compare.optima(SSPACE, TRUEMODEL$bn, PORDER)
	dev.copy2pdf(file = paste(file, '/All_Optima_without_preorder.pdf', sep =''))

	plot.SSPACE(layout = 'twopi')
	title('Landscape (Hill Climbing, no preorder)')
	dev.copy2pdf(file = paste(file, '/Landscape_without_preorder.pdf', sep =''))
}

# Example data
load('Xni4Sg0ahvW8.Rdata')
NAMES = paste("n", 1:N, sep = "")
plotter(file = 'Xni4Sg0ahvW8')

# Stats (ensemble)
files = list.files(pattern = ".Rdata")
length(files)

hclo = list()
boolo = list()

for(i in 1:length(files))
{
	load(files[i])

	hclo = c(hclo, nloptima(SSPACE))
	boolo = c(boolo, nloptima(OSSPACE))
		
	# true_models = append(true_models, list(TRUEMODEL$bn))
}
hclo = unlist(hclo)
boolo = unlist(boolo)

par(mfrow = c(1,2))

boxplot(hclo, main="Hill Climbing", 
   xlab="n = 4", ylab="Number of Local Optima")

boxplot(boolo, main="Bootstrap (poset)", 
   xlab="n = 4", ylab="Number of Local Optima")
dev.copy2pdf(file = 'Stat_number_of_optima.pdf')

hist(hclo)
wilcox.test(hclo, boolo)

for(i in 1:length(files))
{
	load(files[i])
	print(nloptima(OSSPACE))
	if(nloptima(OSSPACE) > 1)
		plotter(substr(files[i], 1, nchar(files[i]) - 6))
}

load('6YwCeQbdq6QM.Rdata')
subgraph.plot(SSPACE, TRUEMODEL$bn, ORDER = PORDER, file = '6YwCeQbdq6QM')
nloptima(SSPACE)
nloptima(OSSPACE)
subgraph.plot(OSSPACE, TRUEMODEL$bn, ORDER = PORDER, file = '6YwCeQbdq6QM')

# # # plot.SSPACE()
# K = 9
# par(mfrow = c(K,K))
# sapply(1:(K*K), function(i) plot.model(amat(true_models[[i]])))
# # title("Some true models used for the test", outer=TRUE)
# dev.copy2pdf(file = 'TRUEMODELS.pdf')



