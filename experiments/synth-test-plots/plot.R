GIULIO = TRUE
DANIELE = !GIULIO

if(GIULIO)
{
	# Giulio
	setwd('/Volumes/Data/Github/Bootstrap-based-Learning/experiments/synth-test-plots')
	git = '/Volumes/Data/Github/Bootstrap-based-Learning/'
} 
if(DANIELE)
{
	setwd('~/Desktop/')
	git = '....'
}

# load the required R packages
library(bnlearn)
library(ggplot2)
library(igraph)

# data
load(paste(getwd(), '/consensus_poset_100_continuous.Rdata', sep =''))
names(consensus_poset_100_binary)

get.Res = function(results, type = 'discrete', nnodes = '10', density = '0.4', reg = 'results.bic' )
{
	x = results[[score]]
	x = x[[type]]
	x = x[[nnodes]]
	x = x[[density]]
	x = x[[reg]]
	return(x)
}

get.TSeries = function(x, algo, samples, noise)
{
	x = x[[algo]]	
	d = unlist(x[samples, noise])
	return(data.frame(score = d, algorithm = algo, samples = samples))
}

score = 'recall'
results = consensus_poset_100_continuous
reg = 'results.bic'
x = get.Res(results, type = 'continuous', nnodes = '15', density = '0.8', reg = reg)

algos = c('agony_adj_matrix_p', 'agony_adj_matrix_q_fdr', 'agony_adj_matrix_q_holm',
	'confidence_adj_matrix_p', 'confidence_adj_matrix_q_fdr', 'confidence_adj_matrix_q_holm',
	'hc_no_restart', 'hc_with_restart')

noise = '0.2'
points = NULL
for(i in algos) points = rbind(points, get.TSeries(x, i, samples='10', noise = noise))
for(i in algos) points = rbind(points, get.TSeries(x, i, samples='50', noise = noise))

# points[points$algorithm == 'agony_adj_matrix_p', ] = 'Agony without MHC'
 
pl = ggplot(points, aes(x = algorithm, y = score))  +
		geom_boxplot(aes(color = algorithm)) +
		facet_wrap(~ samples)  +
		scale_color_brewer(palette="Dark2") +
		xlab("") +
		ylab(score) +
		theme(
			# axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8, face = "italic"), 
			# legend = FALSE,
			legend.position = "bottom",
			legend.title = element_text(face = "italic"),
			axis.text.x=element_blank())+
		labs(title = paste(" Synthetic data (as Figure 3) with continuous variables; noise",noise)) 
print(pl)						
		
#######
get.TSeries = function(x, algo, samples, noise, density)
{
	x = x[[algo]]	
	d = unlist(x[samples, noise])
	return(data.frame(score = d, algorithm = algo, samples = samples, density = density))
}

score = 'density'
results = consensus_poset_100_binary
x = get.Res(results, nnodes = '15', density = '0.4')
y = get.Res(results, nnodes = '15', density = '0.6')
z = get.Res(results, nnodes = '15', density = '0.8')

algos = c('agony_adj_matrix_p', 'agony_adj_matrix_q_fdr', 'agony_adj_matrix_q_holm',
	'confidence_adj_matrix_p', 'confidence_adj_matrix_q_fdr', 'confidence_adj_matrix_q_holm',
	'hc_no_restart', 'hc_with_restart')

noise = '0.2'
points = NULL
for(i in algos) points = rbind(points, get.TSeries(x, i, samples='50', noise = noise, density = '0.4'))
for(i in algos) points = rbind(points, get.TSeries(y, i, samples='50', noise = noise, density = '0.6'))
for(i in algos) points = rbind(points, get.TSeries(z, i, samples='50', noise = noise, density = '0.8'))


pl = ggplot(points, aes(x = algorithm, y = score))  +
		geom_boxplot(aes(color = algorithm)) +
		facet_wrap(~ density)  +
		scale_color_brewer(palette="Dark2") +
		xlab("") +
		ylab(score) +
		theme(
			# axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8, face = "italic"), 
			legend.position = "bottom",
			legend.title = element_text(face = "italic"),
			axis.text.x=element_blank())+
		labs(title = paste(" Synthetic data: 100+100 npb, p<0.05, 750 samples, noise",noise)) 
print(pl)						
		
###########
get.TSeries = function(x, algo, samples, noise)
{
	x = x[[algo]]	
	d = unlist(x[samples, noise])
	return(data.frame(score = d, algorithm = algo, samples = samples))
}

load(paste(getwd(), '/consensus_poset_100_binary_v2.RData', sep =''))
names(consensus_poset_100_binary_v2)
names(consensus_poset_100_binary_v2$poset_comparison$discrete$'15'$'0.8'$'results.bic')
results = consensus_poset_100_binary_v2

score = 'poset_comparison'
x = get.Res(results, nnodes = '15', density = '0.8')

algos = c('in_agony_not_confidence', 'in_confidence_not_agony')
points = NULL
for(i in algos) points = rbind(points, get.TSeries(x, i, samples='10', noise = '0.2'))
for(i in algos) points = rbind(points, get.TSeries(x, i, samples='50', noise = '0.2'))
points

pl = ggplot(points, aes(x = algorithm, y = score))  +
		geom_boxplot(aes(color = algorithm)) +
		facet_wrap(~ samples)  +
		# scale_fill_manual(values=c("black", "black", "#66CC99")) +
		scale_color_brewer(palette="Set1") +
		xlab("") +
		ylab(score) +
		theme(
			# axis.text.x=element_text(color = "black", size=7, angle=30, vjust=.8, hjust=0.8, face = "italic"), 
			legend.position = "bottom",
			legend.title = element_text(face = "italic"),
			axis.text.x=element_blank())+
		labs(title = paste(" Synthetic data: 100+100 npb, p<0.05, 750 samples, noise", noise)) 
print(pl)						





