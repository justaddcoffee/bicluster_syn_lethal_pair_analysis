library(igraph)

source("chords.R")

links <- data.frame(pair.data$gene1, pair.data$gene2)
colnames(links) <- c("source", "target")

network <- graph_from_data_frame(d=links, directed=F)

heatmap(tdata)

rev(sort(degree(network)))
hist(rev(sort(degree(network))), breaks=1000, xlim=c(0,10))
