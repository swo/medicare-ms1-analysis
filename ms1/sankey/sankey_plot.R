edges = read_tsv('edge.tsv') %>% as.data.frame
nodes = read_tsv('node.tsv') %>% as.data.frame

pdf(file='sankey.pdf')
riverplot::makeRiver(nodes=nodes, edges=edges) %>% plot
dev.off()