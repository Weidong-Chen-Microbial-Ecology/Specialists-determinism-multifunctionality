
#不同的微生物网络鲁棒性评估
library(igraph)
library(ggplot2)

#读取不同微生物网络邻接矩阵，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
adj1 <- read.delim('1.txt', row.names = 1, sep = '\t')
adj2 <- read.delim('2.txt', row.names = 1, sep = '\t')
adj3 <- read.delim('3.txt', row.names = 1, sep = '\t')
adj4 <- read.delim('4.txt', row.names = 1, sep = '\t')
g1 <- graph_from_adjacency_matrix(as.matrix(adj1), mode = 'undirected', diag = FALSE)
g2 <- graph_from_adjacency_matrix(as.matrix(adj2), mode = 'undirected', diag = FALSE)
g3 <- graph_from_adjacency_matrix(as.matrix(adj3), mode = 'undirected', diag = FALSE)
g4 <- graph_from_adjacency_matrix(as.matrix(adj4), mode = 'undirected', diag = FALSE)
#计算自然连通度
g1_nc <- nc(g1)
g2_nc <- nc(g2)
g3_nc <- nc(g3)
g4_nc <- nc(g4)
dat <- data.frame(
  network = c(rep('g1', length(g1_nc)), rep('g2', length(g2_nc)), rep('g3', length(g3_nc)),rep('g4', length(g4_nc))), 
  'Proportion of removes nodes' = c((1:length(g1_nc))/length(g1_nc), (1:length(g2_nc))/length(g2_nc), (1:length(g3_nc))/length(g3_nc),(1:length(g4_nc))/length(g4_nc)),
  'Matural Connectivity' = c(g1_nc, g2_nc, g3_nc,g4_nc), check.names = FALSE
)
dat <- subset(dat, `Proportion of removes nodes` <= 1)

#作图
library(ggplot2)

library(ggplot2)

library(ggplot2)

ggplot(dat, aes(`Proportion of removes nodes`, `Matural Connectivity`, color = network, shape = network, size = network)) +
  geom_point() +
  theme_bw() +
  labs(x = 'Proportion of removes nodes', y = 'Matural Connectivity') +
  scale_shape_manual(values = c('g1' = 19, 'g2' = 19, 'g3' = 1, 'g4' = 1)) +  # g1 和 g2 设置为实心圆（shape = 19），g3 和 g4 设置为空心圆（shape = 1）
  scale_size_manual(values = c('g1' = 1, 'g2' = 1, 'g3' = 1, 'g4' = 1)) +
  scale_color_manual(values = c('g1' = "#005a9f", 'g2' = "#9eb340", 'g3' = "#005a9f", 'g4' = "#9eb340"))
