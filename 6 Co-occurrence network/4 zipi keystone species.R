rm(list = ls())
library(microeco)
library(igraph)

#读取 OTU 丰度表，构建 microeco 对象
otu <- read.csv('OTU.csv', row.names = 1,as.is=FALSE)
dataset <- microtable$new(otu_table = otu)

#使用 microeco 包构建一个网络
t1 <- trans_network$new(dataset = dataset, cor_method = 'spearman', filter_thres = 0.001)
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.65)

#然后用我们已有的网络的邻接矩阵将上述网络替换掉，这样就把我们的网络数据也添加至 microeco 对象里面了
adjacency_unweight <- read.delim('adjacency_otu.txt', row.names = 1, sep = '\t', check.names = FALSE)  #读取本地网络数据，邻接矩阵
g <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)  #邻接矩阵 -> igraph 的邻接列表，获得非含权的无向网络
t1$res_network <- g

#计算网络属性，并划分网络模块
t1$cal_network_attr()
t1$cal_module()

#提取各节点的拓扑属性，包括 zi 和 pi 值等
t1$get_node_table(node_roles = TRUE)
result <- t1$res_node_table
result

write.table(result, 'zi_pi.xls', sep = '\t', quote = FALSE, row.names = FALSE)

#简单作图展示
t1$plot_taxa_roles(use_type = 1)
