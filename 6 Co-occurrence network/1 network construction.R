
rm(list = ls())
library(psych)
otu_table <- read.csv('OTU.csv', row.names = 1,as.is=FALSE)
# 建议使用WGCNA里的corAndPvalue来计算相关性
library(WGCNA)
library(igraph)
#也可以过滤一下数据，如将asv丰度低于0.05的过滤掉
rel_abundance = apply(otu_table, 2, function(x) x/sum(x))  

rel_abundance

mean_rel_abundance = rowMeans(rel_abundance)    

low_rel_abundance_otu = rownames(otu_table)[mean_rel_abundance < 0.0005] 

otu_table_filtered = otu_table[!(rownames(otu_table) %in% low_rel_abundance_otu), ]  
#过滤出现频率小于80%的OTU
freq = apply(otu_table_filtered, 1, function(x) sum(x > 0)/length(x))

keep = freq >= 4/5 

otu = otu_table_filtered[keep, ] 
nasv <- t(otu)
write.csv(otu, 'otu.csv')
#再进行相关系数计算

MM_cor <- corAndPvalue(nasv,method ="spearman")

#阈值筛选

r <- MM_cor$cor
p <- MM_cor$p


#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- p.adjust(p, method = 'fdr')    #可选 p 值校正，这里使用 BH 法校正 p 值

r[p>0.05|abs(r)<0.60] = 0

diag(r) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
#构建含权的无向网络，权重代表了asv间的 spearman 相关系数
g <- graph.adjacency(r, weighted = TRUE, mode = 'undirected')
#去除自相关
g <- simplify(g)

#删除的孤立节点（即度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#将相关系数复制一列并取其绝对值为权重
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

#简单绘制网络图
plot(g)

#导出成graphml 格式，使用 gephi 软件打开并进行可视化编辑
write.graph(g, 'net.graphml', format = 'graphml')
#获得邻接矩阵
adjacency_otu <- data.frame(get.adjacency(g,sparse=FALSE))#获得邻接矩阵
#保存邻近矩阵
write.csv(adjacency_otu, 'adjacency_otu.csv')
######计算网络常用的几种拓扑系数#####
nodes_num = length(V(g))                   #节点数
nodes_num

edges_num = length(E(g))                   #边数
edges_num

positive.cor_num = sum(E(g)$corr>0)        #正相关的数量
positive.cor_num

negative.cor_num = sum(E(g)$corr<0)        #负相关的数量
negative.cor_num

average_degree = mean(degree(g))           #平均度
average_degree

average_path_length = average.path.length(g, directed = FALSE)     #平均路径长度
average_path_length

network_diameter = diameter(g, directed = FALSE)                   #网络直径
network_diameter

network_density = graph.density(g)                                 #网络密度
network_density

clustering_coefficient = transitivity(g)                           #聚类系数
clustering_coefficient

network_parameter = data.frame(nodes_num, 
                               edges_num, 
                               positive.cor_num, 
                               negative.cor_num, 
                               average_degree,
                               average_path_length,
                               network_diameter, 
                               network_density,
                               clustering_coefficient                               
)

network_parameter
write.csv(network_parameter, 'network_parameter.csv')                                  

otu1 = otu
otu1[otu1>0] = 1

write.csv(otu1, 'adjacent_matrix.csv')