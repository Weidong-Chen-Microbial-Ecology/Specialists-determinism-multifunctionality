
#install.packages("NST")
library(NST)
comm = read.table("OTU.txt",header=T,row.names=1)
comm<-t(comm)
comm[1:6,1:6]

#group 记录了样本的分组信息
group = read.table("group.txt",header=T,row.names=1)
group

#计算 NST/MST，以下仅展示了一些重要的参数项，更多详情 ?tNST
#指定 OTU 丰度表 comm 和样本分组 group
#dist.method 用于指定相异度量，这里以 jaccard 相异指数为例
#null.model 用于选择零模型算法，这里使用默认值 PF
#rand = 1000，随机化次数，这里使用默认值 1000
#nworker 可用于指定多线程，以加快运算效率
set.seed(123)
tnst <- tNST(comm = comm, group = group, dist.method = 'jaccard', null.model = 'PF', 
             rand = 1000, nworker = 6)

#输出结果 tnst 以列表类型存储，包含了计算的两两样本对之间的 NST/MST
#names(tnst)
#tnst$ndex.pair  #所有样本对之间的 NST/MST
#tnst$index.pair.grp  #仅组内样本对之间的 NST/MST

#例如查看两组群落中，组内样本对的 NST/MST
nst_group <- tnst$index.pair.grp
nst_group

#输出主要的统计结果
write.table(nst_group, 'nst_group.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#在结果中可以看到，个别NST出现了大于1或者小于0的情况，但是MST没有，MST严格分布在[0,1]范围内，因此相较而言更推荐使用MST。

#为了更好地评估确定性和随机过程在两组群落中的相对重要性
#不妨简单绘制个箱线图，将两组群落内 MST 的分布可视化
#并通过统计检验（如这里以非参数的 wilcox test 为例），比较两组群落内 MST 数值是否存在显著差异
#library(ggplot2)
nst_group<-read.delim('D-TU.txt')
library(ggpubr)


p = ggboxplot(data = nst_group, x = 'group', y = 'NST.ij.ruzicka', color = 'group') +
  stat_compare_means(method = 'wilcox.test', comparisons = list(c('Control', 'Treat'))) +
  labs(y = 'Modified Stochasticity Ratio (MST)') +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75), limits = c(0, 0.75)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue")

print(p)
