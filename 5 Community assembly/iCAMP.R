rm(list = ls())
library(ape)
library(picante)
#读取 OTU 丰度表
comm = read.table("OTU.txt",header=T,row.names=1)
comm<-t(comm)
#读取 OTU 的进化树文件
tree <- read.tree('OTUs.tre')
tree$tip.label = gsub("'","",tree$tip.label)
#修剪系统发育树以仅包含群落数据集中存在的物种或具有非缺失性状数据的物种
tree <- prune.sample(comm, tree)
library(iCAMP)

#基于 Stegen et al (2013 and 2015) 的方法的群落构建分析
#更多详情 ?qpen，本示例默认以 βNTI=2 和 RC=0.95 作为不同过程的阈值，随机化 1000 次构建零分布，4 线程执行
set.seed(123)
qpen.out <- qpen(comm = comm, pd = pd, sig.bNTI = 2, sig.rc = 0.95, rand.time = 1000, nworker = 4)

qpen.out$ratio  #异质选择（Heterogeneous.Selection）、同质选择（Homogeneous.Selection）、扩散限制（Dispersal.Limitation）、同质扩散（Homogenizing.Dispersal）以及漂变（Undominated）的样本对占比
head(qpen.out$result)  #各样本对的 bMNTD、BC、bNTI、RC 值以及最重要的生态过程

#输出
write.csv(qpen.out, 'Gd-qpen.out.csv', row.names = FALSE)
