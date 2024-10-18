
rm(list = ls())
#安装 plspm 包

#加载 plspm 包
library(plspm)

#读取数据
dat <- read.csv('data.csv', row.names = 1,as.is=FALSE)
dat
#指定潜变量，在 R 中以列表（list）存储变量和潜变量的关系
#您可以直接指定列名称
dat_blocks <- list(
  Environment = c('TC','TN','NO2.N','NH4.N'), 
  Space='PCNM1',
  NST = 'NST',
  EMF = 'EMF',
  Sbeta = 'Sbeta',
  Gbeta = 'Gbeta'
)

dat_blocks
#通过 0-1 矩阵描述潜变量之间的关联，其中 0 代表变量间没有关联，1 代表有关联
Space <- c( 0, 0, 0, 0, 0,0)
Environment<- c(1, 0, 0, 0, 0,0)
NST<- c(1,1,0,0,0,0)
Gbeta <- c(1, 1, 1,0, 0,0)
Sbeta <- c(1, 1,1,0, 0,0)
EMF <- c(1,1,1,1,1,0)

dat_path <- rbind( Space,Environment,NST,Gbeta, Sbeta,EMF)
colnames(dat_path) <- rownames(dat_path)
dat_path

#指定因果关系，可选 A（代表列是行的因） 或 B（代表行是列的因）
dat_modes <- rep('A',6)
dat_modes
#一个简单的 PLS-PM，更多参数详情 ?plspm
dat_pls <- plspm(dat, dat_path, dat_blocks, modes = dat_modes)
dat_pls
summary(dat_pls)

#查看路径系数的参数估计值，以及相关的统计信息
dat_pls$path_coefs
dat_pls$inner_model

#查看因果关系的路径图，详情 ?innerplot
innerplot(dat_pls, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray', box.lwd = 0)
#查看变量间的直接或间接影响状态
dat_pls$effects
#goodness-of-fit 值可以帮助评估模型优度
dat_pls$gof
#查看路径显著性p值,可在ppt或AI进行美化和显著性标记
dat_pls$inner_model
