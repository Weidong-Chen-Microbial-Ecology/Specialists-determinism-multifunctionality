#清空工作区
rm(list = ls())
#读取网络邻接矩阵，物种相对丰度表
g <- read.delim('OTU.txt', row.names = 1, sep = '\t', check.names = FALSE)
spe <- read.delim('OTU_abundance.txt', row.names = 1, sep = '\t', check.names = FALSE)

#计算各物种的正、负连通性
r_pos <- c()
r_neg <- c()
spe_i <- c()

for (i in names(g)) {
  co <- na.omit(g[[i]])
  r_pos <- c(r_pos, mean(co[co>0]))
  r_neg <- c(r_neg, mean(co[co<0]))
  spe_i <- c(spe_i, i)
}
r <- data.frame(spe_i, r_pos, r_neg)
r
# 将NaN值赋值为0
r$r_neg[is.nan(r$r_neg)] <- 0

#计算各样本的正、负凝聚力
C_pos <- c()
C_neg <- c()
sample <- c()

for (j in names(spe)) {
  C_pos <- c(C_pos, sum(spe[[j]]*r$r_pos))
  C_neg <- c(C_neg, sum(spe[[j]]*r$r_neg))
  sample <- c(sample, j)
}
C <- data.frame(sample, C_pos, C_neg)
C

#输出
write.table(r, '物种连通性.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.csv(C, '凝聚力.csv')
