library(vegan)
rm(list = ls())
#读入物种数据
otu <- read.delim('OTU.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))

#读取环境数据
env <- read.delim('environment.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

##根据原理一步步计算 db-RDA
#计算样方距离，以 Bray-curtis 距离为例，详情 ?vegdist
dis_bray <- vegdist(otu, method = 'bray')


#PCoA 排序，这里通过 add = TRUE校正负特征值，详情 ?cmdscale
pcoa <- cmdscale(dis_bray, k = nrow(otu) - 1, eig = TRUE, add = TRUE)

#提取 PCoA 样方得分（坐标）
pcoa_site <- pcoa$point

#db-RDA，环境变量与 PCoA 轴的多元回归
#通过 vegan 包的 RDA 函数 rda() 执行，详情 ?rda
db_rda <- rda(pcoa_site, env, scale = FALSE)

#被动拟合物种得分
v.eig <- t(otu) %*% db_rda$CCA$u/sqrt(nrow(otu) - 1)
db_rda$CCA$v <- decostand(v.eig, 'normalize', MARGIN = 2)
v.eig <- t(otu) %*% db_rda$CA$u/sqrt(nrow(otu) - 1)
db_rda$CA$v <- decostand(v.eig, 'normalize', MARGIN = 2)

#先作图展示下，详情 ?plot.cca
#样方展示为点，物种暂且展示为“+”，环境变量为向量
par(mfrow = c(1, 1))
# 调整绘图的长宽比例
par(mar = c(5, 5, 4, 2), pin = c(4, 3))

par(lwd = 1.5)  # 设置线条宽度为2
plot(db_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1, main = 'Generalist')
points(db_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('blue', 9), rep('blue3', 9), rep('blue4', 9)), cex = 0.6)
text(db_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'black', cex = 0.8)


#db-RDA 结果解读，以 capscale() 函数结果为例，简介
db_rda <- capscale(otu~., env, distance = 'bray', add = TRUE)

#查看统计结果信息，以 I 型标尺为例
db_rda.scaling1 <- summary(db_rda, scaling = 1)
db_rda.scaling1

#scores() 提取排序得分（坐标），以 I 型标尺为例，前四轴为例
#使用物种加权和计算的样方得分
db_rda_site.scaling1 <- scores(db_rda, choices = 1:4, scaling = 1, display = 'wa')        
#物种变量（响应变量）得分
db_rda_sp.scaling1 <- scores(db_rda, choices = 1:4, scaling = 1, display = 'sp')
#环境变量（解释变量）得分
db_rda_env.scaling1 <- scores(db_rda, choices = 1:4, scaling = 1, display = 'bp')

#RsquareAdj() 提取 R2，详情 ?RsquareAdj() 
r2 <- RsquareAdj(db_rda)
db_rda_noadj <- r2$r.squared #原始 R2
db_rda_adj <- r2$adj.r.squared       #校正后的 R2
#概念如上所述，如下为代码实现
#关于约束轴承载的特征值或解释率，应当在 R2 校正后重新计算
db_rda_exp_adj <- db_rda_adj * db_rda$CCA$eig/sum(db_rda$CCA$eig)
db_rda_eig_adj <- db_rda_exp_adj * db_rda$tot.chi
#置换检验
#所有约束轴的置换检验，即全局检验，基于 999 次置换，详情 ?anova.cca
db_rda_test <- anova.cca(db_rda, permutations = 999)

#各约束轴逐一检验，基于 999 次置换
db_rda_test_axis <- anova.cca(db_rda, by = 'axis', permutations = 999)

#p 值校正（Bonferroni 为例）
db_rda_test_axis$`Pr(>F)` <- p.adjust(db_rda_test_axis$`Pr(>F)`, method = 'bonferroni')

#计算方差膨胀因子，详情 ?vif.cca
vif.cca(db_rda)
#前向选择，以 ordiR2step() 的方法为例，基于 999 次置换检验，详情 ?ordiR2step
db_rda_forward_pr <- ordiR2step(capscale(otu~1, env, distance = 'bray', add = TRUE), scope = formula(db_rda), R2scope = TRUE, direction = 'forward', permutations = 999)
anova(db_rda_forward_pr)#检验模型显著性
db_rda_forward_pr$anova#检验各环境因子显著性，p<0.05进入RDA分析
#以 db_rda 和 db_rda_forward_pr 为例，简要绘制双序图比较变量选择前后结果
par(mfrow = c(1, 2))
dev.new(width = 10, height = 5)
plot(db_rda, scaling = 1, main = '原始模型，I 型标尺', display = c('wa', 'cn'))
plot(db_rda_forward_pr, scaling = 1, main = '前向选择后，I 型标尺', display = c('wa', 'cn'))
#比较选择前后校正后 R2 的差异，详情 ?RsquareAdj
#可以看到变量选择后，尽管去除了很多环境变量，但总 R2 并未损失很多
RsquareAdj(db_rda)$adj.r.squared
RsquareAdj(db_rda_forward_pr)$adj.r.squared

#所有约束轴的全局检验，999 次置换，详情 ?anova.cca
db_rda_forward_pr_test <- anova.cca(db_rda_forward_pr, permutations = 999)

#各约束轴逐一检验，999 次置换
db_rda_forward_pr_test_axis <- anova.cca(db_rda_forward_pr, by = 'axis', permutations = 999)

#p 值校正（Bonferroni 为例）
db_rda_forward_pr_test_axis$`Pr(>F)` <- p.adjust(db_rda_forward_pr_test_axis$`Pr(>F)`, method = 'bonferroni')
