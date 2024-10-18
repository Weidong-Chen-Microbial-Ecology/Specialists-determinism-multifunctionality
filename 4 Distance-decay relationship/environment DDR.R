rm(list = ls())
##计算群落间物种组成相似度，根据原文描述，使用 Bray-curtis 相似度
spe <- read.csv('OTU.csv', row.names = 1)
spe <- data.frame(t(spe))

library(vegan)
#vegan 包 vegdist() 计算群落间物种组成 Bray-curtis 相异度矩阵
#并通过 1-Bray-curtis 相异度获得 Bray-curtis 相似度
comm_sim <- 1 - as.matrix(vegan::vegdist(spe, method = 'bray'))

#将矩阵转换为两两群落对应数值的数据框结构
diag(comm_sim) <- 0  #去除群落相似度矩阵中的对角线值，它们是样本的自相似度
comm_sim[upper.tri(comm_sim)] <- 0  #群落相似度矩阵是对称的，因此只选择半三角（如下三角）区域的数值即可
comm_sim <- reshape2::melt(comm_sim)
comm_sim <- subset(comm_sim, value != 0)
#读取环境组成数据
site_env <- read.csv('environment.csv',row.names = 1)
#对于各环境变量，均计算为欧氏距离测度
env_dis <- vegdist(site_env['env'], method = 'euclidean')
#将矩阵转换为两两群落对应数值的数据框结构
env_dis <- reshape2::melt(env_dis)
env_dis <- subset(env_dis, value != 0)
write.csv(comm_sim,'comm_sim.csv')
write.csv(env_dis,'env_dis.csv')
comm_sim <- read.csv('comm_sim.csv', row.names = 1)
env_dis <- read.csv('env_dis.csv', row.names = 1)
##采样点距离和群落相似度数据合并
comm_dis <- merge(comm_sim, env_dis, by = c('Var1', 'Var2'))
names(comm_dis) <- c('site1', 'site2', 'comm_sim', 'env_dis')
#lm() 拟合群落组成相似度（comm_sim）地理距离（site_dis_km）的一元线性关系
fit <- lm(comm_sim~env_dis, data = comm_dis)
summary(fit)  #展示拟合模型的简单统计

# 绘制散点图和拟合线
p <- ggplot(comm_dis, aes(env_dis, comm_sim)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE, color = 'blue') + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.position = 'none', 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Distance (km)', y = 'Bray-curtis similarity', title = 'Linear Regression')
print(p)
