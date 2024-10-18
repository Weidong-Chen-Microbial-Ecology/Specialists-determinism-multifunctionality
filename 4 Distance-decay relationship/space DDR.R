rm(list = ls())
##根据经纬度计算采样点的地理距离
#读取采样点的地理位置数据
site <- read.delim('latlong.txt', sep = '\t', row.names = 1, check.names = FALSE)
#计算采样点间的地理距离
library(geosphere)
#geosphere 包 distm() 根据经纬度计算地理距离（默认距离单位，米）
#distm() 要求两列数据，第一列是经度，第二列是纬度
site_dis <- geosphere::distm(site[c("Lon", "Lat")])
rownames(site_dis) <- rownames(site)
colnames(site_dis) <- rownames(site)

#将采样点地理距离矩阵转换为两两对应数值的数据框结构
site_dis <- reshape2::melt(site_dis)
site_dis <- subset(site_dis, value != 0)
head(site_dis)

##计算群落间物种组成相似度，使用 Bray-curtis 相似度
spe <- read.delim('OTU.txt', sep = '\t', row.names = 1, check.names = FALSE)
spe <- data.frame(t(spe))
spe <- spe[rownames(site), ]

#vegan 包 vegdist() 计算群落间物种组成 Bray-curtis 相异度矩阵
#并通过 1-Bray-curtis 相异度获得 Bray-curtis 相似度
comm_sim <- 1 - as.matrix(vegan::vegdist(spe, method = 'bray'))

#将矩阵转换为两两群落对应数值的数据框结构
diag(comm_sim) <- 0  #去除群落相似度矩阵中的对角线值，它们是样本的自相似度
comm_sim[upper.tri(comm_sim)] <- 0  #群落相似度矩阵是对称的，因此只选择半三角（如下三角）区域的数值即可
comm_sim <- reshape2::melt(comm_sim)
comm_sim <- subset(comm_sim, value != 0)
head(comm_sim)

##采样点距离和群落相似度数据合并
comm_dis <- merge(comm_sim, site_dis, by = c('Var1', 'Var2'))
names(comm_dis) <- c('site1', 'site2', 'comm_sim', 'site_dis')

comm_dis[which(! comm_dis$edge %in% '0'),'edge'] <- '1'
head(comm_dis)


#注：site_dis 列的单位是米（m），在大尺度范围下，以千米（km）作为度量会更好一些
comm_dis$site_dis_km <- comm_dis$site_dis/1000

#lm() 拟合群落组成相似度（comm_sim）地理距离（site_dis_km）的一元线性关系
fit <- lm(comm_sim~site_dis_km, data = comm_dis)
summary(fit)  #展示拟合模型的简单统计
#例如
coefficients(fit)[1]  #获取截距
coefficients(fit)[2]  #获取 site_dis_km 的斜率

#也可以 names(summary(fit)) 后查看主要的内容项，然后从中提取，例如
summary(fit)$adj.r.squared  #校正后 R2
#ggplot2 绘制带线性拟合线的散点图
library(ggplot2)

p <- ggplot(comm_dis, aes(site_dis_km, comm_sim)) +
  geom_point(aes(color = as.character(edge))) +
  scale_color_manual(values = c('red', 'blue'), limits = c(0, 1)) +
  geom_smooth(method = 'lm', se = FALSE) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.position = 'none', plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limit = c(0, 1)) +
  labs(x = 'Distance (km)', y = 'Bray-curtis similarity', title = 'North')

# 绘制散点图和拟合线
p <- ggplot(comm_dis, aes(site_dis, comm_sim)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE, color = 'blue') + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.position = 'none', 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Distance (km)', y = 'Bray-curtis similarity', title = 'Linear Regression')
print(p)
