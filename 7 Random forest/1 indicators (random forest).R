rm(list = ls())
#读取 OTUs 丰度表
otu <- read.table('OTU.txt', sep = '\t', row.names = 1, header = TRUE, fill = TRUE)

#过滤低丰度 OTUs 类群，它们对分类贡献度低，且影响计算效率
#30个样本，就按 OTUs 丰度的行和不小于30为准吧
otu <- otu[which(rowSums(otu) >= 30), ]

#将总数据集分为训练集（占 70%）和测试集（占 30%）
set.seed(123)
select_train <- sample(30, 30*0.7)
otu_train <- otu[select_train, ]
otu_test <- otu[-select_train, ]
#randomForest 包的随机森林
library(randomForest)

#随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
set.seed(123)
otu_train.forest <- randomForest(EMF ~ ., data = otu_train, importance = TRUE)
otu_train.forest
##关键 OTUs 识别
#查看表示每个变量（OTUs）重要性的得分
#summary(otu_train.forest)
importance_otu <- otu_train.forest$importance
head(importance_otu)

#或者使用函数 importance()
importance_otu <- data.frame(importance(otu_train.forest))
head(importance_otu)

#可以根据某种重要性的高低排个序，例如根据“Mean Decrease Accuracy”指标或者“Mean Decrease Gini”，按需选择
importance_otu <- importance_otu[order(importance_otu$X.IncMSE, decreasing = TRUE), ]
head(importance_otu)

#输出表格
write.csv(importance_otu, 'importance_otu.csv', sep = '\t', col.names = NA, quote = FALSE)
#作图展示 top30 重要的 OTUs
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)), main = 'Top 30 - variable importance')

# 选出重要性排名前30的OTU
top_otu <- head(importance_otu, 30)
top_otu
# 绘制柱状图
# 加载所需的包
library(ggplot2)
top_otu<-read.csv('M.csv')

# 加载所需的包
library(ggplot2)

# 根据 X.IncMSE 对数据进行降序排序
top_otu_sorted <- top_otu[order(-top_otu$IncNodePurity), ]

# 将 OTU 列转换为因子，并按照 X.IncMSE 的降序排列它们的水平
top_otu_sorted$OTU <- factor(top_otu_sorted$OTU, levels = unique(top_otu_sorted$OTU)[order(-top_otu_sorted$IncNodePurity)])

# 绘制纵向柱状图，最大值在顶部，数值依次递减
ggplot(top_otu_sorted, aes(x = OTU, y = IncNodePurity)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(x = "OTU", y = "X.IncMSE", title = "Top OTUs by X.IncMSE")
