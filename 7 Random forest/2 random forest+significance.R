
library(tidyverse)
library(randomForest)
library(rfUtilities)
library(rfPermute)
Rt2<-read.csv('PLSPM.csv',row.names = 1)
set.seed(123)# 设立一个数，为了复现结果
PE6_rf <- randomForest(EMF ~ ., data= Rt2,
                       importance=TRUE,proximity=TRUE)
PE6_rf
PE6_dat <- importance(PE6_rfP, sort.by = NULL, decreasing = TRUE)
PE6_dat
PE6_dat %>%
  as_tibble(rownames = "names") %>%
  data.frame() %>%
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**",
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>%
  arrange(X.IncMSE) %>%
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>%
  ggplot(aes(x = names, y = X.IncMSE))+
  geom_bar(aes(fill = group),stat = "identity")+
  geom_text(aes(y = X.IncMSE + 1,label = label))+
  labs(x = "", y = "%IncMSE")+
  coord_flip()# 转置坐标轴
