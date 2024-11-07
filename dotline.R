
#X vs Y，获取回归线上下SD的数据
setwd("D:\\writing\\20240301_CRISPR\\GSE126486\\dotline")
library(data.table)
library(dplyr)

tmp=fread("gene_effects.csv") %>% 
  tidyr::pivot_longer(c(2:ncol(.)),names_to = 'Gene',values_to = 'Effect')
#交叉比较分析
tmp %>% 
  select(c(cell_line_name,Gene,Effect)) %>% 
  tidyr::pivot_wider(names_from = cell_line_name, values_from = c(Effect)) %>% 
  setnames(c("Gene","y","x")) -> regression_data

#计算回归方程的系数
# alphas=seq(0.1,10,by=0.1)
# absolute_erros = lapply(
#   alphas,
#   FUN = function(alpha) {
#     sum_absolute_error = regression_data %>%
#       #计算垂直距离
#       mutate(y_dist = abs(y - x * alpha)*(1/sqrt(1+alpha^2))) %>%
#       select(y_dist) %>% sum()
#   }
# )
# candidate_alpha=alphas[which.min(absolute_erros)]
# |y-kx-b|/sqrt(1+k^2)
#结果可视化，采用6倍SD，可以根据项目改动
# sigma_sd = regression_data %>% 
#   mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) ) %>% 
#   pull(y_dist) %>% sd()*6

#计算截距和斜率
linear_regression=lm(y~x,data=regression_data)
candidate_alpha=linear_regression$coefficients[2]
candidate_intercept=linear_regression$coefficients[1]
sigma_sd=regression_data %>% mutate(
  y_dist = abs(candidate_alpha*x + candidate_intercept - y)/sqrt(1+candidate_alpha^2)
) %>% pull(y_dist) %>% sd()*8

library(ggrepel)
p1=regression_data %>% 
  ggplot(aes(x=x,y=y))+
  geom_point()+
  geom_point(aes(x=x,y=y),colour='grey',
             data=regression_data %>% 
               mutate(y_dist = abs(candidate_alpha*x + candidate_intercept - y)/sqrt(1+candidate_alpha^2) )%>% 
               filter(y_dist < sigma_sd)
  )+
  geom_text_repel(aes(x=x,y=y,label=Gene),max.overlaps = 15,
                  data=regression_data %>% 
                    mutate(y_dist = abs(candidate_alpha*x + candidate_intercept - y)/sqrt(1+candidate_alpha^2) )%>% 
                    filter(y_dist > sigma_sd)
  )+
  geom_abline(slope = candidate_alpha,
              intercept = candidate_intercept,colour='red',lwd=1.2)+
  #计算上部截距
  geom_abline(slope = candidate_alpha,
              intercept = c(candidate_intercept + sigma_sd/(1/sqrt(1+candidate_alpha^2))), #sina=cos(1-a)=1/sqrt(1+k^2), tana=k,通过三角函数的关系，获取截距值
              colour='grey',lty=2,lwd=1.2)+
  #计算下部截距
  geom_abline(slope = candidate_alpha,
              intercept = c(candidate_intercept - sigma_sd/(1/sqrt(1+candidate_alpha^2))),
              colour='grey',lty=2,lwd=1.2)+
  labs(x="DMSO",y="CompB")+
  theme_bw(
    base_size = 15
  )+xlim(c(-3,1))+ylim(c(-3,1))
ggsave("comparisions_dotplot.jpg",p1,width = 8,height = 6)

#一些筛选结果
regression_data %>% 
  mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
  filter(abs(y_dist) > sigma_sd) %>% arrange(desc(abs(y_dist))) %>% write.csv(.,"results/dependency_genes.csv")
regression_data %>% 
  mutate(y_dist = (y - x * candidate_alpha)*(1/sqrt(1+candidate_alpha^2)) )%>% 
  filter(abs(y_dist) > sigma_sd) %>% 
  filter(abs(x)<0.3 & abs(y)>0.8)