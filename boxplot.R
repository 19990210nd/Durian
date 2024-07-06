install.packages(c("data.table", "ggplot2"))
install.packages("RColorBrewer")
# 加载 RColorBrewer 包
library(RColorBrewer)
library("ggplot2")
library("data.table")
library("eoffice")

#toc-rank
setwd("E:\\tmp\\202312\\Durian\\ongoing")
data=fread("toc.dat",header = T)
#names(data) = c("Method","Rg","Rg_x","Rg_y","Rg_z")
p = ggplot(data = data, mapping = aes(y = Druglikeness, x = Similarity, group = Method,color=Method,shape=Method)) +
  geom_point(size =5) +
  labs(y = "Druglikeness", x = "Similarity") +
  xlim(0,4) +
  ylim(4,10) +
  theme_bw() +
  #theme(element_blank())+
  theme(axis.text.x = element_text(size = 14)) +
  # 调整 y 轴标签字号为 14
  theme(axis.text.y = element_text(size = 14)) +
  # 调整 x 轴标题字号为 16
  theme(axis.title.x = element_text(size = 16)) +
  # 调整 y 轴标题字号为 16
  theme(axis.title.y = element_text(size = 16))
p
ggsave("Rank_S-D.png",width = 4,height = 4)



p = ggplot(data = data, mapping = aes(y = Druglikeness, x = Affinity, group = Method,color=Method,shape=Method)) +
  geom_point(size =5) +
  labs(y = "Druglikeness", x = "Affinity") +
  xlim(0,6) +
  ylim(4,10) +
  theme_bw() +
  #theme(element_blank())+
  theme(axis.text.x = element_text(size = 14)) +
  # 调整 y 轴标签字号为 14
  theme(axis.text.y = element_text(size = 14)) +
  # 调整 x 轴标题字号为 16
  theme(axis.title.x = element_text(size = 16)) +
  # 调整 y 轴标题字号为 16
  theme(axis.title.y = element_text(size = 16))
p
ggsave("Rank_A-D.png",width = 4,height = 4)

p = ggplot(data = data, mapping = aes(y = Affinity, x = Similarity, group = Method,color=Method,shape=Method)) +
  geom_point(size =5) +
  labs(y = "Affinity", x = "Similarity") +
  xlim(0,4) +
  ylim(0,6) +
  theme_bw() +
  #theme(element_blank())+
  theme(axis.text.x = element_text(size = 14)) +
  # 调整 y 轴标签字号为 14
  theme(axis.text.y = element_text(size = 14)) +
  # 调整 x 轴标题字号为 16
  theme(axis.title.x = element_text(size = 16)) +
  # 调整 y 轴标题字号为 16
  theme(axis.title.y = element_text(size = 16))
p
ggsave("Rank_S-A.png",width = 4,height = 4)



# withdraw! see below(final)
#drug-likeness
##violin chart
#/Volumes/tcsnnd/writing/durain/figure/abl1/druglikeness/rw
setwd("/Volumes/tcsnnd/writing/durain/figure/abl1/druglikeness/rw/png")
for (i in dir("/Volumes/tcsnnd/writing/durain/figure/abl1/druglikeness/rw/png")){
  name=sub(".csv","",i)     ###将文件名中.csv前面的名字提出来，也为了后面批量生成文件名
  cc <- paste("/Volumes/tcsnnd/writing/durain/figure/abl1/druglikeness/rw/png","/",i,sep="")   ###正确的找到文件的路径及名称
  #ori_data=fread("abl1-dock.csv",header = T)
  ori_data <- fread(cc, header = TRUE)
  data <- melt(ori_data, id.vars = names(ori_data)[0])
  names(data) = c("Method","Value")
  methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
  custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                     "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
  p<-ggplot(data, aes(x = Method, y = Value, color = Method)) +
    geom_violin(trim = FALSE, outlier.shape = NA) +
    geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
    #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
    geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
    #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
    #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    theme_bw() +
    theme(element_blank())
  ggsave(p,file=paste(name,".png",sep=""),width = 6,height = 4) 
  #p<-ggplot(data, aes(x = BindingEnergy)) + geom_density(aes(color = BindingEnergy_type)) + theme_bw()+ geom_density(aes(fill = BindingEnergy_type), alpha=0.4) +xlim(-15,-5) +ylim(0,0.7)
}
#/Volumes/tcsnnd/writing/durain/figure/abl1/druglikeness/cd/png
setwd("/Volumes/tcsnnd/writing/durain/figure/abl1/druglikeness/cd/png")
for (i in dir("/Volumes/tcsnnd/writing/durain/figure/abl1/druglikeness/cd/png")){
  name=sub(".csv","",i)     ###将文件名中.csv前面的名字提出来，也为了后面批量生成文件名
  cc <- paste("/Volumes/tcsnnd/writing/durain/figure/abl1/druglikeness/cd/png","/",i,sep="")   ###正确的找到文件的路径及名称
  #ori_data=fread("abl1-dock.csv",header = T)
  ori_data <- fread(cc, header = TRUE)
  data <- melt(ori_data, id.vars = names(ori_data)[0])
  names(data) = c("Method","Value")
  methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
  custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                     "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
  p<-ggplot(data, aes(x = Method, y = Value, color = Method)) +
    geom_violin(trim = FALSE, outlier.shape = NA) +
    geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
    #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
    geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
    #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
    #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    theme_bw() +
    theme(element_blank())
  ggsave(p,file=paste(name,".png",sep=""),width = 6,height = 4) 
  #p<-ggplot(data, aes(x = BindingEnergy)) + geom_density(aes(color = BindingEnergy_type)) + theme_bw()+ geom_density(aes(fill = BindingEnergy_type), alpha=0.4) +xlim(-15,-5) +ylim(0,0.7)
}





#similarity
setwd("/Volumes/tcsnnd/writing/durain/figure/similarity")
setwd("E:\\tmp\\202312\\Durian\\ongoing")
list.files()
library(tidyr)
# rw
data <- data.frame(
  Method = rep(c("LiGAN", "Pkt2Mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen"), each = 4),
  Measurement = rep(c("Act. Sim.", "Shape. Sim.", "Esp. Sim.", "Pharm. Sim."), times = 6),
  Value = c(0.8013,0.748,0.9275,0.7389,
            0.6953,0.313,0.919,0.2839,
            0.3653,0.3648,-0.0398,0.3226,
            0.7562,0.378,0.9471,0.3542,
            0.5984,0.0115,0.754,0.0107,
            0.6713,0.4542,0.9448,0.393),
  Std = c(0.0473,0.0242,0.0023,0.0433,
          0.1258,0.0658,0.0074,0.0638,
          0.1750,0.0822,0.073,0.0831,
          0.1152,0.0714,0.0074,0.0637,
          0.1665,0.0287,0.1791,0.0283,
          0.1471,0.0736,0.0082,0.0653)
)
method_order <- c("LiGAN", "Pkt2Mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
# 将 Method 转换为有序因子变量，并指定顺序
data$Method <- factor(data$Method, levels = method_order, ordered = TRUE)
custom_colors <- c("Act. Sim." = '#99cbeb', "Shape. Sim." = '#40826D', "Esp. Sim." = '#db6968',"Pharm. Sim." = '#606f8a')
# 示例代码，绘制均值和标准差点图
ggplot(data, aes(x = Method, y = Value, color = Measurement)) +
  geom_point(position = position_dodge(width = 0.2), size = 1) +
  geom_errorbar(aes(ymin = Value - Std, ymax = Value + Std), position = position_dodge(width = 0.2), width = 0.25) +
  labs(x = "Method", y = "Similarity") +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  theme(element_blank())
ggsave("Similarity_rw.png",width = 6,height = 4)


# cd
data <- data.frame(
  Method = rep(c("LiGAN", "Pkt2Mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen"), each = 5),
  Measurement = rep(c("Sim. Train.", "Sim. Test.", "Shape. Sim.", "Esp. Sim.", "Pharm. Sim."), times = 6),
  Value = c(
    0.1622,0.1333,0.812,0.9705,0.8177,
    0.168,0.1539,0.4085,0.9861,0.3517,
    0.1623,0.1404,0.3311,-0.0415,0.2891,
    0.2478,0.2103,0.3248,0.9847,0.3035,
    0.0833,0.0836,0.0161,0.8357,0.0162,
    0.1033,0.0961,0.5398,0.9906,0.4615
  ),
  Std = c(
    0.1166,0.1333,0.1045,0.0673,0.0452,
    0.1048,0.1539,0.0995,0.0092,0.0688,
    0.1104,0.1404,0.1032,0.0802,0.1046,
    0.1198,0.2103,0.1257,0.0153,0.0705,
    0.1265,0.0836,0.0833,0.2084,0.0405,
    0.0839,0.0961,0.0792,0.0074,0.0833
    )
)
method_order <- c("LiGAN", "Pkt2Mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
# 将 Method 转换为有序因子变量，并指定顺序
data$Method <- factor(data$Method, levels = method_order, ordered = TRUE)

custom_colors <- c("Sim. Train." = '#99cbeb', "Sim. Test." = '#ea9c9d', "Shape. Sim." = '#40826D', "Esp. Sim." = '#db6968',
                   "Pharm. Sim." = '#606f8a')
ggplot(data, aes(x = Method, y = Value, color = Measurement)) +
  geom_point(position = position_dodge(width = 0.2), size = 1) +
  geom_errorbar(aes(ymin = Value - Std, ymax = Value + Std), position = position_dodge(width = 0.2), width = 0.25) +
  labs(x = "Method", y = "Similarity") +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  theme(element_blank())
ggsave("Similarityg_cd.png",width = 6,height = 4)



#affinity
#violin chart
setwd("/Volumes/tcsnnd/writing/durain/figure/abl1")
setwd("/Volumes/ND/durian/groundtruth")
list.files()
# 绘制violin plot
# abl1-docking
ori_data=fread("abl1-dock.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Energy")
groundtruth_value = -13.1
# abl1-scoring
#ori_data2=fread("abl1-score.csv",header = T)
#data2 <- melt(ori_data2, id.vars = names(ori_data)[0])
#names(data2) = c("Method","Energy")
#groundtruth_value2 = -12.48291
#methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")

methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
ggplot(data, aes(x = Method, y = Energy, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  ylim(-20,0) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),  # 修改 x 轴字体大小
    axis.text.y = element_text(size = 12),  # 修改 y 轴字体大小
    axis.title.x = element_text(size = 14),  # 修改 x 轴标题字体大小
    axis.title.y = element_text(size = 14),   # 修改 y 轴标题字体大小
    legend.text = element_text(size = 12),  # 修改图例字体大小
    legend.title = element_text(size = 14) 
  )+
  theme(element_blank())
ggsave("Docking_abl1.png",width = 5,height = 4)

#gria2
ori_data=fread("gria2-dock.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Energy")
groundtruth_value = -7.7
methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
ggplot(data, aes(x = Method, y = Energy, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  ylim(-20,0) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),  # 修改 x 轴字体大小
    axis.text.y = element_text(size = 12),  # 修改 y 轴字体大小
    axis.title.x = element_text(size = 14),  # 修改 x 轴标题字体大小
    axis.title.y = element_text(size = 14),   # 修改 y 轴标题字体大小
    legend.text = element_text(size = 12),  # 修改图例字体大小
    legend.title = element_text(size = 14) 
  )+
  theme(element_blank())
ggsave("Docking_gria2.png",width = 5,height = 4)

#cdk2
ori_data=fread("cdk2-dock.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Energy")
groundtruth_value = -8.1
methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
ggplot(data, aes(x = Method, y = Energy, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  ylim(-20,0) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),  # 修改 x 轴字体大小
    axis.text.y = element_text(size = 12),  # 修改 y 轴字体大小
    axis.title.x = element_text(size = 14),  # 修改 x 轴标题字体大小
    axis.title.y = element_text(size = 14),   # 修改 y 轴标题字体大小
    legend.text = element_text(size = 12),  # 修改图例字体大小
    legend.title = element_text(size = 14) 
  )+
  theme(element_blank())
ggsave("Docking_cdk2.png",width = 5,height = 4)

#parp1
ori_data=fread("parp1-dock.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Energy")
groundtruth_value = -9.8
methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
ggplot(data, aes(x = Method, y = Energy, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  ylim(-20,0) +
  theme_bw() +
  theme(element_blank())
ggsave("Docking_parp1.png",width = 5,height = 4)


#docking
ggplot(data, aes(x = Method, y = Energy, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  ylim(-15,0) +
  theme_bw() +
  theme(element_blank())
ggsave("Docking_abl1.png",width = 5,height = 4)
#scoring
ggplot(data2, aes(x = Method, y = Energy, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_hline(yintercept = groundtruth_value2, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用相同的 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  ylim(-15,0) +
  theme_bw() +
  theme(element_blank())
ggsave("Scoring_abl1.png",width = 5,height = 4)

#Affinity-2
#distribution-curve
ori_data=fread("abl1-dock.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Energy")
groundtruth_value = -13.1
# abl1-scoring
ori_data2=fread("abl1-score.csv",header = T)
data2 <- melt(ori_data2, id.vars = names(ori_data)[0])
names(data2) = c("Method","Energy")
groundtruth_value2 = -12.48291
#methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
#docking
ggplot(data, aes(x = Energy, color = Method)) +
  geom_density(aes(color = Method)) +
  geom_density(aes(fill = Method, color = Method), alpha=0.2) +
  xlim(-18,0) +
  ylim(0,0.6) +
  geom_vline(xintercept = groundtruth_value, linetype = "dashed", color = "red") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(element_blank())
ggsave("Docking_abl1_distribution.png",width = 5,height = 4)
#scoring
ggplot(data2, aes(x = Energy, color = Method)) +
  geom_density(aes(color = Method)) +
  geom_density(aes(fill = Method, color = Method), alpha=0.2) +
  xlim(-18,0) +
  ylim(0,0.6) +
  geom_vline(xintercept = groundtruth_value2, linetype = "dashed", color = "red") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(element_blank())
ggsave("Scoring_abl1_distribution.png",width = 5,height = 4)





# druglikeness new
setwd("/Volumes/tcsnnd/writing/durain/figure/druglikeness")
list.files()
# 绘制violin plot
# logps
ori_data=fread("logps.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Value")
# abl1-scoring
ori_data2=fread("logps-rw.csv",header = T)
data2 <- melt(ori_data2, id.vars = names(ori_data)[0])
names(data2) = c("Method","Value")
#methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
#cd
ggplot(data, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(-10,15) +
  theme_bw() +
  theme(element_blank())
ggsave("logps_cd.png",width = 6,height = 4)
#rw
ggplot(data2, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(-10,15) +
  theme_bw() +
  theme(element_blank())
ggsave("logps_rw.png",width = 6,height = 4)

# lipinskis
list.files()
ori_data=fread("lipinskis.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Value")
# abl1-scoring
ori_data2=fread("lipinskis-rw.csv",header = T)
data2 <- melt(ori_data2, id.vars = names(ori_data)[0])
names(data2) = c("Method","Value")
#methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
#cd
ggplot(data, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(0,5) +
  theme_bw() +
  theme(element_blank())
ggsave("lipinskis_cd.png",width = 6,height = 4)
#rw
ggplot(data2, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(0,5) +
  theme_bw() +
  theme(element_blank())
ggsave("lipinskis_rw.png",width = 6,height = 4)

# mws
list.files()
ori_data=fread("mws.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Value")
# abl1-scoring
ori_data2=fread("mws-rw.csv",header = T)
data2 <- melt(ori_data2, id.vars = names(ori_data)[0])
names(data2) = c("Method","Value")
#methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
#cd
ggplot(data, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(0,800) +
  theme_bw() +
  theme(element_blank())
ggsave("mws_cd.png",width = 6,height = 4)
#rw
ggplot(data2, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(0,800) +
  theme_bw() +
  theme(element_blank())
ggsave("mws_rw.png",width = 6,height = 4)

# qeds
list.files()
ori_data=fread("qeds.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Value")
# abl1-scoring
ori_data2=fread("qeds-rw.csv",header = T)
data2 <- melt(ori_data2, id.vars = names(ori_data)[0])
names(data2) = c("Method","Value")
#methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
#cd
ggplot(data, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(0,1) +
  theme_bw() +
  theme(element_blank())
ggsave("qeds_cd.png",width = 6,height = 4)
#rw
ggplot(data2, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(0,1) +
  theme_bw() +
  theme(element_blank())
ggsave("qeds_rw.png",width = 6,height = 4)

# sas
list.files()
ori_data=fread("sas.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Value")
# abl1-scoring
ori_data2=fread("sas-rw.csv",header = T)
data2 <- melt(ori_data2, id.vars = names(ori_data)[0])
names(data2) = c("Method","Value")
#methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
#cd
ggplot(data, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(0,1) +
  theme_bw() +
  theme(element_blank())
ggsave("sas_cd.png",width = 6,height = 4)
#rw
ggplot(data2, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(0,1) +
  theme_bw() +
  theme(element_blank())
ggsave("sas_rw.png",width = 6,height = 4)


# tpsas
list.files()
ori_data=fread("tpsas.csv",header = T)
data <- melt(ori_data, id.vars = names(ori_data)[0])
names(data) = c("Method","Value")
# abl1-scoring
ori_data2=fread("tpsas-rw.csv",header = T)
data2 <- melt(ori_data2, id.vars = names(ori_data)[0])
names(data2) = c("Method","Value")
#methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
methods <- c("LiGAN", "Pkt2mol", "DiffSBDD", "SBDD", "GraphBP", "SurfGen")
custom_colors <- c("LiGAN" = '#99cbeb', "Pkt2mol" = '#40826D', "DiffSBDD" = '#db6968',
                   "SBDD" = '#4d97cd', "GraphBP" = '#606f8a', "SurfGen" = "#ea9c9d")
#cd
ggplot(data, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(0,300) +
  theme_bw() +
  theme(element_blank())
ggsave("tpsas_cd.png",width = 6,height = 4)
#rw
ggplot(data2, aes(x = Method, y = Value, color = Method)) +
  geom_violin(trim = FALSE, outlier.shape = NA) +
  geom_violin(aes(fill = Method, color = Method), alpha=0.2) +
  #geom_hline(yintercept = groundtruth_value, linetype = "dashed", color = "red") +
  geom_boxplot(width = 0.1, outlier.shape = NA, aes(fill = Method, color = Method), alpha=0.2) +
  #geom_boxplot(aes(fill = Method, color = Method), alpha=0.2) +
  #scale_color_manual(values = brewer.pal(length(methods), "Set1")) +  # 使用 Set1 配色方案
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  ylim(0,1) +
  theme_bw() +
  theme(element_blank())
ggsave("tpsas_rw.png",width = 6,height = 4)










#2bindingdb_rescore
setwd("/Volumes/ND_Solid_2/task/GPCR_classA/figure/2bindingdb_rescore/bindingdb_rescore/png")
for (i in dir("/Volumes/ND_Solid_2/task/GPCR_classA/figure/2bindingdb_rescore/bindingdb_rescore/csv")){
  name=sub(".csv","",i)     ###将文件名中.csv前面的名字提出来，也为了后面批量生成文件名
  cc <- paste("/Volumes/ND_Solid_2/task/GPCR_classA/figure/2bindingdb_rescore/bindingdb_rescore/csv","/",i,sep="")   ###正确的找到文件的路径及名称
  data <- fread(cc, header = TRUE)
  p<-ggplot(data, aes(x = BindingEnergy)) + geom_density(aes(color = BindingEnergy_type)) + theme_bw()+ geom_density(aes(fill = BindingEnergy_type), alpha=0.4) +xlim(-15,-5) +ylim(0,0.7)
  ggsave(p,file=paste(name,".png",sep=""),width = 5,height = 4) 
}



#grountruth
setwd("/Volumes/tcsnnd/task/Durian/analysis/real_world/groundtruth/durain_sbdd")
list.files()
data=fread("vina_groundtruth.csv",header = T)
p <- ggplot(data, aes(x = drd3_dock)) + 
  geom_density() +
  geom_vline(aes(xintercept = as.numeric(data[[1,2]])), color = "red", linetype = "dashed") +
  xlim(-15, 0) + ylim(0, 0.8) +
  theme_bw()
p
ggsave("sbdd_vina_drd3_dock.png",width = 5,height = 4)