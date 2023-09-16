
# 分析ICI治疗(Pembrolizumab-PD1单抗)-STAT4，处理见文件夹
#https://pubmed.ncbi.nlm.nih.gov/35623341/

rm(list = ls() )
setwd("D:\\研究生科研\\生物信息学\\STAT4\\immune_response\\GSE194040_Pembrolizumab\\2.analysed_STAT4")
#加载rdata
load("D:/研究生科研/生物信息学/STAT4/immune_response/GSE194040_Pembrolizumab/2.analysed_STAT4/GSE194040_STAT4.RData")

#0. 加载数据----
#setwd("D:\\研究生科研\\生物信息学\\STAT4\\immune_response\\GSE194040_Pembrolizumab\\2.analysed_STAT4")

#load("D:/研究生科研/生物信息学/STAT4/immune_response/GSE194040_Pembrolizumab/1.refined_data/GSE194040_Pembrolizumab_refined.RData")

#1.0 提取合并数据

boxplot(exp_drug[1:100,]) #查看表达4-18
boxplot(log2(exp_drug[1:100,]) )#2-4

#1.1 提取表达数据
#提取目标基因表达矩阵
genes <-c("STAT4",
          #"CD3D","CD3E","CD3G", # T markers
          "SAMD3","TESPA1","GPR171","GVINP1","LCK","THEMIS","ARHGAP15","TRAV16","TRBV29.1","GFI1","TRAV12.3","CCR7","KLRB1"  #riskgene
          #"CCR7","IL7R","FOXP3","CXCL13","ZFP36","GZMK","IFIT1","IFNG","LAG3","KI67","AREG","FCGR3A" #
) #https://www.nature.com/articles/s41588-021-00911-1/figures/13


gene_exp <- exp_drug[rownames(exp_drug) %in% genes,] #1-69
gene_exp<- as.data.frame(t(gene_exp))
gene_exp$id <-rownames(gene_exp)#11-69

#1.2 处理临床数据
gene_meta_exp <- meta_drug
                                                     # R Responser       NR  Non-responser
gene_meta_exp$response <- ifelse(gene_meta_exp$pcr=="1","R","NR") #31 -38

gene_meta_exp <-merge(gene_meta_exp,gene_exp,by="id")#merge 69-14
str(gene_meta_exp)
table(gene_meta_exp$arm)
#分组信息
#STAT4
gene_meta_exp$STAT4_group <-ifelse(gene_meta_exp$STAT4 > median(gene_meta_exp$STAT4),"High","Low")

#rISKSCORE-匹配到9个基因
gene_meta_exp$riskscore <-gene_meta_exp$CCR7*0.2786+gene_meta_exp$GFI1*-0.333+gene_meta_exp$GPR171*-0.1014+gene_meta_exp$GVINP1*0.6031+
  gene_meta_exp$KLRB1*-0.7562+gene_meta_exp$LCK*0.0163+gene_meta_exp$SAMD3*0.0842+gene_meta_exp$TESPA1*-0.3684+
  gene_meta_exp$ARHGAP15*0.5291 #gene_meta_exp$THEMIS*0.085+

gene_meta_exp$riskscore_group <-ifelse(gene_meta_exp$riskscore > median(gene_meta_exp$riskscore),"High","Low")
table(gene_meta_exp$riskscore_group )

#2.0 绘图----
#2.1 boxplot
library(ggplot2);library(ggpubr)
ggplot(gene_meta_exp, aes(x = response, y = riskscore ) )+ # STAT4 #  取LOG2不影响显著性
  labs(y="Riskscore",x= NULL)+  
  geom_boxplot(aes(fill = response),outlier.size = 0.02,size=0.02,outlier.alpha = 0)+ 
  #geom_bar(stat = "identity")+
  stat_n_text()+ #library(EnvStats) # 显示样本数
  scale_fill_manual(values = c("skyblue","pink") )+ #"#1CB4B8", "#EB7369" "skyblue","pink"
  guides(fill = guide_legend(title = 'STAT4'))+
  theme_bw()+ mytheme +
  stat_compare_means()
# ylim(-0.4,0.65)+
# stat_compare_means(aes(group =  STAT4_group),
#                    label = "p.signif",
#                    label.y = 8.2,
#                    size=3,#0.6
#                    method = "wilcox.test",
#                    hide.ns = F)


mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.title.x = element_blank(),
                 #axis.text.x = element_text(angle = 60, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_blank() )#element_text(size= 12) 
#riskscore p=0.19; stat4 p=0.1 !!!!!!!!!!!无显著性放弃
#ggsave("riskscore_response_1.pdf",width = 2.5,height = 4.5,family="serif")

#3.0 热图----
library(dplyr)
exp_marker1 <- arrange(gene_meta_exp, response)#排序 69-26
rownames(exp_marker1)<-exp_marker1$id

library(pheatmap)
pheatmap::pheatmap(exp_marker1[c(14:23)],cluster_rows = F, 
                   cellwidth=8,cellheight = 8,
                   annotation_row=exp_marker1[c(13,24,26)]
)

p1 <-pheatmap(t(exp_marker1[c(14:23)]),cluster_cols  = F, show_colnames = F,
         border_color = "white",cellwidth=10,cellheight = 10,gaps_col= 38, #table(exp_marker1$response)
         color=colorRampPalette( c("navy", "white", "firebrick3"))(10),
         annotation_col=exp_marker1[c(13,24,26)]
)
p1
#保存热图
save_pheatmap_pdf <- function(x, filename, width, height) {
  library(grid)
  x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.05)#修改聚类树线条宽度
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height,family = "serif")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}#保存热图的函数
save_pheatmap_pdf(p1, "pheatmap_3.pdf",14,4)

#4.0 计算相关性
# code: D:\研究生科研\生物信息学\STAT4\immune_response\GSE173839_Durvalumab_her2-\1processed\GSE173839_process.R

#5.0 score STAT3/STAT4/PDL1 immune response----
#5.1 ssgsea打分----
#5.1.1 制作gmt文件

STAT4_genes <- c("JAK2","STAT3","STAT4","CD274","IL12RB1","IL12RB2")
STAT4_gene_exp <- exp_drug[rownames(exp_drug) %in% STAT4_genes,] #6

library(dbplyr)
gset <- c("Hallmark_my_geneset","NA",STAT4_genes)
gset <- gset%>% 
  as.data.frame() %>% 
  t()
#write.table(gset,file = "my_geneset_STAT4s.gmt",sep = "\t",row.names = F,col.names = F,quote = F)

#读取gmt文件
gmt <- GSEABase::getGmt( "my_geneset_STAT4s.gmt" ) #GeneSetCollection 对象
gmt

#5.1.2 ssgsea富集
library(GSVA)
library(GSEABase)


STAT4_gene_exp_matrix <-as.matrix(STAT4_gene_exp)
gsva_matrix<- gsva(
  expr = as.matrix(exp_drug),
  gset.idx.list = gmt,
  method='ssgsea',
  kcdf='Gaussian'#,abs.ranking=TRUE
  ) #结果全部NA，使用全部基因表达

#5.1.3 GSVA
gsva_matrix_1 <- gsva(as.matrix(exp_drug), gmt,
               mx.diff=FALSE, verbose=FALSE,
               parallel.sz=1)#线程数

#富集结束。整理合并数据
gsva_matrix<-as.data.frame(gsva_matrix);gsva_matrix<-as.data.frame(t(gsva_matrix ))
colnames(gsva_matrix)<-"ssGSEA";gsva_matrix$id <-rownames(gsva_matrix)

gsva_matrix_1 <-as.data.frame(gsva_matrix_1 );gsva_matrix_1 <-as.data.frame(t(gsva_matrix_1  ))
colnames(gsva_matrix_1 )<-"GSVA";gsva_matrix_1 $id <-rownames(gsva_matrix_1 )

exp_marker2 <- merge(exp_marker1,gsva_matrix,by="id")
exp_marker2 <- merge(exp_marker2,gsva_matrix_1,by="id")

#5.2 免疫响应预测----
#5.2.1 boxplot值分别柱状图
library(ggplot2);library(ggpubr)
# ggplot(exp_marker2, aes(x = response, y =  ssGSEA) )+ # GSVA
#   labs(y="ssGSEA score",x= NULL)+  
#   geom_boxplot(aes(fill = response),outlier.size = 0.02,size=0.02,outlier.alpha = 0)+ 
#   #geom_bar(stat = "identity")+
#   stat_n_text()+ #library(EnvStats) # 显示样本数
#   scale_fill_manual(values = c("skyblue","pink") )+ #"#1CB4B8", "#EB7369" "skyblue","pink"
#   guides(fill = guide_legend(title = 'STAT4'))+
#   theme_bw()+ mytheme +
#   stat_compare_means()
#ggsave("./out/boxplot_STAT4s_ssgsea.png",width = 4,height = 4)
#ssgaea值较好分开

ggplot(exp_marker2,aes(x = response, y =ssGSEA,color =response))+ #,palette = c("firebrick3","skyblue")
  labs(y="ssGSEA score",x= NULL)+
  # geom_point(alpha=0.5,size=3,
  #            position=position_jitterdodge(jitter.width = 0.45,
  #                                          jitter.height = 0,
  #                                          dodge.width = 0.8))+
  geom_boxplot(alpha=0.5,width=0.55,
               position=position_dodge(width=0.8),
               size=0.05,outlier.colour = NA)+
  # geom_violin(alpha=0.2,width=0.9,
  #             position=position_dodge(width=0.8),
  #             size=0.25)+
  #stat_n_text()+ #library(EnvStats) # 显示样本数NR 38  N 31
  scale_color_manual(values = c("#CC79A7", "#56B4E9") )+ #color
  stat_compare_means(#aes(group = Group) ,
    comparisons=list(c("NR","R")),#my_comparisons,
    label = "p.format",#"p.format p.signif 1.1e-06
    method = "wilcox",
    show.legend= F,#删除图例中的"a"
    label.x=1.5,bracket.size=0.1,vjust=0.1,
    #label.y = max(log2(GDSC2$Camptothecin_1003)),
    hide.ns = T,size=4)+
  theme_bw()+
  theme(axis.title  = element_text(size=12), axis.text = element_text(size=12),legend.position ="none",  
        panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() )+ 
  scale_x_discrete(labels= c("NR (n=38)","R (n=31)"))
#ggsave("./out/boxplot_STAT4s_ssgsea_3.pdf",width = 4,height = 4,family="serif")

#5.2.2 ROC
library(pROC)
roc_data <- roc(exp_marker2$response,exp_marker2$ssGSEA,ci=TRUE) 
roc_data[["ci"]]
roc_data[["auc"]] #AUC = 0.68
#reportROC::reportROC(gold = exp_marker2$response,predictor = exp_marker2$ssGSEA,important = "se",plot=T)

library(ggplot2)
g <- ggroc(roc_data,legacy.axes = TRUE, alpha = 1, colour = "pink",  size = 1)
g
g + #ggtitle("ROC curve") +
  #geom_ribbon(aes(x=1 - roc_data$specificities, ymin=0, ymax=roc_data$sensitivities),alpha=0.5,fill="skyblue") +#填充
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey50", linetype="dashed")+
  xlab("1-Specificity (FPR)") + ylab("Sensitivity (TPR)")+
  theme_bw()+
  annotate("text",x=0.6,y=0.2,label=paste(" AUC: ",round(roc_data[["auc"]],1) ,"\n ","95%CI :",round(roc_data[["ci"]],2)[1],"-",round(roc_data[["ci"]],2)[3] ),size=6,colour = "red",family="serif")+
  scale_x_continuous(expand = c(0,0.01)) +scale_y_continuous(expand = c(0,0.01))+ #X轴数据两边c(0,0.02),否则横坐标显示不全
  theme(text=element_text(size=16,family="serif")) +#字体TNR
  theme(plot.margin=unit(rep(1,5),'lines') )
#ggsave("./out/ROC_STAT4s_SSGSEA.pdf",family="serif",width = 4,height=4)

#5.2.3 热图
#合并表达STATs
STAT4_gene_exp_1 <-as.data.frame(t(STAT4_gene_exp))
STAT4_gene_exp_1$id <-rownames(STAT4_gene_exp_1)
exp_marker2 <-merge(STAT4_gene_exp_1,exp_marker2,by="id")
rownames(exp_marker2)<-exp_marker2$id
colnames(exp_marker2)[7]<-"STAT4"

library(dplyr)
exp_marker3 <- arrange(exp_marker2, response)#排序 69-26

library(pheatmap)
ann_colors = list(
  response=c(NR="skyblue",R="pink")
)
p2 <-pheatmap(t(scale(exp_marker3[c(2:7)])),cluster_cols  = F, show_colnames = F,
              cellwidth=10,cellheight = 10,
              gaps_col= 38, border_color = NA,#table(exp_marker1$response)
              color=colorRampPalette( c("navy", "white", "firebrick3"))(10) ,#
              annotation_col=exp_marker3[c(19)],annotation_colors = ann_colors,
)
p2
save_pheatmap_pdf(p2, "./out/pheatmap_zscore_STATs_1.pdf",14,4)#SCALE

#5.2.3 比较文章biomarker参数与ssgsea score 对免疫预测响应的差异----
#读取文章biomarker
biomarker <-openxlsx::read.xlsx("D:\\研究生科研\\生物信息学\\STAT4\\immune_response\\GSE194040_Pembrolizumab\\0.RAW_DATA\\Table S2 Patient-level biomarker scores.xlsx",rows = c( 2:1000))
colnames(biomarker)[1] <-"id"
#987-54
#匹配出exp_marker3中有的样本 #69-34
exp_marker4 <-merge(exp_marker3,subset(biomarker,biomarker$id %in% exp_marker3$id),by="id")#69-87
rownames(exp_marker4)<-exp_marker4$id
#table s1 Immune pathway mrna 7
#Module5_TcellBcell ICS5 B_cells Dendritic_cells Mast_cells STAT1_sig Chemokine12 Module3_IFN
#c(60,58,64,63,65,61,59,62) #58:64

#绘制带有注释的ssgsea得分表达热图！
exp_marker4 <- arrange(exp_marker4, response)#排序
p4 <-pheatmap(t(scale(exp_marker4[c(2:7,58:64)])),  #t(scale(exp_marker4[c(2:7)])) #STATs p3
              cluster_cols  = F, show_colnames = F,cluster_rows  = F,
              cellwidth=10,cellheight = 10,
              gaps_col= 38,gaps_row= 6, border_color = NA,#table(exp_marker1$response)
              color=colorRampPalette( c("navy", "white", "firebrick3"))(10),#alpha(rev(RColorBrewer::brewer.pal(11,"PiYG") ), 0.7),   #"BrBG" 褐浅蓝  "PiYG"紫绿
              #alpha(colorRampPalette(colors = c("orange", "white", "purple"))(10),0.6) ,#colorRampPalette( c("navy", "white", "firebrick3"))(10) ,#
              annotation_col=exp_marker4[c(19,44,85)],annotation_colors = ann_colors,
)
p4
#save_pheatmap_pdf(p3, "./out/pheatmap_zscore_STATs_anno_1.pdf",14,4)
#save_pheatmap_pdf(p4, "./out/pheatmap_zscore_STATs_biomarker_1.pdf",14,6)
#合并热图！zscore  STATs+biomarker

##A. 热图-值
pheatmap::pheatmap(exp_marker4[,58:64]) #-2 - 3 代码见上

##B. ROC 森林图
#B.1计算响应预测ROC
library(pROC)
# roc_data <- roc(exp_marker2$response,exp_marker2$ssGSEA,ci=TRUE) 
# roc_data[["ci"]];roc_data[["auc"]] 
# round(roc_data[["auc"]],1);round(roc_data[["ci"]],2)[c(1,3)] #round(roc_data[["ci"]],2)#L95 AUC H95

df_1 <-data.frame()
for (i in c(colnames(exp_marker4[c(2:7,33,58:64)] )) ){
  roc_res<- roc(exp_marker4$response,exp_marker4[,i],ci=TRUE) 
  print(i );print(round(roc_res[["ci"]],2) )
  df <-data.frame(name=i,
        L95=round(roc_res[["ci"]],2)[1],
        H95=round(roc_res[["ci"]],2)[3],
        AUC=round(roc_res[["ci"]],2)[2])
  df_1<-rbind(df_1,df)
}
df_1 #8 - 4


#B.2 绘制ROC森林图
library(ggplot2)
#df_1 <- arrange(df_1, AUC)#排序

df_1$group_col <- c( rep("#e7a40e", 6),"skyblue",#e7a40e"
                     rep("#1c6891", 7)#, rep("#a59d70", 5), rep("#4f4a30", 3)
                    )
df_1$group_col <-factor(df_1$group_col,levels = unique(df_1$group_col),ordered = T)
df_1$name <-factor(df_1$name,levels = unique(df_1$name),ordered = T)

ggplot(df_1)+
  # 0轴竖线：
  geom_hline(yintercept = 0, linewidth = 0.3)+
  # 线条：
  #geom_linerange(aes(name, ymin = L95, ymax = H95, color = name), show.legend = F)+
  geom_linerange(aes(name, ymin = 0.4, ymax = AUC), show.legend = F,color = df_1$group_col)+ #, color = group_col
  #geom_label()+
  
  # 散点：
  geom_point(aes(name, AUC, size=AUC),color = df_1$group_col)+ #color = group_col,
  # annotate("rect",
  #          xmin = c(0.5,7.5),  #7.5 8.5 top ssGSEA
  #          xmax = c(7.5,8.5),
  #          ymin = 0.4, ymax = 1, alpha = 0.2, fill = rev(unique(df_1$group_col))) + #背景填充
  annotate("text", label = df_1$AUC, x = df_1$name, y = df_1$AUC+0.1,size=3.0 )+ #,  colour =df_1$group_col
  annotate("rect",
           xmin = c(0.5,6.5,7.5),  #6.5 7.5 top ssGSEA
           xmax = c(6.5,7.5,14.5),
           ymin = 0.4, ymax = 1, alpha = 0.1, fill =unique(df_1$group_col)) + #rev(unique(df_1$group_col))
  scale_y_continuous(expand = c(0,0),limits = c(0.4,1))+
  xlab("")+
  ylab("AUC")+
  theme_bw()+
  theme(legend.position = "none",axis.text.y = element_text(color =df_1$group_col) )+ #,axis.text.y = element_text(color =df_1$name )
  coord_flip()
  
#ggsave("./out/plots_forest_auc_all_1.pdf", height = 3.5, width = 4,family="serif")
#更换名称  ssGSEA 

#6.0 单细胞解聚 免疫响应----
#将单细胞分群数据应用到此矩阵


#7.0 20230820 STAT4/CD274亚型响应曲线？----
#柱状图比例
#https://geek-docs.com/r-language/r-tutorials/g_change-y-axis-to-percentage-points-in-ggplot2-barplot-in-r.html
##7.1 STAT4-NR/R
ggplot(exp_marker4, aes(x = response, fill = STAT4_group)) +
  geom_bar(width = 0.5, position = "fill")+ # 百分比柱状图
  scale_fill_brewer(palette = "Blues") + # 调色板{RColorBrewer}
  scale_y_continuous(labels = scales::percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "STAT4")) +
  labs(title = "", x = "Response",y = "Percent") +
  theme_minimal() 
  coord_flip() # 倒转x与y轴
#ggsave("./out/percent_STAT4_response.pdf", height = 3, width = 2.5,family="serif")

#添加分组CD274-HIGH/LOW 
exp_marker5  <-exp_marker4
exp_marker5$CD274_group <- ifelse(exp_marker5$CD274 > median(exp_marker5$CD274),"High","Low")
table(exp_marker5$CD274_group)#34 35
#合并分组STAT4+CD274
exp_marker5$Group <- paste0("STAT4_",exp_marker5$STAT4_group," ","CD274_",exp_marker5$CD274_group)
table(exp_marker5$Group)

ggplot(exp_marker5, aes(x = response, fill = Group )) + #CD274_group
  geom_bar(width = 0.5, position = "fill")+ # 百分比柱状图
  scale_fill_brewer(palette = "Blues") + #Reds Blues调色板{RColorBrewer} ,direction = -1
  scale_y_continuous(labels = scales::percent) + # 修改y轴刻度为百分比
  guides(fill=guide_legend(title = "Subgroup")) + #CD274
  labs(title = "", x = "Response",y = "Percent") +
  theme_minimal() 
#ggsave("./out/percent_CD274_STAT4_response_1.pdf", height = 3, width = 3.5,family="serif")

#显示百分比ggstatsplot-边框太粗
#https://blog.csdn.net/weixin_42933967/article/details/96200319
#diamonds2 <- dplyr::filter(diamonds, color %in% c('J', 'H')) #选取diamonds数据集中color为J和H的数据，保存为diamonds2
library(ggstatsplot);library(dplyr) #加载ggstatsplot包
#使用ggstatsplot的ggbarstats函数
#diamonds2 %>% 
  ggbarstats(#diamonds2,x = color, y = clarity,
    exp_marker5,x = Group, y = response,xlab="",ylab = "Percent (%)",lwd=0.1,
             bar.proptest = T, #不显示p值的显著与否
             palette = 'Blues',#设置颜色板
             perc.k=0, #小数点
             results.subtitle = F #副标题不显示统计结果  #y区间？
  )+scale_y_continuous(labels = c(0,25,50,75,100) )  +
    theme(axis.ticks.x = element_blank(),panel.grid = element_blank()  ) 
    #scale_size_manual (values= c(0.1,2,2,2))+
   
    
    #scale_size_manual(values = c(0.1,1)) #无效 修改框线宽


  #coord_flip() #旋转坐标轴
#ggsave("./out/percent_CD274_STAT4_response_4.pdf", height = 3.5, width = 3.5,family="serif")



#20230908 STAT4 SRPS TNBC response----
table(exp_marker5$Receptor.Subtype )#TN 29
exp_marker5_tnbc <-subset(exp_marker5,exp_marker5$Receptor.Subtype=="TN")#29

plot(density(exp_marker5_tnbc$ssGSEA))

ggplot(exp_marker5_tnbc,aes(x = response, y =ssGSEA,color =response))+ #,palette = c("firebrick3","skyblue")
  labs(y="Srps",x= NULL,title = "TNBC")+ #Srps  STAT4
  geom_point(alpha=0.5,size=3,
             position=position_jitterdodge(jitter.width = 0.45,
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  # geom_boxplot(alpha=0.5,width=0.55,
  #              position=position_dodge(width=0.8),
  #              size=0.05,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.25)+
  stat_n_text()+ #library(EnvStats) # 显示样本数NR 38  N 31
  scale_color_manual(values = c("#CC79A7", "#56B4E9") )+ #color
  stat_compare_means(#aes(group = Group) ,
    comparisons=list(c("NR","R")),#my_comparisons,
    label = "p.format",#"p.format p.signif 1.1e-06
    method = "wilcox",
    show.legend= F,#删除图例中的"a"
    label.x=1.5,bracket.size=0.1,vjust=0.1,
    #label.y = max(log2(GDSC2$Camptothecin_1003)),
    hide.ns = T,size=4)+
  theme_bw()+
  theme(axis.title  = element_text(size=12), axis.text = element_text(size=12),legend.position ="none",  
        panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() ) +
  #scale_x_discrete(labels= c("NR (n=38)","R (n=31)"))
  scale_y_continuous(limits = c(-0.8,0.3))
#ggsave("./out/20230908_tnbc_Srps_violin.pdf",width = 3,height = 3,family="serif")

ggplot(exp_marker5_tnbc,aes(x = response, y =STAT4,color =response))+ #,palette = c("firebrick3","skyblue")
  labs(y="STAT4",x= NULL,title = "TNBC")+ #Srps  STAT4
  geom_point(alpha=0.5,size=3,
             position=position_jitterdodge(jitter.width = 0.45,
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  # geom_boxplot(alpha=0.5,width=0.55,
  #              position=position_dodge(width=0.8),
  #              size=0.05,outlier.colour = NA)+
  # geom_violin(alpha=0.2,width=0.9,
  #             position=position_dodge(width=0.8),
  #             size=0.25)+
  stat_n_text()+ #library(EnvStats) # 显示样本数NR 38  N 31
  scale_color_manual(values = c("#CC79A7", "#56B4E9") )+ #color
  stat_compare_means(#aes(group = Group) ,
    comparisons=list(c("NR","R")),#my_comparisons,
    label = "p.format",#"p.format p.signif 1.1e-06
    method = "wilcox",
    show.legend= F,#删除图例中的"a"
    label.x=1.5,bracket.size=0.1,vjust=0.1,
    #label.y = max(log2(GDSC2$Camptothecin_1003)),
    hide.ns = T,size=4)+
  theme_bw()+
  theme(axis.title  = element_text(size=12), axis.text = element_text(size=12),legend.position ="none",  
        panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() ) 
  #scale_x_discrete(labels= c("NR (n=38)","R (n=31)"))
  scale_y_continuous(limits = c(-0.8,0.3))
#ggsave("./out/20230908_tnbc_stat4.pdf",width = 3,height = 3,family="serif")

#roc
library(pROC)

df_2 <-data.frame()
for (i in c(colnames(exp_marker5_tnbc[c(2:7,33,58:64)] )) ){
  roc_res<- roc(exp_marker5_tnbc$response,exp_marker5_tnbc[,i],ci=TRUE) 
  print(i );print(round(roc_res[["ci"]],2) )
  df <-data.frame(name=i,
        L95=round(roc_res[["ci"]],2)[1],
        H95=round(roc_res[["ci"]],2)[3],
        AUC=round(roc_res[["ci"]],2)[2])
  df_2<-rbind(df_2,df)
}
df_2 #8 - 4  
#表现差
# STAT4 0.39 0.85 0.62
# ssGSEA 0.58 0.99 0.78
