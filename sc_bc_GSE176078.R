
#!Sc_BC_GSE176078_seurat.Rdata 后续处理
#rm(list=ls() )
rm(list=ls() )
setwd("/home/data/t100333/R/sc_bc")



#0. 加载数据----
#setwd("/home/data/t100333/R/sc_bc")
#load("~/R/sc_bc/sc_bc_GSE176078.RData")#分析数据
#load("~/R/sc_bc/data/Sc_BC_GSE176078_seurat.Rdata")



#0.1 小提琴violinplot----
#https://www.maimengkong.com/?post=1150
#install.packages( "remotes") #已经安装
remotes::install_github( "lyc-1995/MySeuratWrappers") #通过链接安装包 
library(MySeuratWrappers) 
library(ggplot2)
#需要展示的基因 features
my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
markers <- c( 'CD3D', 'S100A8', 'S100A9', 'CD79A', 'CCL5', 'NKG7')
VlnPlot(mydata, features = markers, 
        stacked=T,pt.size= 0, 
        cols = my36colors, #颜色 
        direction = "horizontal", #水平作图 
        x.lab = '', y.lab = '')+ #横纵轴不标记任何东西 
  theme(axis.text.x = element_blank, 
        axis.ticks.x = element_blank) #不显示坐标刻度 #失败

VlnPlot(object = mydata, features = 'STAT4', 
          split.by = 'celltype_major')

VlnPlot(mydata,group.by = 'celltype_major', features = "STAT4",pt.size=0)+NoLegend()
#STAT4样本表达异质性----
library(Seurat);library(ggplot2);library(ggpubr)
VlnPlot(mydata, group.by  = 'subtype',#split.by  = 'subtype',
        features = "STAT4",pt.size=0)+ # +NoLegend()
  #scale_y_log10()+
  # scale_fill_manual(values = paletteer::paletteer_d('ggsci::category20c_d3'))+
  # stat_compare_means(comparisons = list(c("TNBC","ER+"),
  #                                       c("TNBC", "HER2+"),
  #                                       c("HER2+","ER+") ),
  #                    label = "p.signif", label.x = 1.5)
  stat_compare_means()#label.y  = 5.0
#ggsave("./out_1/VIOLIN_stat4_sample_nolog.pdf",family="serif",width = 4,height = 6)# 4-5

#山脊图
RidgePlot(mydata, features = c("STAT4", "CD3D"), group.by  = 'orig.ident')


DotPlot(mydata , features =c("STAT4"),group.by = "subtype"#,
        #cols = c("skyblue", "pink")
        )+ #celltype_major orig.ident
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 

#ggsave("./out_1/dotplot_stat4__sample_3.pdf",family="serif",width = 2.5,height = 3.6)

#,"CD274" 放弃表达低，无异质性？

#0.2 寻找高变基因----
mydata <- FindVariableFeatures(object = mydata, selection.method = "vst", nfeatures = 1500)
# 10个表达变化最为剧烈的基因
top100 <- head(VariableFeatures(mydata), 100) #head(mydata$RNA@var.features,10)
top100
plot1 <- VariableFeaturePlot(mydata) # 画出表达变化的基因，从而观察其分布
plot2 <- LabelPoints(plot = plot1, points = top100, repel = TRUE)# 画出表达变化的基因，标记前10个基因
#plot1+plot2 

#归一化数据
mydata <- ScaleData(mydata)

#0.3线性降维分析----
#A. PCA：对缩放后的数据进行PCA分析，默认使用前面鉴定表达变化大的基因。使用features参数可以重新定义数据集。
mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata)) #慢
print(mydata[["pca"]], dims = 1:5, nfeatures = 5)#查看对每个主成分影响比较大的基因集
VizDimLoadings(mydata, dims = 1:2, reduction = "pca") #可视化对每个主成分影响比较大的基因集/#绘制每个PCA成分的相关基因
DimPlot(mydata, reduction = "pca",group.by = "subtype")#,split.by = 'ident' #PCA
DimHeatmap(mydata, dims = 1, cells = 500, balanced = TRUE) #主成分的热图
DimHeatmap(object = mydata, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2) #1：4是热图分成几个，ncol是分为几列。

#B.#确定数据集的"主成分个数"
#详见教程：https://cloud.tencent.com/developer/article/1841480  #https://www.jianshu.com/p/fef17a1babc2
mydata <- JackStraw(mydata, num.replicate = 100) #37min  建议跳过！
mydata <- ScoreJackStraw(mydata, dims = 1:20)
JackStrawPlot(mydata, dims = 1:20)#图像异常
ElbowPlot(mydata)

#0.4 细胞聚类----
#Seurat v3应用基于图形的聚类方法，例如KNN方法。具有相似基因表达模式的细胞之间绘制边缘，然后将他们划分为一个内联群体。
mydata <- FindNeighbors(mydata, dims = 1:10)
mydata <- FindClusters(mydata, resolution = 0.8)#默认0.8,约0.5min
#head(Idents(mydata), 5)
table(mydata$RNA_snn_res.0.8)#31

##运行非线性降维（UMAP/tSNE）
mydata <- RunUMAP(mydata, dims = 1:10)#3min
DimPlot(mydata, reduction = "umap") #28个  31 

mydata@meta.data$celltype_major <-factor(mydata@meta.data$celltype_major) #设置为factor
DimPlot(mydata, reduction = "umap")
DimPlot(mydata, reduction = "umap",group.by = "subtype",label = TRUE) #TNBC,ER2+,HER2+
DimPlot(mydata, reduction = "umap",group.by = "celltype_major",label = TRUE) #9个
DimPlot(mydata, reduction = "umap",group.by = "celltype_subset",label = TRUE)#太多
DimPlot(mydata, reduction = "umap",group.by = "celltype_minor",label = TRUE)#也太多
DimPlot(mydata, reduction = "umap",group.by = "orig.ident",label = F) #Extended Data Fig. 1
#LabelClusters(DimPlot(mydata, reduction = "umap"),id = 'ident')#同理

# umap <- data.frame(mydata@reductions$umap@cell.embeddings)
# library(ggplot2)
# ggplot(umap,aes(UMAP_2,UMAP_1))+geom_point() +coord_flip() #与文献UMAP不一致

mydata <- RunTSNE(object = mydata, dims = 1:30)       #TSNE聚类3min-失败
TSNEPlot(object = mydata, label = F, pt.size = 2)    #TSNE可视化
DimPlot(mydata, reduction = "tsne",label = TRUE) #同上,区分不开，建议umap




#0.5 寻找差异基因----
mydata.markers <- FindAllMarkers(mydata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #久！5min*30
#mydata.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#save.image("~/R/sc_bc/sc_bc_GSE176078.RData")
head(mydata.markers)

VlnPlot(object = mydata, features = c("MS4A1", "CD79A"))

#0.5.1 其他代码美化-差异基因-小提琴图----
#:https://zhuanlan.zhihu.com/p/536774748
marker_selected_1 <- mydata.markers  %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 >= 0.5 & pct.2 <= 0.6) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = 2) #64-7

source('SeuratPlot.R')
library(patchwork)
#ViolinPlot(object = mydata, groupBy = 'cluster', MarkerSelected = marker_selected_1[1:10] )

table(duplicated(marker_selected_1$gene) ) #有重复基因

#解析SeuratPlot.R代码
# (1)获取绘图数据1
MarkerSelected <-marker_selected_1[1:10,]
groupBy= "celltype_major"#"seurat_clusters" #"celltype_major"

  
plot_data = FetchData(object = mydata,#object, 
                      vars = c(MarkerSelected$gene, groupBy), 
                      slot = 'data') %>% 
  dplyr::rename(group = as.name(groupBy)) %>% 
  tidyr::pivot_longer(cols = -group, names_to = 'Feat', values_to = 'Expr')

# (2)获取绘图数据2
ident_plot = MarkerSelected %>% 
  dplyr::select(cluster, gene)

# (3)绘图
library(forcats)
figure_1 = ggplot(data = plot_data, mapping = aes(x = Expr,
                                                  y = fct_relevel(factor(x = Feat, 
                                                                     levels = MarkerSelected$gene)), 
                                                  fill = group, 
                                                  label = group)) +
  geom_violin(scale = 'width', adjust = 1, trim = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = '', times = length(x) - 2), x[length(x) - 1], '')) +
  facet_grid(cols = vars(group), scales = 'free') +
  cowplot::theme_cowplot(font_family = 'Arial') +
  scale_fill_manual(values = paletteer::paletteer_d('ggsci::category20c_d3')) + # BiocManager::install("paletteer")
  xlab('Expression Level') + 
  ylab('') +
  theme(legend.position = 'none', 
        panel.spacing = unit(x = 0, units = 'lines'),
        axis.line = element_blank(), #去除x和y轴坐标线(不包括axis tick)；
        panel.background = element_rect(fill = NA, color = 'black'),
        strip.background = element_blank(), #去除分页题头背景；
        strip.text = element_text(color = 'black', size = 10, family = 'Arial', face = 'bold'),
        axis.text.x = element_text(color = 'black', family = 'Arial', size = 11),
        axis.text.y = element_blank(),
        axis.title.x = element_text(color = 'black', family = 'Arial', size = 15),
        axis.ticks.x = element_line(color = 'black', lineend = 'round'),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(x = 0.1, units = 'cm'))

figure_2 = ggplot(data = ident_plot, aes(x = 1,
                                         y = fct_rev(factor(x = gene, levels = MarkerSelected$gene)),
                                         fill = cluster)) +
  geom_tile() +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = paletteer::paletteer_d('ggsci::category20c_d3') ) + #
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_legend(direction = 'vertical',
                             label.position = 'right',
                             title.theme = element_blank(),
                             keyheight = 0.5,
                             nrow = 2)) +
  xlab('Feature') +
  theme(legend.text = element_text(family = 'Arial', color = 'black', size = 11),
        legend.position = 'bottom',
        legend.justification = 'left',
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,05,0,0),
        panel.spacing = unit(0, 'lines'),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        plot.margin = unit(x = c(0,0,0,0), units = 'cm'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, color = 'black', family = 'Arial'),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

figure_2 + figure_1 + patchwork::plot_layout(nrow = 1, widths = c(0.03, 0.97))


#1.0 查看基因----
#1.1 DotPlot
features = c("STAT4","JAK2","CD274","IL12RB1","IL12RB2","P35") #CD212=IL12RB1  IL-12A=NKSF1/CLMF/P35

DotPlot(mydata , features = features,group.by = "celltype_minor" )+ #celltype_major
#coord_flip()+ #theme_bw()+ #去除背景，旋转图片
scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) #文字90度呈现

#ggsave("./out_1/dotplot_stat4__celltype_major.pdf",family="serif",width = 4,height = 6)#保存3,4 

#2.0 计算相关性----
#2.1 加载数据
rm(list=ls() )
setwd("/home/data/t100333/R/sc_bc")
load("~/R/sc_bc/sc_bc_GSE176078.RData")

#2.2.1 相关性-FeatureScatter
library(Seurat)
table(mydata$celltype_major)
p1=FeatureScatter(mydata,group.by="celltype_major",shuffle = T, feature1 = "STAT4", feature2 = "STAT3") 
p1
p2=FeatureScatter(mydata,group.by="celltype_major",feature1 = "STAT4", feature2 = "CD274") #0.01
p1|p2 #有相关性无P值：https://www.jianshu.com/p/9850c4d88a67

#A.取T细胞计算相关性-失败FeatureScatter
mydata[, mydata$celltype_major %in% c( "T-cells")] #Cancer Epithelial

table(Idents(mydata))#1-31
table(mydata$celltype_major)
table(mydata[, mydata$celltype_major %in% c( "T-cells")]$celltype_major)# T 35001

FeatureScatter(mydata[, mydata$celltype_major %in% c( "T-cells")]$celltype_major  ,
               #group.by="celltype_major",
               feature1 = "STAT4", feature2 = "CD274")


#B.修改代码FeatureScatter放弃太难:https://zhuanlan.zhihu.com/p/474703308
FeatureScatter()
getAnywhere(FeatureScatter)#查看代码

FeatureScatter(mydata, group.by=c("celltype_major"),
               #span=0.1, #越大越平滑，也越失去细节 #无效
               jitter=T,
               feature1 = "CD274", feature2 = "STAT4")

中
#C.取个细胞亚群计算相关性-
table(mydata$celltype_major)
mydata_t <- mydata[, mydata$celltype_major %in% c( "T-cells")] #Cancer Epithelial

#table(Idents(mydata))#1-31
#table(mydata[, mydata$celltype_major %in% c( "T-cells")]$celltype_major)# T 35001

FeatureScatter(mydata_t, group.by=c("celltype_major"),
               jitter=T,
               feature1 = "CD274", feature2 = "STAT4") 

#0:  T B Normal Epithelial; Plasmablasts 0.09;PVL 0.05; CAFS 0.03; Cancer Epithelial,Myeloid 0.01;Endothelial -0.01
#相关性太差-T相关性太低

#D.计算与各细胞群marker相关性----
library(ggplot2)

features = c("STAT4","CD274") #,"JAK2","STAT3","IL12RB1","IL12RB2"
marker0=c("EPCAM",# epithelial cells , 
          "MKI67",#proliferating cells , 
          "CD3D",#"CD3E","CD3G",#T cells (CD3D),  全是负相关~-0.03
          "CD68",#myeloid cells (CD68), 
          "MS4A1",#B cells (MS4A1), 
          "JCHAIN",#plasmablasts (JCHAIN), 
          "PECAM1",#endothelial cells (PECAM1) and 
          "PDGFRB")#mesenchymal cells (fibroblasts/perivascular-like cells; PDGFRB))  #marker？

library(ggplot2)
DotPlot(mydata , features = features,group.by = "seurat_clusters" )+ #celltype_major  celltype_major
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5))

table(mydata$orig.ident) #样本名
table(mydata$seurat_clusters)#31 clusters

FeatureScatter(mydata_t, group.by=c("celltype_major"),
               jitter=T,
               feature1 = "CD3G", feature2 = "STAT4")
#BULK RNA计算STAT4与marker相关性




#2.2.2 提取表达数据计算相关性----
#A.提取表达数据
exprSet <- mydata@assays[["RNA"]]@data #scale.data #1500-99304所有基因表达矩阵
exprSet<-as.data.frame(exprSet) #27719-99304 20G

features = c("STAT4","JAK2","STAT3","CD274","IL12RB1","IL12RB2") #提取目标基因
gene_exp <- exprSet[rownames(exprSet) %in% features,]
gene_exp <-as.data.frame(gene_exp)#5-99304
rm(exprSet)
gene_exp <-as.data.frame(t(gene_exp)) #"IL12RB2" "STAT4"   "JAK2"    "CD274"   "IL12RB1"
gene_exp$id <-rownames(gene_exp)

#匹配META数据
met <- mydata@meta.data;colnames(met) #99304-11
met$id <- rownames(met)
  #合并！

#B. ggstatsplot
boxplot(as.numeric(t(gene_exp)$STAT4))

library(ggstatsplot) #失败
ggscatterstats(data = t(gene_exp), y = STAT4, x = CD274, centrality.para = "mean", 
               margins = "both", xfill = "#CC79A7", yfill = "#009E73", 
               marginal.type = "densigram", # #类型可以换成density,boxplot,violin,densigram 
               title = "Relationship between S100A8 and S100A9")
#C. ggplot2
df1<-merge(df1,ccle_ann[,c("cell","type_refined","tcga_code","Pathology")],by="cell")#匹配mRNA的meta数据 312-6 #匹配后细胞数量减少



#3.0 计算样本中基因阳性细胞比例----
#参考研究购买教程：https://mp.weixin.qq.com/s/hR15bo8a91AhkW-cVnGV-A

genes <-c("STAT4", "CD274","JAK2","STAT3","IL12RB1","IL12RB2")

#3.0.0 PercentageFeatureSet 计算细胞中基因/基因集表达表达比例！----
#参考：http://www.360doc.com/content/22/1017/11/66626516_1052039459.shtml

library(Seurat)
p.percent <- PercentageFeatureSet(mydata, "STAT4") #FeaturePlot(mydata,"STAT4")
#colnames(mydata@meta.data)
#median(p.percent$nCount_RNA)#0值太多中位值意义？ #table(p.percent$nCount_RNA > 0)

mydata[["p.percent_STAT4"]]<- PercentageFeatureSet(mydata, "STAT4")*100#% of STAT4 #增加到seurat中
FeaturePlot(mydata,"p.percent_STAT4")                       #genes

#3.0.1 AddModuleScore基因/基因集细胞表达打分----
mydata <-AddModuleScore(mydata,features = list(c("STAT4","CD274","JAK2","STAT3","IL12RB1","IL12RB2")),
                        name = "Srps")#as.list否则变为单独每个基因在seurat中
FeaturePlot(mydata,"Srps1")  #colnames(mydata@meta.data)[12]<-"Srps"
DotPlot(mydata , features = "Srps",group.by = "celltype_major" ) #20230325 查看stats srps################
DotPlot(mydata , features = "Srps",group.by = "celltype_minor" ) 

DotPlot(mydata , features =  c("STAT4","CD274","JAK2","STAT3","IL12RB1","IL12RB2","Srps"),group.by = "celltype_minor" )+ #celltype_major  
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  scale_color_viridis(8, option = "G",alpha = 0.7)+ #library(viridis) #A G 
  #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5))
#ggsave("./out_1/dotplot_stats_srps_minor.pdf",family="serif",width = 4.5,height = 4)#保存4,3.5 


#3.1 计算样本中的表达比值----
##A.计算样本中的比值-设置矩阵----
#矩阵仅有表达值，亚组细胞请修改C.代码！
# ap<-AverageExpression(mydata,features = genes,group.by = "celltype_major",slot="data")[[1]] #6-9
# table(mydata$celltype_major)
# ap<-as.data.frame(ap)

ap<-AverageExpression(mydata,features = genes,group.by = "orig.ident",slot="data")[[1]] #6-9 
table(mydata$orig.ident)#26
ap<-as.data.frame(ap)#各样本表达数据，接396

##20230326计算样本srps，每个样本srps平均值
Srps <-data.frame(srps=mydata$Srps,sample=mydata$orig.ident )#99301-1 如何计算每个样本比值？
table(Srps$sample)

mean_srps<-data.frame()
for (i in unique(Srps$sample)){
  mean_1 <- data.frame(Srps=mean(Srps[Srps$sample==i,]$srps),id=i)
  mean_srps<-rbind(mean_srps,mean_1)
}
#合并细胞比例到mean_srps
Cellratio_srps <- merge(mean_srps,cellper,by="id")  #26-11
#相关性绘图507 line
##结束计算srps

##B计算样本中的比值-循环取值----
rowgroup = genes
colgroup = colnames(ap)
##先取每种细胞，再取每个基因，计算该细胞种类中阳性（表达>0）的细胞占所有细胞的比例百分比
#提取自己的sobj,在SOBJ中提取S4对象中的RNA矩阵
#稍等3min
for( i in rowgroup){
  for (j in colgroup){
    sobj <-subset(mydata,orig.ident==j) #orig.ident  #celltype_major
          ap[i,j] <-as.numeric(table(sobj@assays$RNA@data[i,] >0) [2]  ) / length(sobj@assays$RNA@data[i,])  #*100 #不乘以100
                               
  }
}

#NA是0，替换
ap[is.na(ap)] <-0 

##热图
pheatmap::pheatmap(ap[rownames(ap) %in% c("STAT4","CD274"),] ,
                   color = c(colorRampPalette(colors = c("skyblue","firebrick3"))(10))
                   )
#仅在细胞亚群中表达值

#美化
library(pheatmap)
p<-pheatmap(ap[rownames(ap) %in% c("STAT4","CD274"),], 
            cluster_col = T,cluster_rows  = F,
            #legend_breaks = c( -8,seq(from=-0.2, to=max(a), by=2)),
            legend_labels = c("NA",seq(from=0, to=8, by=2)),
            border = T,border_color ="black", #gaps_col=3,
            display_numbers = F,na_col = "#DDDDDD", #
            color = c(colorRampPalette(colors = c("skyblue","firebrick3"))(10)), #colorRampPalette(c("navy", "white", "firebrick3"))(50)
            cellwidth = 15,cellheight = 15,
            legend = TRUE,fontsize = 12,angle_col = 45#,filename = "1.pdf" 
)

p
save_pheatmap_pdf <- function(x, filename, width=8, height=3) {
  library(grid)
  x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.05)#修改聚类树线条宽度
  #x$gtable$grobs[[2]]$vp <- vpar(lwd = 0.05)
  
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height,family = "serif")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}#保存热图的函数
#save_pheatmap_pdf(p, "./out_1/pheatmap_ratio_sample.pdf")




##C.查看聚类细胞比例----
#table(mydata$orig.ident)#样本
#prop.table(table(mydata$celltype_major) )#亚群比例
#table(mydata$celltype_major, mydata$orig.ident)#样本-亚群比例

Cellratio <- prop.table(table(mydata$celltype_major ,#celltype_minor1206细亚群  #celltype_major 20221204
                              mydata$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <-data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
colnames(cellper)[1]<-"id"

#查看各细胞比例

colourCount = length(unique(Cellratio$Var1)) #celltype_major 9
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"))

##D.合并细胞比例和阳性比例----
ap1 <-as.data.frame(t(ap))#26-6
ap1$id <- rownames(ap1)
cellper1<-merge(cellper,ap1,by="id")#26-36  #16
#rownames(cellper1)<-cellper1$id;cellper1<-cellper1[,-1] 
#添加注释信息
cellper1$Subtype <-ifelse(cellper1$id %in% c("CID44041","CID4465","CID4495","CID44971","CID44991","CID4513","CID4515","CID4523","CID3946","CID3963"),"TNBC",
                          ifelse(cellper1$id %in% c("CID4461","CID4463","CID4471","CID4530N","CID4535","CID4040","CID3941","CID3948","CID4067","CID4290A","CID4398"),"ER+","HER2+")
                          )

# c("CID44041","CID4465","CID4495","CID44971","CID44991","CID4513","CID4515","CID4523","CID3946","CID3963")#TNBC-10
# c("CID4461","CID4463","CID4471","CID4530N","CID4535","CID4040","CID3941","CID3948","CID4067","CID4290A","CID4398")#ER+ 11
# c("CID3586","CID3921","CID45171","CID3838","CID4066") #HER2+ 5


##E.计算相关性及绘制散点图----
# 相关性分析
colnames(cellper1)
y1=Cellratio_srps$Srps#cellper1$'Cancer Epithelial' #`T-cells` CD274
x1=Cellratio_srps$`T-cells`#cellper1$STAT4

z1=lm(y1~x1)
corT1=cor.test(x1,y1)
cor1=corT1$estimate
cor1=round(cor1,3)

pvalue1=corT1$p.value
if(pvalue1<0.001){
  pval1=signif(pvalue1,4)
  pval1=format(pval1,scientific = TRUE)
}else(
  pval1=round(pvalue1,3)
)

#ggplot2绘图
library(ggplot2);library(ggrepel);library(ggExtra)
p2 <-ggplot(Cellratio_srps,#cellper1, #df1
            #aes(STAT4,cellper1$`Cancer Epithelial`,#CD274,#`T-cells`,
            aes(Srps,`T-cells`,
                color=id) #id subtype
)+
  geom_point(size = 2,alpha=0.5)+
  geom_smooth(method="lm",formula = y~x,#y ~ x,
              linetype=2,color="red",fill="#D3D3D3",size=0.5)+#拟合线 color="#6495ED"
  theme_bw(base_size = 14,base_family = "serif")+
  xlab("Srps")+#xlab("STAT4+ cell Ratio")+
  ylab("T cell Ratio")+ #"T-cell Ratio"
  #guides(color = guide_legend(title = "Patients"))+ #修改legend名称 0.1 -0.9
  annotate("text", x=-0.1, y=0.65, label=paste0("Cor = ",round(cor1,3),"\n","p = ", signif(as.numeric(pval1), 3)) , family='serif', fontface='italic',colour='red', size=5)+
  #x=0.1, y=0.12
  #scale_color_manual(values = c("#CC79A7", "#56B4E9"))+ #"#93cc82","#88c4e8" 浅蓝，浅绿 
  theme(panel.grid =element_blank() ) #,legend.position = "none"
p2 
#ggMarginal(p2,size=10,type = "histogram", xparams = list(colour = "skyblue", size = 0.2,fill = "#66CC66"),yparams = list(colour = "skyblue", size = 0.2,fill = "#FFCC66"))
#ggsave("./out_1/cor_Srps_Tratio_patients.pdf",width = 7,height = 4)#4.5-3 #74 #44  #6.5/9 5

