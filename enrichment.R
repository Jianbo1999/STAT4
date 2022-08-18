#对riskscore得分进行差异分析-高表达基因LASSO-COX模型
#rm(list=ls())

#1.读取数据+处理----
setwd("D:\\研究生科研\\生物信息学\\STAT4\\Subgroup_STAT4\\3.Enrichment")

clinical <- read.csv("clinical_sup3.csv")[,-1] #来自D:\研究生科研\生物信息学\STAT4\Subgroup_STAT4\2.Prognosis  #clinical_sup3.csv
clinical <-clinical[,-c(139:336)]#删除139-336 之前的差异基因，1001-174  ;142-154 lassogene
clinical <-clinical[,-c(96:138)]#删除之前的GSVA #1001-131
#后面加入数据后再保存！

exp1 <-read.csv("D:\\研究生科研\\生物信息学\\STAT4\\Subgroup_STAT4\\1.Riskscore\\exp1.csv",row.names = 1)
colnames(exp1) <-gsub("\\.","-",colnames(exp1)) #1067
colnames(exp1)
#取1001样本表达矩阵
exp1 <-exp1[,colnames(exp1) %in% clinical$sample]#1001

#未保存矩阵exp1

#2. 差异分析-limma----
##2.0 检测数据
par(mfrow=c(1,2)) #一式两图
boxplot(exp1,outline=F)
boxplot(log2(exp1+1),outline=F) #0-25   #log2   0-8
#过滤在所有样本都是零表达的基因~2000个
exprSet=exp1[apply(exp1,1, function(x) sum(x>1) > 1),] #17193 1001
dim(exp1);dim(exprSet)
#标准化expr
exprSet <-log2(exprSet+1);dim(exprSet) #无差异基因
#boxplot(exprSet)#0-15

##2.1 limma----
group_list <- ifelse(colnames(exp1) %in% subset(clinical,clinical$Group =="High")$sample,"High","Low") 

draw_h_v <- function(exprSet,need_DEG,n='DEseq2',group_list,logFC_cutoff){
  ## we only need two columns of DEG, which are log2FoldChange and pvalue
  ## heatmap
  need_DEG=nrDEG#降序前50  #nrDEG
  library(pheatmap) #nrDEG$log2FoldChange 降序 need_DEG[order(need_DEG$log2FoldChange),]
  choose_gene= head(rownames(need_DEG),50)#head(rownames(need_DEG[order(need_DEG$log2FoldChange,decreasing = T),]),50) ## 50 maybe better need_DEG/need_DEG[order(need_DEG$log2FoldChange),]
  choose_matrix=exprSet[choose_gene,]#提取前50的表达矩阵
  choose_matrix[1:4,1:4]
  choose_matrix=t(scale(t(log2(choose_matrix+1)))) #画热图要先归一化
  ## http://www.bio-info-trainee.com/1980.html
  annotation_col = data.frame( group_list=group_list  )#变成一个表
  rownames(annotation_col)=colnames(exprSet)#再取个名字
  pheatmap(choose_matrix,show_colnames = F,cluster_cols = T,cluster_rows = T,annotation_col = annotation_col,
           filename = paste0(1,'_need_DEG_top50_heatmap.png'))
  
  library(ggfortify)
  df=as.data.frame(t(choose_matrix))#行列转换
  choose_matrix
  df$group=group_list
  png(paste0(2,'_DEG_top50_pca.png'),res=120)#存为什么文件
  p=autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = "group")+theme_bw()#画图
  print(p)
  dev.off()
  
  #画火山图
  if(! logFC_cutoff){
    logFC_cutoff <- with(need_DEG,mean(abs( log2FoldChange)) + 2*sd(abs( log2FoldChange)) )
    
  }#如果没有给定阈值（差异倍数认为多少是上下调），则用平均值加上两倍的方差
  logFC_cutoff=0.5 #
  
  need_DEG$change = as.factor(ifelse(need_DEG$pvalue < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                     ifelse(need_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(need_DEG[need_DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(need_DEG[need_DEG$change =='DOWN',])
  )
  library(ggplot2)
  g = ggplot(data=need_DEG, 
             aes(x=log2FoldChange, y=-log10(pvalue), 
                 color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
  print(g)
  ggsave(g,filename = paste0(3,'_volcano.png'))
  dev.off()
}


library(edgeR)
if(T){ 
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)#注意limma包的表达数据是需要经过log的
  v <- voom(dge,design,plot=TRUE, normalize="quantile")#说明书需要 normalize
  fit <- lmFit(v, design)
  group_list
  cont.matrix=makeContrasts(contrasts=c('High-Low'),levels = design)#注意这里不要反了，反了就是相反的结果
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  tempOutput = topTable(fit2, coef='High-Low', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)  
  head(DEG_limma_voom)  
  nrDEG=DEG_limma_voom[,c(1,4)] 
  colnames(nrDEG)=c('log2FoldChange','pvalue') 
  draw_h_v(exprSet,nrDEG,'limma',group_list,0.5)
}
#log2+1标准化后，logFC=1,无差异；改小logFC=0.5 ，up 1 down 97
#不标准化，logFC=1,up 3, down 53 =56   #热图差异不明显？

#write.csv(DEG_limma_voom,"./output/DEG_limma_voom_fpkm.csv")#保存结果-所有基因
nrDEG1 <- subset(nrDEG,(nrDEG$log2FoldChange >0.5 | nrDEG$log2FoldChange < -0.5)& nrDEG$pvalue <0.05) #56   #0.5 98

##2.2 热图----
#：D:\研究生科研\生物信息学\STAT4\Subgroup_STAT4\1.Riskscore\Riskscore_STAT4.R
library(pheatmap) #nrDEG$log2FoldChange 降序 need_DEG[order(need_DEG$log2FoldChange),]
choose_gene= rownames(nrDEG1)#head(rownames(need_DEG[order(need_DEG$log2FoldChange,decreasing = T),]),50) ## 50 maybe better need_DEG/need_DEG[order(need_DEG$log2FoldChange),]
choose_matrix=exprSet[choose_gene,]#提取差异基因的表达矩阵
choose_matrix[1:4,1:4]
#colnames(choose_matrix) <- group_list#列名改为Low High
choose_matrix=t(scale(t(log2(choose_matrix+1)))) #画热图要先归一化
annotation_col = data.frame( Group=group_list  )#变成一个表
rownames(annotation_col)=colnames(exprSet)#再取个名字


pheatmap(choose_matrix,show_colnames = F,cluster_cols = T,cluster_rows = T,annotation_col = annotation_col,
         #cutree_col = 2,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(25))

#排列表达矩阵顺序-列完全按照annotation_col分布
choose_matrix_1 <- choose_matrix[,ID_high_low]
ID_high_low <- c(rownames(subset(annotation_col,annotation_col$Group =="Low")),#501
                 rownames(subset(annotation_col,annotation_col$Group =="High")) )

ann_colors <- list(
  STAT4 = c(High = "#FF6699", Low = "skyblue"))#分组颜色

pheat<- pheatmap(choose_matrix_1,annotation_legend=T,
             gaps_col=500,fontsize = 8,fontsize_col=1,#fontsize_row = 6,
             color<- colorRampPalette(c("navy", "white", "firebrick3"))(25),#colorRampPalette(c('#436eee','white','#EE0000'))(25),
             annotation_colors = ann_colors,
             #filename = "pheatmap_limma_197_fpkm.pdf",
             show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = T,annotation_col = annotation_col)#,
#cutree_cols = 2 )
pheat

save_pheatmap_pdf <- function(x, filename, width, height) {
  library(grid)
  x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.02)#修改聚类树线条宽度
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height,family = "serif")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}#保存热图的函数
#save_pheatmap_pdf(pheat, "./output/pheatmap_limma_98_LOGFC0.5_1.pdf",4,3) #10,7

##2.3 火山图volcano----
library("ggthemes")
library("ggplot2")

logFC_cutoff <- 0.5  #所有基因nrDEG
nrDEG$change = as.factor(
  ifelse(nrDEG$pvalue < 0.05 & abs(nrDEG$log2FoldChange) > logFC_cutoff, 
         ifelse(nrDEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')   #判断基因变化
) 
head(nrDEG) ;table(nrDEG$change) 


p <- ggplot(data = nrDEG, aes(x = -log10(as.numeric(pvalue)), y =log2FoldChange) ) +  #注意加载data
  geom_point(size = 2, aes(color= change), show.legend = T,na.rm =TRUE,alpha = 0.5) + #,alpha = 0.5
  scale_color_manual(values = c('#436eee','grey70','#EE0000')) + # 'skyblue', 'grey70', 'pink'
  labs(x = '-log10(pvalue)', y = 'Log2(fold change)')+  #设置坐标轴名称
  geom_hline(yintercept = c(-0.5, 0.5), linetype = 'dotdash', color = 'grey30') +  #设置Y轴虚线
  geom_vline(xintercept = -log10(0.05), linetype = 'dotdash', color = 'grey30')+    #设置X轴虚线
  theme_bw()+
  theme(text=element_text(size=16,family="serif"),#字体改为Times New Roman
        axis.text = element_text(size = 15, color = "black",family="serif"),#坐标轴字体
        legend.title = element_blank(),legend.background =element_rect(color='grey',linetype = 1),legend.text = element_text(size=15) ,legend.position=c(0.9,0.6))
p
#标出高表达基因
#interest_genes  <- nrDEG[which(abs(limma_voom_DEG$logFC) > 5 & limma_voom_DEG$P.Value < 0.05),]
#interest_genes_char <- rownames(nrDEG1) #c(rownames(interest_genes))
#interest_genes  <- rownames(subset(nrDEG1,nrDEG1$log2FoldChange >2 | nrDEG1$log2FoldChange < -1 ))  #FC >2 < -1 
#添加图层来标记兴趣基因-有边框版:注意边框是一个新的数据框!
nrDEG1$color <- ifelse(nrDEG1$change =="UP","pink","skyblue")  #pink skyblue
nrDEG2 <- subset(nrDEG1,nrDEG1$log2FoldChange >1.5 | nrDEG1$log2FoldChange < -1 ) #过滤
write.csv(nrDEG1,"./220719/Figure2A_DGEs.csv")

#对所有差异基因标出
nrDEG1$change = as.factor(
  ifelse(nrDEG1$pvalue < 0.05 & abs(nrDEG1$log2FoldChange) > logFC_cutoff, 
         ifelse(nrDEG1$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')   #判断基因变化
) 
nrDEG1$color <- ifelse(nrDEG1$change =="UP","pink","skyblue")

p+#geom_point(size=4,data=nrDEG,aes(color= change),alpha = 0.8)+
  ggrepel::geom_label_repel(aes(label=rownames(nrDEG1)),data = nrDEG1,color=nrDEG1$color,family="serif",max.overlaps=20)

#ggsave("./output/volcano_limma_98_LOGFC0.5.pdf.pdf",width = 8,height = 6)

#3. 富集分析----
library(clusterProfiler)
eg = bitr(rownames(nrDEG1), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(eg)#差异基因ID转化-98

##3.1 GO富集分析----
go_all <- enrichGO(eg$ENTREZID, "org.Hs.eg.db", 
                   keyType = "ENTREZID",
                   ont = 'ALL',
                   pvalueCutoff  = 0.05,
                   pAdjustMethod = "BH",  
                   qvalueCutoff  = 0.1, readable=T)  #一步到位

#结果处理及GO绘图
go_all_result <- as.data.frame(go_all@result)#所有GO结果-353
#write.csv(go_all_result,"./output/go_all_result.csv")

library(ggplot2)
library(ggrepel)
mf <-subset(go_all_result,go_all_result$ONTOLOGY=="MF")#[1:20,]  #mf 总42
mf <-mf[mf$Count > 6,] #10

p <-ggplot(mf,aes(x= -log10(p.adjust),y=reorder(Description,p.adjust) ))+#画图指定x轴与y轴
  geom_point(aes(size=Count,color=Description),alpha=0.5)+#指定气泡大小和颜色填充
  scale_size(range=c(5,15))+#图的大小设置
  theme_bw()+ #设置背景为白色
  geom_text_repel(aes(label= Description,color=Description),size= 4,segment.color= "grey", show.legend= T)+
  #geom_text( show_guide  = F )+#删除图例中的"a"
  theme(legend.position= c(0.9,0.65),# "right",#
        axis.text.y = element_blank(),
        text=element_text(size=16,family="serif"),#字体改为Times New Roman
        axis.text = element_text(size = 16, color = "black",family="serif"),
        #legend.key.size = unit(0.01, "inches"),#legend之间间距
        legend.background = element_rect(fill = NA),
        legend.box.background   = element_rect(fill = "white",size = 0.05),#填充及边框线宽
        legend.text  = element_text(color = "blue",size = 10),legend.title = element_text(color = "blue",size = 10),
  )+labs(y="Terms")+guides(color = FALSE, #删除color的legend
                           size = guide_legend(title="Counts",#修改图例标题
                                               override.aes = list(colour = "grey",alpha = 0.5,label = "")) #删除有颜色的legend ,
  )+coord_cartesian(xlim = c(2,10))#+
  #coord_cartesian(ylim = c(-0.1,11))#拓展Y轴显示范围
p  
#保存图片ggsave("./output/GO_MF.pdf",family="serif",width = 6,height = 5)
##3.2 KEGG----
#install.packages('R.utils') 
R.utils::setOption("clusterProfiler.download.method","auto")
library(DO.db)
library(clusterProfiler)
library(org.Hs.eg.db)
kegg_all <- enrichKEGG(gene = eg$ENTREZID, 
                       organism ='hsa', #"hsa" #Homo sapiens (human)
                       pAdjustMethod = "BH",keyType = "kegg",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1,
                       use_internal_data =FALSE)
#失败
# kegg_all <- enrichKEGG(gene= eg$ENTREZID,organism  = 'hsa', qvalueCutoff = 0.05)
# kegg_all <- enrichKEGG(eg$ENTREZID, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 3, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)

# write.table(eg,"eg.txt")


# fail to download KEGG data...
# Error in download.KEGG.Path(species) : 
#   'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'..

kegg_all_results <- as.data.frame(kegg_all@result)#148
#write.csv(kegg_all_results,"./output/kegg_all_results.csv")#保存未过滤的KEGG结果
##3.2.1 KEGG可视化----
#a.服务器运行KEGG,本地可视化
kegg_all_results <-read.table("./output/kegg_all_results.txt")#98
kegg <- subset(kegg_all_results,kegg_all_results$pvalue <0.05) #26   全部统一为P值
#b.数据处理：手动挑选15个 immune and cancer-related pathways
#kegg <-kegg[c("hsa04514","hsa04060","hsa04658","hsa04659","hsa04612","hsa04062","hsa04672","hsa05340","hsa04660","hsa04650","hsa0462","hsa05235","hsa04064","hsa04210","hsa04668"),]
kegg <- subset(kegg,kegg$Count > 6)#8 个
kegg <-  kegg[-c(3,4),]#6,删除Hematopoietic cell lineage Viral protein interaction with cytokine and cytokine receptor 
library(DOSE)
kegg$GeneRatio <- c(parse_ratio(kegg$GeneRatio) )
#c. 可视化
#气泡图bubble
library(ggplot2) #pvalue
ggplot(kegg,aes(GeneRatio,reorder(Description,Count)))+ #reorder(Description,Count),Description
  geom_point(aes(size=Count,color=-1*log10(pvalue)))+
  scale_color_gradient(low="green",high="red")+
  labs(color=expression(-log10(pvalue) ),size="Count",x="GeneRatio",y="",title="KEGG")+#添加图注
  theme_bw(base_family = "serif") +
  theme(text=element_text(size=15,family="serif"),axis.text=element_text(size=15))
#ggsave("./output/KEGG_bubble.pdf",width = 8,height = 4)
#条形图barplot
ggplot(kegg,aes(y=reorder(Description,Count),x=Count,fill=-1*log10(pvalue)))+
  geom_bar(stat = "identity")+
  scale_fill_gradient(low = "skyblue",high = "pink")+ #默认pvalue,勿加kegg_results_p
  labs(title = "KEGG",
       x = "Counts", 
       y = "",
       fill="-log10(pvalue)",family = "serif")+#fill修改图例名称
  theme_bw(base_family = "serif")+
  annotate("text",x=15.5,y=3.7,label="-log10(pvalue)",size=5,family = "serif")+ #手动添加注释legend
  theme(text=element_text(size=16,family="serif"),axis.text=element_text(size=16),
        panel.grid = element_blank(),legend.position = c(0.9,0.3),legend.title = element_blank(),
        legend.background = element_blank() ) #legend.text.align=1 legend对齐,   太长无效
#ggsave("./output/KEGG_barplot.pdf",width = 8,height = 4)
#write.csv(kegg,"./220719/Figure2B_kegg.csv")

##3.3 GSEA----
##a,处理数据：在差异结果中注释一行ENTREZID :https://cloud.tencent.com/developer/article/1838918
eg  #98
nrDEG1$SYMBOL <- rownames(nrDEG1)
library(dplyr)
data_all <- nrDEG1 %>% 
  inner_join(eg,by="SYMBOL")
dim(data_all)#98 6
head(data_all)
##b, 排序：将基因按照logFC进行从高到低排序，只需要基因列和logFC即可
data_all_sort <- data_all %>% 
  arrange(desc(log2FoldChange))
head(data_all_sort)

geneList = data_all_sort$log2FoldChange #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList) #"numeric"

##c,特定参考基因集GSEA分析
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
gsea_gmt <- read.gmt("./output/c7.all.v7.5.1.entrez.gmt") #读gmt文件 c7.all.v7.5.1.entrez.gmt富集失败
gsea <- GSEA(geneList,                                                    #    h.all.v7.5.1.entrez.gmt  富集2个   
             TERM2GENE = gsea_gmt) #GSEA分析
head(gsea) #no term enriched under specific pvalueCutoff...
##d,所有基因集富集分析
library(msigdbr) 
msigdbr_species()
Hs_msigdbr <- msigdbr(species="Homo sapiens")
colnames(Hs_msigdbr)
Hs_df <- as.data.frame(Hs_msigdbr[,c('gs_name','entrez_gene','gene_symbol')])
head(Hs_df)

# MSigDb数据富集分析:https://blog.csdn.net/qq_27390023/article/details/121301371
em_msig <- GSEA(geneList,pvalueCutoff = 0.1,TERM2GENE=Hs_df[,c(1,2)]) #约1分钟

head(em_msig,20)  # sorted by pvalue   #no term enriched under specific pvalueCutoff.
em_msig@result$Description #32
gsea_results <- em_msig@result  #未富集到GSEA

##3.4 GSVA----
##3.4.1.准备及运行
library(GSVA)
library(GSEABase)
library(clusterProfiler)
#加载注释gmt
gsea_gmt <- read.gmt("D:/研究生科研/生物信息学/STAT4/Subgroup_STAT4/1.Riskscore/fig_tab/enrichment/h.all.v7.5.1.symbols.gmt") #读gmt文件 

gsea_gmt_list = split(gsea_gmt$gene, gsea_gmt$term)   ##分隔成list  
gsva_gsea_gmt<- gsva(as.matrix (log2(exp1+1) ), #标准化矩阵log2(exp1+1)
                     gset.idx.list=gsea_gmt_list, #基因集
                     kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                     method = "gsva",
                     parallel.sz=12) # 并行线程数目  5000+ 30min， 50 3min
dim(gsva_gsea_gmt)  #50!
gsva_gsea_gmt <- as.data.frame(gsva_gsea_gmt)
#注意保存Rdata

#3.4.2绘制热图
library(pheatmap)
pheatmap(gsva_gsea_gmt,cluster_rows = T,scale = "row",
         show_colnames=F,
         cutree_cols=2) #不要运行太大的热图数据！！！

##3.4.3.差异的GSVA基因集-heatmap-volcano
library(limma)
head(design) #分组
compare <- makeContrasts(High - Low, levels = design)
fit <- lmFit(gsva_gsea_gmt, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff_gsva <- topTable(fit3, coef = 1, number = 10000)
head(Diff_gsva )  #https://blog.csdn.net/weixin_41368414/article/details/123967524
#hallmark基因集，差异分析-|FC|小于0.2,自定义cutoff
#write.csv(Diff_gsva,"./output/Diff_gsva_hallmark_limma.csv")

##3.4.4 可视化差异GSVA
Diff_gsva_padjust <- subset(Diff_gsva,Diff_gsva$adj.P.Val < 0.05) #43
gsva_gsea_gmt_padjust <- subset(gsva_gsea_gmt,rownames(gsva_gsea_gmt) %in% rownames(Diff_gsva_padjust) )#43
#删除去掉"HALLMARK_"
library(stringr)
rownames(gsva_gsea_gmt_padjust) <- str_replace(rownames(gsva_gsea_gmt_padjust), "HALLMARK_","")
rownames(Diff_gsva_padjust) <- str_replace(rownames(Diff_gsva_padjust), "HALLMARK_","")

#c.2.1 差异热图
pheatmap(gsva_gsea_gmt_padjust,cluster_rows = T,scale = "row",
         show_colnames=F,cutree_cols=2,
         annotation_col = annotation_col)
#排列表达矩阵顺序-列完全按照annotation_col分布
gsva_gsea_gmt_padjust_1 <- gsva_gsea_gmt_padjust[,ID_high_low]
# ID_high_low <- c(rownames(subset(annotation_col,annotation_col$STAT4 =="Low")),#534
#                  rownames(subset(annotation_col,annotation_col$STAT4 =="High")) )
#ann_colors <- list(STAT4 = c(High = "#FF6699", Low = "skyblue"))#分组颜色

p0 <- pheatmap(gsva_gsea_gmt_padjust_1 ,annotation_legend=T,gaps_col=500,#cutree_cols=2,#gaps_col=533,
               fontsize = 7,fontsize_col=1,
               color<- colorRampPalette(c('skyblue','white','firebrick'))(25),#colorRampPalette(c('#436eee','white','#EE0000'))(25),
               annotation_colors = ann_colors,
               #filename = "pheatmap_limma_197_fpkm.pdf",
               show_rownames = T,show_colnames = F,cluster_cols = F,cluster_rows = T,annotation_col = annotation_col)#,
p0 #gsva 分组热图

save_pheatmap_pdf_1 <- function(x, filename, width=6, height=4.5) { 
  stopifnot(!missing(x)) 
  stopifnot(!missing(filename)) 
  pdf(filename, width=width, height=height,family = "serif") 
  grid::grid.newpage()
  grid::grid.draw(x$gtable) 
  dev.off() 
}
save_pheatmap_pdf_1(p0, "./output/gsva_diff_heatmap.pdf")
write.csv(gsva_gsea_gmt_padjust_1,"./220719/Figure2C_gsva.csv")


##GSVA结果解度与展示？score  gsva_gsea_gmt_padjust  cell cycle  EMT and metabtasis
#合并43个显著差异到临床数据中
gsva_gsea_gmt_padjust_1 <- as.data.frame(t(gsva_gsea_gmt_padjust))
gsva_gsea_gmt_padjust_1$sample <- rownames(gsva_gsea_gmt_padjust_1)
clinical_1 <- merge(clinical,gsva_gsea_gmt_padjust_1,by="sample")#1001-174

#EMT通路与riskscore
cor(clinical_1$lasso.risk.score,clinical_1$EPITHELIAL_MESENCHYMAL_TRANSITION) #-0.19
plot(clinical_1$lasso.risk.score,clinical_1$EPITHELIAL_MESENCHYMAL_TRANSITION)
boxplot(factor(clinical_1$Group),clinical_1$EPITHELIAL_MESENCHYMAL_TRANSITION)

library(ggplot2);library(ggpubr)
clinical_1$STAT4_group
ggplot(data=clinical_1,aes(x=STAT4_group,y=EPITHELIAL_MESENCHYMAL_TRANSITION,color=STAT4_group))+  #fill
  # geom_point(alpha=0.5,size=1.5,
  #            position=position_jitterdodge(jitter.width = 0.35,
  #                                          jitter.height = 0,
  #                                          dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
             position=position_dodge(width=0.8),
             size=0.05,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.25)+
  scale_color_manual(values = c("#CC79A7", "#56B4E9") )+ #color
  
  theme_bw() +ylab("EMT gasv score")+xlab("")+
  #theme(legend.position="none") 
  theme(axis.title  = element_text(size=16), axis.text = element_text(size=12),legend.position ="top",  
        panel.grid = element_blank(),legend.key = element_blank() ) +
  labs(color="STAT4")+
  stat_compare_means(#aes(group = Group) ,
    comparisons=list(c("High","Low")),#my_comparisons,
    
    label = "p.signif",#"p.format p.signif 9.4e-12
    method = "wilcox",
    show.legend= F,#删除图例中的"a"
    label.x=1.5,bracket.size=0.1,vjust=0.2,
    #label.y = max(log2(GDSC2$Camptothecin_1003)),
    hide.ns = T,size=5)

##全部GSVA  BOXPLOT
##cox gsva


##COX GSVA boxplot

#A.循环分组！！
group_df<-list()

for (i in  133:175){
  gene <- colnames(clinical_1)[i]  #
  group_name <- paste0(gene,"_gsva_group")
  group=ifelse(clinical_1[,gene] >= median(clinical_1[,gene] ),"High","Low") 
  #print(group_name)
  #print(group)
  group_df[[group_name]] <- group #list()
}
View(group_df)
group_df <- as.data.frame(group_df)
rownames(group_df) <-clinical_1$sample


#添加临床信息
group_df$sample <- rownames(group_df)
group_df<- merge(group_df,clinical_1[,c("OS.time.y","OS.x","sample")],by="sample")
rownames(group_df) <-group_df$sample

#B ##批量OS 
pFilter=0.05                        #定义单因素显著性
library(survival)

sigGenes=c("OS.time.y","OS.x","sample")#注意修改此处及和下面,需要取出的列
outTab_TF=data.frame()
for(i in colnames(group_df[,2:44]) ){   #注意此处为参与COX的名称
  cox <- coxph(Surv(OS.time.y, OS.x) ~ group_df[,i], data = group_df) #as.numeric
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab_TF=rbind(outTab_TF,
                  cbind(id=i,
                        HR=coxSummary$conf.int[,"exp(coef)"],
                        HR.95L=coxSummary$conf.int[,"lower .95"],
                        HR.95H=coxSummary$conf.int[,"upper .95"],
                        pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
  }
}

View(outTab_TF)

#write.table(outTab_TF,file="./output/uniCox_outTab_GSVA.xls",sep="\t",row.names=F,quote=F)
outTab_TF_sig <- subset(outTab_TF,outTab_TF$pvalue <0.05)#6 gsva
outTab_TF_sig$id
gsub("_gsva_group","",outTab_TF_sig$id)

##C. OS显著KM绘图  :tnfa_nfkb[1] apoptosis[2]  勉强mycv1 [5]  angiogenesis[6]
library(survminer)
for (i in outTab_TF_sig$id[6] ){
  group_name <- paste0(i,"_TF_group")
  fit<-survfit(Surv(OS.time.y/365, OS.x) ~ group_df[[i]], data =group_df)#group_df[[group_name]]
  #print(fit)#查看median风险50%生产率对应的时间！
  p<-ggsurvplot(fit,#palette=c("#FC8D62","#66C2A5"),
                data=group_df,
                risk.table = F,
                risk.table.col = 'strata',
                risk.table.y.text = T,
                risk.table.title = "",#group_name,#paste0('Number at risk',i),
                pval = T,pval.size=6,pval.method = T,pval.method.size=6,
                conf.int = T, #置信区间
                xlab = 'Time (years)',
                ggtheme = theme_light(),
                font.x=c(16),font.y=c(16),font.lenged=c(16),font.tickslab = c(16), #Y轴数字大小
                size=0.5,censor.size=3,censor.shape=NA,#NA不要+号
                surv.median.line = 'none',#不要中位值竖线
                title=gsub("_gsva_group","",i),#paste0("Riskscore",' Survival Curve'),
                legend.title = "Group",
                legend.labs = c("High", "Low")#, fontsize=20,font.family="serif"
                )
  
  
}
p #theme(text = element_text(size=16)) #+guides(color = FALSE)
#ggsave("./output/gsva_OS_angiogenesis.pdf",family="serif")

##D. os cox forest 
#准备森林图数据
forest_data <- outTab_TF_sig
colnames(forest_data)<- c("id","HR","Low_95CI","High_95CI","Pvalue")#,"Coeffcient")
forest_data$Pvalue <- round(as.numeric(forest_data$Pvalue),4)
#forest_data$Coeffcient <- round(as.numeric(forest_data$Coeffcient),4)
#forest_data$Pvalue <- ifelse(as.numeric(forest_data$Pvalue)>0.05,"ns",ifelse(as.numeric(forest_data$Pvalue)>0.01,"*",ifelse(as.numeric(forest_data$Pvalue)>0.001,"**","***")))
forest_data$HR <-round(as.numeric(forest_data$HR),3)
forest_data$id <- gsub("_gsva_group","",forest_data$id)

forest_data$yend1 <-nrow(forest_data):1
forest_data[,3:7] <- as.data.frame(sapply(forest_data[,3:7] ,as.numeric))
forest_data <- forest_data[order(forest_data$HR),] #HR排序
rownames(forest_data) <- 1:6 #排序有效
forest_data$yend <-1:nrow(forest_data)

library(ggplot2)
ggplot(data = forest_data,aes(x = HR, y = reorder(id,HR),group = id)
) +
  geom_segment( aes(x =Low_95CI, xend = High_95CI, y = id, yend = id),color = "black",size=0.1 )+ # 主要横线
  #geom_segment( aes(x =as.numeric(Low_95CI), xend = as.numeric(Low_95CI), y = yend, yend = yend+0.2),color = "black" ,size=0.1)+
  #geom_segment( aes(x =as.numeric(Low_95CI), xend = as.numeric(Low_95CI), y = yend, yend = yend-0.2),color = "black" ,size=0.1)+
  #geom_segment( aes(x =as.numeric(High_95CI), xend = as.numeric(High_95CI), y = yend, yend = yend+0.2),color = "black",size=0.1 )+
  #geom_segment( aes(x =as.numeric(High_95CI), xend = as.numeric(High_95CI), y = id, yend =yend-0.2 ),color = "black" ,size=0.1)+
  geom_segment( aes(x=1,xend=1,y=0.5,yend=6.5),color="grey",linetype='dashed',size=0.1)+ #HR=1
  
  geom_point(size=4,aes(color = 100*Pvalue),pch=15)+
  #scale_fill_manual(values = c("#93cc82","#88c4e8"))+ # c("#93cc82","#88c4e8") #"#CC79A7", "#56B4E9" 天蓝，浅紫   #"#238b45","#2171b5"深蓝，绿
  labs(x = 'Hazard Ratio', y = ''#,color="Sign."
       )+theme_bw()+
  theme(text = element_text(size=16,family = "serif"),#legend.position = c(0.8,0.8),
        panel.grid = element_blank() )

ggsave("./output/Forestplot_GSVA_2.pdf",width=8, height=3.5) 

##类似SSGSVA画法--挑OS通路
gsva_boxplot <- merge(gsva_gsea_gmt_padjust_1,clinical_1[,c("sample","STAT4_group")],by="sample")
rownames(gsva_boxplot) <- gsva_boxplot[,1]
gsva_boxplot <- gsva_boxplot[,-1] #1-43
library(reshape2)
gsva_New = melt(gsva_boxplot[,c(forest_data$id,"STAT4_group")])# 融合数据:宽数据变成长数据43043-3   取OS GSVA!!!!!   6006-3
colnames(gsva_New )=c("group_stat4","Celltype","Score")  #设置行名 #无sample?
head(gsva_New )
# 按免疫细胞占比中位数排序绘图（可选）
library(dplyr)
plot_order_1 = gsva_New [gsva_New $group_stat4=="High",] %>% 
  group_by(Celltype) %>% 
  summarise(m = median(Score)) %>% #
  arrange(desc(m)) %>% 
  pull(Celltype)

gsva_New$Celltype = factor(gsva_New $Celltype,levels = plot_order_1)

if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 0, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
library(ggplot2);library(ggpubr)
box_TME <- ggplot(gsva_New , aes(x = Celltype, y = Score))+ ##Celltype/T cells follicular /helper Neutrophils
  labs(y="ssGSEA score",x= NULL)+  
  geom_boxplot(aes(fill = group_stat4),outlier.size = 0.02,size=0.02,outlier.alpha = 0)+ 
  #geom_bar(stat = "identity")+
  scale_fill_manual(values = c("skyblue","pink") )+ #"#1CB4B8", "#EB7369" "skyblue","pink"
  guides(fill = guide_legend(title = 'STAT4'))+
  theme_bw()+ mytheme + 
  #ylim(-0.4,0.65)+
  stat_compare_means(aes(group =  group_stat4),
                     label = "p.signif",
                     #label.y = 0.6,
                     size=3,#0.6
                     method = "wilcox.test",
                     hide.ns = T)+
  coord_flip()
# geom_signif(annotations = c(""),color="black",size = 0.1, #手动加括号
#             y_position = c(0.12,0.12),
#             xmin = c(8.81,21.81),
#             xmax = c(9.188,22.188),
#             tip_length = c(0.01,0.02,0.1,0.1))


box_TME
ggsave("./immune/ssGSEA_box_OS6_2.pdf",box_TME,height=5,width=6,family="serif") 
write.csv(gsva_New,"./220719/Figure2E_gsva_boxplot.csv")

###
plotList_3 <- list()
corPlotNum <- 50
if(nrow(outTab)<corPlotNum){
  corPlotNum=nrow(gsva_boxplot[,1:43])  
}  #同样，定义一个空的列表plotList_2用于保存输出结果
for(i in 1:corPlotNum){ #1:corPlotNum
  Gene <- Gene
  Drug <- gsva_boxplot[,i]#outTab_GDSC2[abs(as.numeric(outTab_GDSC2$cor)) > 0.4 ,][i,2] #3 个#outTab_GDSC2[i,2]
  
  #x <- as.numeric(exp[Gene,])
  #y <- as.numeric(drug[Drug,])
  df1 <- as.data.frame(cbind(Gene,gsva_boxplot[,Drug]))
  colnames(df1)[2] <- "IC50" #GSVA
  df1$IC50 <- as.numeric(df1$IC50)
  #df1$group <- ifelse(df1$x > median(df1$x), "high", "low")
  df1$group <- gsva_boxplot$STAT4_group #subset(GDSC2,colnames(GDSC2)== Drug)$STAT4_group
  compaired <- list(c("Low", "High"))
  p1 <-ggplot(df1,aes(x =group, y =IC50,color =group))+ #,palette = c("firebrick3","skyblue")
    labs(title= gsub("\\_.*","",Drug),ylab="Log2(IC50)")+
    geom_point(alpha=0.2,size=1.5,
               position=position_jitterdodge(jitter.width = 0.45,
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
    geom_boxplot(alpha=0.5,width=0.55,
                 position=position_dodge(width=0.8),
                 size=0.05,outlier.colour = NA)+
    # geom_violin(alpha=0.2,width=0.9,
    #             position=position_dodge(width=0.8),
    #             size=0.25)+
    scale_color_manual(values = c("#CC79A7", "#56B4E9") )+ #color
    
    # stat_compare_means(comparisons = compaired,vjust=0.4,
    #                    method = "wilcox.test",size=4,   #设置统计方法
    #                    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
    #                                     symbols = c("***", "**", "*", "ns")))+
    #coord_cartesian(ylim = c( NA,max(log2(df1$IC50))+1.5 )  )+   #否则*显示不全,log2(0)=inf
    ylim(c( NA,max(log2(df1$IC50))+1.5 ))+
    stat_compare_means(#aes(group = Group) ,
      comparisons=list(c("High","Low")),#my_comparisons,
      label = "p.signif",#"p.format p.signif 9.4e-12
      method = "wilcox",
      show.legend= F,#删除图例中的"a"
      label.x=1.5,bracket.size=0.1,vjust=0.1,
      #label.y = max(log2(GDSC2$Camptothecin_1003)),
      hide.ns = T,size=4)+
    theme_bw()+
    theme(axis.title  = element_text(size=12), axis.text = element_text(size=12),legend.position ="none",  
          panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() )
  
  plotList_3[[i]]=p1
}
ggarrange(plotlist=plotList_3,nrow=5,ncol=5) #STAT4_group



#ggsave("./output/GSVA_boxplot_1.pdf",family="serif",width = 6,height = 2)



#放弃：high riskscore should be more higher EMT gsva score, but the fact is opposite.
#3.5 EMT 评分----
#3.5.1 log2 Z scores https://zhuanlan.zhihu.com/p/366181769  ；#https://www.sciencedirect.com/science/article/pii/S0169500219306932#sec0095
#EMT genes:
Epithelial_genes <- c("CDH1","CDH3","CLDN4","EPCAM","ST14","MAL2") #6
Mesencyhmal_genes<-c("VIM","SNAI2","ZEB2","FN1","MMP2","AGER","SNAI1","TWIST1","ZEB1","FOXC1","TWIST2","CDH2")#12
EMT_genes_logZ <- data.frame(genes=Epithelial_genes,
                             group="Epithelial")
EMT_genes_logZ <- rbind(EMT_genes_logZ, data.frame(genes=Mesencyhmal_genes,
                                                   group="Mesencyhmal")  )
#a. 选取高表达六个基因计算
EMT_exp_logZ <- subset(exp1,rownames(exp1) %in% EMT_genes_logZ$genes) #18 1001
boxplot(log2(EMT_exp_logZ[1,] ))#查看基因表达
for (i in 1:18){
  print(rownames(EMT_exp_logZ[i,]))
  print(mean(as.numeric(EMT_exp_logZ[i,]) ))
} #查看均值
 

#计算EMT Zscore
zscore <- matrix(nrow = 18,ncol = 1001)
zscore <- as.data.frame(zscore)
for (i in 1:1001){  #1001
  for (j in 1:18){
    varname = rownames(EMT_exp_logZ)[j]
    zscore0 = (EMT_exp_logZ[varname,i] - mean(as.numeric(EMT_exp_logZ[varname,])) )/ sd(as.numeric(EMT_exp_logZ[varname,]) ) #！注意公式！
    zscore[j,i] <- zscore0
    
    rownames(zscore)[j]<-varname
    colnames(zscore)[i]<-colnames(EMT_exp_logZ)[i]
  }
  
} #时间较久~10min

boxplot(zscore[1,] )  #查看zscore与log2 zscore
boxplot(log2(zscore[1,]) ) 

zscore_mean_1 <-data.frame()
for (i in 1:18){
  
  # ifelse(rownames(zscore[i,]) %in% Epithelial_genes,
  # print(rownames(zscore[i,])+ mean(as.numeric(zscore[i,]) )),
  # "")
  # #print(mean(as.numeric(zscore[i,]) ))
  
  zscore_mean <- data.frame(mean=mean(as.numeric(zscore[i,])),gene= rownames(zscore[i,]) )
  zscore_mean_1 <- rbind(zscore_mean_1,zscore_mean)
  } #查看均值
zscore_mean_1 <- subset(zscore_mean_1,zscore_mean_1$gene %in%  Mesencyhmal_genes) #12 #选取高表达六个Mesencyhmal_genes基因
 #score EMT

zscore1<- subset(zscore,rownames(zscore)  %in%c("VIM","FN1","MMP2","SNAI2","FOXC1","TWIST1",
                                                "CDH1","CDH3","CLDN4","EPCAM","ST14","MAL2")  )  #最终EMT 基因 12-1001
#计算每个样本的EMT score； log2 Z
#zscore1<-log2(zscore1) #负数取log2 会有NA？？
zscore1<- as.data.frame(t(zscore1))
zscore1$EMT_Zscore <- rowSums(zscore1[,c("VIM","FN1","MMP2","SNAI2","FOXC1","TWIST1")])-rowSums(zscore1[,c("CDH1","CDH3","CLDN4","EPCAM","ST14","MAL2")])
  
#合并结果到clinical


#3.5.2 https://zhuanlan.zhihu.com/p/442875983  
#a.读取基因
emt315 <- readxl::read_xlsx("./emt/EMT315.xlsx",col_names = NA)
colnames(emt315) <- c("gene","group")
Epithelial_genes145 <- subset(emt315,emt315$group=="Epi")$gene #145
Mesencyhmal_genes170 <- subset(emt315,emt315$group=="Mes")$gene #170

#b. 计算得分
emt315_exp <- subset(exp1,rownames(exp1) %in% c(emt315$gene,"GZMA","PRF1") ) #取EMT+"GZMA","PRF1"基因矩阵:305-1001
emt315_exp <- as.data.frame(t(emt315_exp))#转置矩阵后计算得分-1001-305

Mesencyhmal_genes163<- colnames(emt315_exp)[colnames(emt315_exp) %in% Mesencyhmal_genes170] #163   -7
Epithelial_genes140<- colnames(emt315_exp)[colnames(emt315_exp) %in% Epithelial_genes145]  #140   -5

emt315_exp$EMT_score <-(rowSums(emt315_exp[,Mesencyhmal_genes163])/1001-rowSums(emt315_exp[,Epithelial_genes140])/1001 ) #/ (sd(emt315_exp[,Mesencyhmal_genes163])-sd(emt315_exp[,Epithelial_genes140]) )#计算EMT315
sd(emt315_exp[2,Mesencyhmal_genes163])

EMT_score_df <- data.frame()
for (i in 1:1001){ #i 为样本数
  mes_exp <- emt315_exp[,colnames(emt315_exp) %in%  Mesencyhmal_genes163][i,]#1001-163
  epi_exp <- emt315_exp[,colnames(emt315_exp) %in%  Epithelial_genes140][i,]
  EMT_score <- (mean(as.numeric(mes_exp)) -mean(as.numeric(epi_exp)) )/(sd(mes_exp)-sd(epi_exp) )##-mean(emt315_exp[i,Epithelial_genes140]) )#/()
  CYT=sqrt(emt315_exp[i,]$GZMA * emt315_exp[i,]$PRF1)
  EMT_score_df <- rbind(data.frame(sample=rownames(emt315_exp)[i],EMT_score=EMT_score,CYT=CYT),EMT_score_df) 
}
emt <- (EMT_score_df$EMT_score - mean(EMT_score_df$EMT_score) )/sd(EMT_score_df$EMT_score)
cyt <- (EMT_score_df$CYT - mean(EMT_score_df$CYT) )/sd(EMT_score_df$CYT)

EMT_score_df$ECI <- log( exp(emt) / exp(cyt)) #CALCULATE ECI
##合并2个EMT到临床数据中：EMT_score_df，zscore1，clinical_1
#clinical_1 <- merge(clinical,gsva_gsea_gmt_padjust_1,by="sample")#1001-174

clinical_1 <-merge(clinical_1,EMT_score_df,by="sample") #1001-177
clinical_1 <- cbind(clinical_1,as.data.frame(zscore1$EMT_Zscore)) #1001-178
colnames(clinical_1)[178] <-  "EMT_Zscore"
#write.csv(clinical_1,"clinical_1.csv") #保存clinical_1 ，以后直接读取！！！
clinical_1 <-read.csv("clinical_1.csv")

##3.5.3 EMT score 与临床

library(ggplot2);library(ggpubr)                 #ECI  CYT   Group
ggplot(data=clinical_1,aes(x=STAT4_group,y=log2(EMT_Zscore),color=STAT4_group))+  #fill
  geom_point(alpha=0.5,size=1.5,
             position=position_jitterdodge(jitter.width = 0.35,
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.5,
               position=position_dodge(width=0.8),
               size=0.05,outlier.colour = NA)   +
  # geom_violin(alpha=0.2,width=0.9,
  #             position=position_dodge(width=0.8),
  #             size=0.25)+
  #scale_color_manual(values = c("#CC79A7", "#56B4E9") )+ #color
  
  theme_bw() +#ylab("EMT gasv score")+xlab("")+
  #theme(legend.position="none") 
  theme(axis.title  = element_text(size=16), axis.text = element_text(size=12),legend.position ="top",  
        panel.grid = element_blank(),legend.key = element_blank() ) +
  stat_compare_means(#aes(group = Group) ,
    comparisons=list(c("High","Low")),#my_comparisons,
    
    label = "p.signif",#"p.format p.signif 9.4e-12
    method = "wilcox",
    show.legend= F,#删除图例中的"a"
    label.x=1.5,bracket.size=0.1,vjust=0.2,
    #label.y = max(log2(GDSC2$Camptothecin_1003)),
    hide.ns = T,size=5)
#RISKSCORE GROUP ECI 合理；EMT RIKSCORE/ EMT zscore 无显著性;STAT4 group趋势相反
# ECI,CYT 结果不受log2影响；

#CYT 预后：https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8332854/#!po=20.4545
#解释CYT,EMT,ECI与riskscore  ?或者放弃

#4. 免疫评分----
#4.1 estimate ----
##安装包：下次使用，请设置工作目录
# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)#引用包
##输入表达矩阵，基因在行，sample在列。检查文件格式

#write.table(exp1,"exp1.txt",sep = "\t",quote=F)#输入矩阵：https://www.jianshu.com/p/e72914fbcfaf
in.file2 <- "D:/\u7814\u7a76\u751f\u79d1\u7814/\u751f\u7269\u4fe1\u606f\u5b66/STAT4/Subgroup_STAT4/3.Enrichment/exp1.txt"
filterCommonGenes(input.f=in.file2,output.f="./immune/breast_19551genes.gct")#9876 genes (536 mismatched)
##计算评分
estimateScore("./immune/breast_19551genes.gct", "./immune/brca_estimate_score.gct", platform="illumina") #illumina
##读取结果及可视化
estimatescores=read.table("./immune/brca_estimate_score.gct",skip = 2,header = T)
rownames(estimatescores)=estimatescores[,1]
estimatescores.final=t(estimatescores[,3:ncol(estimatescores)])
#save(estimatescores.final, file="estimatescores.final.Rda")
#write.table(estimatescores.final, file = "./immune/estimatescores.final.txt",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖score
#D:\研究生科研\生物信息学\STAT4\Subgroup_STAT4\3.Enrichment\immune
##合并结果：https://www.jianshu.com/p/800ab6347d7e
rownames(estimatescores.final)<- gsub("\\.","-",rownames(estimatescores.final) )#处理结果
estimatescores.final <- as.data.frame(estimatescores.final)
estimatescores.final$sample <- rownames(estimatescores.final) #1003-4
clinical_1 <- merge(clinical_1,estimatescores.final,bvy="sample")#合并到临床数据中-178->181
##简要可视化
#画图
library(ggpubr)
ggboxplot(clinical_1, x = "STAT4_group", y = "StromalScore",#"ESTIMATEScore",# ImmuneScore   StromalScore
          fill = "STAT4_group", palette = "lancet")+
  stat_compare_means(aes(group = STAT4_group),
                     method = "wilcox.test",
                     label = "p.signif",label.x = 1.5,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=0, hjust=1)) #重新绘制
#ESTIMATEScore= ImmuneScore +  StromalScore   #合理:STAT4 LOW均低于STAT HIGH ！！！
##可视化二：宽数据转为长数据，绘图！https://www.jianshu.com/p/800ab6347d7e
estimate_test <-clinical_1[,c("sample","Group","STAT4_group","StromalScore","ESTIMATEScore","ImmuneScore")]  #estimatescores.final
estimate_test <-tidyr::gather(estimate_test,key=category,value = score,ESTIMATEScore,ImmuneScore,StromalScore,-STAT4_group)
#画图
ggboxplot(estimate_test, x = "category", y = "score",
          fill = "STAT4_group", palette = "lancet",outlier.shape=NA,size=0.1,alpha=0.5 )+
  #labs(ylab="Estimate Score",xlab="")+
  ylab("Scores")+xlab("")+
  scale_fill_manual(values = c("#CC79A7", "#56B4E9") )+
  guides(fill = guide_legend(title = 'STAT4'))+
  stat_compare_means(aes(group = STAT4_group),#
                     method = "wilcox.test",
                     label = "p.signif",size=5,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme_bw()+
  theme(text = element_text(size=14)#,axis.text.x = element_text(angle=30, hjust=1)
        ) 
#ggsave("./immune/estimate_box.pdf",family="serif",width = 6,height = 4)
#write.csv(estimate_test,"./220719/Figure3B_estimate_box.csv")

##图丑，重新画
library(ggplot2);library(ggpubr)
ggplot(estimate_test,aes(x =category, y =score,color =STAT4_group))+ #,palette = c("firebrick3","skyblue")
  ylab("Score")+
  geom_point(alpha=0.1,size=1.0,
             position=position_jitterdodge(jitter.width = 0.2,
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  geom_boxplot(alpha=0.5,width=0.45,
               position=position_dodge(width=0.8),
               size=0.05,outlier.colour = NA)+
  # geom_violin(alpha=0.2,width=0.9,
  #             position=position_dodge(width=0.8),
  #             size=0.25)+
  scale_color_manual(values = c("#56B4E9","#CC79A7") )+ #color "#56B4E9","#CC79A7" 天蓝 朱红
  #coord_cartesian(ylim = c(NA,max(estimate_test$score)+100 ))+
  stat_compare_means(aes(group = STAT4_group),label = "p.signif",size=3,show.legend = F,label.y = c(5200,3500,3000))+
  geom_signif(annotations = c("","",""),color="black",size = 0.1, #手动加括号
              y_position = c(5200,3500,3000),
              xmin = c(0.8,1.8,2.8),
              xmax = c(1.2,2.2,3.2),
              tip_length = c(0.02,0.2,0.02,0.18,0.04,0.1))+
  guides(color = guide_legend(title = 'STAT4'))+ ##/ labs(fill/ color="")+ # legend https://www.jianshu.com/p/b50c8393e3ed
  theme_bw()+
  theme(axis.title  = element_text(size=12), #axis.text = element_text(size=12),#legend.position ="none",  
        axis.text.x = element_text(size=10),#angle=0, hjust=0,
        legend.position = "top",
        panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() )
#list(c("High","Low")) 需要TRUE/FALSE值的地方不可以用缺少值 可能与长数据有关
ggsave("./immune/estimate_boxplot_1.pdf",family="serif",width = 5,height = 3)


##4.2 cibersort----
#a.安装包：https://blog.csdn.net/m0_58549466/article/details/124255582
devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)
if(!require(CIBERSORT))devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT) #version 0.1.0
#b. 加载数据
data(LM22) ;LM22
data(mixed_expr);mixed_expr#TCGA的演示数据，正式情况下就用自己的数据:matrix


# 分别定义signature矩阵LM22和我的数据（演示）矩阵mixed_expr
results <- cibersort(sig_matrix = LM22, mixture_file = exp1_matrix)#表达矩阵-1001时间较久~10min

#按行（样本内部）标准化可以看出在各类样本内部，M2浸润程度（占比）最高
rowscale <- results[,1:ncol(LM22)]#只是相当于备份了一下results
rowscale <- rowscale[,apply(rowscale, 2, function(x){sum(x)>0})]#删除全是0的列
rowscale <- as.data.frame(rowscale)  #cibersort结果
rowscale1 <- cbind(estimatescores.final[1:1001,1:3],rowscale) #cibersort+合并estimate结果[,1:3]  [,3:25]
rowscale1[,1:3] <- as.data.frame(sapply(rowscale1[,1:3],as.numeric))

library(pheatmap)

annotation_col1 <- data.frame( Group= clinical_1$Group,STAT4=clinical_1$STAT4_group,PAM50=clinical_1$PAM50)#1001
row.names(annotation_col1) <- rownames(rowscale1)

pheatmap::pheatmap(t(rowscale1[,1:3]),
         #scale = 'row',#按行标准化，不标准化就会按绝对值显示，很诡异
         #cluster_col=T,#是否对列聚类，不聚类，坐标轴就按照原来的顺序显示
         #cluster_row=F,#是否对行聚类
         #angle_col = "315",#调整X轴坐标的倾斜角度
         show_colnames = F,
         annotation_col = annotation_col1) #注意方向！Error in check.length("fill") : 'gpar' element 'fill' must not be length 0

# 各类样本之间也具有自己占比高的特异性免疫细胞
pheatmap(rowscale,
         scale = 'column',
         cluster_col=F,
         cluster_row=T,
         angle_col = "315")
# 堆积比例图
# results <- as.data.frame(results)
# results$sample <- rownames(results)
# results<- merge(results,clinical_1[,c("sample","Group","STAT4_group")],by="sample")

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175'
)
cellnum <- results[,1:22] #ncol(LM22)
cell.prop<- apply(cellnum, 1, function(x){x/sum(x)})
data4plot <- data.frame()
for (i in 1:ncol(cell.prop)) {
  data4plot <- rbind(
    data4plot,
    cbind(cell.prop[,i],rownames(cell.prop),
          rep(colnames(cell.prop)[i],nrow(cell.prop)
          )
    )
  )
}
colnames(data4plot)<-c('proportion','celltype','sample')
data4plot$proportion <- as.numeric(data4plot$proportion)

library(ggplot2)
ggplot(data4plot,aes(sample,proportion,fill=celltype) )+
  geom_bar(stat="identity",position="fill",alpha=0.8)+
  scale_fill_manual(values=my36colors)+#自定义fill的颜色
  ggtitle("cell portation")+ylab("Cell propotion")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.title.x=element_text(size=1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+#把x坐标轴横过来
  guides(fill=guide_legend(title=NULL)) #样本太大，分组？堆叠图
#绘制热图？

##可视化版本二：添加分组信息：https://blog.csdn.net/qazplm12_3/article/details/116773501
#1.分组信息准备
group_riskscore <-clinical_1$Group %>% 
  factor(.,levels = c("High","Low"))
group_stat4<-clinical_1$STAT4_group %>% 
  factor(.,levels = c("High","Low"))
#2.合并到结果中
TME_data <- as.data.frame(results[,1:22])

TME_data$Group <- group_riskscore
TME_data$group_stat4 <- group_stat4
TME_data$sample <- row.names(TME_data)
library(reshape2)
TME_New = melt(TME_data)# 融合数据:宽数据变成长数据
colnames(TME_New)=c("Group","group_stat4","Sample","Celltype","Composition")  #设置行名
head(TME_New)
#3. 按免疫细胞占比中位数排序绘图（可选）
plot_order = TME_New[TME_New$group_stat4=="High",] %>% 
  group_by(Celltype) %>% 
  summarise(m = median(Composition)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltype)

TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)
#4. 绘制箱线图

if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 90, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
library(ggplot2);library(ggpubr)
box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition))+ ##Celltype/T cells follicular /helper Neutrophils
  labs(y="Cell composition",x= NULL)+  
  geom_boxplot(aes(fill = group_stat4),outlier.size = 0.05,size=0.02,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("skyblue","pink") )+ #"#1CB4B8", "#EB7369"
  guides(fill = guide_legend(title = 'STAT4'))+
  theme_bw() + mytheme + 
  #ylim(0,0.009)+
  stat_compare_means(aes(group =  group_stat4),
                     label = "p.signif",
                     label.y = 0.12,size=5,#0.6
                     method = "wilcox.test",
                     hide.ns = T)+
  geom_signif(annotations = c(""),color="black",size = 0.1, #手动加括号
              y_position = c(0.12,0.12),
              xmin = c(8.81,21.81),
              xmax = c(9.188,22.188),
              tip_length = c(0.01,0.02,0.1,0.1))


box_TME#cibersort
#ggsave("./immune/cibersort_box_1.pdf",box_TME,height=5,width=10,family="serif") #Neutrophils <0.5%不讨论
#write.csv(TME_New,"./220719/Figure3D_cibersort.csv")
#write.csv(t(results[,1:22]),"./220719/Figure4A3_cibersort.csv")

#BOXPLOT 
colnames(TME_data)[22] #Neutrophils
mean(TME_data$Neutrophils)
ggplot(subset(TME_data,TME_data$Neutrophils <0.1), aes(x = group_stat4, y = Neutrophils))+ ##Celltype/T cells follicular /helper Neutrophils
  labs(y="Cell composition",x= NULL)+  
  #geom_point(aes(fill = group_stat4))+
  geom_boxplot(aes(fill = group_stat4),outlier.size = 0.05,size=0.02,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("skyblue","pink") )+ #"#1CB4B8", "#EB7369"
  guides(fill = guide_legend(title = 'STAT4'))+
  theme_bw() + mytheme + 
  #ylim(0,0.1)+
  stat_compare_means(aes(group =  group_stat4),
                     label = "p.signif",
                     label.y = 0.03,
                     size=5,#0.6
                     method = "wilcox.test",
                     hide.ns = T)


##4.3 cibersortX----
##4.4 ssGSEA 免疫浸润----
#读取基因集:http://events.jianshu.io/p/da1dbbf793ab

gene_set<-read.csv("./immune/mmc3.csv",skip = 2)[, 1:2]
head(gene_set)
list<- split(as.matrix(gene_set)[,1], gene_set[,2])        

library(GSVA)
ssgsea<-gsva(exp1_matrix, list, method = "ssgsea",  #FPKM  exp1_matrix1是否需要LOG2标准化？-0.4-0.5
             min.sz = 10)
ssgsea[1:3, 1:2]#matrix  # 标准化到0-1 ？


annotation_col2 = data.frame( STAT4=clinical_1$STAT4_group,Riskscore=clinical_1$Group)#变成一个表
rownames(annotation_col2)=clinical_1$sample#colnames(exprSet)#再取个名字
ID_high_low_STAT4 <- c(rownames(subset(annotation_col2,annotation_col2$STAT4 =="Low")),#502
                 rownames(subset(annotation_col2,annotation_col2$STAT4 =="High")) )#499  排序取出

library(pheatmap)
ssgsea_pheat <-pheatmap(ssgsea[,ID_high_low_STAT4],
         show_colnames = F,show_rownames = T,
         cluster_rows = T,cluster_cols = F,
         annotation_col = annotation_col2, #1
         #cellwidth=1,cellheight=3,
         fontsize=6,#,
         gaps_col = 502,color = colorRampPalette(c("skyblue", "white", "pink") )(10) ,#c("navy", "white", "firebrick3")
         #,filename = 'ssgsea.pdf',width = 8,height = 4,family="serif"
         )
ssgsea_pheat
save_pheatmap_pdf(ssgsea_pheat, "./output/ssgsea_pheat_1.pdf",8,2.5)
#dev.off()
#write.csv(ssgsea,"./220719/Figure4A1_ssgaes.csv")

#合并免疫浸润
cibor_ssgsea_estim_df <- merge(ssgsea1,TME_data_1,by="sample")#1001-53
#estimatescores.final[,1:3] <- as.data.frame(sapply(estimatescores.final[,1:3], as.numeric))#
cibor_ssgsea_estim_df <-merge(estimatescores.final,cibor_ssgsea_estim_df,by="sample")#1001-56

rownames(cibor_ssgsea_estim_df) <-cibor_ssgsea_estim_df$sample
cibor_ssgsea_estim_df <-cibor_ssgsea_estim_df[,-1]#1001-52

#cibor_ssgsea_estim heatmap
library(pheatmap)
pheat<-pheatmap(t(estimatescores.final[ID_high_low_STAT4,][1:3]),#t(TME_data[ID_high_low_STAT4,][1:22]),#ssgsea[,ID_high_low_STAT4],# #t(estimatescores.final[ID_high_low_STAT4,][1:3]),
         show_colnames = F,show_rownames = T,
         cluster_rows = T,cluster_cols = F,
         annotation_col = annotation_col2[1], #[1]
         ##cellwidth=1,cellheight=3,
         fontsize=6,#,
         gaps_col = 502,          #"#db6989"朱红   #"#f47720"橙色
         color=colorRampPalette( c("navy", "white", "firebrick3"))(10)  # CIBORT colorRampPalette(brewer.pal(1, "PiYG"))(10)#color =colorRampPalette( c("#459943","#db6989") )(10) #colorRampPalette( c( "#f47720" , "white", "#459943") )(10) #c("navy", "white", "firebrick3") #c("skyblue", "white", "pink")
         #,filename = 'ssgsea.pdf',width = 8,height = 4,family="serif"
)  #数据差异大，分开画热图，手动拼,每个配色一种方案

#save_pheatmap_pdf(pheat, "./output/2.pdf",8,0.6)#2.5
#write.csv(t(estimatescores.final[ID_high_low_STAT4,][1:3]),"./220719/Figure4A2_estimate.csv")

##ssGSEA  boxplot----
##1. 添加CIBORSORT 名字 
TME_data_1 <- TME_data
colnames(TME_data_1)[1:22]<- paste0("CIBORSORT: ",colnames(TME_data_1)[1:22])
#添加ssGSEA 名字 
ssgsea1 <-as.data.frame(t(ssgsea))#1001-28
#colnames(ssgsea1)<-paste0("ssGSEA: ",colnames(ssgsea1) )
ssgsea1$sample <- rownames(ssgsea1)
##BOXLPOT
ssgsea1<-merge(ssgsea1,clinical_1[,c(1,96,129)],by="sample") #合并临床数据到SSGSEA
rownames(ssgsea1)<-ssgsea1[,1]
ssgsea1 <-ssgsea1[,-1]
#变为长数据
library(reshape2)
ssgsea1_New = melt(ssgsea1)# 融合数据:宽数据变成长数据
colnames(ssgsea1_New)=c("group_stat4","Group","Celltype","Score")  #设置行名 #无sample?
head(ssgsea1_New)
# 按免疫细胞占比中位数排序绘图（可选）
library(dplyr)
plot_order_1 = ssgsea1_New[ssgsea1_New$group_stat4=="High",] %>% 
  group_by(Celltype) %>% 
  summarise(m = median(Score)) %>% #
  arrange(desc(m)) %>% 
  pull(Celltype)

ssgsea1_New$Celltype = factor(ssgsea1_New$Celltype,levels = plot_order_1)

if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 90, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
library(ggplot2);library(ggpubr)
box_TME <- ggplot(ssgsea1_New, aes(x = Celltype, y = Score))+ ##Celltype/T cells follicular /helper Neutrophils
  labs(y="ssGSEA score",x= NULL)+  
  geom_boxplot(aes(fill = group_stat4),outlier.size = 0.02,size=0.02,outlier.alpha = 0)+ 
  #geom_bar(stat = "identity")+
  scale_fill_manual(values = c("skyblue","pink") )+ #"#1CB4B8", "#EB7369" "skyblue","pink"
  guides(fill = guide_legend(title = 'STAT4'))+
  theme_bw()+ mytheme + 
  ylim(-0.4,0.65)+
  stat_compare_means(aes(group =  group_stat4),
                     label = "p.signif",
                     label.y = 0.6,size=3,#0.6
                     method = "wilcox.test",
                     hide.ns = T)#+
# geom_signif(annotations = c(""),color="black",size = 0.1, #手动加括号
#             y_position = c(0.12,0.12),
#             xmin = c(8.81,21.81),
#             xmax = c(9.188,22.188),
#             tip_length = c(0.01,0.02,0.1,0.1))


box_TME

ggsave("./immune/ssGSEA_box_1.pdf",box_TME,height=5,width=10,family="serif") 
#write.csv(ssgsea1_New,"./220719/Figure3C_ssgsea_boxplot.csv")

#分组柱状图
library(ggplot2);library(ggpubr)
box_TME <- ggplot(ssgsea1_New, aes(x = Celltype, y = Score,fill= group_stat4))+ #,fill= group_stat4# Celltype/T cells follicular /helper Neutrophils
  labs(y="ssGSEA score",x= NULL)+  
  geom_boxplot(aes(fill = group_stat4),outlier.size = 0.02,size=0.02,outlier.alpha = 0)+ 
  #geom_point(aes(group = group_stat4) )+
  #geom_bar(stat = "identity")+
  #stat_summary(mapping=aes(group=group_stat4),geom = 'point',fun = 'mean',cex=1.3,width=.6)+
  scale_fill_manual(values = c("skyblue","pink") )+ #"#1CB4B8", "#EB7369" "skyblue","pink"
  guides(fill = guide_legend(title = 'STAT4'))+
  theme_bw()+ mytheme + 
  ylim(-0.4,0.65)+
  stat_compare_means(aes(group =  group_stat4),
                     label = "p.signif",
                     label.y = 0.6,size=3,#0.6
                     method = "wilcox.test",
                     hide.ns = T)
box_TME
##2. 重做数据为MEAN SD 画分组柱状图
summary(ssgsea1[1:28])#类似此
#https://blog.csdn.net/qq_42804713/article/details/124548178
library(dplyr)

ssgsea2 <- data.frame()
for (i in colnames(ssgsea1[1:28])) {
  ssgsea_group <-subset(ssgsea1,ssgsea1$STAT4_group=="High")
  ssgsea2 <-rbind(ssgsea2,
                  cbind(i,summarise(ssgsea_group ,M=mean(ssgsea_group[,i] ),S=sd(abs(ssgsea_group[,i]) )) )
  )
}
colnames(ssgsea2) <-c("Celltype","MEAN","SD")
ssgsea2$group <- "High"
#高低组分开计算再合并
ssgsea3 <- data.frame()
for (i in colnames(ssgsea1[1:28])) {
  ssgsea_group <-subset(ssgsea1,ssgsea1$STAT4_group=="Low")
  ssgsea3 <-rbind(ssgsea3,
                  cbind(i,summarise(ssgsea_group ,M=mean(ssgsea_group[,i] ),S=sd(ssgsea_group[,i] ) ) )
  )
}
colnames(ssgsea3) <-c("Celltype","MEAN","SD")
ssgsea3$group <- "Low"

ssgsea_df <- rbind(ssgsea2,ssgsea3)#56-4
#降序？

#绘图柱状图
library(ggplot2)
ggplot(ssgsea_df, aes(x=Celltype, y=MEAN, fill=group)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=MEAN, ymax=MEAN+SD),size=0.75,width=0.08,position=position_dodge(0.5))+
  #geom_col(aes(fill=group), position = position_dodge2(preserve = 'single')) +
  labs(y="ssGSEA Score",title=NULL)+
  theme_bw()+
  mytheme  #放弃


#DF and plot
#https://blog.csdn.net/qq_50522851/article/details/122061089
#https://zhuanlan.zhihu.com/p/487192059

#5. 药物预测----
#5.1 oncoPredict GDSC2----
#教程https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247507376&idx=1&sn=6993a13200be452c9dfa527e9fdc8ae4&scene=21#wechat_redirect
##a.安装包#包指导手册：https://cran.r-project.org/web/packages/oncoPredict/vignettes/calcPhenotype.html
#install.packages("oncoPredict")
# install.packages('./output/oncoPredict_0.2.tar.gz',repos=NULL, type="source")
# install.packages('./output/oncoPredict_0.2/oncoPredict',repos=NULL, type="source")
# lapply(c( 'preprocessCore', 'TxDb.Hsapiens.UCSC.hg19.knownGene'), install.packages)#安装依赖
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
##b.
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
#dir='./DataFiles/Training Data/'
dir="./output/oncoPredict_0.2/oncoPredict/vignettes/"
# GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr_short.rds'))
# GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
# GDSC2_Res <- exp(GDSC2_Res) #
trainingExprData = readRDS(file=file.path(dir,"CTRP2_Expr.rds") )
trainingExprData = readRDS(file = "CTRP2_Expr.rds")
trainingPtype = readRDS(file=file.path(dir,"CTRP2_Res.rds"))
##c. 输入测试数据
# testExpr<- GDSC2_Expr[,sample(1:ncol(GDSC2_Expr),10)]
# testExpr[1:4,1:4]  
# colnames(testExpr)=paste0('test',colnames(testExpr))
# dim(testExpr)  
##d.运行函数
exp1_matrix <-as.matrix(exp1)
calcPhenotype(trainingExprData = trainingExprData,#GDSC2_Expr,
              trainingPtype = trainingPtype,#GDSC2_Res,
              testExprData = exp1_matrix,#testExpr, #测试集！输入exp1进行计算
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' ) #3min
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = exp1_matrix,#testExpr, #测试集！输入exp1进行计算
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              cc=TRUE,
              removeLowVaringGenesFrom = 'rawData' )

##预测结果
library(data.table)
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
testPtype[1:4, 1:4]

#读取结果
GDSC2_cor  <-read.table("./calcPhenotype_Output/cors.txt")#338 gene 198 durg
GDSC2_pvalue  <-read.table("./calcPhenotype_Output/pvalues.txt")#338 gene 198 durg
GDSC2 <-testPtype
#与分组合并
colnames(GDSC2)[1] <-"sample" #1001-199
GDSC2 <-merge(GDSC2,clinical[,c(1,95,128)],by="sample") #201
GDSC2 <-GDSC2[-c(155,565,592,700),]#201-删除异常值 #test <- subset(GDSC2,GDSC2$Camptothecin_1003 > 6000)
rownames(GDSC2) <- GDSC2$sample;GDSC2<-GDSC2[,-1]

GDSC2 <- as.data.frame(zoo::na.fill(GDSC2,0) ) #NA替换为0       #log2(0)=inf
GDSC2[,1:198] <-as.data.frame(sapply(GDSC2[,1:198],as.numeric))#997-200
rownames(GDSC2) <- testPtype[-c(155,565,592,700),]$V1
  
##绘图
boxplot(GDSC2$sample,log2(GDSC2$Camptothecin_1003) )
ggplot2::ggplot(aes(GDSC2$Group,log2(GDSC2$Camptothecin_1003),fill=Group),data=GDSC2)+
  geom_boxplot( ) 
#小提琴图
cbPalette <- c("#E69F00", "#CC79A7", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#999999","#0072B2","#D55E00")
library(ggplot2);library(ggpubr)
ggplot(data=GDSC2,aes(x=Group,y=log2(GDSC2$Cisplatin_1005),color=Group))+  #fill
  geom_point(alpha=0.5,size=1.5,
             position=position_jitterdodge(jitter.width = 0.35,
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.05,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.25)+
  scale_color_manual(values = c("#CC79A7", "#56B4E9") )+ #color
  
  theme_bw() +ylab("log2(IC50)")+xlab("")+
  #theme(legend.position="none") 
  theme(axis.title  = element_text(size=16), axis.text = element_text(size=12),legend.position ="top",  
        panel.grid = element_blank(),legend.key = element_blank() ) +
  stat_compare_means(#aes(group = Group) ,
                     comparisons=list(c("High","Low")),#my_comparisons,
                     
                     label = "p.signif",#"p.format p.signif 9.4e-12
                     method = "wilcox",
                     show.legend= F,#删除图例中的"a"
                     label.x=1.5,bracket.size=0.1,vjust=0.2,
                     #label.y = max(log2(GDSC2$Camptothecin_1003)),
                     hide.ns = T,size=5)
#ggsave("./output/drug_test.pdf",family="serif",width = 5,height =5 )


#批量计算分组显著性差异
wilcox <- wilcox.test(Camptothecin_1003~Group,data=GDSC2)
wilcox$p.value

GDSC2_pvalue_Group <- data.frame()
 for (i in colnames(GDSC2)[2:199] ){   #
   wilcox <- wilcox.test(as.numeric(GDSC2[,i])~Group,data=GDSC2) #(x=factor(GDSC2$Group),y=GDSC2[,i])# 
   #pvalue<-cbind(i,wilcox$p.value)
   GDSC2_pvalue_Group<-rbind(GDSC2_pvalue_Group,cbind(i,wilcox$p.value))
 }
colnames(GDSC2_pvalue_Group)<- c("drug","pvalue") #198
#https://www.jianshu.com/p/a68ce9459edd
GDSC2_pvalue_Group$p_adj <- round(p.adjust(GDSC2_pvalue_Group$pvalue,"BH"),3)#校正Pvalue  BH
GDSC2_pvalue_Group <- subset(GDSC2_pvalue_Group,GDSC2_pvalue_Group$pvalue <0.001) #p_adj <0.001=14,0.01=28!, 0.05=40 ##挑选Pvalue<0.05 的药物
GDSC2_pvalue_Group$drug<- gsub("\\_.*","",GDSC2_pvalue_Group$drug)


##提取STAT相关性/手动计算
GDSC2_cor_stat4 <- subset(GDSC2_cor,rownames(GDSC2_cor)=="STAT4")#338个基因均不在其中

outTab_GDSC2 <- data.frame() #首先，新建一个空的数据框，用于保存后续分析结果
  Gene <- "STAT4"
  x <- as.numeric(exp1_matrix[Gene,]) 
  #对药物循环 
  for (Drug in colnames(testPtype[,2:199]) ){ 
    y <- as.numeric(testPtype[,2:199][,Drug]) 
    corT <- cor.test(x,y,method= "pearson") 
    cor <- corT$estimate 
    pvalue <- corT$p.value 
    if( pvalue < 0.05) { 
      outVector <- cbind(Gene,Drug,cor,pvalue) #as.data.frame( 
      #print(outVector)
      #outVector <- as.data.frame(outVector) 
      outTab_GDSC2<- rbind(outTab_GDSC2,outVector)
    } 
  } 
dim(outTab_GDSC2)#STAT4与药物相关性pearson #75-4   p均小于0.01
outTab_GDSC2$p_adj <- round(p.adjust(outTab_GDSC2$pvalue,"BH"),4)# cor >0.4  ~6   75-5
outTab_GDSC2 <- outTab_GDSC2[abs(as.numeric(outTab_GDSC2$cor)) >0.3,]#13 cor >0.3
#outTab_GDSC2$Drug<-gsub("\\_.*","",outTab_GDSC2$Drug)

#批量绘图函数-cor >0.3
GDSC2[,c("MN-64_1854","AZD8186_1918","AZD6482_2169","STAT4_group")]

##方法二，箱线图
plotList_3 <- list()
corPlotNum <- 16
if(nrow(outTab)<corPlotNum){
  corPlotNum=nrow(outTab_GDSC2)
}  #同样，定义一个空的列表plotList_2用于保存输出结果
for(i in 1:3){ #1:corPlotNum
  Gene <- Gene
  Drug <- outTab_GDSC2[abs(as.numeric(outTab_GDSC2$cor)) > 0.4 ,][i,2] #3 个#outTab_GDSC2[i,2]
  
  #x <- as.numeric(exp[Gene,])
  #y <- as.numeric(drug[Drug,])
  df1 <- as.data.frame(cbind(Gene,GDSC2[,Drug]))
  colnames(df1)[2] <- "IC50"
  df1$IC50 <- as.numeric(df1$IC50)
  #df1$group <- ifelse(df1$x > median(df1$x), "high", "low")
  df1$group <- GDSC2$STAT4_group #subset(GDSC2,colnames(GDSC2)== Drug)$STAT4_group
  compaired <- list(c("Low", "High"))
  p1 <-ggplot(df1,aes(x =group, y =log2(IC50),color =group))+ #,palette = c("firebrick3","skyblue")
    labs(title= gsub("\\_.*","",Drug),ylab="Log2(IC50)")+
    geom_point(alpha=0.2,size=1.5,
               position=position_jitterdodge(jitter.width = 0.45,
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
    geom_boxplot(alpha=0.5,width=0.55,
                 position=position_dodge(width=0.8),
                 size=0.05,outlier.colour = NA)+
    # geom_violin(alpha=0.2,width=0.9,
    #             position=position_dodge(width=0.8),
    #             size=0.25)+
    scale_color_manual(values = c("#CC79A7", "#56B4E9") )+ #color

    # stat_compare_means(comparisons = compaired,vjust=0.4,
    #                    method = "wilcox.test",size=4,   #设置统计方法
    #                    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
    #                                     symbols = c("***", "**", "*", "ns")))+
    #coord_cartesian(ylim = c( NA,max(log2(df1$IC50))+1.5 )  )+   #否则*显示不全,log2(0)=inf
    ylim(c( NA,max(log2(df1$IC50))+1.5 ))+
    stat_compare_means(#aes(group = Group) ,
      comparisons=list(c("High","Low")),#my_comparisons,
      label = "p.signif",#"p.format p.signif 9.4e-12
      method = "wilcox",
      show.legend= F,#删除图例中的"a"
      label.x=1.5,bracket.size=0.1,vjust=0.1,
      #label.y = max(log2(GDSC2$Camptothecin_1003)),
      hide.ns = T,size=4)+
    theme_bw()+
    theme(axis.title  = element_text(size=12), axis.text = element_text(size=12),legend.position ="none",  
          panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() )
   
  plotList_3[[i]]=p1
}
ggarrange(plotlist=plotList_3,nrow=1,ncol=3) #STAT4_group


#ggsave("./output/GDSC2_boxplot_cor0.3_13.pdf",family="serif",width = 8,height = 9)
#ggsave("./output/GDSC2_boxplot_cor0.4_3.pdf",family="serif",width = 6,height = 2)

##棒棒糖图lollipop
library(ggpubr)
outTab_GDSC2[,3:4] <- as.data.frame(sapply(outTab_GDSC2[,3:4],as.numeric ))
outTab_GDSC2$'-log10(pvalue)' <- -log10(outTab_GDSC2$pvalue)

outTab_GDSC2_1 <-outTab_GDSC2
outTab_GDSC2_1$Drug <- gsub("\\_.*","",outTab_GDSC2_1$Drug)
ggdotchart(outTab_GDSC2_1, x = "Drug", y = "cor",
           color ="-log10(pvalue)",
           sorting = "descending",#排序“ascending”, “descending”, “none”
           add.params = list(color = "lightgray", size = 1.0), #画棒棒
           add = "segments",#画棒棒
           dot.size = "cor",#-log10(pvalue), #"cor"
           xlab="",ylab="Correlation Coefficient",                                
           ggtheme = theme_bw()) + 
  coord_flip()+  #翻转坐标轴
  #scale_color_gradient(low = "skyblue",high="orange")+  #颜色范围
  font("x.text", size = 13, vjust = 0.1,angle = 0)+  #X坐标轴
  font("y.text",size = 13,face = "bold")+  #坐标加粗
  geom_hline(yintercept=c(-0.4,0.4), linetype="dashed",size=0.1)+ #添加虚线
  scale_color_continuous(low = "skyblue",high="orange"
                      #range  = c(0.1, 10), 
                      #limits = c(0, 5e-22), 
                      #breaks = c( 5e-60, 5e-23, 3e-22,5e-22) 
                      )+
  geom_hline(yintercept=c(0),size=0.1)+
  theme(panel.grid = element_blank(),axis.text = element_text( ) )
#ggsave("./output/GDSC2_lollipop_13_1.pdf",family="serif",width = 6,height = 4)
#write.csv(outTab_GDSC2_1,"./220719/Figure4D_oncopredict.csv")
#  # 保存三个显著性 #write.csv(GDSC2[,c("MN-64_1854","AZD8186_1918","AZD6482_2169","STAT4_group")],"./220719/Figure4E_oncopredict_boxplot.csv")
#改pvalue刻度范围显示

#分组STAT4!!!!



#5.2 CellMiner https://mp.weixin.qq.com/s/re0Rw2c6h6mPlFXkB7TDCQ----
##A.数据收集及处理：D:\研究生科研\生物信息学\STAT4\Subgroup_STAT4\3.Enrichment\cellminer
library(readxl) 
rt1 <- read_excel( path= "D:/研究生科研/生物信息学/STAT4/Subgroup_STAT4/3.Enrichment/cellminer/DTP_NCI60_ZSCORE.xlsx", skip= 7)
colnames(rt1) <-rt1[1,]
rt1<-rt1[-1,-c(67,68)]
table(rt1$`FDA status`)  #筛选药物标准-FDA approved-218

rt1 <- rt1[rt1$`FDA status` %in% c("FDA approved"),]  #218-66
rt1 <- rt1[,-c( 1, 3: 6)] #删除无关列 218-61
#write.table(rt1, file = "./cellminer/drug.txt",sep = "\t",row.names = F,quote = F)

##B. 基因表达数据的准备-示例数据
rt2<- read_excel(path = "D:/研究生科研/生物信息学/STAT4/Subgroup_STAT4/3.Enrichment/cellminer/RNA__RNA_seq_composite_expression.xls", skip = 9) 
colnames(rt2) <- rt2[ 1,] 
rt2 <- rt2[- 1,-c( 2: 6)] #230808-61
#write.table(rt2, file = "./cellminer/geneExp.txt",sep = "\t",row.names = F,quote = F)

##重新读取药物输入文件
library(impute) 
library(limma)
#rt <- read.table( "./cellminer/drug.txt",sep= "\t",header=T,check.names=F) 
rt <- as.matrix(rt1) 
rownames(rt) <- rt[, 1] 
drug <- rt[, 2:ncol(rt)] 
dimnames <- list(rownames(drug),colnames(drug)) 
data <- matrix( as.numeric( as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)#218-60

mat<- impute.knn(data);drug<- mat$data;drug <- avereps(drug) #补齐NA值

##重新读取表达输入文件
exp<- read.table( "./cellminer/geneExp.txt", sep= "\t", header=T, row.names = 1, check.names=F) 
dim(exp) #23808-60   19551-1001
exp[1: 4, 1: 4]
#表达谱只能是NCI60？

##提取特定基因表达
#准备基因名称
genelist <- c("STAT4")#c("FANCD2","BRCA1","ABCC1","TP53","EGFR")
exp<- exp[genelist,] #5-60 #注意矩阵名称exp  与exp1内容！ 1-60
##药物敏感性计算
outTab <- data.frame() #首先，新建一个空的数据框，用于保存后续分析结果

for (Gene in row.names(exp) ){ 
  x <- as.numeric(exp[Gene,]) 
  #对药物循环 
  for (Drug in row.names(drug)){ 
    y <- as.numeric(drug[Drug,]) 
    corT <- cor.test(x,y,method= "pearson") 
    cor <- corT$estimate 
    pvalue <- corT$p.value 
    if( pvalue < 0.05) { 
      outVector <- cbind(Gene,Drug,cor,pvalue) #as.data.frame( 
  #print(outVector)
  #outVector <- as.data.frame(outVector) 
      outTab <- rbind(outTab,outVector)
    } 
  } 
} #55-4
outTab <- outTab[order( as.numeric( as.vector(outTab$pvalue))),]   #p <0.05 17  0.01-4
outTab$p_adj <- round(p.adjust(outTab$pvalue,"BH"),4)#校正P value 0.01-0.05
cellminer <- subset(outTab,abs(as.numeric(outTab$cor))>0.3 )#p <0.05 abs(cor) >0.3  - 7个
#write.table(outTab, file= "./cellminer/drugCor.txt", sep= "\t", row.names=F, quote=F)
#write.table(cellminer, file= "./cellminer/cellminer_7.txt", sep= "\t", row.names=F, quote=F)
#write.csv(cellminer,"./220719/Figure4B_cellminer.csv")

##棒棒糖图:cor pvalue
library(ggpubr)
outTab[,3:4] <- as.data.frame(sapply(outTab[,3:4],as.numeric ))
ggdotchart(outTab, x = "Drug", y = "cor",
           color = "pvalue",
           sorting = "descending",#排序“ascending”, “descending”, “none”
           add.params = list(color = "lightgray", size = 1.0), #画棒棒
           add = "segments",#画棒棒
           dot.size = "cor", 
           xlab="",ylab="Coefficient",  
           #font.label = list(size = 14, face = "bold", color ="red"),
           ggtheme = theme_bw()) + 
  coord_flip()+  #翻转坐标轴
  scale_color_gradient(low = "skyblue",high="orange")+  #颜色范围
  font("x.text", size = 13, vjust = 0.1,angle = 0)+  #X坐标轴
  font("y.text",size = 13,face = "bold")+  #坐标加粗
  geom_hline(yintercept=c(-0.3,0.3), linetype="dashed",size=0.1)+ #添加虚线
  geom_hline(yintercept=c(0),size=0.1)+
  theme(panel.grid = element_blank(),axis.text = element_text( ) )
#ggsave("./cellminer/cellminer_lollipop_7_2.pdf",family="serif",width = 5,height = 4)

##可视化
library(ggplot2) ;library(ggpubr)
#定义一个空的列表plotList_1用于保存输出结果。提取分析结果中最显著的前16个结果；当然，如果结果outTab的行数少于corPlotNum的话，则将其行数赋值给corPlotNum。
plotList_1 <- list()
corPlotNum <- 8
outTab <-cellminer #选出7个展示

if(nrow(outTab)<corPlotNum){
  corPlotNum=nrow(outTab)
}
#使用ggplot()函数结合for循环，逐个绘制散点图，进行可视化展示。

for(i in 1:corPlotNum){
  Gene <- outTab[i,1]
  Drug <- outTab[i,2] #
  x <- as.numeric(exp[Gene,])
  y <- as.numeric(drug[Drug,])
  cor <- sprintf("%.03f",as.numeric(outTab[i,3]))
  pvalue=0
  if(as.numeric(outTab[i,4])<0.001){
    pvalue="p<0.001"
  }else{
    pvalue=paste0("p=",sprintf("%.03f",as.numeric(outTab[i,4])))
  }
  df1 <- as.data.frame(cbind(x,y))
  p1=ggplot(data = df1, aes(x = x, y = y))+
    geom_point(size=1.5,color="skyblue",alpha=0.5)+
    geom_smooth(method = lm, se = T,size=0.1,color="firebrick3",fill="grey90")+
    #stat_smooth(method="lm",se=FALSE, formula=y~x,span = 0.1)+
    labs(x="",y="",title = paste0(Drug),#paste0(Gene,", ",Drug),
         subtitle = paste0("Cor=",cor,", ",pvalue))+
    theme(axis.ticks = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank())+
    theme_bw()
  plotList_1[[i]]=p1
}
#通过for循环，结合使用ggpubr包中的ggboxplot()函数，逐个绘制箱线图，并使用stat_compare_means()函数进行两组间统计分析，从而进行可视化展示。
nrow <- ceiling(sqrt(corPlotNum))
ncol <- ceiling(corPlotNum/nrow)
ggarrange(plotlist=plotList_1,nrow=nrow,ncol=ncol)
#ggsave("./cellminer/cellminer_point.pdf",family="serif",width = 10,height = 10)
#ggsave("./cellminer/cellminer_point_7_1.pdf",family="serif",width = 8,height = 8)

##方法二，箱线图
plotList_2 <- list()
corPlotNum <- 2#8
if(nrow(outTab[c(5,7),])<corPlotNum){  #outTab[c(5,7),] outTab
  corPlotNum=nrow(outTab[c(5,7),])
}  #同样，定义一个空的列表plotList_2用于保存输出结果
for(i in 1:corPlotNum){
  Gene <- outTab[i,1]
  Drug <- outTab[c(5,7),][i,2]#Vinblastine TYROTHRICIN
  x <- as.numeric(exp[Gene,])
  y <- as.numeric(drug[Drug,])
  df1 <- as.data.frame(cbind(x,y)) #df2 <- as.data.frame(cbind(as.numeric(exp["STAT4",]),as.numeric(drug["Vinblastine",]),as.numeric(drug["TYROTHRICIN",]) ) )
  colnames(df1)[2] <- "IC50"
  df1$group <- ifelse(df1$x > median(df1$x), "high", "low")
  compaired <- list(c("low", "high"))
  p1 <- 
    # ggboxplot(df1,
    #               x = "group", y = "IC50",
    #               fill = "group", palette = c("firebrick3","skyblue"),#c("#00AFBB", "#E7B800"),
    #               size = 0.5,alpha=0.5,
    #               #xlab = paste0("The_expression_of_", Gene),
    #               ylab = paste0("IC50_of_", Drug)) +
    ggplot(df1,aes(x =group, y =IC50,color =group))+ #,palette = c("firebrick3","skyblue")
    labs(title= Drug)+
    geom_point(alpha=0.5,size=1.5,
               position=position_jitterdodge(jitter.width = 0.45,
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
    geom_boxplot(alpha=0.2,width=0.55,
                 position=position_dodge(width=0.8),
                 size=0.05,outlier.colour = NA)+
    # geom_violin(alpha=0.2,width=0.9,
    #             position=position_dodge(width=0.8),
    #             size=0.25)+
    scale_color_manual(values = c("#CC79A7", "#56B4E9") )+ #color
    coord_cartesian(ylim = c(min(df1$IC50),max(df1$IC50)+0.5 ))+
    stat_compare_means(comparisons = compaired,
                       method = "wilcox.test",size=5,   #设置统计方法
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")))+
      theme_bw()+
      theme(axis.title  = element_text(size=12), axis.text = element_text(size=12),legend.position ="none",  
                       panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() )
  plotList_2[[i]]=p1
}

ggarrange(plotlist=plotList_2,nrow=1,ncol=2)#STAT4！分组！趋势符合+显著性-2个
ggsave("./cellminer/cellminer_box_7_3.pdf",family="serif",width = 4,height = 2)# 8 8 12 2

df2 <- as.data.frame(cbind(as.numeric(exp["STAT4",]),as.numeric(drug["Vinblastine",]),as.numeric(drug["TYROTHRICIN",]) ) )
colnames(df2) <- c("STAT4","Vinblastine","TYROTHRICIN")
df2$group <- ifelse(df2$STAT4 > median(df2$STAT4), "high", "low")
#plot(df2$Vinblastine)plot(log2(df2$Vinblastine))
write.csv(df2,"./220719/Figure4C_cellminer_boxplot.csv")

#5.3 RNAactDrug----
#RNAactDrug online database: http://bio-bigdata.hrbmu.edu.cn/RNAactDrug/index.jsp
#D:\研究生科研\生物信息学\STAT4\Subgroup_STAT4\3.Enrichment\output\RNAactDrug.csv



#6.蛋白质组富集分析----
##6.1 读取数据及合并临床数据
RPPA_data_BRCA <- read.csv("D:\\研究生科研\\生物信息学\\STAT4\\Subgroup_STAT4\\2.Prognosis\\RPPA\\RPPA_data_BRCA.csv")[,-1]
rownames(RPPA_data_BRCA) <- RPPA_data_BRCA[,1]
RPPA_data_BRCA <- RPPA_data_BRCA[,-1]
colnames(RPPA_data_BRCA)<- gsub("\\.","-",colnames(RPPA_data_BRCA))  #226 937sample
##6.2 填充NA
library(impute)
RPPA_dimnames <- list(rownames(RPPA_data_BRCA),colnames(RPPA_data_BRCA)) 
RPPA_matrix <- matrix( as.numeric( as.matrix(RPPA_data_BRCA)),nrow=nrow(RPPA_data_BRCA),dimnames=RPPA_dimnames)
RPPA_data_BRCA1 <- impute.knn(as.matrix(RPPA_data_BRCA))
RPPA_data_BRCA1 <- as.data.frame(RPPA_data_BRCA1$data) #RPPA_data_BRCA1 为补齐NA后的数据框
##6.3 转置合并临床数据-RPPA_data_clinical
RPPA_data_clinical <- as.data.frame(t(RPPA_data_BRCA1))
RPPA_data_clinical$sample <- rownames(RPPA_data_clinical)
RPPA_data_clinical <- merge(RPPA_data_clinical,clinical_1[,1:132],by="sample")#801 sample-358   1：226
rownames(RPPA_data_clinical)<-RPPA_data_clinical[,1]
RPPA_data_clinical <-RPPA_data_clinical[,-1]
##6.4 绘制分组热图
library(pheatmap)

annotation_col_RPPA <- data.frame( Riskscore= RPPA_data_clinical$Group,STAT4=RPPA_data_clinical$STAT4_group)#,PAM50=RPPA_data_clinical$PAM50) #
row.names(annotation_col_RPPA) <- rownames(RPPA_data_clinical)
annotation_colors_RPPA <- list(Riskscore=c(High="#D5B32B", Low="#2BD54D"),STAT4=c(High="pink", Low="skyblue"))#调色板 https://www.bejson.com/ui/getcolor/

#cols = c("#999999","#FF0099", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#990000","#9900cc","#66FF66","#663300","#0000FF","#CC0033","#FF0000","#000099","#660066","#333333","#FFCCCC","#993366","#33CC33","#000099","#CC9900")#24
#cols1 = c("#57b4e9","#f0e442","#999999","#019e74","#e69f00","#0073b0","#d55e00","#cc7aa5")#文献配色
  #           蓝       橙色    灰          绿
#热图顺序排
pheat1 <-pheatmap(t(RPPA_data_clinical[,1:226]),show_rownames = T,show_colnames = F,cluster_cols = T,cluster_rows = T,annotation_col  = annotation_col_RPPA,
         #cutree_col = 2,
         annotation_colors= annotation_colors_RPPA,fontsize = 8,fontsize_col=1,
         color = colorRampPalette(c("navy", "white", "firebrick3") )(10) ) # c("blue", "white", "red"))(10)
pheat1
#save_pheatmap_pdf(pheat1, "./output/pheatmap_RPPA_1.pdf",16,20)
#其他绘图教程：https://zhuanlan.zhihu.com/p/475394286

##6.5 RPPA富集分析
#https://blog.csdn.net/tu__zi/article/details/105288555
#6.5.1 RPPA limma差异分析
#a.设置分组-STAT4
group_list_1 <-RPPA_data_clinical$STAT4_group

design_1 <- model.matrix(~0+factor(group_list_1))#801 2
colnames(design_1)=levels(factor(group_list_1))
rownames(design_1)=rownames(RPPA_data_clinical)
head(design_1) #分组  801 High 394; Low   407

#b. limma analysis
library(limma)
compare <- makeContrasts(High - Low, levels = design_1)
fit <- lmFit(t(RPPA_data_clinical[,1:226] ), design_1)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff_RPPA <- topTable(fit3, coef = 1, number = 10000)
head(Diff_RPPA ) #226 6

Diff_RPPA_pvalue <- subset(Diff_RPPA,Diff_RPPA$P.Value <0.05)#96 6
#是否筛选LOGFC？或者直接T/wilcox检验？

#准备GSEA数据
#排序：将基因按照logFC进行从高到低排序，只需要基因列和logFC即可
Diff_RPPA_pvalue_sort <- Diff_RPPA_pvalue %>% 
  arrange(desc(logFC))
head(Diff_RPPA_pvalue_sort)#96
Diff_RPPA_pvalue_sort$protein <-toupper(rownames(Diff_RPPA_pvalue_sort))
#fix(Diff_RPPA_pvalue_sort)

#
geneList_RPPA = data.frame(logFC =Diff_RPPA_pvalue_sort$logFC,protein=Diff_RPPA_pvalue_sort$protein)
geneList_RPPA_a<-bitr(geneList_RPPA$protein,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")#62.5 % fail to map
geneList_RPPA_1 <-merge(geneList_RPPA_a,geneList_RPPA,by="SYMBOL") #取出匹配的部分-34
geneList_RPPA_1 <-geneList_RPPA_1[order(geneList_RPPA_1$logFC,decreasing = T),]#排序
#kegg
geneList_RPPA_2 <- as.numeric(geneList_RPPA_1$logFC);names(geneList_RPPA_2)<-geneList_RPPA_1$ENTREZID

GSEA_RPPA <- GSEA(geneList_RPPA_2,pvalueCutoff = 0.1,TERM2GENE=Hs_df[,c(1,2)])#见上
GSEA_RPPA_df <-GSEA_RPPA@result #5个太少，放弃

##6.5.2 注释及富集分析
#A.注释蛋白名称
## 取得蛋白名
Gene <- rownames(Diff_RPPA_pvalue)  #Diff_RPPA_pvalue
Gene <- gsub("\\_.*","",Gene)
Gene <-toupper(Gene) #全部大写
#fix(Gene)#手动修改


## 将蛋白名转化为基因ID
library(clusterProfiler)
library(org.Hs.eg.db)
#将SYMBOL转化为ENTREZID
a <- bitr(Gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")#59.79% FAIL  #35

# B通路富集
library(clusterProfiler)#重新安装：4.4.4 https://www.jianshu.com/p/ef79aa155800
#detach("package:clusterProfiler", unload = TRUE)

# install.packages("DOSE")
# devtools::install_github('GuangchuangYu/clusterProfiler')
# R.utils::setOption( "clusterProfiler.download.method",'curl')

KEGGresult_RPPA <- enrichKEGG(gene = a$ENTREZID,organism = 'human',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
## 画气泡图
dotplot(KEGGresult_RPPA,font.size=8,showCategory=10) ## 最多显示10个通路
## 在网页中显示出通路的信息
browseKEGG(KEGGresult_RPPA,'hsa01521')
#结果保存
KEGGresult_RPPA_df <- KEGGresult_RPPA@result #148 
KEGGresult_RPPA_df$Description
##C.绘制KEGG富集结果 挑选结果，重新绘制图片

#蛋白富集意义？

#7. 转录因子活性分析----
##7.1转录活性分析准备及demo----
library(dorothea)## We load the required packages
library(bcellViper)
library(dplyr)
library(viper)
# accessing expression data from bcellViper
data(bcellViper, package = "bcellViper")

# acessing (human) dorothea regulons
# for mouse regulons: data(dorothea_mm, package = "dorothea")
data(dorothea_hs, package = "dorothea")

regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))

tf_activities <- run_viper(dset, regulons, 
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))
tf_activities@assayData[["exprs"]]
pheatmap::pheatmap(tf_activities@assayData[["exprs"]][1:10,],
                   cluster_cols = T,cluster_rows = T)
##7.2 BC计算Dorothea----
tf_activities_bc <- run_viper(exp1_matrix, regulons, 
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))


tf_activities_bc <- as.data.frame(tf_activities_bc)#118-1001
#保存数据write.csv(tf_activities_bc,"./output/tf_activities_bc.csv")

##7.3 可视化----
#A.热图
library(pheatmap)
my_color = mycolor<-colorRampPalette(c("royalblue","white","red"))(100)
pheatmap(tf_activities_bc,  #[1:20,]/interest_tf
         display_numbers = F, number_color="darkblue",number_format = "%.2f",#"%.1e",#显示两位小数
         fontsize=9, 
         fontsize_row = 10, 
         color=my_color, #breaks = my_breaks, 
         main = "DoRothEA", angle_col = 0,
         #gaps_col=c(3), 
         cluster_col=F,cluster_row = TRUE, #gaps_col应该禁止聚类
         treeheight_col = 0,  border_color ="white") 
#B.分组热图:  group_df FROM 7.4
library(pheatmap)

annotation_col_TF <- data.frame( Riskscore= group_df$Group,STAT4=group_df$STAT4_group)#,PAM50=RPPA_data_clinical$PAM50) #
row.names(annotation_col_TF) <- rownames(group_df)
annotation_colors_RPPA <- list(Riskscore=c(High="#D5B32B", Low="#2BD54D"),STAT4=c(High="pink", Low="skyblue"))#调色板 https://www.bejson.com/ui/getcolor/

pheat_TF <-pheatmap(t(group_df[,120:237]),show_rownames = T,show_colnames = F,cluster_cols = T,cluster_rows = T,
                  annotation_col  = annotation_col_TF,
                  #cutree_col = 2,
                  annotation_colors= annotation_colors_RPPA,fontsize = 8,fontsize_col=1,
                  color = colorRampPalette(c("navy", "white", "firebrick3") )(10) ) # c("blue", "white", "red"))(10)
pheat_TF
save_pheatmap_pdf(pheat_TF, "./output/pheatmap_TF_DOROTHEA_1.pdf",16,10)
#热图顺序排
#排列表达矩阵顺序-列完全按照annotation_col分布

ID_high_low_TF <- c(rownames(subset(annotation_col_TF,annotation_col_TF$STAT4 =="Low")),#501
                 rownames(subset(annotation_col_TF,annotation_col_TF$STAT4 =="High")) )
TF_matrix <-as.data.frame( t(group_df[,120:237]))
TF_matrix <-TF_matrix[,ID_high_low_TF] #118
library(pheatmap)
pheat_TF <- pheatmap(TF_matrix,show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = T,
                    annotation_col  = annotation_col_TF,
                    gaps_col=502 ,
                    annotation_colors= annotation_colors_RPPA,fontsize = 8,fontsize_col=1,
                    color = colorRampPalette(c("navy", "white", "firebrick3") )(10) ) # c("blue", "white", "red"))(10)
pheat_TF

#标基因在旁边-ComplexHeatmap
genes_TF <- as.data.frame(c(outTab_TF_sig$id,"STAT4.x")) #rownames(TF_matrix)
colnames(genes_TF)<-"genes"
library(ComplexHeatmap)#BiocManager::install("ComplexHeatmap")
samples <- rep(c('Low', 'High'), c(502, 499)) #定义样本分组信息 

B <- Heatmap(as.matrix(TF_matrix),#表达矩阵
        col = colorRampPalette(c("navy","white","firebrick3"))(100),#颜色定义
        show_row_names = T,show_column_names = F,#不展示行名
        cluster_columns  = F,column_split = samples,#份开
        row_names_gp = gpar(fontsize = 6),#行名称大小
        #row_dend_width = unit(1, "cm"),#聚类树宽度
        #column_dend_height = unit(10, "mm"),
        row_dend_gp = gpar(lwd=0.1), #行树宽度
        heatmap_legend_param = list(title = "",title_position = "topleft"),#matrix legend
        top_annotation = HeatmapAnnotation(Group = samples, name = "Exp",
                                           simple_anno_size = unit(2, 'mm'), 
                                           #annotation_legend_param = list(foo = list(title = "foo_top_anno")),
                                           col = list(Group = c('Low' = '#00DAE0', 'High' = '#FF9289') ),
                                           show_annotation_name = F) )#分组注释
B

 pdf("./output/pheatmap_TF_DOROTHEA_6.pdf",family="serif",width = 5,height = 5)
 B
 dev.off()
        
B+ rowAnnotation(link = anno_mark(at = which(rownames(TF_matrix) %in% genes_TF$genes), 
                                labels = genes_TF$genes, labels_gp = gpar(fontsize = 10) ))
#save_pheatmap_pdf(pheat_TF, "./output/pheatmap_TF_DOROTHEA_.pdf",16,10)
write.csv(TF_matrix,"./220719/Figure2F_TF.csv")

#C.
#C.boxplot-STAT4
#所有样本矩阵分析-正常与肿瘤-活性差异：及STAT4分组
#GTEX匹配为对照？
 

##7.4 dorothea预后TF浅探----
##A.合并临床数据
tf_activities_bc_clinical <- as.data.frame(t(tf_activities_bc))
tf_activities_bc_clinical$sample <- rownames(tf_activities_bc_clinical)#118
tf_activities_bc_clinical <- merge(tf_activities_bc_clinical,clinical_1,by="sample")

##KM
library(survminer)
library(survival) 

tf_activities_bc_clinical$STAT4_TF_GROUP <- ifelse(tf_activities_bc_clinical$STAT4.x >= median(tf_activities_bc_clinical$STAT4.x ),"High","Low")

fit<-survfit(Surv(OS.time.y/365, OS.x) ~ tf_activities_bc_clinical$STAT4_TF_GROUP, data =tf_activities_bc_clinical)
print(fit)#查看median风险50%生产率对应的时间！
ggsurvplot(fit,#palette=c("#FC8D62","#66C2A5"),
           data=tf_activities_bc_clinical,
           risk.table = T,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = T,#TRUE
           conf.int = F, #置信区间
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste0("Riskscore",' Survival Curve'),
           legend.title = "Group"#,
           #legend.labs = c("High", "Low")
)

#B.循环分组！！
group_df<-list()#data.frame()
for (i in  2:119){
  gene <- colnames(tf_activities_bc_clinical)[i]
  group_name <- paste0(gene,"_TF_group")
  group=ifelse(tf_activities_bc_clinical[,gene] >= median(tf_activities_bc_clinical[,gene] ),"High","Low") 
  group_df[[i-1]] <- group
  names(group_df)[i-1] <-group_name
  group_df <-cbind(data.frame(group_name<-group_df) )
  #colnames(group_df[[i-1]]) <-group_name
  #group_df[group_name] <-group
  #tf_activities_bc_1 <-list(group_name=ifelse(tf_activities_bc_clinical[,i] >= median(tf_activities_bc_clinical[,i] ),"High","Low") ) 
  #tf_activities_bc_1= cbind(group_df,tf_activities_bc_1) 
  #print(group)
}
rownames(group_df) <-tf_activities_bc_clinical$sample


#添加临床信息
group_df$sample <- rownames(group_df)
group_df<- merge(group_df,tf_activities_bc_clinical,by="sample")
rownames(group_df) <-group_df$sample

#C.单因素OS分析后KM分析
##批量OS 

pFilter=0.05                        #定义单因素显著性
library(survival)


sigGenes=c("OS.time.y","OS.x","sample")#注意修改此处及和下面,需要取出的列
outTab_TF=data.frame()
for(i in colnames(group_df[,120:237]) ){   #注意此处为参与COX的名称
  cox <- coxph(Surv(OS.time.y, OS.x) ~ as.numeric(group_df[,i]), data = group_df)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab_TF=rbind(outTab_TF,
                    cbind(id=i,
                          HR=coxSummary$conf.int[,"exp(coef)"],
                          HR.95L=coxSummary$conf.int[,"lower .95"],
                          HR.95H=coxSummary$conf.int[,"upper .95"],
                          pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
  }
}

View(outTab_TF)

#write.table(outTab_TF,file="./RPPA/uniCox_216_outTab_TF.xls",sep="\t",row.names=F,quote=F)
outTab_TF_sig <- subset(outTab_TF,outTab_TF$pvalue <0.05)#15个显著性TF  #无STAT4！！！ outTab_TF_sig$id

##OS显著KM绘图  :#略交叉NFATC2  REL 不交叉NFKB2 [4]
for (i in outTab_TF_sig$id[4] ){
  group_name <- paste0(i,"_TF_group")
  fit<-survfit(Surv(OS.time.y/365, OS.x) ~ group_df[[group_name]], data =group_df)
  #print(fit)#查看median风险50%生产率对应的时间！
  p<-ggsurvplot(fit,#palette=c("#FC8D62","#66C2A5"),
             data=group_df,
             risk.table = T,
             risk.table.col = 'strata',
             risk.table.y.text = T,
             risk.table.title = "",#group_name,#paste0('Number at risk',i),
             pval = T,#TRUE
             conf.int = T, #置信区间
             xlab = 'Time (years)',
             ggtheme = theme_light(),
             font.tickslab = c(12), #Y轴数字大小
             size=0.5,censor.size=3,censor.shape=NA,#NA不要+号
             surv.median.line = 'none',#不要中位值竖线
             title=i,#paste0("Riskscore",' Survival Curve'),
             legend.title = "Group",
             legend.labs = c("High", "Low") ) 
  

}
p #+guides(color = FALSE)
#ggsave("./output/dorothea_TF_OS_NFKB2.pdf",family="serif")


pdf("./output/dorothea_TF_OS_NFKB2_0.pdf",family="serif",width = 4,height = 5)
p
dev.off()

p$plot +
  theme(text = element_text(family = "serif",size = 16),title= element_text(family = "serif",size = 18),# )+ #,
        legend.title = element_text(size = 14),legend.text = element_text(),legend.key = element_blank(),legend.position = c(0.8,0.85))+#改字体,face = "bold",c(0.8,0.8) legend.title = element_blank(),
  scale_x_continuous(expand = c(0.05,0)) +scale_y_continuous(expand = c(0.02,0)) +
  guides(color = FALSE)+
  #scale_fill_brewer(palette = c("#FC8D62","#66C2A5"), name = "Riskscore", labels = c("High", "Low") )#对fill重新填充,再之前基础上叠加
  #scale_fill_discrete(name = "Riskscore", labels = c("High", "Low ")) #new legend  颜色无法改变
  scale_color_manual( values  = c("#FC8D62","#66C2A5"),#name = "Riskscore", 
                     labels = c("High", "Low") ) #改变颜色legend。改线条颜色scale_color_manual() fill

##D.单因素COX后森林图

#准备森林图数据

forest_data <- outTab_TF_sig
colnames(forest_data)<- c("id","HR","Low_95%CI","High_95%CI","Pvalue")#,"Coeffcient")
forest_data$Pvalue <- round(as.numeric(forest_data$Pvalue),4)
#forest_data$Coeffcient <- round(as.numeric(forest_data$Coeffcient),4)
forest_data$Pvalue <- ifelse(as.numeric(forest_data$Pvalue)>0.05,"ns",ifelse(as.numeric(forest_data$Pvalue)>0.01,"*",ifelse(as.numeric(forest_data$Pvalue)>0.001,"**","***")))
forest_data$HR <-round(as.numeric(forest_data$HR),3)
#forest_data$id[14] <- "Riskscore"
forest_table <- cbind(c("Signature", forest_data$id),
                      c("HR", forest_data$HR),
                      c("Pvalue", forest_data$Pvalue) )#,
#c("Coeffcient", forest_data$Coeffcient))##设置标签
forest_table
csize <- data.frame(mean=c(NA, as.numeric(forest_data$HR)),
                    lower=c(NA, as.numeric(forest_data$`Low_95%CI`)),
                    upper=c(NA, as.numeric(forest_data$`High_95%CI`)))##设置图形数据
csize
#https://www.jianshu.com/p/b460e3cd3bc5
library("forestplot")
#png(filename = "Forestplot-2.png",width=960, height=640)
pdf("./output/Forestplot_TF.pdf",width=6, height=4)
#write.csv(forest_table,"./220719/Figure2G_TF_forest.csv")
forestplot(labeltext = forest_table, 
           csize,
           graph.pos = 3,
           txt_gp=fpTxtGp(label=gpar(fontfamily = "serif",cex=1.25,fontsize=8),
                          ticks=gpar(cex=1.2),
                          xlab=gpar(cex = 0.2),
                          title=gpar(cex = 1.2,fontface=2)), # 修改label，ticks 字体和大小
           xticks=c(0.5,0.75,1,1.25,1.5),#xticks.digits=1,
           col=fpColors(box="deeppink", lines="skyblue", zero = "gray50"),##线条颜色设置 ##1c61b6
           zero=1, cex=0.9, lineheight = "auto", boxsize=0.2, colgap=unit(3,"mm"),
           lwd.ci=1, ci.vertices=TRUE, ci.vertices.height = 0.15,  #ci置信区间线条
           #is.summary = c(TRUE,rep(FALSE,14)), #加粗行
           xlog = F)#坐标轴log压缩
dev.off() 

##D.1 重新绘制森林图geom_point+geom_segment----
View(forest_table)#
View(forest_data)
#colnames(forest_data)[3:4] <-c("Low_95CI","High_95CI")
#forest_data$yend <-c(1:15)
#write.csv(forest_data,"./220719/Figure2D_forest_plot.csv")

library(ggplot2)
ggplot(data = forest_data#,aes(x = HR, y = id)
       ) +
  geom_segment( aes(x =as.numeric(Low_95CI), xend = as.numeric(High_95CI), y = id, yend = id),color = "black",size=0.1 )+ # 主要横线
  geom_segment( aes(x =as.numeric(Low_95CI), xend = as.numeric(Low_95CI), y = yend, yend = yend+0.2),color = "black" ,size=0.1)+
  geom_segment( aes(x =as.numeric(Low_95CI), xend = as.numeric(Low_95CI), y = yend, yend = yend-0.2),color = "black" ,size=0.1)+
  geom_segment( aes(x =as.numeric(High_95CI), xend = as.numeric(High_95CI), y = yend, yend = yend+0.2),color = "black",size=0.1 )+
  geom_segment( aes(x =as.numeric(High_95CI), xend = as.numeric(High_95CI), y = yend, yend = yend-0.2),color = "black" ,size=0.1)+
  geom_segment( aes(x=1,xend=1,y=0.5,yend=15.5),color="grey",linetype='dashed',size=0.1)+ #HR=1
  
  geom_point(size=4,aes(x = HR, y = id, group = id, color = Pvalue),pch=15)+
  scale_color_manual(values = c("#93cc82","#88c4e8") )+   #"#CC79A7", "#56B4E9" 天蓝，浅紫   #"#238b45","#2171b5"深蓝，绿
  labs(x = 'Hazard Ratio', y = '',color="Sign.")+theme_bw()+
  theme(text = element_text(size=16,family = "serif"),legend.position = c(0.8,0.8),panel.grid = element_blank() )

#ggsave("./output/Forestplot_TF_3.pdf",width=4, height=3.5) 
#D.2 COX TF boxplot
outTab_TF_sig$id
tf_boxplot <- tf_activities_bc_clinical[,c(outTab_TF_sig$id,"STAT4_group")]
tf_boxplot <-melt(tf_boxplot)
colnames(tf_boxplot )=c("group_stat4","TF","Score")
#write.csv(tf_boxplot,"./220719/Figure2H_tf_boxplot.csv")

library(dplyr)
plot_order_2 = tf_boxplot [gsva_New $group_stat4=="High",] %>% 
  group_by(TF) %>% 
  summarise(m = median(Score)) %>% #
  arrange(desc(m)) %>% 
  pull(TF)
tf_boxplot$TF = factor(tf_boxplot$TF ,levels = plot_order_2)

if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 60, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
library(ggplot2);library(ggpubr)
box_TME <- ggplot(tf_boxplot, aes(x = TF, y = Score))+ ##Celltype/T cells follicular /helper Neutrophils
  labs(y="TF Score",x= NULL)+  
  geom_boxplot(aes(fill = group_stat4),outlier.size = 0.02,size=0.02,outlier.alpha = 0)+ 
  #geom_bar(stat = "identity")+
  scale_fill_manual(values = c("skyblue","pink") )+ #"#1CB4B8", "#EB7369" "skyblue","pink"
  guides(fill = guide_legend(title = 'STAT4'))+
  theme_bw()+ mytheme + 
  #ylim(-0.4,0.65)+
  stat_compare_means(aes(group =  group_stat4),
                     label = "p.signif",
                     #label.y = 0.6,
                     size=3,#0.6
                     method = "wilcox.test",
                     hide.ns = T)
  #coord_flip()
box_TME
ggsave("./output/TF_boxplot_15_1.pdf",family="serif",width=8,height = 4)
##D.3 DE TF KEGG?


##E.multi-cox构建 DE_TF_riskscore----
library(survival)
multiCox=coxph(Surv(OS.time.y, OS.x) ~ CTCF+ESR2+NFATC2+NFKB2+PPARA+RARA+REL+RELB+RFX5+RUNX1+STAT3+STAT5A+STAT5B+STAT6+TCF7L2, data = tf_activities_bc_clinical)
multiCox=step(multiCox, direction="both")
multiCoxSum=summary(multiCox) #7 signature
#输出模型参数
outTab1=data.frame()
outTab1=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab1=cbind(id=row.names(outTab1),outTab1)
outTab1=gsub("`","",outTab1)
#write.table(outTab1,file="./outcomes/multiCox_forest.txt",sep="\t",row.names=F,quote=F)
#riskscore_TF
multiCox$coefficients

riskScore=predict(multiCox, type="risk", newdata=tf_activities_bc_clinical)#计算得分-
tf_activities_bc_clinical$riskscore_tf <- riskScore #合并
tf_activities_bc_clinical$riskscore_tf_1 <- tf_activities_bc_clinical$CTCF*(0.1564397)+tf_activities_bc_clinical$ESR2*(-0.2089292)+tf_activities_bc_clinical$PPARA*(-0.2462478)+
  tf_activities_bc_clinical$REL*(-0.2091362)+tf_activities_bc_clinical$STAT3*(0.3006313)+tf_activities_bc_clinical$STAT5A*(-0.2253238)+tf_activities_bc_clinical$TCF7L2*(-0.1595570) #不相同
##KM
tf_activities_bc_clinical$riskscore_tf_group <- ifelse(tf_activities_bc_clinical$riskscore_tf >= median(tf_activities_bc_clinical$riskscore_tf ),"High","Low")#中部略交叉 riskscore1 与riskscore 结果相同
#为何predict()函数计算的Riskscore不等于基因的表达量与其系数的乘积的加权呢? https://www.136.la/jingpin/show-205406.html

fit<-survfit(Surv(OS.time.y/365, OS.x) ~ riskscore_tf_group, data =tf_activities_bc_clinical)
print(fit)#查看median风险50%生产率对应的时间！
library(survminer)
ggsurvplot(fit,palette=c("#FC8D62","#66C2A5"),
           data=tf_activities_bc_clinical, #
           risk.table = F,
           risk.table.col = 'strata',
           risk.table.y.text = F,
           risk.table.title = 'Number at risk',
           pval = T,#TRUE
           conf.int = F, #置信区间
           xlab = 'Time (years)',
           ggtheme = theme_light(),
           surv.median.line = 'hv',
           title=paste0("Riskscore",' Survival Curve'),
           legend.title = "Group"#,
           #legend.labs = c("High", "Low")
) #分得开，但是ROC太差-放弃riskscore_tf

##ROC
library(pROC)
roc_data <- roc(tf_activities_bc_clinical$OS.x,tf_activities_bc_clinical$riskscore_tf,ci=TRUE) 
roc_data[["ci"]]
roc_data[["auc"]] #NFKB2  DFI 0.56  OS 0.56  #ROC太差了-建立RISKSCORE模型？ riskscore_tf / _1 同 0.60 
library(ggplot2)
g <- ggroc(roc_data,legacy.axes = TRUE, alpha = 1, colour = "red",  size = 1)
g
g + ggtitle("ROC curve") +
  geom_ribbon(aes(x=1 - roc_data$specificities, ymin=0, ymax=roc_data$sensitivities),alpha=0.5,fill="skyblue") +#填充
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey50", linetype="dashed")+
  xlab("1-Specificity (FPR)") + ylab("Sensitivity (TPR)")+
  theme_bw()+
  annotate("text",x=0.6,y=0.2,label=paste(" AUC: ",round(roc_data[["auc"]],3) ,"\n ","95%CI :",round(roc_data[["ci"]],3)[1],"-",round(roc_data[["ci"]],3)[3] ),size=8,colour = "red",family="serif")+
  scale_x_continuous(expand = c(0,0.01)) +scale_y_continuous(expand = c(0,0.01))+ #X轴数据两边c(0,0.02),否则横坐标显示不全
  theme(text=element_text(size=16,family="serif")) +#字体TNR
  theme(plot.margin=unit(rep(1,5),'lines') )
##time-ROC :任然较差 1，3，5年分别为0.6 0.7 0.68
library(timeROC)
library(survival)
#写为函数？
with(tf_activities_bc_clinical,
     ROC <<- timeROC(T=OS.Time,#结局时间 
                     delta=OS.x,#生存结局 
                     marker=riskscore_tf ,#预测变量  #lasso_rs   
                     cause=1,#阳性结局赋值，比如死亡与否
                     weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
                     times=c(365,365*3,365*5),#时间点，选取5年(60个月)和8年生存率
                     ROC = TRUE,
                     iid = TRUE)
)

#优化曲线-ggplot2
time_ROC_df <- data.frame(
  TP_1year = ROC$TP[, 1],
  FP_1year = ROC$FP[, 1],
  TP_3year = ROC$TP[, 2],
  FP_3year = ROC$FP[, 2],
  TP_5year = ROC$TP[, 3],
  FP_5year = ROC$FP[, 3]
)
library(ggplot2)
ggplot(data = time_ROC_df)+ #禁止填充否则颜色改变！
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +#geom_ribbon(aes(x=time_ROC_df[,2],ymin=0, ymax=time_ROC_df[,1]),alpha=0.4,fill="#BC3C29FF")+
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#0072B5FF") +#geom_ribbon(aes(x=time_ROC_df[,4],ymin=0, ymax=time_ROC_df[,3]),alpha=0.4,fill="#0072B5FF")+
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +#geom_ribbon(aes(x=time_ROC_df[,6],ymin=0, ymax=time_ROC_df[,5]),alpha=0.4,fill="#E18727FF")+
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2)+
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 6,
           label = paste0("AUC at 1 years = ", sprintf("%.2f",ROC$AUC[[1]])), color = "#BC3C29FF")+ #95CI未标注   #"%.2f"保留2位置数据好看
  annotate("text",
           x = 0.75, y = 0.15, size = 6,
           label = paste0("AUC at 3 years = ", sprintf("%.2f", ROC$AUC[[2]])), color = "#0072B5FF") +
  annotate("text",
           x = 0.75, y = 0.05, size = 6,
           label = paste0("AUC at 5 years = ", round(ROC$AUC[[3]],2) ), color = "#E18727FF") + #sprintf("%.2f", ROC$AUC[[3]])  signif(ROC$AUC[[3]],2)
  labs(x = "1−Specificity (FPR)", y = "Sensitivity (TPR)",size=16,family="serif")+
  theme(
    axis.text = element_text(size = 15, color = "black",family="serif"),#坐标轴字体
    axis.title.x = element_text( size = 20, color = "black",family="serif"),
    axis.title.y = element_text( size = 20, color = "black",family="serif"),panel.grid = element_blank()#删除网格
  )+
  scale_x_continuous(expand = c(0,0.03)) +scale_y_continuous(expand = c(0,0.01))
#ggsave("./output/timeROC_lasso_rs_DGEs_135years_OS.pdf",family="serif",width = 5.4)

#补充RISKSCORE KM PAM50 220719----
library(survival)
library(survminer)
table(clinical_1$IHC_sub) #IHC_sub
fit<-survfit(Surv(OS.time.y/365, OS.x) ~ STAT4_group, data =subset(clinical_1,clinical_1$IHC_sub=="HER2"))##clinical_1$PAM50=="Normal" Basal 165  #group_df[[group_name]]
#print(fit)#查看median风险50%生产率对应的时间！
p<-ggsurvplot(fit,#palette=c("#FC8D62","#66C2A5"),
              data=subset(clinical_1,clinical_1$IHC_sub=="HER2"),#修改此次，否则P值一致 #"T1","T2","T3","T4"
              risk.table = F,
              risk.table.col = 'strata',
              risk.table.y.text = T,
              risk.table.title = "",#group_name,#paste0('Number at risk',i),
              pval = T,pval.size=6,pval.method = F,pval.method.size=6,
              conf.int = F, #置信区间
              xlab = "",#'Time (years)',
              ylab="",
              font.title=c(16),font.x=c(16),font.y=c(16),font.lenged=c(20),font.tickslab = c(16), #Y轴数字大小
              font.family="serif",
              size=1,censor.size=3,censor.shape=NA,#NA不要+号 SIZE= 0.5
              surv.median.line = 'none',#不要中位值竖线
              #title="T1 & T2",#gsub("_gsva_group","",i),#paste0("Riskscore",' Survival Curve'),
              legend.title = "STAT4",#PAM50
              legend.labs =c("High", "Low"),# c("T1","T2","T3","T4"),#c("HER2","Luminal","TNBC"),#c("Basal","Her2","LumA","LumB","Normal"),#c("High", "Low"),#, fontsize=20,font.family="serif"
              ggtheme = theme_bw()+theme(text = element_text(family="serif"),legend.title = element_text(size=14),legend.text = element_text(size=14),panel.grid = element_blank() )#theme_light(),
)
p
#STAT4_IHC_HER2_44


##补ROC  basal  ROC TIMEROC太差
pam50 <-clinical_1# subset(clinical_1,clinical_1$PAM50 =="Normal")#165 clinical_1$PAM50=="Basal"
library(pROC)
roc_data <- roc(pam50$OS.x,pam50$lasso.risk.score,ci=TRUE) #lasso.risk.score
roc_data[["ci"]]
roc_data[["auc"]] #lasso.risk.score 0.64   basal 0.64  TNBC 0.62  / STAT4LOG2   OS ROC 0.57  luminal 0.6
#time roc  135 0.53 0.67 0.69 差
library(timeROC)
library(survival)
#写为函数？
with(pam50,
     ROC <<- timeROC(T=OS.Time,#结局时间 
                     delta=OS.x,#生存结局 
                     marker=lasso.risk.score,#STAT4_log2,#,#预测变量  #lasso_rs   
                     cause=1,#阳性结局赋值，比如死亡与否
                     weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
                     times=c(365,365*3,365*5),#时间点，选取5年(60个月)和8年生存率
                     ROC = TRUE,
                     iid = TRUE)
)

#优化曲线-ggplot2
time_ROC_df <- data.frame(
  TP_1year = ROC$TP[, 1],
  FP_1year = ROC$FP[, 1],
  TP_3year = ROC$TP[, 2],
  FP_3year = ROC$FP[, 2],
  TP_5year = ROC$TP[, 3],
  FP_5year = ROC$FP[, 3]
)
library(ggplot2)
ggplot(data = time_ROC_df)+ #禁止填充否则颜色改变！
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +#geom_ribbon(aes(x=time_ROC_df[,2],ymin=0, ymax=time_ROC_df[,1]),alpha=0.4,fill="#BC3C29FF")+
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#0072B5FF") +#geom_ribbon(aes(x=time_ROC_df[,4],ymin=0, ymax=time_ROC_df[,3]),alpha=0.4,fill="#0072B5FF")+
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +#geom_ribbon(aes(x=time_ROC_df[,6],ymin=0, ymax=time_ROC_df[,5]),alpha=0.4,fill="#E18727FF")+
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2)+
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 6,
           label = paste0("AUC at 1 years = ", sprintf("%.2f",ROC$AUC[[1]])), color = "#BC3C29FF")+ #95CI未标注   #"%.2f"保留2位置数据好看
  annotate("text",
           x = 0.75, y = 0.15, size = 6,
           label = paste0("AUC at 3 years = ", sprintf("%.2f", ROC$AUC[[2]])), color = "#0072B5FF") +
  annotate("text",
           x = 0.75, y = 0.05, size = 6,
           label = paste0("AUC at 5 years = ", round(ROC$AUC[[3]],2) ), color = "#E18727FF") + #sprintf("%.2f", ROC$AUC[[3]])  signif(ROC$AUC[[3]],2)
  labs(x = "1−Specificity (FPR)", y = "Sensitivity (TPR)",size=16,family="serif")+
  theme(
    axis.text = element_text(size = 15, color = "black",family="serif"),#坐标轴字体
    axis.title.x = element_text( size = 20, color = "black",family="serif"),
    axis.title.y = element_text( size = 20, color = "black",family="serif"),panel.grid = element_blank()#删除网格
  )+
  scale_x_continuous(expand = c(0,0.03)) +scale_y_continuous(expand = c(0,0.01))
#ggsave("./output/timeROC_lasso_rs_basal_OS.pdf",family="serif",width = 5.4)
write.csv(time_ROC_df,"./220719/Figure5C_timeroc.csv")

##220719 重新做STAT4 clinical  ----
write.csv(clinical_1,"./220719/clinical_1.csv")


#age TNM pathologic stage PAM50 boxplot
median(clinical_1$age)#58  clinical_1$STAT4_log2
clinical_1$age_group <- ifelse(clinical_1$age >= 58,"> 58","< 58")#high 512  low 489
table(clinical_1$IHC_sub)
clinical_1_pam50 <-clinical_1 #PAM50 排序
clinical_1_pam50$PAM50 <- with(clinical_1_pam50,reorder(PAM50,STAT4_log2,median))

library(ggplot2);library(ggpubr)
ggplot(clinical_1$PAM50,#subset(clinical_1,clinical_1$N_stage %in% c("N0","N1","N2","N3")), #"N0","N1","N2","N3" "T1","T2","T3","T4" "Stage I","Stage II","Stage III","Stage IV"
       aes(N_stage,STAT4_log2,color=N_stage))+  #age_group
  geom_point(alpha=0.2,size=1.5,
             position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  geom_boxplot(alpha=1,width=0.45,fill=NA,
               position=position_dodge(width=0.8),
               size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
               outlier.stroke = 0.5)+
  # geom_violin(alpha=0.2,width=0.9,
  #             position=position_dodge(width=0.8),
  #             size=0.25)+
  #scale_color_manual(values = my36colors[4:7])+ #c("#CC79A7", "#56B4E9")
  #geom_jitter(width = 0.25)+
  theme_bw() +ylab("STAT4 expression" )+xlab("")+labs(color="")+
  #ylim(0,3.5)+#theme(legend.position="none") 
  theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14),legend.position ="top",  
        panel.grid = element_blank(),legend.key = element_blank() ) +
  #stat_compare_means(#aes(group = Group) ,
    # comparisons=#list(c("Luminal","TNBC") ),
    # # list(c("LumB", "LumA"),c("LumB","Normal"),c("LumB","Basal"),c("LumB","Her2"),
    # #                   c("LumA", "Normal"),c("LumA","Basal"),c("LumA","Her2") #,# c("Normal", "Basal"),c("Normal","Her2"),c("Basal","Her2")
    # #                   ),
    #   #list(c("Stage I","Stage II"),c("Stage I","Stage III"),c("Stage I","Stage IV"),c("Stage II","Stage III"),c("Stage II","Stage IV"),c("Stage III","Stage IV")),#list(c("M0","M1")),#
    #   #list(c("T1","T2"),c("T1","T3"),c("T1","T4"),c("T2","T3"),c("T2","T4"),c("T3","T4")),#list(c("> 58","< 58")),#"High","Low"
    #   #list(c("T1","T4"),c("T2","T4"),c("T3","T4")),
    #   label = "p.signif",#"p.format p.signif 9.4e-12
    # method = "wilcox",
    # show.legend= F,#删除图例中的"a"
    # label.x=1.5,bracket.size=0.1,#vjust=0.5,
    # #label.y = max(clinical_1_msi_tmb$tmb)-1,
    # hide.ns = F,size=4)+
  stat_compare_means(method = "kruskal.test",show.legend= F,label.x = 1)# ,label.y = 4 #非参数检验kruskal.test not anova

#ggsave("./220719/stat4_boxplot_Nstage_1.pdf",width = 4,height = 4,family="serif")

##补差异表达STAT4 normal tumor
STAT4_mRNA <- readxl::read_xlsx("D:\\研究生科研\\生物信息学\\STAT4\\Subgroup_STAT4\\1.Riskscore\\STAT4_mRNA_1221.xlsx")
STAT4_mRNA$sample <-substr(STAT4_mRNA$ID,1,16)
STAT4_mRNA$Type <-substr(STAT4_mRNA$ID,14,16)#1222
table(STAT4_mRNA$Type)#01A-1077  01B-24  01C-1  06A-7  11A-99  11B-14   tumor 1077+7=1084  normal 99
STAT4_mRNA <-subset(STAT4_mRNA,STAT4_mRNA$Type %in% c("01A","06A","11A"))#1183
STAT4_mRNA$class <-ifelse(STAT4_mRNA$Type %in% c("01A","06A"),"Tumor","Normal")
table(STAT4_mRNA$class)#Normal 99  Tumor 1084
plot(STAT4_mRNA$STAT4)#0-15
plot(log2(STAT4_mRNA$STAT4) )#-6 - 4
plot(log2(STAT4_mRNA$STAT4+1) ) #0-4
STAT4_mRNA$STAT4_log2 <- log2(STAT4_mRNA$STAT4+1)
STAT4_mRNA$patient <- substr(STAT4_mRNA$ID,1,12)#相同的数据进行分组配对

#boxplot
library(ggplot2);library(ggpubr)
ggplot(STAT4_mRNA,aes(class,STAT4_log2,color=class))+  #age_group
  geom_point(alpha=0.2,size=1.5,
             # position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
             #                               jitter.height = 0,
             #                               dodge.width = 0.8)
             )+
  geom_boxplot(alpha=1,width=0.45,fill=NA,
               position=position_dodge(width=0.8),
               size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
               outlier.stroke = 0.5)+
  # geom_violin(alpha=0.2,width=0.9,
  #             position=position_dodge(width=0.8),
  #             size=0.25)+
  #scale_color_manual(values = my36colors[4:7])+ #c("#CC79A7", "#56B4E9")
  geom_line(aes(group = patient), color = 'gray70', lwd = 0.1,linetype = 1) + #配对线 linetype = 2虚线
  theme_bw() +ylab("STAT4 expression" )+xlab("")+labs(color="")+
  #ylim(0,3.5)+#theme(legend.position="none") 
  theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14),legend.position ="top",  
        panel.grid = element_blank(),legend.key = element_blank() ) +
  stat_compare_means(#aes(group = Group) ,
    comparisons=list(c("Tumor","Normal") ),
    label = "p.signif",#"p.format p.signif 9.4e-12
    method = "wilcox",
    show.legend= F,#删除图例中的"a"
    label.x=1.5,bracket.size=0.1,vjust=0.5,
    #label.y = max(clinical_1_msi_tmb$tmb)-1,
    hide.ns = T,size=4)

ggsave("./220719/stat4_boxplot_tumorvsnormal_paired_0.pdf",width = 4,height = 4,family="serif")
#save write.csv(STAT4_mRNA,"./220719/Figure1A_STAT4_normalvstumor.csv")

##补STAT4  OS DSS PFI
plot(clinical_1$STAT4_log2) #0-3
table(clinical_1$IHC_sub)  clinical_1$PFI.time
fit<-survfit(Surv(PFI.time/365, PFI.x) ~ STAT4_group, data =subset(clinical_1,clinical_1$PFI.time <= 10*365) ) #clinical_1$PAM50=="Normal" Basal 165  #group_df[[group_name]]
#print(fit)#查看median风险50%生产率对应的时间！
survminer::ggsurvplot(fit,palette=c("#66C2A5","#FC8D62"),
              data=subset(clinical_1,clinical_1$PFI.time <= 10*365) ,
              risk.table = F,
              risk.table.col = 'strata',
              risk.table.y.text = T,
              risk.table.title = "",#group_name,#paste0('Number at risk',i),
              pval = T,pval.size=6,pval.method = F,pval.method.size=6,
              conf.int = T, #置信区间
              xlab = "",#'Time (years)',
              ylab="",#ylim=c(0.22,1),
              font.title=c(16),font.x=c(16),font.y=c(16),font.lenged=c(20),font.tickslab = c(16), #Y轴数字大小
              font.family="serif",
              size=0.5,censor.size=3,censor.shape=NA,#NA不要+号
              surv.median.line = 'none',#不要中位值竖线
              #title="HER2",#gsub("_gsva_group","",i),#paste0("Riskscore",' Survival Curve'),
              legend.title = "STAT4",
              legend.labs = c("High", "Low"),#, fontsize=20,font.family="serif"
              ggtheme = theme_bw()+theme(text = element_text(family="serif"),legend.title = element_text(size=14),legend.text = element_text(size=14),panel.grid = element_blank() )#theme_light(),
)

#保存数据
OS_10 <-subset(clinical_1,clinical_1$OS.Time <= 10*365)#960
write.csv(OS_10[,c("sample","Barcode","OS.Time","OS.x","STAT4_group")],"./220719/Figure1L_OS_10.csv")
DSS_10 <-subset(clinical_1,clinical_1$DSS.time <= 10*365)#960
write.csv(DSS_10[,c("sample","Barcode","DSS.time","DSS","STAT4_group")],"./220719/Figure1L_DSS_10.csv")
PFI_10 <-subset(clinical_1,clinical_1$PFI.time <= 10*365)#969
write.csv(PFI_10[,c("sample","Barcode","PFI.time","PFI.x","STAT4_group")],"./220719/Figure1N_PFI_10.csv")

#8 TMB MSI  MATH----
#https://mp.weixin.qq.com/s?__biz=MzA4NDAzODkzMA==&mid=2651277130&idx=1&sn=715f2054e4f89898f48c245d89fc7847&chksm=841ea337b3692a21d70abf4600839781f7706941c85e328db7cd3779a7c66dc3cab0f24d9b39&scene=21#wechat_redirect



#8.2 TMB MATH----
##A. 数据准备
library(maftools)
library(tidyverse)
mafFilePath = dir(path = "./TMB_MSI/snv_gdc_download_20220709_104103",#"./0000_all_maf",#注意文件夹
                  pattern = "masked.maf.gz$",full.names = T,recursive=T)#下载：https://cloud.tencent.com/developer/article/1983809
for (extracted_maf in mafFilePath) { #length(file_extracted_maf)
  file_appended <- read.delim(extracted_maf, header = T, sep = '\t', comment.char = '#',stringsAsFactors = F)
  first_file <- rbind(first_file,file_appended)
}  #较久
BRCA_maf <- first_file#89258-140
BRCA_maf <- read.maf(BRCA_maf)#232.7MB
# 将工作地址设置为之前创建的目录-改回来
#setwd("D:/\u7814\u7a76\u751f\u79d1\u7814/\u751f\u7269\u4fe1\u606f\u5b66/STAT4/Subgroup_STAT4/3.Enrichment")

##计算TMB
tmb <- tmb(maf = BRCA_maf,
           captureSize = 50,  #Default 50MB
           logScale = T)#987
head(tmb) #结果与cbioportal下载数据的临床数据中TMB不同; 计算值接近Breast Invasive Carcinoma (TCGA, Firehose Legacy);而不是Breast Invasive Carcinoma (TCGA, PanCancer Atlas),但有MSI

tcgaCompare(
  maf = BRCA_maf, cohortName = 'BRCA_maf', 
  logscale = TRUE, capture_size = 50
)#样本数不同,版本不同

#8.2.1 TMB  MATH
##A. cbioportal 下载 :TCGA-BRCA pan-cancer data_mutations.txt --> TMB MATH 
file_cnv<- read.delim("./TMB_MSI/CBIOPORTAL_brca_tcga_pan_can_atlas_2018/brca_tcga_pan_can_atlas_2018/data_mutations.txt", header = T, sep = '\t', comment.char = '#',stringsAsFactors = F)#130495-111
BRCA_maf_bioportal <- read.maf(file_cnv) #232.3 MB
##B. calculate TMB
#B.1 比较评价
tcgaCompare(
  maf = BRCA_maf_bioportal, cohortName = 'BRCA_maf', 
  logscale = TRUE, capture_size = 50
)
tmb <- tmb(maf = BRCA_maf_bioportal,
           captureSize = 50,  #Default 50MB
           logScale = T)#1009 median 0.78/MB
head(tmb)
tmb <- tmb[-1009,] # 删除异常值106.80   #1008

#plot(tmb$Tumor_Sample_Barcode,log2(tmb$total_perMB)) #-4-6
plot(tmb$Tumor_Sample_Barcode,log2(tmb$total_perMB+1)) #0-6 !

##C.计算MATH：https://mp.weixin.qq.com/s?__biz=MzA4NDAzODkzMA==&mid=2651277130&idx=1&sn=715f2054e4f89898f48c245d89fc7847&chksm=841ea337b3692a21d70abf4600839781f7706941c85e328db7cd3779a7c66dc3cab0f24d9b39&scene=21#wechat_redirect
#计算mutant-allele tumor heterogeneity
barcode <- unique(BRCA_maf_bioportal@data$Tumor_Sample_Barcode)
head(barcode)#1009
MATH <- data.frame()
for (i in barcode){
  out.math = inferHeterogeneity(maf = BRCA_maf_bioportal, tsb = i)
  Tumor_Sample_Barcode=unique(out.math$clusterData$Tumor_Sample_Barcode)
  m = unique(out.math$clusterData$MATH)
  out = data.frame(Tumor_Sample_Barcode, m)
  MATH = rbind(MATH, out)
}#安装mclust
head(MATH) #1007 
plot(MATH$Tumor_Sample_Barcode,log2(MATH$m))#3-6  #all >1 ,no log(x+1)

##D 合并TMB MATH ? 放弃，直接使用cbioportal TMB MSI,不讨论MATH



#8.3MSI,TMB数据直接下载cbioportal----
##A.直接从cbioportal 下载 TCGA-BRCA pan-cancer 临床数据中的MSI 
msi_tmb <- read.table("./TMB_MSI/CBIOPORTAL_brca_tcga_pan_can_atlas_2018/brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt",skip=4,sep="\t",header=T) #1084 -18
msi_tmb <-msi_tmb[,c(1,15,17)] #1084-3
colnames(msi_tmb)[1] <- "sample"
##B.匹配到临床数据中clinical_1
clinical_1_msi_tmb <- merge(clinical_1,msi_tmb,by="sample")#997-184


##查看是否需要log2
plot(msi_tmb$sample,log2(msi_tmb$TMB_NONSYNONYMOUS+1) )
plot(msi_tmb$sample,msi_tmb$TMB_NONSYNONYMOUS)

plot(clinical_1_msi_tmb$MSI_SENSOR_SCORE)#0-30,多数在0-5
plot(log2(clinical_1_msi_tmb$MSI_SENSOR_SCORE+1) ) #0-5
plot(clinical_1_msi_tmb$TMB_NONSYNONYMOUS)#0-150 多在 0-50
plot(log2(clinical_1_msi_tmb$TMB_NONSYNONYMOUS+1) ) #0-6
##C.TMB MSI log2(x+1) 处理
clinical_1_msi_tmb$tmb <- log2(clinical_1_msi_tmb$TMB_NONSYNONYMOUS+1)
clinical_1_msi_tmb$msi <- log2(clinical_1_msi_tmb$MSI_SENSOR_SCORE+1)
plot(clinical_1_msi_tmb$msi)

##D. 相关性及分组
#D.1 计算相关性; NA value error COR NA
cor(clinical_1_msi_tmb$msi,clinical_1_msi_tmb$STAT4_log2, use = "complete.obs")#NA -0.17  p 2.409962e-08
cor(clinical_1_msi_tmb$tmb,clinical_1_msi_tmb$STAT4_log2, use = "complete.obs")#NA  0.1064
cor(clinical_1_msi_tmb$MSI_SENSOR_SCORE,clinical_1_msi_tmb$STAT4, use = "complete.obs")#NA -0.06
cor(clinical_1_msi_tmb$TMB_NONSYNONYMOUS,clinical_1_msi_tmb$STAT4, use = "complete.obs")#NA , use = "complete.obs"0.05
#太差了，不显示COR
#D.2 cor 散点图
library(ggplot2)#引用包
library(ggpubr)
library(ggExtra)
gene1="CD274"#"IDO1"#"msi"#"CD274"#"tmb"             #第一个基因名字
gene2="STAT4_log2"              #第二个基因名字
y=as.numeric(clinical_1_msi_tmb[,gene1])
x=as.numeric(clinical_1_msi_tmb[,gene2])#STAT4
df1=as.data.frame(cbind(x,y))
head(df1)

#计算相关系数及P值：
corT=cor.test(x,y,method="pearson") # spearman
cor=corT$estimate
pValue=corT$p.value
cor;pValue #-0.45 4.43926e-52

## ggplot的传统做法：TMB  MSI 均log2(x+1)
p1=ggplot(df1, aes(x, y) ) +#ggtitle("BRCA")+ 
  xlab("STAT4 expression")+ylab( paste0(toupper(gene1)," Score"))+                        #signif(pValue, 3)  round(pValue,3)
  annotate("text", x=2.5, y=4.5, label=paste0("Cor = ",round(cor,3),"\n","pvalue = ", signif(pValue, 3)) , family='serif', fontface='italic',colour='red', size=6)+
  geom_point(shape = 21, colour = "#87CEFA", fill = "#87CEFA", size = 3, stroke = .5,alpha=0.5)+ geom_smooth(method="lm",formula = y ~ x,linetype=2,color="#6495ED",fill="#D3D3D3") + 
  theme_bw(base_size = 16,base_family = "serif")+
  theme(axis.text  = element_text(family = "serif",size = 16),panel.grid=element_blank())
  

#+stat_cor(method = 'pearson', aes(x =x, y =y))#spearman 字太小
p1  #"#FC8D62","#66C2A5"  color="#6495ED",fill="#D3D3D3"

p1+#theme(text = element_text(family = "serif",face = "bold",size = 16))+#改字体
  scale_x_continuous(expand = c(0.01,0.01)) +scale_y_continuous(expand = c(0,0.05))

ggMarginal(p1,size=10,type = "histogram", xparams = list(colour = "skyblue", size = 0.2,fill = "#66CC66"),yparams = list(colour = "skyblue", size = 0.2,fill = "#FFCC66"))
ggsave("./TMB_MSI/1.outcomes/cor_msi_stst4log.pdf",width = 4,height = 4)

#colnames(df1) <- c("STAT4","CD274")
#write.csv(df1,"Figure6J_cor_STAT4_CD274.csv")

#D.3 分组TMB MSI boxplot
library(ggplot2);library(ggpubr)
#gene1 <-"tmb"                           #  "IDO1","PDCD1"(PD1),"FOXP3","CTLA4"  CD274(PDL1)
ggplot(clinical_1_msi_tmb, aes( STAT4_group,log2(IDO1) ,group=STAT4_group,color=STAT4_group))+  #fill
  # geom_point(alpha=0.5,size=1.5,
  #            position=position_jitterdodge(jitter.width = 0.35,
  #                                          jitter.height = 0,
  #                                          dodge.width = 0.8))+
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.05,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.25)+
  scale_color_manual(values = c("#CC79A7", "#56B4E9") )+ #color
  
  theme_bw() +ylab(paste0(toupper(gene1)," ") )+xlab("")+labs(color="STAT4")+
  #theme(legend.position="none") 
  theme(axis.title  = element_text(size=16), axis.text = element_text(size=12),legend.position ="top",  
        panel.grid = element_blank(),legend.key = element_blank() ) +
  stat_compare_means(#aes(group = Group) ,
    comparisons=list(c("High","Low")),#my_comparisons,
    
    label = "p.signif",#"p.format p.signif 9.4e-12
    method = "wilcox",
    show.legend= F,#删除图例中的"a"
    label.x=1.5,bracket.size=0.1,vjust=0.5,
    #label.y = max(clinical_1_msi_tmb$tmb)-1,
    hide.ns = T,size=5)
#ggsave(paste0("./TMB_MSI/1.outcomes/boxplot_",gene1,"_stst4log.pdf"),family="serif",width=4,height = 4)
#趋势符合！ 低表达预后差，TMB低，MSI高，合理


#8.3.1 A mRNA  cor---
##a.提取CD274表达
cd274 <- exp1[rownames(exp1) %in% "CD274" ,] #exp1是否标准化log2(x+1)?
cd274 <- as.data.frame(t(cd274))
cd274$sample <- rownames(cd274)

IC <- exp1[rownames(exp1) %in% c("IDO1","PDCD1","FOXP3","CTLA4") ,] #      #PD1  PDL1不存在
IC <- as.data.frame(t(IC))
IC$sample <- rownames(IC)
##b. 合并及可视化
clinical_1_msi_tmb <- merge(clinical_1_msi_tmb,IC,by="sample")#997  191

clinical_1_msi_tmb$PDCD1
cor(clinical_1_msi_tmb$PDCD1,clinical_1_msi_tmb$STAT4_log2) #0.69

#8.3.1 B Protein  cor
#蛋白质组定量106个样,无PDL1/CD274  /RPPA无STAT4，PDL1

#8.3.2 提取多个免疫检查点，计算相关性？
IC_genes <- c("IDO1","PDCD1","FOXP3","CTLA4","CD274")
for (j in IC_genes) {
  
  corT=cor.test(clinical_1_msi_tmb[,"STAT4_log2"],clinical_1_msi_tmb[,j],method="pearson") # spearman
  cor=corT$estimate
  pValue=corT$p.value
  print(c(j,cor,pValue))  #全是正相关
}
##boxplot
IC_df <- exp1[c(IC_genes),]#1001-5-  clinical_1_msi_tmb
IC_df <-as.data.frame(t(IC_df)) #1001-5
IC_df$sample <- rownames(IC_df)
IC_df <- merge(IC_df,clinical_1[,c("sample","STAT4_group")],by="sample")#1001 -7 
rownames(IC_df) <- IC_df$sample
IC_df <-IC_df[,-1]

library(reshape2)
IC_df <-melt(IC_df)
colnames(IC_df)<-c("STAT4","Gene","Expression")  #取log2????
IC_df$Expression <-log2(IC_df$Expression)
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   #axis.text.x = element_text(angle = 0, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
ggplot(IC_df, aes(x = Gene, y = Expression) )+ ##Celltype/T cells follicular /helper Neutrophils
  labs(y="Expression",x= NULL)+  
  
  geom_point(aes(color = STAT4),alpha=0.3,size=1.5,
             position=position_jitterdodge(jitter.width = 0.3,
                                           jitter.height = 0,
                                           dodge.width = 0.8))+
  geom_boxplot(aes(fill = STAT4),alpha=0.5,outlier.size = 0.02,size=0.02,outlier.alpha = 0)+ 
 
  scale_color_manual(values = c("skyblue","pink") )+ #"#1CB4B8", "#EB7369" "skyblue","pink"
  scale_fill_manual(values = c("skyblue","pink") )+
  guides(fill = guide_legend(title = 'STAT4'))+
  theme_bw()+ mytheme + 
  #ylim(-0.4,0.65)+
  stat_compare_means(aes(group =  STAT4),
                     label = "p.signif",
                     #label.y = 0.6,size=3,#0.6
                     method = "wilcox.test",
                     hide.ns = T)

ggsave("./immune/IC5_1.pdf",family="serif",width = 6,height = 4)#免疫抑制因子，趋势相反：参考https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-020-02267-2
#write.csv(IC_df,"./220719/Figure3A_immune_checkpoint.csv")

#8.4 TIDE and submap immune predict----
##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8861524/

#TIDE注册不上  SUBMAP不会设置

##寻找BC PDL1 cohort --
#cor 未找到

#STAT4_group pCR?



