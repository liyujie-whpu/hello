getwd()#启动工作目录 #更改工作目录：setwd("D:/../../..")
#读取基因表达矩阵
library(DESeq2)#一般运行都放在最前面
cs<-read.csv("E:/R Language/cs.csv",row.names = 1)#文件夹名+文件名把数据读进来
head(cs)#查看数据集前六行
dim(cs)#查看数据集维度（几行几列）
cs_1<-cs[rowSums(cs) != 0,]#将所有表达量是0的基因都去掉
dim(cs_1)
meta<-read.csv("E:/R Language/meta.csv",stringsAsFactors = T)#读入样本的分组文件
meta
colnames(cs_1) == meta$id#判断列名是否匹配，TRUE则都匹配
#DESeq2差异表达分析的函数，是分组数据集的列名
dds <-DESeqDataSetFromMatrix(countData = cs_1,#countData用于矩阵输入：一个非负整数的矩阵
                             colData = meta,#colData用于矩阵输入：一个至少有一列的DataFrame或data.frame。 colData的行与countData的列相对应
                             design = ~sample)#有多个变量的设计“~+组+条件“
dds <- DESeq(dds)#差异表达分析数据存储到dds
res <- results(dds)#结果提取到res
head(res)#查看前六行的结果
class(res)
res_1<-data.frame(res)#res转化为数据框存储为新的变量res_1
#library(dplyr)增加新的列
#res_1 %>%
#列名group，等于条件判断 mutate(group = case_when(
#                        log2FoldChange >= 2 & padj <= 0.05 ~ "UP",
#                        log2FoldChange <= -2 & padj <= 0.05 ~"DOWN",
#                        TRUE ~ "NOT_CHANGE"
#)) -> res_2
#统计情况 table(res_2$group)
#输出结果到文件 write.csv(res_2,file="DESeq2/diff_expr_result.csv",
#quote = F,row.names = TRUE)

class(res_1)
head(res_1)
res_1 <- res_1[order(res_1$padj, res_1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

#log2FC≥1 & padj<0.01 标识 up，代表显著上调的基因
#log2FC≤-1 & padj<0.01 标识 down，代表显著下调的基因
#其余标识 none，代表非差异的基因
res_1[which(res_1$log2FoldChange >= 1 & res_1$padj < 0.01),'sig'] <- 'up'
res_1[which(res_1$log2FoldChange <= -1 & res_1$padj < 0.01),'sig'] <- 'down'
res_1[which(abs(res_1$log2FoldChange) <= 1 | res_1$padj >= 0.01),'sig'] <- 'none'

#输出选择的差异基因总表
res_1_select <- subset(res_1, sig %in% c('up', 'down'))
write.table(res_1_select, file = 'D1_D180.DESeq2.select.txt', sep = '\t', col.names = NA, quote = FALSE)

#根据 up 和 down 分开输出
res_1_up <- subset(res_1, sig == 'up')
res_1_down <- subset(res_1, sig == 'down')

write.table(res_1_up, file = 'D1_D180.DESeq2.up.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res_1_down, file = 'D1_D180.DESeq2.down.txt', sep = '\t', col.names = NA, quote = FALSE)
##ggplot2 差异火山图
library(ggplot2)

#默认情况下，横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 padj
p <- ggplot(data = res_1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 1) +  #绘制散点图
  scale_color_manual(values = c('red', 'gray', 'green'), limits = c('up', 'none', 'down')) +  #自定义点的颜色
  labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'D1 vs D180', color = '') +  #坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), #背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-12, 12) + ylim(0, 35)  #定义刻度边界

p

#差异分析（知网上那个）实验一下，接上面代码的匹配成功。
dds <-  DESeqDataSetFromMatrix(countData = cs_1,colData = meta,design = ~ sample)
dim(dds)
#过滤
dds <- dds[rowSums(cs_1(dds)) > 1,]
nrow(dds) 
## 差异比较
dep <- DESeq(dds)
res <- results(dep)
diff = res
diff <- na.omit(diff)  ## 去除缺失值NA
dim(diff)
write.csv(diff,"all_diff.csv")#diff就是差异分析的总分析结果，保存
#进一步筛选差异基因，使用Padj值和log2FC进行筛选。Padj是P值矫正之后的数值，一般
#选取小于等于0.05（显著差异）的基因；同时log2FC是基因表达量的差异倍数。
#例如log2FC为1，证明这个基因在两种不同处理中的表达量相差了一倍，
#通常以大于1或小于-1为标准，大于1的为上调表达，少于-1的为下调表达。
foldChange = 1
padj = 0.05
diffsig <- diff[(diff$pvalue < padj & abs(diff$log2FoldChange) > foldChange),]
dim(diffsig)
write.csv(diffsig, "All_diffsig.csv")#diffsig是筛选出具有差异表达的基因，保存
#①PCA主成分分析初图01
vsd=vst(dds, blind = F)     #vst()函数效果和rlog（）一样，且速度更快。
plotPCA(vsd, intgroup=c("sample"))
#美化02
plotPCA(vsd, intgroup=c("sample"),returnData = TRUE)#返回坐标数据
PCAreturn <- "
             PC1        PC2 group  sample   name
D1.1   -31.16520   4.958110    D1   D1   D1.1
D1.2   -28.65333  -1.118512    D1   D1   D1.2
D1.3   -25.45342  -6.308166    D1   D1   D1.3
D180.1  22.45407  12.814777  D180 D180 D180.1
D180.2  31.35207   2.507937  D180 D180 D180.2
D180.3  31.46579 -12.854146  D180 D180 D180.3
"        
PCAreturn <-read.table(header = TRUE,text = PCAreturn)#将PCAreturn转换为表格形式
#header=T表示第一行的数据为列名
#将样品列改为因子型并排序
PCAreturn$sample=factor(PCAreturn$sample,levels = c("D1","D180"))
#绘制PCA图：1.aes()设置美学映射，横纵坐标为PC1,PC2，颜色按sample排序，画点大小为2；
#2.labs()设置横纵坐标题目；face="bold"设置加粗；
#3.theme_light(base_size=16)设置散点图风格，变成白底灰色网格；
ggplot(PCAreturn)+
  geom_point(aes(PC1,PC2,color=sample),size = 2)+
  labs(x="PC1:85%variance",y="PC2:7%variance",face="bold")+
  theme_light(base_size = 10)
#处理数据：1.读取基因计数表：
counts_fpkm<-read.csv("E:/R Language/counts_fpkm.csv",row.names = 1)
#for循环，依次提取cs_fpkm的1：6列的列名，
#将col_fpkm以此命名为：clm所代表的列名_fpkm
#将total=cs_fpkm文件中clm所代表的那一列的所有基因reads数的总和
#对样本进行fpkm标准化，并新建出样本_fpkm列
#for(clm in colnames(counts_fpkm)[1:6]){
#  col_fpkm=paste0(clm,"_FPKM")
#  total=sum(counts_fpkm[clm])
#  counts_fpkm[col_fpkm] = (counts_fpkm[clm] * 10^6)/(counts_fpkm$Length * as.numeric(total)
#/ 1000)
#}
#计算D1相对于D180上调或下调(这里我用已有的cs_fpkm表来接着做)
counts_fpkm$logFC_D1_D180_fpkm=
  log2((rowSums(counts_fpkm[c("D1_1_fpkm","D1_2_fpkm","D1_3_fpkm")])+1)/(rowSums(counts_fpkm[c("D180_1_fpkm","D180_2_fpkm","D180_3_fpkm")]+1)))
#提取dds中的DESeq后生成的“sample_D1_vs_D180"结果并生成D1_vs_D180数据：
D1_vs_D180=results(dds, name = "meta_D1_vs_D180")#没成功，不知道是不就是前面的res
#将res？转换成数据框格式;加后缀#
colnames(res)=paste0(colnames(res),"_D1_D180")
res=as.data.frame(res)
#检测缺失是否存在，存在NA则输出T
res$padj_D1_D180[is.na(res$padj_D1_D180)] = 1
counts_fpkm=merge(counts_fpkm,res,by="row.names")#按行名合并
rownames(counts_fpkm)=counts_fpkm$Row.names#以counts_fpkm的ROW.NAMES列作为列名
counts_fpkm=counts_fpkm[,-1]#将counts_fpkm的第一列删除
#筛选差异，从counts_fpkm中取满足条件的子集:log2Fo1dchange_D1_D180>2并且padj_D1_D180 <0.01
#subset(取子集，abs(取绝对值，#生成满足条件的res(就是组间对比)
res = subset(counts_fpkm,abs(counts_fpkm$log2FoldChange_D1_D180) >2 & counts_fpkm$padj_D1_D180 <0.01)
#pheatmap(绘制热图:%in%取前者在后者中出现的元素: unique()去除重复元素: union()取并集
#paste0不带分隔的字符串连接,用来命名热图中的列名
#cluster_rows = TRUE按行进行聚集或谱系聚类；
#cluster_cols =FALSE不按列进行聚集或谱系聚类；
#scale ="row"在行方向上居中并且缩放；
#show_rownames = FALSE，show_colnames = TRUE只显示列名不显示行名；
#cutree_rows = NA不将行划分为集群的数量；
#color = color按 color编辑热图中的颜色矢量；
#annotation_col=NA不按列显示热图的颜色注样；
#treeheight_row = 0行高度
library(pheatmap)
#初步热图
#设置color
color=colorRampPalette(c('blue4', "white","red"))(256)
htmap_D1_D180 = pheatmap(counts_fpkm[
  rownames(counts_fpkm) %in% unique(union(rownames(res),rownames(res))),
  paste0(rep(c("D1","D180"),each = 3),"_",rep(1:3,2),"_fpkm")],
  cluster_cols = FALSE,
  clustering_method = "single",
scale = "row",
show_rownames = FALSE,show_colnames = TRUE,
cutree_rows = NA,
color = color,
annotation_col = NA,
treeheight_row = 0)

#GO，筛选基因
D1_only = intersect(rownames(res[res$logFC_counts_fpkm <-2,]),rownames(res[res$logFC_counts_fpkm < -2,]))
GO_D1= enrichGO(D1_only,keyType = "SYMBOL",ont = "BP",OrgDb = org.Ss.eg.db)
dotplot(GO_D1,showCategory=20,orderBy = "X")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")  #下载包
BiocManager::install('og.Ss.eg.db')





