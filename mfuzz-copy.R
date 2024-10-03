library(performance)
library(tidyverse)
library(Seurat)
library(readxl)
library(rstatix)
library(patchwork)
library(fs)
library(Mfuzz)
library(RColorBrewer)
library(Matrix)
#BiocManager::install("Mfuzz")

# 函数
file = "./Dropbox/中山/FL/DATA/merge/Progenitor_cell_merge.rds"

check_outlier <- function(df, .groupvar, .checkvar){
  .groupvar <- sym(.groupvar)
  .checkvar <- sym(.checkvar)
  df_outlier_check <- df %>%
    dplyr::group_by(!! .groupvar) %>%
    dplyr::summarise(cv = sd(!! .checkvar, na.rm = TRUE)/mean(!! .checkvar, na.rm = TRUE),)
  return(df_outlier_check)
}

check_outlier_mean <- function(df, .groupvar, .checkvar){
  .groupvar <- sym(.groupvar)
  .checkvar <- sym(.checkvar)
  df_outlier_check <- df %>%
    dplyr::group_by(!! .groupvar) %>%
    dplyr::summarise(cv = mean(!! .checkvar, na.rm = TRUE),)
  return(df_outlier_check)
}
  
  #mutate(sampleid = paste(Donor_id, Age, sep = "_"))
  # glue到一起的单细胞，不同condition/sample直接做成pseudobulk
  sr=readRDS(file)
  all_patient=as.data.frame(unique(cbind(sr@meta.data$patient,sr@meta.data$age)))
  rownames(all_patient)=all_patient[,1]
  colnames(all_patient)=c('patient','age_group')
  all_patient$age_group=factor(all_patient$age_group,levels =c("7-8w","9-11w",'12-14w','15-17w','P2w','P4w','P18w','P50w','1-2y','2-3y','>3y') )
  #sr=Progenitor_cell_merge
  Pseudobulk <- AggregateExpression(sr, group.by = c("patient"), return.seurat = FALSE, slot = "data")
  df <- Pseudobulk[[1]] %>% as.data.frame()
  df <- df[,unique(sr@meta.data$patient)]
  
  ##去离群值患者样本的2种方案：
  # 1、用所有测得基因中的高变基因，检查样本是否有离群值
  # 制作bulk seuratobject
  seurat_object <- CreateSeuratObject(counts = df, project = "PseudobulkProject")
  meta_df <- data.frame(age_group = all_patient[colnames(df),2], 
                        row.names = colnames(df))  # 确保行名与df中的列名一致
  seurat_object <- AddMetaData(seurat_object, metadata = meta_df)
  # 如果还没有做归一化和寻找高变基因
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object,npcs=22,seed.use=123)
  seurat_object <- RunTSNE(seurat_object, dims = 1:5,perplexity =5,seed.use = 1)
  #seurat_object <- RunTSNE(seurat_object, dims = 1:7,perplexity =7)
  #seurat_object <- RunUMAP(seurat_object, dims = 1:20,n_neighbors=0.1)
  
  # 提取 t-SNE 结果和 metadata 信息
  tsne_results <- Embeddings(seurat_object, reduction = "tsne")
  age_groups <- seurat_object@meta.data$age
  patient=rownames(seurat_object@meta.data)
  # 合并 t-SNE 结果和 age 信息到一个数据框中
  tsne_df <- data.frame(tsne_results, age_group = age_groups,patient=patient)
  morandi_colors <- c(
    "#E6A5A1",  # 鲜明的粉色
    "#D4B996",  # 明亮的米色
    "#B8D4B8",  # 显眼的浅绿色
    "#D2B8A3",  # 显眼的浅棕色
    "#9E9FBD",  # 鲜明的淡蓝色
    "#BDA7C3",  # 明亮的淡紫色
    "#8E8D8F",  # 经典的灰色
    "#A3C1E1",  # 鲜明的蓝色
    "#B0A8A0",  # 明亮的灰棕色
    "#D0A8A0",  # 温和的浅橙色
    "#A7A9B9"   # 柔和的灰蓝色
  )
  gradient_colors <- colorRampPalette(c("darkred",'yellow', "blue"))(11)
  
  # 使用 ggplot2 绘制散点图
  library(ggplot2)
  p=#ggplot(tsne_df, aes(x = PC_1, y = PC_2, color = age_group)) +
    ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, color = age_group)) +
    geom_point(size = 3) +
    geom_text(aes(label = patient), size = 3, vjust=-1,family = "Arial",show.legend = FALSE) +
    labs(title = "t-SNE Visualization by Age Group", x = "tSNE_1", y = "tSNE_2") +
    #labs(title = "PCA Visualization by Age Group", x = "PCA_1", y = "PCA_2") +
    scale_color_manual(values = gradient_colors)+theme_minimal()  # 可根据需要调整颜色
  plot(p)
  ggsave('./Dropbox/中山/FL/DATA/merge/年龄组上渐变色.png',bg='white')
  
  #2、检查基因是否有outlier：保留（50%的人表达+方差前25%）基因均值/变异系数,剔除同一年龄组内差异大的基因，减少聚类目标
  # 50%的样本都有表达的基因【去掉】
  ## df =df[apply(df,1,function(x) sum(x>0)>=0.5*ncol(df)),]
  #2（v2） 对每个年龄组进行筛选，并找到在所有年龄组中都表达大于 0 的基因
  genes_to_keep <- apply(df, 1, function(gene_expression) {
    all(sapply(split(gene_expression, all_patient$age), function(group) any(group > 0)))
  })
  df <- df[genes_to_keep, ]
  m.vars=apply(df,1,var)
  # 取方差前25%的基因
  df.upper <- df[which(m.vars > quantile(m.vars, probs = seq(0,1,0.25))[3]),]
  df_t <- t(df.upper)
  df_t = df_t %>% as.data.frame()
  df_t$group = all_patient[rownames(df_t),2]
  # 行为sample，列为基因（50%的人表达+方差前25%）,group为年龄组
  
  
  df_t$group <-factor(df_t$group,levels =c("7-8w","9-11w",'12-14w','15-17w','P2w','P4w','P18w','P50w','1-2y','2-3y','>3y') )
  
  dfuse <- check_outlier(df_t, "group", colnames(df_t)[1])
  colnames(dfuse)[2] =colnames(df_t)[1]
  for (q in 2:(ncol(df_t)-1)) {
    name <- colnames(df_t)[q]
    m <- check_outlier(df_t, "group", name)
    colnames(m)[2] <-name
    dfuse <- cbind(dfuse,m[,2])
    print(q)
  }
  
  write.csv(dfuse,"./Dropbox/中山/FL/DATA/merge/all_gene_cv_byage.csv")
  dfuse = read.csv("./Dropbox/中山/FL/DATA/merge/all_gene_cv_byage.csv")
  #dfuse = read.csv("./Dropbox/中山/FL/DATA/merge/all_gene_mean_byage.csv")
  # 把na换成0
  dfuse <-  dfuse%>% mutate_all(~replace(., is.na(.), 0))
  # 储存样本对应年龄组和基因
  df_group=dfuse$group
  df_genes=colnames(dfuse)[3:ncol(dfuse)]
  # 展示各个年龄组CV的分布
  long_df <- melt(dfuse[,2:ncol(dfuse)], id.vars = "group", measure.vars = df_genes)
  colnames(long_df)=c('Age_group','Gene','CV')
  long_df$Age_group=factor(long_df$Age_group,levels =c("7-8w","9-11w",'12-14w','15-17w','P2w','P4w','P18w','P50w','1-2y','2-3y','>3y') )
  p<-ggplot(long_df, aes(x = Age_group, y = CV))
  p=p+geom_boxplot()+geom_violin(aes(fill = Age_group))+ylim(0,1000)+theme_classic()
  plot(p)
  ggsave(p,'./Dropbox/中山/FL/DATA/merge/年龄组内Mean.png')
  # 年龄组和样本对应表格
  patient_N <- aggregate(. ~ age_group, data = as.data.frame(all_patient), FUN = function(x) paste(x, collapse = ","))
  write.csv(patient_N,'./Dropbox/中山/FL/DATA/merge/patient_agegroup.csv')
  # 探索去掉任何年龄组有CV大于1的基因
  gene_delet<- apply(dfuse[,2:ncol(dfuse)], 2, function(x) sum(x > 1) > 0)
  
##重新按照时间点和一个bulk， 只保留1\2步骤过滤的gene_use
  Pseudobulk <- AggregateExpression(sr, group.by = c("age"), return.seurat = FALSE, slot = "data")
  df <- Pseudobulk[[1]] %>% as.data.frame()
  df_delet <- df[gene_delet,]
  df_use=df[-gene_delet,]
  # 指定时间顺序
  df_use=df_use[,c("g7-8w","g9-11w",'g12-14w','g15-17w','P2w','P4w','P18w','P50w','g1-2y','g2-3y','g>3y')]
  
  #准备Mfuzz
  fpkm_raw <- df_use
  fpkm <- fpkm_raw
  fpkm <- as.matrix(fpkm)
  fpkm <- apply(fpkm,2,as.numeric)
  rownames(fpkm)=rownames(fpkm_raw)
  mfuzz_class <- new('ExpressionSet',exprs = fpkm)
  #mfuzz_class <- filter.NA(mfuzz_class, thres = 0.2)
  #mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
  #mfuzz_class <- standardise(mfuzz_class)

  set.seed(123)
  cluster_num <- 20
  mfrow <- c(4, 5)
  mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
  pdf("./Dropbox/中山/FL/DATA/merge/plot_test.pdf")
  mfuzz.plot(mfuzz_class,
             cl = mfuzz_cluster,
             mfrow = mfrow ,
             new.window = FALSE,
             time.labels = colnames(fpkm))
  dev.off()
  clu <- mfuzz_cluster$cluster %>%as.data.frame()
  write.csv(clu,"./Dropbox/中山/FL/DATA/merge/cluster_test.csv")
  #确定最佳聚类数
  library(vegan)
  vegan <- cascadeKM(fpkm, inf.gr = 2, sup.gr = cluster_num, iter = 1000)
  result = vegan$results %>% as.data.frame() %>% t()
  png("./Dropbox/中山/FL/DATA/merge/vegan_test.png", width = 800, height = 600, res = 100)
  plot(vegan, sortg=TRUE,grpmts.plot = TRUE)
  dev.off()

##检查adult gene表达 08.24
  ##重新按照时间点和一个bulk， 只保留adult_gene
  library(xlsx)
  # 读取xlsx文件
  adult_gene <- read.xlsx("/Users/liyaqi/Dropbox/中山/FL/DATA/mmc3.xlsx", sheetIndex = 2,header = FALSE)
  Pseudobulk <- AggregateExpression(sr, group.by = c("age"), return.seurat = FALSE, slot = "data",normalization.method = "LogNormalize")
  df <- Pseudobulk[[1]] %>% as.data.frame()
  adult_gene_remain=intersect(toupper(adult_gene$X1),toupper(rownames(df)))
  rownames(df)=toupper(rownames(df))
  df_use=df[adult_gene_remain,]
  # 指定时间顺序
  df_use=df_use[,c("g7-8w","g9-11w",'g12-14w','g15-17w','P2w','P4w','P18w','P50w','g1-2y','g2-3y','g>3y')]
  
  #准备Mfuzz
  fpkm_raw <- df_use
  fpkm <- fpkm_raw
  fpkm <- as.matrix(fpkm)
  #fpkm <- apply(fpkm,2,as.numeric)
  rownames(fpkm)=rownames(fpkm_raw)
  mfuzz_class <- new('ExpressionSet',exprs = fpkm)
  mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
  mfuzz_class <- standardise(mfuzz_class)
  mfuzz_class <- filter.NA(mfuzz_class, thres = 0.01)
  
  
  set.seed(123)
  cluster_num <- 15
  mfrow <- c(4, 5)
  mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
  pdf("./Dropbox/中山/FL/DATA/merge/plot_adult_gene_mfuzz.pdf")
  mfuzz.plot(mfuzz_class,
             cl = mfuzz_cluster,
             mfrow = mfrow ,
             new.window = FALSE,
             time.labels = colnames(fpkm))
  dev.off()
  clu <- mfuzz_cluster$cluster %>%as.data.frame()
  write.csv(clu,"./Dropbox/中山/FL/DATA/merge/cluster_adult_gene_mfuzz.csv")
  #确定最佳聚类数
  library(vegan)
  vegan <- cascadeKM(fpkm, inf.gr = 2, sup.gr = cluster_num, iter = 1000)
  result = vegan$results %>% as.data.frame() %>% t()
  png("./Dropbox/中山/FL/DATA/merge/vegan_adult_gene_mfuzz.png", width = 800, height = 600, res = 100)
  plot(vegan, sortg=TRUE,grpmts.plot = TRUE)
  dev.off()
  
##检查分布
  Idents(sr) <- "age"
  markers <- FindMarkers(object = sr, slot='data',
                         ident.1 = "15-17w", 
                         ident.2 = "P2w",
                         min.pct = 0.01,            # 设置最小表达百分比
                         logfc.threshold = 0.25,   # 设置log-fold变化阈值
                         test.use = "MAST", 
                         # 使用 MAST 统计测试方法
                         min.cells.group = 3)      # 设置最小细胞数量
  
  # 查看结果
  markers$avg_log2FC=as.numeric(markers$avg_log2FC)
  markers_sig=markers[markers$p_val_adj<0.05,]
  markers_sig
  write.csv(markers_sig,'./Dropbox/中山/FL/DATA/merge/All_HSC_HPC_markers.csv')
  
  marker_adult=markers[markers$avg_log2FC<0 & markers$p_val_adj<0.05,]
  
  marker_adult_intersect_gene=intersect(toupper(adult_gene$X1),toupper(rownames(marker_adult)))
  marker_adult_intersect=marker_adult[marker_adult_intersect_gene,]
  marker_adult_intersect=marker_adult_intersect[order(marker_adult_intersect$avg_log2FC),]
  marker_adult_intersect
  length(marker_adult_intersect[,1])
  
  genes_to_plot= rownames(markers_sig)[1:9]
  sr@meta.data$age=factor(sr@meta.data$age,levels=c("7-8w","9-11w",'12-14w','15-17w','P2w','P4w','P18w','P50w','1-2y','2-3y','>3y'))
  p=VlnPlot(object = sr, layer='data',
          features = genes_to_plot, 
          group.by = "age", 
          pt.size = 0.001)
 plot(p)  
 ggsave('./Dropbox/中山/FL/DATA/merge/adult_intersect_HSC_HPC_vln.png',height=25,width=40,units ='cm')
 ##
 table(sr@meta.data$age)
 write.csv(marker_adult_intersect,'./Dropbox/中山/FL/DATA/merge/adult_intersect_HSC_HPC_markers.csv')
 
 
 ##检查findmarker找到的 gene表达 08.24
 ##重新按照时间点和一个bulk， 只保留adult_gene
 library(xlsx)
 # 读取xlsx文件
 Pseudobulk <- AggregateExpression(sr, group.by = c("age"), return.seurat = FALSE, slot = "data",normalization.method = "LogNormalize")
 df <- Pseudobulk[[1]] %>% as.data.frame()
 adult_gene_remain=intersect(toupper(rownames(markers_sig)),toupper(rownames(df)))
 rownames(df)=toupper(rownames(df))
 df_use=df[adult_gene_remain,]
 # 指定时间顺序
 df_use=df_use[,c("g7-8w","g9-11w",'g12-14w','g15-17w','P2w','P4w','P18w','P50w','g1-2y','g2-3y','g>3y')]
 
 #准备Mfuzz
 fpkm_raw <- df_use
 fpkm <- fpkm_raw
 fpkm <- as.matrix(fpkm)
 #fpkm <- apply(fpkm,2,as.numeric)
 rownames(fpkm)=rownames(fpkm_raw)
 mfuzz_class <- new('ExpressionSet',exprs = fpkm)
 mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
 mfuzz_class <- standardise(mfuzz_class)
 mfuzz_class <- filter.NA(mfuzz_class, thres = 0.01)
 
 
 set.seed(123)
 cluster_num <- 15
 mfrow <- c(4, 5)
 mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
 pdf("./Dropbox/中山/FL/DATA/merge/plot_markers_mfuzz.pdf")
 mfuzz.plot(mfuzz_class,
            cl = mfuzz_cluster,
            mfrow = mfrow ,
            new.window = FALSE,
            time.labels = colnames(fpkm))
 dev.off()
 clu <- mfuzz_cluster$cluster %>%as.data.frame()
 write.csv(clu,"./Dropbox/中山/FL/DATA/merge/cluster_findmarkers_mfuzz.csv")
 #确定最佳聚类数
 library(vegan)
 vegan <- cascadeKM(fpkm, inf.gr = 2, sup.gr = cluster_num, iter = 1000)
 result = vegan$results %>% as.data.frame() %>% t()
 png("./Dropbox/中山/FL/DATA/merge/vegan_findmarkers_mfuzz.png", width = 800, height = 600, res = 100)
 plot(vegan, sortg=TRUE,grpmts.plot = TRUE)
 dev.off()
 
 ## 按年龄组检查 n_feature\n_count\cell count
 all_sc=readRDS('/Users/liyaqi/Dropbox/中山/FL/DATA/merge/merge_HBFL.rds')
 all_sc@meta.data$age=factor(all_sc@meta.data$age,levels=c("7-8w","9-11w",'12-14w','15-17w','P2w','P4w','P18w','P50w','1-2y','2-3y','>3y'))
 library(dplyr)
 library(dplyr)
 
 # 假设你的Seurat对象叫做`seurat_obj`
 meta_data <- all_sc@meta.data
 
 # 生成细胞类型和年龄的细胞数目表
 cell_count_table <- meta_data %>%
   group_by(Celltype_V2, age) %>%
   summarise(Cell_Count = n()) %>%
   spread(key = age, value = Cell_Count, fill = 0)
 
 # 计算每个年龄段的细胞总数
 total_counts <- colSums(cell_count_table[, -1])
 
 # 创建一个新的数据框用于存储总和行，并确保列名对齐
 total_row <- as.data.frame(t(total_counts))
 total_row <- cbind(Celltype_V2 = "Total", total_row)
 
 # 合并总和行到细胞数目表
 cell_count_table <- bind_rows(cell_count_table, total_row)
 
 # 检查结果
 write.csv(cell_count_table,'./Dropbox/中山/FL/DATA/merge/cell_count_table.csv')
 
 ###按照data slot统计,后检查feature count
 library(Seurat)
 library(ggplot2)
 library(dplyr)
 ###08.26待运行
 all_sc <- all_sc %>%
   AddMetaData(
     metadata = Matrix::colSums(all_sc@assays$RNA@counts),
     col.name = "n_Count"
   ) %>%
   AddMetaData(
     metadata = Matrix::colSums(all_sc@assays$RNA@counts > 0),
     col.name = "n_Feature"
   )
 
 meta_data <- all_sc@meta.data
 
 # 创建细胞类型列表
 cell_types <- unique(meta_data$Celltype_V2)
 
 # 初始化plotlist
 plotlist <- list()
 
 # 为每个细胞类型绘制nFeature和nCount的小提琴图并存储在plotlist中
 for (cell_type in cell_types) {
   # 筛选特定细胞类型的数据
   cell_data <- meta_data %>% filter(Celltype_V2 == cell_type)
   
   # 绘制nFeature的小提琴图
   p1 <- ggplot(cell_data, aes(x = age, y = nFeature_RNA)) +
     geom_violin(trim = FALSE) +
     labs(title = paste(cell_type, "- nFeature by Age"), x = "Age", y = "nFeature") +
     theme_minimal()
   
   # 绘制nCount的小提琴图
   p2 <- ggplot(cell_data, aes(x = age, y = nCount_RNA)) +
     geom_violin(trim = FALSE) +
     labs(title = paste(cell_type, "- nCount by Age"), x = "Age", y = "nCount") +
     theme_minimal()
   
   # 将图形存入plotlist
   plotlist[[paste(cell_type, "nFeature", sep = "_")]] <- p1
   plotlist[[paste(cell_type, "nCount", sep = "_")]] <- p2
 }
 
 # 输出结果
 plot(plotlist[[12]])
 
 
 # 文件分割传输
 
 # 设置每个文件最大行数或分割数目
 chunk_size <- 2000  # 每个切片的行数 (可以调整这个数字)
 
 # 获取细胞名列表
 cells <- Cells(all_sc)
 num_cells <- length(cells)
 num_chunks <- ceiling(num_cells / chunk_size)
 
 # 创建保存分割数据的文件夹
 dir.create("seurat_chunks_2")
 
 # 分割并保存
 for (i in 1:num_chunks) {
   start_idx <- ((i - 1) * chunk_size) + 1
   end_idx <- min(i * chunk_size, num_cells)
   
   # 获取子集
   cell_subset <- cells[start_idx:end_idx]
   seurat_subset <- subset(all_sc, cells = cell_subset)
   
   # 保存子集为 RDS 文件
   saveRDS(seurat_subset, file = paste0("seurat_chunks/seurat_chunk_", i, ".rds"))
 }
 
 
 ##检查E16.5类似时期和出生后差异基因
 Idents(sr) <- "age"
 markers <- FindMarkers(object = sr, slot='data',
                        ident.1 = "9-11w", 
                        ident.2 = "P2w",
                        min.pct = 0.01,            # 设置最小表达百分比
                        logfc.threshold = 0.25,   # 设置log-fold变化阈值
                        test.use = "MAST", 
                        # 使用 MAST 统计测试方法
                        min.cells.group = 3)      # 设置最小细胞数量
 
 # 查看结果
 markers$avg_log2FC=as.numeric(markers$avg_log2FC)
 markers_sig=markers[markers$p_val_adj<0.05,]
 markers_sig
 write.csv(markers_sig,'./Dropbox/中山/FL/DATA/merge/All_HSC_HPC_markers.csv')
 
 
## public
 
