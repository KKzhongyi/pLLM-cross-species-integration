rm(list=ls())
library(Seurat)
library(ggplot2)
library(Matrix)
library(patchwork)
library(dplyr)
library(tidyr)
library(clue)
library(cluster)
library(mclust)
library(rhdf5)
library(lisi)
library(MLmetrics)  
library(foreach)
library(doParallel)
library(proxy)       
library(qs)  
library(ggsci)
options(future.globals.maxSize = 1e9)
options(timeout = 3000)
print("UMAPgeneration_Match.R starts")
#获取传递的参数
args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
celltypenumbaseline <- as.numeric(args[2])
resolution <- as.numeric(args[3])
file_dir <- args[4]
h5file <- args[5]
species_list <- strsplit(args[6], ",")[[1]]
species_list <- trimws(species_list)
cpu_number <- as.numeric(args[7])
options(qs.nthreads = 8)
integrationmethod <- args[8]



# output_dir <- "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/merge_new/Cell_Atlas_of_Aqueous_Humor/Results/Cell_Atlas_of_Aqueous_Humor_0.5_1_0.5_0.5_TRUE_KNN/"
# celltypenumbaseline=0
# resolution=0.04
# file_dir <- "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/merge_new/Cell_Atlas_of_Aqueous_Humor"
# h5file <- "Cell_Atlas_of_Aqueous_Humor_0.5_1_0.5_0.5_bycluster"
# bycluster <- TRUE
# species_list <- "crab_eating_macaque, rhesus_macaque,human, mice, pig"
# species_list <- strsplit(species_list, ",")[[1]]
# species_list <- trimws(species_list)

find_resolution <- function(
    obj, target_k,
    reduction = "integrated", dims = 1:30,
    lo = 0, hi = 1, max_iter = 25,
    tol = 0, seed = 1
) {
  set.seed(seed)
  k_of <- function(res) {
    length(unique(
      FindClusters(obj, resolution = res, verbose = FALSE)$seurat_clusters))
  }

  best <- NULL        # 只在 k == target_k 时赋值
  best_above <- NULL  # 记录第一个 k > target_k 的结果

  for (i in seq_len(max_iter)) {
    mid <- (lo + hi) / 2
    k   <- k_of(mid)
    cat(sprintf("iter %d — mid = %.6f, k = %d\n", i, mid, k))

    if (k == target_k) {
      best <- list(resolution = mid, k = k)
      break
    }
    if (k > target_k) {
      hi <- mid
      if (is.null(best_above) ||
          abs(k - target_k) <= abs(best_above$k - target_k)) {
        best_above <- list(resolution = mid, k = k)
      }

    } else {
      lo <- mid
    }
  }

  # 判断返回哪个结果
  if (!is.null(best)) {
    return(best)
  } else if (!is.null(best_above)) {
    return(best_above)
  } else {
    warning("No resolution found with k >= target_k.")
    return(list(resolution = NA, k = NA))
  }
}


obj_path <- paste0(file_dir,"/",h5file,".qs")
processed_obj_path <- file.path(file_dir, paste0(h5file,"processed_obj_", integrationmethod, ".qs"))

if(file.exists(processed_obj_path)){
  obj <- qs::qread(processed_obj_path)
}else{
if (file.exists(obj_path)) {
  # 文件存在，直接读取
  obj <- qs::qread(obj_path)
} else {
    # 打开 HDF5 文件
  file_dir <- paste0(file_dir,"/",h5file,".h5")

  # 读取数值型数据
  numeric_data <- h5read(file_dir, "data/numeric")

  # 转置数值数据
  numeric_data <- t(numeric_data)

  # 读取列名和行名
  numeric_columns <- as.character(h5read(file_dir, "meta/numeric_columns"))
  index <- as.character(h5read(file_dir, "meta/index"))

  # 检查数值数据维度是否匹配
  if (ncol(numeric_data) != length(numeric_columns)) {
    stop("数值数据的列数与 numeric_columns 不匹配！")
  }
  if (nrow(numeric_data) != length(index)) {
    stop("数值数据的行数与 index 不匹配！")
  }

  # 转换为数据框并设置列名和行名
  numeric_data <- as.data.frame(numeric_data)
  colnames(numeric_data) <- numeric_columns
  rownames(numeric_data) <- index

  # 读取字符串型数据
  string_data <- h5read(file_dir, "data/strings")

  # 转置字符串数据
  string_data <- t(string_data)

  # 读取字符串列名
  string_columns <- as.character(h5read(file_dir, "meta/string_columns"))

  # 检查字符串数据维度是否匹配
  if (ncol(string_data) != length(string_columns)) {
    stop("字符串数据的列数与 string_columns 不匹配！")
  }
  if (nrow(string_data) != length(index)) {
    stop("字符串数据的行数与 index 不匹配！")
  }

  # 转换为数据框并设置列名和行名
  string_data <- as.data.frame(string_data, stringsAsFactors = FALSE)
  colnames(string_data) <- string_columns
  rownames(string_data) <- index

  # 打印结果
  print(dim(numeric_data))  
  print(numeric_data[1:12, (ncol(numeric_data) - 4):ncol(numeric_data)])  # 查看最后几列  
  print(dim(string_data))
  # 查看string_data前12行
  head(string_data, 12)
  # ######
  # 创建 Seurat 对象
  # 将 numeric_data 转换为稀疏矩阵
  obj <- CreateSeuratObject(counts = t(numeric_data), meta.data = string_data)#注意这里的counts和meta.data都是来源于python，需要行列转置
  obj
  ## 清理内存
  rm(numeric_data, string_data,numeric_columns, string_columns, index)                                       
  gc()                                     
  # 保存 Seurat 对象为 qs 文件
  qs::qsave(obj, file = obj_path, preset = "balanced")
}
# 计算每个 cell_type 在每个 species 中的细胞数
cell_type_counts <- obj@meta.data %>%
  group_by(species, cell_type) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = species, values_from = count, values_fill = list(count = 0))

print(paste0("head(cell_type_counts)",head(cell_type_counts,5)))

# 只取 species 的细胞数都大于 celltypenumbaseline 的 cell_type
valid_cell_types <- cell_type_counts %>%
  filter(if_all(all_of(species_list), ~ .x >= celltypenumbaseline)) %>%
  pull(cell_type)
print(paste0("length of valid_cell_types:",length(valid_cell_types)))

#清理中间变量
rm(cell_type_counts)                       
gc()

# 构建新的 Seurat 对象 obj
obj <- subset(obj, cells = WhichCells(obj, expression = cell_type %in% valid_cell_types))
print(obj)
rm(valid_cell_types)                      
gc()
print(paste0("unique(obj@meta.data$cell_type)",head(unique(obj@meta.data$cell_type)),5))
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$species)
for (ly in names(obj@assays$RNA@layers)) {
  layer_mat <- obj@assays$RNA@layers[[ly]]
  keep      <- Matrix::rowSums(layer_mat) > 0
  obj       <- subset(obj, features = rownames(layer_mat)[keep])
}
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# 基于 PCA，整合跨物种数据
obj <- IntegrateLayers(
object = obj, 
method = integrationmethod,
orig.reduction = "pca", 
new.reduction = "integrated",
verbose = FALSE
)
obj <- FindNeighbors(obj, reduction = "integrated", dims = 1:30)
qs::qsave(obj, file = processed_obj_path, preset = "balanced")  # ← 保存为 .qs
}
#得到整合后数据的聚类结果，即用integraged.cca替代原始的pca降维
print(paste0("celltype数量:",length(unique(obj@meta.data$cell_type))))
res_out <- find_resolution(obj, target_k = length(unique(obj@meta.data$cell_type)),hi=resolution)
best_res <- res_out$resolution
write(best_res, file = paste0(output_dir,"/best_res.txt"))
message("chosen resolution = ", best_res, "; clusters = ", res_out$k)
obj <- FindClusters(obj, resolution = best_res, cluster.name = "Clusters")
obj <- RunUMAP(obj, reduction = "integrated", dims = 1:30, reduction.name = "umap.cca")


# —— 1）准备原始调色板 —— #
palette_funcs <- list(
  lancet         = function() pal_lancet("lanonc")(9),
  nejm           = function() pal_nejm("default")(8),
  npg            = function() pal_npg("nrc")(10),
  aaas           = function() pal_aaas("default")(8),
  jco            = function() pal_jco()(10),
  jama           = function() pal_jama("default")(7),
  bmj            = function() pal_bmj("default")(5),
  ucscgb         = function() pal_ucscgb("default")(8),
  d3             = function() pal_d3()(10),
  observable     = function() pal_observable()(10),
  locuszoom      = function() pal_locuszoom()(12),
  igv            = function() pal_igv()(8),
  cosmic         = function() pal_cosmic()(10),
  uchicago       = function() pal_uchicago()(10),
  startrek       = function() pal_startrek("uniform")(7),
  tron           = function() pal_tron()(6),
  flatui         = function() pal_flatui()(12),
  futurama       = function() pal_futurama()(9),
  rickandmorty   = function() pal_rickandmorty()(9),
  simpsons       = function() pal_simpsons()(9),
  frontiers      = function() pal_frontiers()(8)
)

# 2) 同时选择多款调色板给 cell_type/Clusters
ctcl_choice    <- c("lancet", "nejm", "npg")  
species_choice <- "startrek"  

# 3) 合并所选原色
base_ctcl_palette <- unlist(
  lapply(ctcl_choice, function(nm) palette_funcs[[nm]]()),
  use.names = FALSE
)

# 单一 species 调色板
base_species_palette <- palette_funcs[[species_choice]]()

# 5) 根据实际类别数切片并绘图
n_ct    <- length(unique(obj@meta.data$cell_type))
ct_cols <- base_ctcl_palette[1:n_ct]

n_sp    <- length(unique(obj@meta.data$species))
sp_cols <- base_species_palette[1:n_sp]


# —— 3）绘图 —— #
p_celltype <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by  = "cell_type",
  pt.size   = 1.2,
  cols      = ct_cols
) + ggtitle("Cell Type")

p_clusters <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by  = "Clusters",
  pt.size   = 1.2,
  cols      = ct_cols
) + ggtitle("Clusters")

p_species  <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by  = "species",
  pt.size   = 1.2,
  cols      = sp_cols
) + ggtitle("Species")


# —— 4）并排展示并保存 —— #
combined_plot <- wrap_plots(p_species, p_celltype, p_clusters, nrow = 1)
ggsave(
  filename = paste0(output_dir, "UMAP_", integrationmethod, ".pdf"),
  plot     = combined_plot,
  width    = 30,
  height   = 8
)
obj <- JoinLayers(obj)
obj


# 提取真实标签和嵌入
species <- obj@meta.data$species
cell_type <- obj@meta.data$cell_type
clusters <- obj@meta.data$Clusters
emb <- Embeddings(obj,reduction="integrated")[,1:30,drop=FALSE]
chunk_size <- 500
ncores <- cpu_number
registerDoParallel(cores=ncores)
cell_type_levels <- sort(unique(cell_type))
cluster_levels   <- sort(unique(clusters)) 
cell_type <- factor(cell_type, levels = cell_type_levels)
clusters  <- factor(clusters , levels = cluster_levels)
# 计算 NMI
nmi_cell_type <- cl_agreement(as.cl_partition(cell_type), as.cl_partition(clusters), method = "NMI")

# 计算 ARI
ari_cell_type <- adjustedRandIndex(cell_type, clusters)

#计算norm_cLISI
cLISI <- compute_lisi(emb,obj@meta.data,c("cell_type"))
print(head(cLISI))
K <- unique(cell_type) %>% length()
norm_cLISI <- mean(1 - (cLISI$cell_type - 1) / (K - 1), na.rm = TRUE)
print(paste0("norm_cLISI:",norm_cLISI))
#释放不需要的变量
rm(obj)
gc()
#计算F1指数
# 匈牙利一对一映射
tab   <- table(cell_type, clusters)   
assignment <- solve_LSAP(tab, maximum = TRUE)
cluster2ct <- rep(NA, ncol(tab))
names(cluster2ct) <- colnames(tab) 
for (i in seq_along(assignment)) {
  cl <- colnames(tab)[ assignment[i] ]
  ct <- rownames(tab)[ i ]
  cluster2ct[cl] <- ct
}
pred <- cluster2ct[ as.character(clusters) ]  
truth <- cell_type

# 生成混淆矩阵 
cm <- table(truth, pred)                                    
classes <- rownames(cm)
support <- rowSums(cm)                                      
n_total <- sum(support)

# 每类别 F1（供 Macro / Weighted）
f1_per_class <- sapply(classes, function(ct){
  MLmetrics::F1_Score(
      y_true = as.integer(truth == ct),
      y_pred = as.integer(pred  == ct),
      positive = 1)
})
#计算weighted_F1和macro_F1
valid_indices <- !is.na(f1_per_class) & !is.na(support)
weighted_F1 <- sum(f1_per_class[valid_indices] * support[valid_indices]) / n_total

macro_F1 <- mean(f1_per_class, na.rm = TRUE)
#计算micro_F1
TP_micro <- sum(diag(cm))
FP_micro <- sum(colSums(cm)) - TP_micro
FN_micro <- sum(rowSums(cm)) - TP_micro

precision_micro <- TP_micro / (TP_micro + FP_micro)
recall_micro    <- TP_micro / (TP_micro + FN_micro)
micro_F1        <- 2 * precision_micro * recall_micro /
                   (precision_micro + recall_micro)


print(paste("NMI (cell_type vs clusters):", nmi_cell_type))
print(paste("ARI (cell_type vs clusters):", ari_cell_type))
#print(paste("cell_type ASW:", cell_type_ASW))
print(paste("norm_cLISI:", norm_cLISI))
print(paste("Macro‑F1:", macro_F1))
print(paste("Micro‑F1:", micro_F1))
print(paste("Weighted‑F1:", weighted_F1))


output  <- c(
  paste0("NMI (cell_type vs clusters):", nmi_cell_type),
  paste0("ARI (cell_type vs clusters):", ari_cell_type),
  paste0("norm_cLISI:", norm_cLISI),
  paste0("Macro-F1:", macro_F1),
  paste0("Micro-F1:", micro_F1),
  paste0("Weighted-F1:", weighted_F1)
)


# 保存结果到文件
writeLines(output, con = paste0(output_dir, "cell_type_cluster_metrics_with_",integrationmethod,".txt"))
