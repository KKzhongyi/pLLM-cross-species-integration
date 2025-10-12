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
options(future.globals.maxSize = 1e9)
options(timeout = 3000)
print("UMAPgeneration_Match.R starts")
#获取传递的参数
args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
file_dir <- args[2]
h5file <- args[3]
cpu_number <- as.numeric(args[4])
options(qs.nthreads = 8)
integrationmethod <- args[5]

best_res <- as.numeric(readLines(file.path(output_dir, "best_res.txt"), n = 1))
processed_obj_path <- file.path(file_dir, paste0(h5file,"processed_obj_", integrationmethod, ".qs"))


obj <- qs::qread(processed_obj_path)
obj <- FindClusters(obj, resolution = best_res, cluster.name = "Clusters")
# 提取真实标签和嵌入
emb <- Embeddings(obj,reduction="integrated")[,1:30,drop=FALSE]
labels <- obj@meta.data$cell_type
rm(obj)
gc()
N <- length(labels)
chunk_size <- 2000
idx_chunks <- split(seq_len(N), ceiling(seq_len(N) / chunk_size))
ncores <- cpu_number
registerDoParallel(cores=ncores)

#计算cell_type ASW
s_sampled <- foreach(i = seq_along(idx_chunks), .combine ="c") %dopar% {
  chunk <- idx_chunks[[i]]
  cat(sprintf("[进度] 正在处理第 %d/%d 块\n", i, length(idx_chunks)))
  D_sub <- as.matrix(proxy::dist(emb[chunk, , drop = FALSE], emb))
  s_local <- numeric(length(chunk))
  for (j in seq_along(chunk)) {
    gi    <- chunk[j]
    dists <- D_sub[j, ]
    same <- which(labels == labels[gi])
    a_i  <- if (length(same) > 1) mean(dists[same[same != gi]]) else 0
    other_types <- setdiff(unique(labels), labels[gi])
    b_vals      <- vapply(other_types,
                          function(ct) mean(dists[labels == ct]),
                          numeric(1))
    b_i <- min(b_vals)
    s_local[j] <- (b_i - a_i) / max(a_i, b_i)
  }
  names(s_local) <- chunk
  s_local
}
# 6. 汇总，并计算最终 ASW
s_vec           <- numeric(N)
s_vec[as.integer(names(s_sampled))] <- s_sampled
cell_type_ASW   <- mean(s_vec, na.rm = TRUE)

print(paste("cell_type ASW:", cell_type_ASW))

#读取已有的结果文件
output_path <- file.path(
  output_dir,
  paste0("cell_type_cluster_metrics_with_", integrationmethod, ".txt")
)

# 1. 读取原有内容
output <- readLines(output_path)

# 2. 生成要追加的新行
ASW <- paste0("cell_type ASW: ", cell_type_ASW)

# 3. 合并到末尾
output <- c(output, ASW)

# 保存结果到文件
writeLines(output, con = paste0(output_dir, "cell_type_cluster_metrics_with_",integrationmethod,".txt"))
