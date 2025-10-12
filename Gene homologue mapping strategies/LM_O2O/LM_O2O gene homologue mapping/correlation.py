from re import T
import torch
import pandas as pd
import os
import numpy as np
import argparse
import h5py


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate adjacency matrix from embeddings."
    )
    parser.add_argument(
        "embeddings_path_1",
        type=str,
        help="Path to the gene symbol to embedding mapping.",
    )
    parser.add_argument(
        "embeddings_path_2",
        type=str,
        help="Path to the gene symbol to embedding mapping.",
    )
    parser.add_argument(
        "NAME_1",
        type=str,
        help="Path to the gene symbol to embedding mapping.",
    )
    parser.add_argument(
        "NAME_2",
        type=str,
        help="Path to the gene symbol to embedding mapping.",
    )
    parser.add_argument(
        "output_path", type=str, help="Path to save the correlation matrix."
    )
    parser.add_argument("cuda", default="cuda:0", type=str, help="CUDA device to use.")
    return parser.parse_args()



args = parse_args()
DEVICE = torch.device(args.cuda)

# 加载 gene_embedding_1
gene_embedding_1 = torch.load(args.embeddings_path_1)
print("Number of embeddings_1:", len(gene_embedding_1))
gene_embedding_1 = pd.DataFrame.from_dict(gene_embedding_1, orient="index")
gene_embedding_1.index.name = "Gene"
gene_embedding_1 = gene_embedding_1.astype(float)
print("Shape of embeddings_1:", gene_embedding_1.shape)
gene_embedding_1_ids = gene_embedding_1.index.tolist()

# 加载 gene_embedding_2
gene_embedding_2 = torch.load(args.embeddings_path_2)
print("Number of embeddings_2:", len(gene_embedding_2))
gene_embedding_2 = pd.DataFrame.from_dict(gene_embedding_2, orient="index")
gene_embedding_2.index.name = "Gene"
gene_embedding_2 = gene_embedding_2.astype(float)
print("Shape of embeddings_2:", gene_embedding_2.shape)
gene_embedding_2_ids = gene_embedding_2.index.tolist()

# 将数据加载到 GPU 上
gene_embedding_1_tensor = torch.tensor(
    gene_embedding_1.values, dtype=torch.float32
).to(DEVICE)
del gene_embedding_1
gene_embedding_2_tensor = torch.tensor(
    gene_embedding_2.values, dtype=torch.float32
).to(DEVICE)
del gene_embedding_2
# 计算相关系数矩阵
x_mean = gene_embedding_1_tensor.mean(dim=1, keepdim=True)
y_mean = gene_embedding_2_tensor.mean(dim=1, keepdim=True)
xm = gene_embedding_1_tensor - x_mean
ym = gene_embedding_2_tensor - y_mean
r_num = torch.mm(xm, ym.t())
r_den = torch.sqrt(torch.sum(xm**2, dim=1, keepdim=True)) * torch.sqrt(
    torch.sum(ym**2, dim=1)
)
r = r_num / (r_den + 1e-8)
del x_mean, y_mean, xm, ym, r_num, r_den  
torch.cuda.empty_cache()

# 初始化 gene_embedding_1 和 gene_embedding_2 的偏好矩阵
gene_embedding_1_matrix = torch.empty(
    (len(gene_embedding_1_ids), len(gene_embedding_2_ids)),
    dtype=torch.long,
    device=DEVICE,
)
gene_embedding_2_matrix = torch.empty(
    (len(gene_embedding_2_ids), len(gene_embedding_1_ids)),
    dtype=torch.long,
    device=DEVICE,
)

# 生成 gene_embedding_1 的偏好矩阵（学生）
for i, gene_id in enumerate(gene_embedding_1_ids):
    correlations = r[i]  # 获取 gene_id 与所有 gene_embedding_2 的相关系数
    _, sorted_indices = torch.topk(
        correlations, k=len(gene_embedding_2_ids), largest=True
    )
    gene_embedding_1_matrix[i] = sorted_indices  # 按排序填入偏好矩阵

# 生成 gene_embedding_2 的偏好矩阵（大学）
for j, gene_id in enumerate(gene_embedding_2_ids):
    correlations = r[:, j]  # 获取 gene_id 与所有 gene_embedding_1 的相关系数
    _, sorted_indices = torch.topk(
        correlations, k=len(gene_embedding_1_ids), largest=True
    )
    gene_embedding_2_matrix[j] = sorted_indices  # 按排序填入偏好矩阵

# 将相关系数矩阵和偏好矩阵从 GPU 转移到 CPU 并转换为 DataFrame
correlations = r.cpu()
del r, gene_embedding_1_tensor, gene_embedding_2_tensor
torch.cuda.empty_cache() 
correlation_df = pd.DataFrame(
    correlations.numpy(), index=gene_embedding_1_ids, columns=gene_embedding_2_ids
) 
gene_embedding_1_matrix_cpu = gene_embedding_1_matrix.cpu().numpy()
del gene_embedding_1_matrix
gene_embedding_2_matrix_cpu = gene_embedding_2_matrix.cpu().numpy()
del gene_embedding_2_matrix
gene_embedding_1_df = pd.DataFrame(
    gene_embedding_1_matrix_cpu,
    index=gene_embedding_1_ids,
    columns=gene_embedding_2_ids,
)
del gene_embedding_1_matrix_cpu
gene_embedding_2_df = pd.DataFrame(
    gene_embedding_2_matrix_cpu,
    index=gene_embedding_2_ids,
    columns=gene_embedding_1_ids,
)
del gene_embedding_2_matrix_cpu

# 将偏好矩阵转变为matchingR可以直接读入的形式
gene_embedding_1_df = gene_embedding_1_df.T
gene_embedding_2_df = gene_embedding_2_df.T

# 保存相关系数矩阵和偏好矩阵为 CSV 文件
correlation_df.to_csv(
    f"{args.output_path}"
    + f"{args.NAME_1}_{args.NAME_2}_gene_embeddings_correlation_matrix.csv",index=True,header=True
)
gene_embedding_1_df.to_csv(
    f"{args.output_path}" + "gene_embedding_1_matrix.csv", index=True
)
gene_embedding_2_df.to_csv(
    f"{args.output_path}" + "gene_embedding_2_matrix.csv", index=True
)

print(
    "已生成并保存相关系数矩阵和偏好矩阵;\n"
    "gene_embedding_1偏好矩阵,行名为gene_embedding_2，列名为gene_embedding_1;\n"
    "每一列是gene_embedding_1的偏好排序,可以作为matchingR的输入;\n"
)
