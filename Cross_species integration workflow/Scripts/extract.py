import pandas as pd
import numpy as np
import h5py
import argparse
import os
from functools import reduce
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Correlation cutoff")
    parser.add_argument(
        "file_dir",
        type=str,
    )
    parser.add_argument(
        "h5filename",
        type=str,
    )
    parser.add_argument(
        "species_list",
        type=str,
    )
    parser.add_argument(
        "dataset",
        type=str,
    )
    return parser.parse_args()


args = parse_args()
os.chdir(args.file_dir.strip())
species_list = [s.strip() for s in args.species_list.split(",")]

# # 整合三个物种的基因


homologue_table = pd.read_csv(
    args.file_dir.strip()+"/"+args.h5filename.strip()+".csv",
    header=0,
    index_col=False,
)
n_species = len(species_list)

homologue_table = homologue_table.dropna(
    subset=homologue_table.columns[:n_species]
)

orig_cols = list(homologue_table.columns[:n_species])
homologue_table[orig_cols] = homologue_table[orig_cols].apply(lambda col: col.str.lower())

for idx, orig in enumerate(orig_cols):
    new_name = f"{species_list[idx]} gene name"
    homologue_table.rename(columns={orig: new_name}, inplace=True)
    
print(homologue_table.shape)
print(homologue_table.iloc[-5:, -5:])

#读取单个数据集
def load_dataset(dataset: str, species: str) -> pd.DataFrame:
    file_path = f"/cluster2/home/zeyu/Projects/Program/cross_species_integration/data/datasets/{dataset.strip()}/raw_data/{species}_df.h5"

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"文件不存在: {file_path}")

    with h5py.File(file_path, "r") as h5f:
        # 读取数值型数据
        numeric_data = h5f["data/numeric"][:]
        numeric_columns = h5f["meta/numeric_columns"][:].astype(str)

        # 读取字符串型数据
        string_data = h5f["data/strings"][:].astype(str)
        string_columns = h5f["meta/string_columns"][:].astype(str)

        # 读取索引
        index = h5f["meta/index"][:].astype(str)

    # 构建 DataFrame
    numeric_df = pd.DataFrame(numeric_data, columns=numeric_columns, index=index)
    string_df = pd.DataFrame(string_data, columns=string_columns, index=index)

    # 合并数据
    df = pd.concat([numeric_df, string_df], axis=1)
    df.columns = df.columns.str.lower()
    # 标准化 cell_type
    if "cell_type" in df.columns:
        df["cell_type"] = df["cell_type"].str.lower()

    return df

#读取多个数据集
def load_multiple_species_datasets(species_list: list, dataset: str) -> dict:
    data_dict = {}
    for species in species_list:
        print(f"Loading dataset for species: {species}")
        df = load_dataset(dataset, species)
        df["species"] = species  # 添加物种列
        data_dict[species] = df
    return data_dict

dataset = args.dataset.strip()
all_species_data = load_multiple_species_datasets(species_list, dataset)



metadata_cols = 2 
expr_data = {}
meta_data = {}
for sp in species_list:
    df = all_species_data[sp]
    expr_data[sp] = df.iloc[:, :-metadata_cols]
    meta_data[sp] = df.iloc[:, -metadata_cols:]
    print(f"{sp} expr shape: {expr_data[sp].shape}")
    print(f"{sp} meta shape: {meta_data[sp].shape}")

for sp in species_list:
    expr_data[sp] = expr_data[sp].loc[:, ~expr_data[sp].columns.duplicated()]


for sp in species_list:
    gene_col = f"{sp} gene name"
    valid = set(homologue_table[gene_col])
    expr_data[sp] = expr_data[sp].loc[:, expr_data[sp].columns.isin(valid)]
    print(f"After intersect with homology ({sp}):", expr_data[sp].shape)

expanded = {sp: {} for sp in species_list}
seen = {sp: set() for sp in species_list}

for _, row in homologue_table.iterrows():
    # new gene name 按 species_list 顺序拼接
    new_gene = "_".join(row[f"{sp} gene name"] for sp in species_list)
    for sp in species_list:
        orig = row[f"{sp} gene name"]
        if orig in expr_data[sp].columns and new_gene not in seen[sp]:
            expanded[sp][new_gene] = expr_data[sp][orig]
            seen[sp].add(new_gene)
            
expanded_df = {sp: pd.DataFrame(expanded[sp]) for sp in species_list}

common_genes = reduce(
    lambda a, b: a & b,
    (set(expanded_df[sp].columns) for sp in species_list)
)
print("Common new genes:", len(common_genes))

final_expr = {
    sp: expanded_df[sp].loc[:, sorted(common_genes)]
    for sp in species_list
}

celltype_sets = {
    sp: set(meta_data[sp]["cell_type"])
    for sp in species_list
}
common_celltypes = set.intersection(*celltype_sets.values())
print("Common cell types:", common_celltypes)

filtered_meta = {
    sp: meta_data[sp][meta_data[sp]["cell_type"].isin(common_celltypes)]
    for sp in species_list
}
filtered_expr = {
    sp: final_expr[sp].loc[filtered_meta[sp].index]
    for sp in species_list
}

# 10. 合并所有物种
combined_expr = pd.concat([filtered_expr[sp] for sp in species_list], axis=0)
combined_meta = pd.concat([filtered_meta[sp] for sp in species_list], axis=0)
combined_meta = combined_meta.loc[combined_expr.index]

# 11. 最终合并
combined_expr = combined_expr[~combined_expr.index.duplicated(keep='first')]
combined_meta = combined_meta[~combined_meta.index.duplicated(keep='first')]
final_df = pd.concat([combined_expr, combined_meta], axis=1)
print("Final dataframe shape:", final_df.shape)
print(final_df.iloc[-5:, -5:])

# 保存到 HDF5 文件
with h5py.File(args.h5filename.strip()+".h5", "w") as h5f:
    # 保存数值型数据
    numeric_data = final_df.iloc[:, :-2].to_numpy()  # 提取数值列
    h5f.create_dataset("data/numeric", data=numeric_data)

    # 保存字符串型数据
    string_data = final_df.iloc[:, -2:].to_numpy().astype("S")  # 提取字符列并转换为字节
    h5f.create_dataset("data/strings", data=string_data)

    # 保存列名
    numeric_columns = final_df.columns[:-2].to_numpy(dtype="S")  # 数值型列名
    string_columns = final_df.columns[-2:].to_numpy(dtype="S")  # 字符串型列名
    h5f.create_dataset("meta/numeric_columns", data=numeric_columns)
    h5f.create_dataset("meta/string_columns", data=string_columns)

    # 保存索引
    index = final_df.index.to_numpy(dtype="S")
    h5f.create_dataset("meta/index", data=index)
