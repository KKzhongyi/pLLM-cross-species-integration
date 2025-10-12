import pandas as pd
from pathlib import Path
from typing import List, Set, Tuple
import os


 
BASE_DIR  = Path("/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods")
O2O_DIR   = BASE_DIR / "HM_O2O/cross_species_homologue_genes"
HOMO_RES  = O2O_DIR / "o2oResults"              # homologue greedy 结果
ONN_RES   = BASE_DIR / "ONN/initial/ONNResults" # correlation greedy 结果
RESULT_DIR = Path(                                   
    "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/HL_O2O/initial/HL_O2O_Results"
)
RESULT_DIR.mkdir(exist_ok=True, parents=True)        
 
PAIRS: List[Tuple[str, str, str]] = [
    #species_1 ,      species_2 ,                                  match_dir
    # ("human",            "mice",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Homo_sapiens.GRCh38.pep.all_Mus_musculus.GRCm39.pep.all"),
    # ("chimpanzee",      "gorilla",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Pan_troglodytes.Pan_tro_3.0.pep.all_Gorilla_gorilla.gorGor4.pep.all"),
    # ("crab_eating_macaque", "mice",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all_Mus_musculus.GRCm39.pep.all"),
    # ("crab_eating_macaque", "rhesus_macaque",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all_Macaca_mulatta.Mmul_10.pep.all"),
    # ("gorilla",          "rhesus_macaque",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Gorilla_gorilla.gorGor4.pep.all_Macaca_mulatta.Mmul_10.pep.all"),
    # ("human",             "chimpanzee",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Homo_sapiens.GRCh38.pep.all_Pan_troglodytes.Pan_tro_3.0.pep.all"),
    ("microcebus",        "human",
     "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Microcebus_murinus.Mmur_3.0.pep.all_Homo_sapiens.GRCh38.pep.all"),
    # ("mice",              "frog",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Mus_musculus.GRCm39.pep.all_Xenopus_tropicalis.UCB_Xtro_10.0.pep.all"),
    # ("mice",              "pig",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Mus_musculus.GRCm39.pep.all_Sus_scrofa.Sscrofa11.1.pep.all"),
    # ("rhesus_macaque",    "human",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Macaca_mulatta.Mmul_10.pep.all_Homo_sapiens.GRCh38.pep.all"),
    # ("rhesus_macaque",    "commonmarmoset",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Macaca_mulatta.Mmul_10.pep.all_Callithrix_jacchus.mCalJac1.pat.X.pep.all"),
    # ("human",    "pig",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Homo_sapiens.GRCh38.pep.all_Sus_scrofa.Sscrofa11.1.pep.all"),
    # ("zebrafish",    "frog",
    #  "/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/ONN_o2o/Species_mapping/Results/Danio_rerio.GRCz11.pep.all_Xenopus_tropicalis.UCB_Xtro_10.0.pep.all"),
     
]

 
def reprocess_raw(txt_path: Path) -> pd.DataFrame:
    raw = pd.read_csv(txt_path, sep="\t")
    raw[[raw.columns[3], raw.columns[4]]] = raw[[raw.columns[3], raw.columns[4]]].apply(
        pd.to_numeric, errors="coerce"
    )
    raw.rename(columns={raw.columns[-1]: "orthology confidence"}, inplace=True)
    raw["orthology confidence"] = pd.to_numeric(
    raw["orthology confidence"], errors="coerce"
    ).fillna(0).astype(int)
    raw["identical_scores"] = raw[[raw.columns[3], raw.columns[4]]].mean(axis=1)
    return raw[[raw.columns[0], raw.columns[1], "identical_scores", "orthology confidence"]]

def greedy_union(df: pd.DataFrame, col1: str, col2: str) -> pd.DataFrame:
    """贪心匹配：先按 orthology confidence 再按组合得分排，保证一对一。
    返回列：col1, col2, correlation, identical_scores, orthology confidence, combo_score"""

    # 1. 缺失值处理
    for col in ["correlation", "identical_scores"]:
        df[col] = df[col].fillna(df[col].median())
    df["orthology confidence"] = df["orthology confidence"].fillna(0)

    # 2. 计算组合得分 combo_score
    for col in ["correlation", "identical_scores"]:
        rng = df[col].max() - df[col].min()
        df[f"norm_{col}"] = 0.5 if rng == 0 else (df[col] - df[col].min()) / rng
    df["combo_score"] = df["norm_correlation"] + df["norm_identical_scores"]

    # 3. 贪心挑选
    used_f, used_s, rows = set(), set(), []
    sorter = df.sort_values(
        by=["orthology confidence", "combo_score"],
        ascending=[False, False],
        kind="mergesort",
    )
    for _, r in sorter.iterrows():
        a, b = r[col1], r[col2]
        if a not in used_f and b not in used_s:
            rows.append(
                (
                    a,
                    b,
                    r["correlation"],           
                    r["identical_scores"],        
                    r["orthology confidence"],    
                    r["combo_score"],            
                )
            )
            used_f.add(a)
            used_s.add(b)

    return pd.DataFrame(
        rows,
        columns=[
            col1,
            col2,
            "correlation",
            "identical_scores",
            "orthology confidence",
            "combo_score",
        ],
    )

def save_diff(base: Set[Tuple[str,str]], new: Set[Tuple[str,str]], path: Path):
    d_new, d_old = new-base, base-new
    if d_new or d_old:
        rows = ([(a,b,"only_in_union")  for a,b in d_new] +
                [(a,b,"only_in_single") for a,b in d_old])
        pd.DataFrame(rows, columns=["A","B","source"]).to_csv(path, sep="\t", index=False)
        print(f"  ⚠ 差异 {len(rows)} 行 → {path.name}")

def save_dup(df: pd.DataFrame, sp1: str, sp2: str):
    dup = df[df[df.columns[0]].duplicated(keep=False) |
             df[df.columns[1]].duplicated(keep=False)]
    if not dup.empty:
        p = RESULT_DIR / f"{sp1}_{sp2}_union_duplicates.txt"   
        dup.to_csv(p, sep="\t", index=False)
        print(f"  ⚠ 重复基因 {dup.shape[0]} 行 → {p.name}")


 
for sp1, sp2, match_dir in PAIRS:                              
    print(f"\n=== 综合配对 {sp1} ↔ {sp2} ===")
    match_dir = Path(match_dir)                                

    raw_txt  = O2O_DIR / f"{sp1}_{sp2}_o2o.txt"
    hom_csv  = HOMO_RES / f"{sp1}_{sp2}_o2o.csv"
    corr_csv = ONN_RES  / f"{sp1}_{sp2}_o2o_all.csv"

    # 2.1 correlation matrix（原始）                          
    try:
        corr_mat_path = next(match_dir.glob("*gene_embeddings_correlation_matrix.csv"))
    except StopIteration:
        print("  ❌ correlation matrix 缺失，跳过")
        continue
    corr_mat = pd.read_csv(corr_mat_path, index_col=0)
    corr_med = corr_mat.stack().median()
    gene_set_A, gene_set_B = set(corr_mat.index), set(corr_mat.columns)

    # 2.2 输入文件检查
    if not (raw_txt.exists() and hom_csv.exists() and corr_csv.exists()):
        print("  ❌ homologue / ONN 文件缺失，跳过")
        continue
    
    colA = f"{sp1} gene name"
    colB = f"{sp2} gene name" 

    # 3) 读取两张单方法结果作为候选
    hom_df  = pd.read_csv(hom_csv,usecols=range(2))
    hom_df.columns  = [colA, colB]
    
    corr_df = pd.read_csv(corr_csv,usecols=range(2))
    corr_df.columns = [colA, colB]

    # 3) homologue 原 txt
    raw_hom = reprocess_raw(raw_txt)
    raw_hom.rename(
        columns={raw_hom.columns[0]: colA, raw_hom.columns[1]: colB},
        inplace=True
    )
    ident_med = raw_hom["identical_scores"].median()

    # 5) 最大候选集并补 identical / confidence
    all_pairs = (
        pd.concat([hom_df, corr_df], ignore_index=True)
          .drop_duplicates()
          .merge(
              raw_hom[[colA, colB, "identical_scores", "orthology confidence"]],
              how="left", on=[colA, colB]
          )                             
    )
    all_pairs["identical_scores"].fillna(ident_med, inplace=True)
    all_pairs["orthology confidence"].fillna(0, inplace=True)

    # 6) 用原始 correlation matrix 补缺值 --------------------- 
    def fetch_corr(r):
        a, b = r[colA], r[colB]
        return corr_mat.at[a, b] if a in gene_set_A and b in gene_set_B else pd.NA
    all_pairs["correlation"] = pd.NA
    all_pairs["correlation"] = all_pairs.apply(fetch_corr, axis=1)
    all_pairs["correlation"].fillna(corr_med, inplace=True)
    print(all_pairs.columns)
    # 7) 联合 greedy
    union_df = greedy_union(all_pairs, colA, colB)
    out_path = RESULT_DIR / f"{sp1}_{sp2}_union_greedy.csv"        
    union_df.to_csv(out_path, index=False)
    print(f"  ✅ 综合配对 {union_df.shape[0]} 对 → {out_path.name}")
    union_df_genes=union_df[[colA, colB]]
    print(union_df_genes.shape)
    # 8) 差异对比
    save_diff(set(hom_df.apply(tuple, axis=1)),
            set(union_df_genes.apply(tuple, axis=1)),
            RESULT_DIR / f"{sp1}_{sp2}_diff_homologue.txt")      
    save_diff(set(corr_df[[colA, colB]].apply(tuple, axis=1)),
            set(union_df_genes.apply(tuple, axis=1)),
            RESULT_DIR / f"{sp1}_{sp2}_diff_correlation.txt")

    # 9) 重复检测
    save_dup(union_df_genes, sp1, sp2)

 



