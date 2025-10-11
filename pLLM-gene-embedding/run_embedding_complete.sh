#function from SATURN https://github.com/snap-stanford/saturn
#!/bin/bash
#SBATCH -J Matching
#SBATCH --nodelist=node3
#SBATCH --cpus-per-task=8
#SBATCH --error=/dev/null                      # 禁用 sbatch 的默认错误文件
#SBATCH --output=/dev/null                     # 禁用 sbatch 的默认输出文件
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1 
#下载物种蛋白表达的fasta数据，清洗数据，提取蛋白质嵌入，将蛋白质嵌入转换为基因嵌入（键值对，基因名为键，嵌入为值）
# 设置变量
NAME="Sus_scrofa.Sscrofa11.1.pep.all"  # 修改为你想要的参考蛋白质组的名称
DATA_PATH="/cluster2/home/zeyu/Projects/Program/cross_species_integration/data"             # 数据目录的路径
SCRIPT_PATH="/cluster2/home/zeyu/Projects/Program/cross_species_integration/SATURN-main/protein_embeddings"     # 脚本目录的路径
ESM_PATH="/cluster2/home/zeyu/Projects/Program/cross_species_integration/SATURN-main/esm-main/scripts"  # ESM 仓库路径 调用ESM模型脚本
TORCH_HOME="/cluster2/home/zeyu/Projects/Program/cross_species_integration/Torch_home"  # Torch Home 调用ESM大语言模型
DEVICE=0  # GPU 编号，根据需要进行修改
#FASTA_URL="http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"  #human(Homo-sapiens)
#FASTA_URL="http://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz" #mice(Mus musculus)
#FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz"  #zebrafish（Danio rerio）
#FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/xenopus_tropicalis/pep/Xenopus_tropicalis.UCB_Xtro_10.0.pep.all.fa.gz"  #	Tropical clawed frog(Xenopus tropicalis)
#FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/microcebus_murinus/pep/Microcebus_murinus.Mmur_3.0.pep.all.fa.gz"  #	Microcebus(Microcebus_murinus)
#FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/pan_troglodytes/pep/Pan_troglodytes.Pan_tro_3.0.pep.all.fa.gz"  #	chimpanzee(Pan_troglodytes)
#FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/gorilla_gorilla/pep/Gorilla_gorilla.gorGor4.pep.all.fa.gz"  #	gorilla(Gorilla gorilla gorilla)
#FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa.gz"  #	rhesus macaques(Macaca mulatta)
#FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/callithrix_jacchus/pep/Callithrix_jacchus.mCalJac1.pat.X.pep.all.fa.gz"  #	White-tufted-ear marmoset(Callithrix jacchus)
FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/sus_scrofa/pep/Sus_scrofa.Sscrofa11.1.pep.all.fa.gz"  #	pig(Sus scrofa)
#FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/macaca_fascicularis/pep/Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all.fa.gz"  #	Crab-eating macaque(Macaca fascicularis)


timestamp=$(date +%Y%m%d_%H%M%S)
log_dir="/cluster2/home/zeyu/Projects/Program/cross_species_integration/logfiles"
log_file="${log_dir}/log_${timestamp}_embeddings_${NAME}.log"
mkdir -p $log_dir

export TORCH_HOME="$TORCH_HOME"

# 设置CUDA_VISIBLE_DEVICES环境变量并运行/cluster2/home/zeyu/miniconda3/envs/torch/bin/python脚本
{
    echo "Start time: $(date)"
    echo "——————Downloading and cleaning the data...——————"
    source  "/cluster/home/zeyu/miniconda3/etc/profile.d/conda.sh"
    conda activate torchI
    if [ ! -f $DATA_PATH/$NAME.fa.gz ]; then
        wget -r $FASTA_URL -O $DATA_PATH/$NAME.fa.gz
    else
        echo "$DATA_PATH/$NAME.fa.gz already exists, skipping download."
    fi

    if [ ! -f $DATA_PATH/$NAME.fa ]; then
        gunzip $DATA_PATH/$NAME.fa.gz
    else
        echo "$DATA_PATH/$NAME.fa already exists, skipping gunzip."
    fi

    if [ ! -f $DATA_PATH/$NAME.clean.fa ]; then
        python $SCRIPT_PATH/clean_fasta.py \
        --data_path=$DATA_PATH/$NAME.fa \
        --save_path=$DATA_PATH/$NAME.clean.fa
    else
        echo "$DATA_PATH/$NAME.clean.fa already exists, skipping clean_fasta.py."
    fi

    echo "——————Protein embeddings are processing...——————" 
    CUDA_VISIBLE_DEVICES="$DEVICE"  \
    python -u $ESM_PATH/extract.py esm2_t48_15B_UR50D \
    "$DATA_PATH/$NAME.clean.fa" \
    "$DATA_PATH/$NAME.clean.fa_esm2_15B" \
    --include mean --truncation_seq_length 1022 --toks_per_batch 2048
    wait
    echo "——————Genes symbols to protein IDs are mapping...——————"
    python -u $SCRIPT_PATH/map_gene_symbol_to_protein_ids.py \
    --fasta_path $DATA_PATH/$NAME.fa \
    --save_path $DATA_PATH/$NAME.gene_symbol_to_protein_ID.json 
    wait
    echo "——————Protein embeddings are converting to gene embeddings...——————"
    python -u $SCRIPT_PATH/convert_protein_embeddings_to_gene_embeddings.py \
    --embedding_dir $DATA_PATH/$NAME.clean.fa_esm2_15B \
    --gene_symbol_to_protein_ids_path $DATA_PATH/$NAME.gene_symbol_to_protein_ID.json \
    --embedding_model ESM2 \
    --save_path $DATA_PATH/$NAME.gene_symbol_to_embedding_ESM2_15B.pt 
wait
echo "$(realpath "$DATA_PATH/$NAME.gene_symbol_to_embedding_ESM2_15B.pt")"
echo "End time: $(date)"
} >> $log_file 2>&1
