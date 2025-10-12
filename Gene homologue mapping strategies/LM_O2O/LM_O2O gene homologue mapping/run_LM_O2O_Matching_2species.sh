#!/bin/bash
#SBATCH -J match
#SBATCH --nodelist=node3
#SBATCH --cpus-per-task=1
#SBATCH --error=/dev/null                      # 禁用 sbatch 的默认错误文件
#SBATCH --output=/dev/null                     # 禁用 sbatch 的默认输出文件
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1 
# 获取当前时间，格式为 YYYYMMDD_HHMMSS
timestamp=$(date +%Y%m%d_%H%M%S)
# 定义日志目录和日志文件名

file_dir="/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/LM_O2O/Species_mapping"




# 将标准输出和标准错误重定向到包含时间戳的日志文件

#####
#human_mice ✅
#chimpanzee_gorilla ✅
#crab_eating_macaque_mice ✅
#crab_eating_macaque_rhesus_macaque ✅
#gorilla_rhesus_macaque ✅
#human_chimpanzee ✅
#human_microcebus ✅
#microcebus_human ✅
#mice_frog ✅
#mice_pig ✅
#rhesus_macaque_human ✅
#rhesus_macaques_commonmarmoset ✅
#human_pig ✅
#zebrafish_frog ✅
#####
NAME_1="Microcebus_murinus.Mmur_3.0.pep.all" #更改为整合的物种1
NAME_2="Homo_sapiens.GRCh38.pep.all" #更改为整合的物种2
DEVICE="cuda:0"
#####
NAME="${NAME_1}_${NAME_2}"
log_dir="$file_dir/logfiles"
log_file="${log_dir}/log_${timestamp}_${NAME}.log"
exec > "${log_file}" 2>&1
mkdir -p "${log_dir}"
RESULTPATH="$file_dir/Results" #保存至指定的跨物种任务的Match方法文件夹下 
mkdir -p "$RESULTPATH"
DATAPATH="/cluster2/home/zeyu/Projects/Program/cross_species_integration/data" 
Embeddingspath_1="${DATAPATH}/${NAME_1}.gene_symbol_to_embedding_ESM2_15B.pt" 
Embeddingspath_2="${DATAPATH}/${NAME_2}.gene_symbol_to_embedding_ESM2_15B.pt" 
outputlocation="${RESULTPATH}/${NAME_1}_${NAME_2}/"
mkdir -p "$outputlocation"

source  "/cluster/home/zeyu/miniconda3/etc/profile.d/conda.sh"
echo "Start time: $(date)"

echo "——————Generatingcorrelationmatrix...——————" 
conda  activate torchI
echo "开始生成相关系数矩阵及偏好矩阵"
python -u /cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/correlation.py \
 $Embeddingspath_1 \
 $Embeddingspath_2 \
 $NAME_1 \
 $NAME_2 \
 $outputlocation \
 $DEVICE
echo "End time: $(date)"  
# conda activate MatchI
# echo "开始匹配"
# cd "$outputlocation"
# Rscript /cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Matching.R "$NAME_1" "$NAME_2"
# echo "End time: $(date)" 

