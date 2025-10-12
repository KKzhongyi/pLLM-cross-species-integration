#!/bin/bash
#SBATCH -J Aqueous_Humor
#SBATCH --nodelist=node3
#SBATCH --cpus-per-task=8
#SBATCH --error=/dev/null                      # 禁用 sbatch 的默认错误文件
#SBATCH --output=/dev/null                     # 禁用 sbatch 的默认输出文件
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1 
# 获取当前时间，格式为 YYYYMMDD_HHMMSS
timestamp=$(date +%Y%m%d_%H%M%S)
# 定义日志目录和日志文件名
NAME="Aqueous_Humor"

celltypenumbaseline=0
resolution=0.1
centroid_number=6
reduction="pca" #pca,integrated.cca
initial_k=20
cpu_number=8
integrationmethod=CCAIntegration #CCAIntegration,HarmonyIntegration
file_dir="/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/Cell_Atlas_of_Aqueous_Humor"
species_list="crab_eating_macaque, rhesus_macaque, human, mice,pig"
dataset="Cell_Atlas_of_Aqueous_Humor"
scripts_dir="/cluster2/home/zeyu/Projects/Program/cross_species_integration/Matching/Methods/scripts"

####
IFS=',' read -r -a arr <<< "$species_list"
for i in "${!arr[@]}"; do
  arr[$i]="$(echo "${arr[$i]}" | xargs)"
done
h5file=$(IFS=_; echo "${arr[*]}")
echo "$h5file"

log_dir="${file_dir}/logfiles"
log_file="${log_dir}/log_${timestamp}_${NAME}_${celltypenumbaseline}_${resolution}_${cutoff}_$integrationmethod.log"
mkdir -p "${log_dir}"
output_dir="${file_dir}/Results/${NAME}_${resolution}_$integrationmethod/"
mkdir -p "$output_dir"
# 将标准输出和标准错误重定向到包含时间戳的日志文件
exec > "${log_file}" 2>&1

source  "/cluster2/home/zeyu/miniconda3/etc/profile.d/conda.sh"
echo "Start time: $(date)"

echo "——————Generating combined datasets——————" 
conda activate torchI
echo "Generate combined datasets"
python -u $scripts_dir/extract.py "$file_dir" "$h5file" "$species_list" "$dataset"
echo "Combined datasets finished"

echo "——————Generating UMAP——————" 
conda activate MatchI
echo "Generate Matching UMAP"
Rscript $scripts_dir/UMAPgeneration.R "$output_dir" "$celltypenumbaseline" "$resolution" "$file_dir" "$h5file"  "$species_list" "$cpu_number" "$integrationmethod"  
echo "Matching UMAP finished"
echo "End time: $(date)" 

echo "——————Calculating ASW——————" 
conda activate MatchI
Rscript ${scripts_dir}/ASW.R "$output_dir" "${file_dir}" "${h5file}" "${cpu_number}" "$integrationmethod"
echo "Matching UMAP finished"
echo "End time: $(date)" 



