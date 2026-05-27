#!/bin/bash
#SBATCH -J lm_o2o_match
#SBATCH --nodelist=node3
#SBATCH --cpus-per-task=1
#SBATCH --error=/dev/null
#SBATCH --output=/dev/null
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LM_O2O_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

NAME_1="${NAME_1:-Homo_sapiens.GRCh38.pep.all}"
NAME_2="${NAME_2:-Mus_musculus.GRCm39.pep.all}"
GENE_EMBEDDING_METHOD="${GENE_EMBEDDING_METHOD:-aggregate}"
DEVICE="${DEVICE:-cuda:0}"
GENE_EMBEDDING_ROOT="${GENE_EMBEDDING_ROOT:-$LM_O2O_DIR/pLLM_gene_embedding/Results}"
LEGACY_DATA_PATH="${LEGACY_DATA_PATH:-/cluster2/home/zeyu/Projects/Program/cross_species_integration/data}"
RESULT_ROOT="${RESULT_ROOT:-$LM_O2O_DIR/Species_mapping/Results}"
LOG_DIR="${LOG_DIR:-$LM_O2O_DIR/Species_mapping/logfiles}"
PYTHON_BIN="${PYTHON_BIN:-/cluster2/home/zeyu/miniconda3/envs/torchI/bin/python}"
CONDA_PROFILE="${CONDA_PROFILE:-/cluster2/home/zeyu/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-torchI}"
EMBEDDINGS_PATH_1="${EMBEDDINGS_PATH_1:-}"
EMBEDDINGS_PATH_2="${EMBEDDINGS_PATH_2:-}"

print_usage() {
    cat <<'EOF'
Usage:
  bash run_LM_O2O_Matching_2species.sh [options]

Options:
  --name-1 NAME                       First species Ensembl pep.all prefix
  --name-2 NAME                       Second species Ensembl pep.all prefix
  --gene-embedding-method METHOD      aggregate | max_pooling | canonical_isoform. Default: aggregate
  --method METHOD                     Alias for --gene-embedding-method
  --gene-embedding-root PATH          Root containing pLLM_gene_embedding/Results outputs
  --legacy-data-path PATH             Fallback root for legacy aggregate .pt files
  --embeddings-path-1 PATH            Optional explicit .pt file for species 1
  --embeddings-path-2 PATH            Optional explicit .pt file for species 2
  --result-root PATH                  Output root for correlation and preference matrices
  --log-dir PATH                      Log directory
  --device DEVICE                     Torch device passed to correlation.py. Default: cuda:0
  --python-bin PATH                   Python executable
  --conda-profile PATH                Optional conda profile script
  --conda-env NAME                    Optional conda environment name
  -h, --help                          Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --name-1)
            NAME_1="$2"; shift 2 ;;
        --name-2)
            NAME_2="$2"; shift 2 ;;
        --gene-embedding-method|--method)
            GENE_EMBEDDING_METHOD="$2"; shift 2 ;;
        --gene-embedding-root)
            GENE_EMBEDDING_ROOT="$2"; shift 2 ;;
        --legacy-data-path)
            LEGACY_DATA_PATH="$2"; shift 2 ;;
        --embeddings-path-1)
            EMBEDDINGS_PATH_1="$2"; shift 2 ;;
        --embeddings-path-2)
            EMBEDDINGS_PATH_2="$2"; shift 2 ;;
        --result-root)
            RESULT_ROOT="$2"; shift 2 ;;
        --log-dir)
            LOG_DIR="$2"; shift 2 ;;
        --device)
            DEVICE="$2"; shift 2 ;;
        --python-bin)
            PYTHON_BIN="$2"; shift 2 ;;
        --conda-profile)
            CONDA_PROFILE="$2"; shift 2 ;;
        --conda-env)
            CONDA_ENV="$2"; shift 2 ;;
        -h|--help)
            print_usage; exit 0 ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            print_usage >&2
            exit 1 ;;
    esac
done

case "$GENE_EMBEDDING_METHOD" in
    aggregate|max_pooling|canonical_isoform) ;;
    *)
        echo "[ERROR] Unsupported --gene-embedding-method: $GENE_EMBEDDING_METHOD" >&2
        echo "Allowed values: aggregate, max_pooling, canonical_isoform" >&2
        exit 1 ;;
esac

resolve_embedding_path() {
    local name="$1"
    local path=""
    case "$GENE_EMBEDDING_METHOD" in
        aggregate)
            path="$GENE_EMBEDDING_ROOT/$name/aggregate/$name.gene_symbol_to_embedding_ESM2_15B.pt"
            if [[ ! -f "$path" && -f "$LEGACY_DATA_PATH/$name.gene_symbol_to_embedding_ESM2_15B.pt" ]]; then
                path="$LEGACY_DATA_PATH/$name.gene_symbol_to_embedding_ESM2_15B.pt"
            fi
            echo "$path"
            ;;
        max_pooling)
            path="$GENE_EMBEDDING_ROOT/$name/max_pooling/$name.max_pooling_gene_symbol_to_embedding_ESM2_15B.pt"
            echo "$path"
            ;;
        canonical_isoform)
            path="$GENE_EMBEDDING_ROOT/$name/canonical_isoform/$name.canonical_gene_symbol_to_embedding_ESM2_15B.pt"
            echo "$path"
            ;;
    esac
}

if [[ -z "$EMBEDDINGS_PATH_1" ]]; then
    EMBEDDINGS_PATH_1="$(resolve_embedding_path "$NAME_1")"
fi
if [[ -z "$EMBEDDINGS_PATH_2" ]]; then
    EMBEDDINGS_PATH_2="$(resolve_embedding_path "$NAME_2")"
fi

if [[ ! -f "$EMBEDDINGS_PATH_1" ]]; then
    echo "[ERROR] Missing species 1 embedding file: $EMBEDDINGS_PATH_1" >&2
    exit 1
fi
if [[ ! -f "$EMBEDDINGS_PATH_2" ]]; then
    echo "[ERROR] Missing species 2 embedding file: $EMBEDDINGS_PATH_2" >&2
    exit 1
fi

PAIR_NAME="${NAME_1}_${NAME_2}"
OUTPUT_DIR="$RESULT_ROOT/$GENE_EMBEDDING_METHOD/$PAIR_NAME"
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/log_$(date +%Y%m%d_%H%M%S)_${GENE_EMBEDDING_METHOD}_${PAIR_NAME}.log"

{
    echo "Start time: $(date)"
    echo "NAME_1: $NAME_1"
    echo "NAME_2: $NAME_2"
    echo "GENE_EMBEDDING_METHOD: $GENE_EMBEDDING_METHOD"
    echo "EMBEDDINGS_PATH_1: $EMBEDDINGS_PATH_1"
    echo "EMBEDDINGS_PATH_2: $EMBEDDINGS_PATH_2"
    echo "OUTPUT_DIR: $OUTPUT_DIR"
    echo "DEVICE: $DEVICE"

    source "$CONDA_PROFILE"
    conda activate "$CONDA_ENV"

    "$PYTHON_BIN" -u "$SCRIPT_DIR/correlation.py" \
        "$EMBEDDINGS_PATH_1" \
        "$EMBEDDINGS_PATH_2" \
        "$NAME_1" \
        "$NAME_2" \
        "$OUTPUT_DIR/" \
        "$DEVICE"

    echo "End time: $(date)"
} >> "$LOG_FILE" 2>&1

echo "Log: $LOG_FILE"
echo "Output: $OUTPUT_DIR"
