#!/bin/bash
#SBATCH -J UnifiedGeneEmbedding
#SBATCH --nodelist=node3
#SBATCH --cpus-per-task=4
#SBATCH --error=/dev/null
#SBATCH --output=/dev/null
#SBATCH --partition=normal
#SBATCH --gres=gpu:2

# Unified gene embedding workflow.
#
# Three strategies are supported through --gene-embedding-method:
#   aggregate         
#   max_pooling       
#   canonical_isoform  
#
# Example:
#   bash run_embedding_complete.sh --name Homo_sapiens.GRCh38.pep.all --gene-embedding-method aggregate
#   bash run_embedding_complete.sh --name Homo_sapiens.GRCh38.pep.all --gene-embedding-method max_pooling
#   bash run_embedding_complete.sh --name Homo_sapiens.GRCh38.pep.all --gene-embedding-method canonical_isoform \
#       --canonical-isoform-file /path/to/gene_to_canonical_protein.json

set -euo pipefail

# -----------------------------
# User-editable defaults
# -----------------------------
# Human:                 NAME="Homo_sapiens.GRCh38.pep.all"
# Mouse:                 NAME="Mus_musculus.GRCm39.pep.all"
# Zebrafish:             NAME="Danio_rerio.GRCz11.pep.all"
# Tropical clawed frog:  NAME="Xenopus_tropicalis.UCB_Xtro_10.0.pep.all"
# Microcebus:            NAME="Microcebus_murinus.Mmur_3.0.pep.all"
# Chimpanzee:            NAME="Pan_troglodytes.Pan_tro_3.0.pep.all"
# Gorilla:               NAME="Gorilla_gorilla.gorGor4.pep.all"
# Rhesus macaque:        NAME="Macaca_mulatta.Mmul_10.pep.all"
# Marmoset:              NAME="Callithrix_jacchus.mCalJac1.pat.X.pep.all"
# Pig:                   NAME="Sus_scrofa.Sscrofa11.1.pep.all"
# Crab-eating macaque:   NAME="Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all"
NAME="${NAME:-Sus_scrofa.Sscrofa11.1.pep.all}"

# Default reviewer-control parameter. aggregate keeps the previous Matching behavior.
GENE_EMBEDDING_METHOD="${GENE_EMBEDDING_METHOD:-aggregate}"

# Optional file for canonical mode. Supports JSON/CSV/TSV gene -> canonical protein ID.
# If empty and method=canonical_isoform, this script will run select_canonical_isoforms.py.
CANONICAL_ISOFORM_FILE="${CANONICAL_ISOFORM_FILE:-}"

DATA_PATH="${DATA_PATH:-/cluster2/home/zeyu/Projects/Program/cross_species_integration/data}"
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMMON_SCRIPT_PATH="$SCRIPT_PATH"
CANONICAL_METHOD_SCRIPT_PATH="$SCRIPT_PATH/canonical_isoform"
ESM_PATH="${ESM_PATH:-/cluster2/home/zeyu/Projects/Program/cross_species_integration/SATURN-main/esm-main/scripts}"
TORCH_HOME="${TORCH_HOME:-/cluster2/home/zeyu/Projects/Program/cross_species_integration/Torch_home}"
PYTHON_BIN="${PYTHON_BIN:-/cluster2/home/zeyu/miniconda3/envs/torchI/bin/python}"
CONDA_PROFILE="${CONDA_PROFILE:-/cluster2/home/zeyu/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-torchI}"
DEVICE="${DEVICE:-1}"
START_DELAY_HOURS="${START_DELAY_HOURS:-0}"
FASTA_URL="${FASTA_URL:-}"
RUN_ESM="${RUN_ESM:-1}"
ESM_TOKS_PER_BATCH="${ESM_TOKS_PER_BATCH:-2048}"
APPRIS_BATCH_SIZE="${APPRIS_BATCH_SIZE:-50}"
APPRIS_SPECIES="${APPRIS_SPECIES:-}"
UNIPROT_ORGANISM="${UNIPROT_ORGANISM:-}"
ALL_ISOFORM_EMBEDDING_DIR="${ALL_ISOFORM_EMBEDDING_DIR:-}"
RESULTS_ROOT="${RESULTS_ROOT:-$SCRIPT_PATH/Results}"
LOG_DIR="${LOG_DIR:-$SCRIPT_PATH/logfile}"

print_usage() {
    cat <<'EOF'
Usage:
  bash run_embedding_complete.sh [options]

Options:
  --name NAME                         Ensembl pep.all prefix, e.g. Homo_sapiens.GRCh38.pep.all
  --gene-embedding-method METHOD      aggregate | max_pooling | canonical_isoform. Default: aggregate
  --method METHOD                     Alias for --gene-embedding-method
  --canonical-isoform-file PATH       Optional JSON/CSV/TSV gene -> canonical protein ID mapping for canonical_isoform mode
  --data-path PATH                    Directory containing NAME.fa or NAME.fa.gz
  --fasta-url URL                     Optional URL used only if NAME.fa and NAME.fa.gz are missing
  --device N                          PyTorch cuda device index passed to ESM extract.py. Default: 1
  --python-bin PATH                   Python executable. Default: torchI environment Python
  --conda-profile PATH                Optional conda profile script
  --conda-env NAME                    Optional conda environment name
  --esm-path PATH                     Directory containing ESM extract.py
  --run-esm 0|1                       If 0, skip ESM extraction and only convert existing .pt files. Default: 1
  --all-isoform-embedding-dir PATH     Optional existing all-isoform ESM .pt directory for aggregate/max_pooling
  --results-root PATH                  Output root directory. Default: pLLM_gene_embedding/Results
  --log-dir PATH                       Log directory. Default: pLLM_gene_embedding/logfile
  --start-delay-hours N               Sleep N hours before loading conda/CUDA. Default: 0
  -h, --help                          Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --name)
            NAME="$2"; shift 2 ;;
        --gene-embedding-method|--method)
            GENE_EMBEDDING_METHOD="$2"; shift 2 ;;
        --canonical-isoform-file)
            CANONICAL_ISOFORM_FILE="$2"; shift 2 ;;
        --data-path)
            DATA_PATH="$2"; shift 2 ;;
        --fasta-url)
            FASTA_URL="$2"; shift 2 ;;
        --device)
            DEVICE="$2"; shift 2 ;;
        --python-bin)
            PYTHON_BIN="$2"; shift 2 ;;
        --conda-profile)
            CONDA_PROFILE="$2"; shift 2 ;;
        --conda-env)
            CONDA_ENV="$2"; shift 2 ;;
        --esm-path)
            ESM_PATH="$2"; shift 2 ;;
        --run-esm)
            RUN_ESM="$2"; shift 2 ;;
        --all-isoform-embedding-dir)
            ALL_ISOFORM_EMBEDDING_DIR="$2"; shift 2 ;;
        --results-root)
            RESULTS_ROOT="$2"; shift 2 ;;
        --log-dir)
            LOG_DIR="$2"; shift 2 ;;
        --start-delay-hours)
            START_DELAY_HOURS="$2"; shift 2 ;;
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

export TORCH_HOME="$TORCH_HOME"

FASTA_PATH="$DATA_PATH/$NAME.fa"
FASTA_GZ_PATH="$DATA_PATH/$NAME.fa.gz"
OUTPUT_BASE="$RESULTS_ROOT/$NAME"
METHOD_OUTPUT_DIR="$OUTPUT_BASE/$GENE_EMBEDDING_METHOD"
ALL_ISOFORMS_DIR="$OUTPUT_BASE/all_isoforms"
mkdir -p "$METHOD_OUTPUT_DIR" "$ALL_ISOFORMS_DIR"

FULL_CLEAN_FASTA="$ALL_ISOFORMS_DIR/$NAME.clean.fa"
FULL_ESM_DIR="$ALL_ISOFORMS_DIR/$NAME.clean.fa_esm2_15B"
GENE_TO_PROTEIN_MAP="$ALL_ISOFORMS_DIR/$NAME.gene_symbol_to_protein_ID.json"

LEGACY_FULL_ESM_DIR="$DATA_PATH/$NAME.clean.fa_esm2_15B"
if [[ -n "$ALL_ISOFORM_EMBEDDING_DIR" ]]; then
    FULL_ESM_DIR="$ALL_ISOFORM_EMBEDDING_DIR"
elif [[ -d "$LEGACY_FULL_ESM_DIR" ]] && [[ -n "$(find "$LEGACY_FULL_ESM_DIR" -maxdepth 1 -name '*.pt' -print -quit 2>/dev/null)" ]]; then
    # Reuse all-isoform embeddings from the legacy aggregate workflow when available.
    FULL_ESM_DIR="$LEGACY_FULL_ESM_DIR"
fi

CANONICAL_FASTA="$METHOD_OUTPUT_DIR/$NAME.canonical_isoform.fa"
CANONICAL_CLEAN_FASTA="$METHOD_OUTPUT_DIR/$NAME.canonical_isoform.clean.fa"
CANONICAL_ESM_DIR="$METHOD_OUTPUT_DIR/$NAME.canonical_isoform.clean.fa_esm2_15B"
CANONICAL_MAPPING="$METHOD_OUTPUT_DIR/$NAME.canonical_gene_symbol_to_protein_ID.json"
CANONICAL_MAPPING_FOR_CONVERSION="$CANONICAL_MAPPING"
CANONICAL_SELECTION_TABLE="$METHOD_OUTPUT_DIR/$NAME.canonical_isoform_selection.tsv"
CANONICAL_API_CACHE="$METHOD_OUTPUT_DIR/$NAME.canonical_isoform_api_cache.json"
CANONICAL_FASTA_BUILD_REPORT="$METHOD_OUTPUT_DIR/$NAME.canonical_fasta_build_report.tsv"

case "$GENE_EMBEDDING_METHOD" in
    aggregate)
        GENE_EMBEDDING_PATH="$METHOD_OUTPUT_DIR/$NAME.gene_symbol_to_embedding_ESM2_15B.pt"
        CONVERSION_REPORT="$METHOD_OUTPUT_DIR/$NAME.aggregate_embedding_conversion_report.tsv"
        ;;
    max_pooling)
        GENE_EMBEDDING_PATH="$METHOD_OUTPUT_DIR/$NAME.max_pooling_gene_symbol_to_embedding_ESM2_15B.pt"
        CONVERSION_REPORT="$METHOD_OUTPUT_DIR/$NAME.max_pooling_embedding_conversion_report.tsv"
        ;;
    canonical_isoform)
        GENE_EMBEDDING_PATH="$METHOD_OUTPUT_DIR/$NAME.canonical_gene_symbol_to_embedding_ESM2_15B.pt"
        CONVERSION_REPORT="$METHOD_OUTPUT_DIR/$NAME.canonical_embedding_conversion_report.tsv"
        ;;
esac

ANNOTATION_GTF_URL=""
case "$NAME" in
    "Homo_sapiens.GRCh38.pep.all")
        ANNOTATION_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz"
        ;;
    "Mus_musculus.GRCm39.pep.all")
        ANNOTATION_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.annotation.gtf.gz"
        ;;
esac
ANNOTATION_GTF=""
if [[ -n "$ANNOTATION_GTF_URL" ]]; then
    ANNOTATION_GTF="$METHOD_OUTPUT_DIR/$(basename "$ANNOTATION_GTF_URL")"
fi

timestamp=$(date +%Y%m%d_%H%M%S)
log_dir="$LOG_DIR"
log_file="$log_dir/log_${timestamp}_${GENE_EMBEDDING_METHOD}_embeddings_${NAME}.log"
mkdir -p "$log_dir"

prepare_fasta() {
    if [[ -f "$FASTA_PATH" ]]; then
        echo "$FASTA_PATH already exists; using it."
        return
    fi
    if [[ -f "$FASTA_GZ_PATH" ]]; then
        echo "Decompressing existing FASTA gzip without deleting it: $FASTA_GZ_PATH"
        gunzip -c "$FASTA_GZ_PATH" > "$FASTA_PATH"
        return
    fi
    if [[ -n "$FASTA_URL" ]]; then
        echo "Downloading FASTA from FASTA_URL: $FASTA_URL"
        wget -c "$FASTA_URL" -O "$FASTA_GZ_PATH"
        gunzip -c "$FASTA_GZ_PATH" > "$FASTA_PATH"
        return
    fi
    echo "[ERROR] Missing FASTA: $FASTA_PATH" >&2
    echo "[ERROR] Provide the file, provide $FASTA_GZ_PATH, or pass --fasta-url." >&2
    exit 1
}

check_cuda() {
    export UNIFIED_CUDA_DEVICE="$DEVICE"
    "$PYTHON_BIN" - <<'PY'
import os
import torch
requested = int(os.environ.get("UNIFIED_CUDA_DEVICE", "0"))
print(f"CUDA_VISIBLE_DEVICES env: {os.environ.get('CUDA_VISIBLE_DEVICES')}")
print(f"torch.cuda.is_available(): {torch.cuda.is_available()}")
print(f"torch.cuda.device_count(): {torch.cuda.device_count()}")
if torch.cuda.is_available():
    if requested >= torch.cuda.device_count():
        raise SystemExit(f"[ERROR] Requested cuda:{requested}, but only {torch.cuda.device_count()} CUDA device(s) are visible.")
    torch.cuda.set_device(requested)
    print(f"torch.cuda.current_device(): {torch.cuda.current_device()}")
    print(f"torch.cuda.get_device_name({requested}): {torch.cuda.get_device_name(requested)}")
PY
}

run_esm_if_needed() {
    local fasta="$1"
    local out_dir="$2"
    local description="$3"
    if [[ "$RUN_ESM" == "0" ]]; then
        echo "RUN_ESM=0; skipping ESM extraction for $description."
        return
    fi
    if [[ -d "$out_dir" ]] && [[ -n "$(find "$out_dir" -maxdepth 1 -name '*.pt' -print -quit 2>/dev/null)" ]]; then
        echo "$out_dir already contains .pt files; skipping ESM extraction for $description."
        return
    fi
    echo "Extracting $description protein embeddings with ESM2-15B..."
    if command -v nvidia-smi >/dev/null 2>&1; then
        nvidia-smi --query-gpu=index,name,memory.total,memory.used,memory.free,utilization.gpu --format=csv,noheader
    fi
    "$PYTHON_BIN" -u "$ESM_PATH/extract.py" esm2_t48_15B_UR50D \
        "$fasta" \
        "$out_dir" \
        --include mean --truncation_seq_length 1022 --toks_per_batch "$ESM_TOKS_PER_BATCH" --cuda_device "$DEVICE"
}

{
    echo "Start time: $(date)"
    echo "NAME: $NAME"
    echo "GENE_EMBEDDING_METHOD: $GENE_EMBEDDING_METHOD"
    echo "CANONICAL_ISOFORM_FILE: ${CANONICAL_ISOFORM_FILE:-NA}"
    echo "DATA_PATH: $DATA_PATH"
    echo "OUTPUT_BASE: $OUTPUT_BASE"
    echo "METHOD_OUTPUT_DIR: $METHOD_OUTPUT_DIR"
    echo "PYTHON_BIN: $PYTHON_BIN"
    echo "DEVICE: $DEVICE"
    echo "RUN_ESM: $RUN_ESM"
    echo "ALL_ISOFORM_EMBEDDING_DIR: ${ALL_ISOFORM_EMBEDDING_DIR:-AUTO}"
    echo "FULL_ESM_DIR: $FULL_ESM_DIR"
    echo "START_DELAY_HOURS: $START_DELAY_HOURS"

    if [[ "$START_DELAY_HOURS" != "0" ]]; then
        delay_seconds=$((START_DELAY_HOURS * 3600))
        echo "Delaying workflow for ${START_DELAY_HOURS} hour(s) before loading conda/CUDA."
        echo "Expected resume time: $(date -d "+${delay_seconds} seconds" 2>/dev/null || echo "after ${delay_seconds} seconds")"
        sleep "$delay_seconds"
        echo "Delay finished. Resume time: $(date)"
    fi

    source "$CONDA_PROFILE"
    conda activate "$CONDA_ENV"
    check_cuda
    prepare_fasta

    if [[ "$GENE_EMBEDDING_METHOD" == "aggregate" || "$GENE_EMBEDDING_METHOD" == "max_pooling" ]]; then
        echo "Preparing all-isoform FASTA and mapping..."
        if [[ ! -f "$FULL_CLEAN_FASTA" ]]; then
            "$PYTHON_BIN" -u "$COMMON_SCRIPT_PATH/clean_fasta.py" \
                --data_path="$FASTA_PATH" \
                --save_path="$FULL_CLEAN_FASTA"
        else
            echo "$FULL_CLEAN_FASTA already exists; skipping clean_fasta.py."
        fi

        if [[ ! -f "$GENE_TO_PROTEIN_MAP" ]]; then
            "$PYTHON_BIN" -u "$COMMON_SCRIPT_PATH/map_gene_symbol_to_protein_ids.py" \
                --fasta_path "$FASTA_PATH" \
                --save_path "$GENE_TO_PROTEIN_MAP"
        else
            echo "$GENE_TO_PROTEIN_MAP already exists; skipping gene-to-protein mapping."
        fi

        echo "Extracting all-isoform protein embeddings..."
        run_esm_if_needed "$FULL_CLEAN_FASTA" "$FULL_ESM_DIR" "all-isoform"

        echo "Converting all-isoform protein embeddings to gene embeddings using $GENE_EMBEDDING_METHOD..."
        "$PYTHON_BIN" -u "$COMMON_SCRIPT_PATH/convert_protein_embeddings_to_gene_embeddings.py" \
            --embedding-dir "$FULL_ESM_DIR" \
            --gene-symbol-to-protein-ids-path "$GENE_TO_PROTEIN_MAP" \
            --method "$GENE_EMBEDDING_METHOD" \
            --embedding-model ESM2 \
            --save-path "$GENE_EMBEDDING_PATH" \
            --save-report "$CONVERSION_REPORT"
    else
        echo "Preparing canonical isoforms..."
        if [[ -n "$CANONICAL_ISOFORM_FILE" ]]; then
            echo "Using user-provided canonical isoform mapping: $CANONICAL_ISOFORM_FILE"
            CANONICAL_MAPPING_FOR_CONVERSION="$METHOD_OUTPUT_DIR/$(basename "$CANONICAL_ISOFORM_FILE")"
            if [[ "$(realpath "$CANONICAL_ISOFORM_FILE")" != "$(realpath -m "$CANONICAL_MAPPING_FOR_CONVERSION")" ]]; then
                cp "$CANONICAL_ISOFORM_FILE" "$CANONICAL_MAPPING_FOR_CONVERSION"
            fi
            "$PYTHON_BIN" -u "$CANONICAL_METHOD_SCRIPT_PATH/build_canonical_fasta_from_mapping.py" \
                --fasta-path "$FASTA_PATH" \
                --canonical-map "$CANONICAL_MAPPING_FOR_CONVERSION" \
                --save-fasta "$CANONICAL_FASTA" \
                --save-report "$CANONICAL_FASTA_BUILD_REPORT"
        else
            echo "[WARN] --gene-embedding-method canonical_isoform was selected, but no --canonical-isoform-file was provided."
            echo "[WARN] The workflow will run APPRIS/UniProt/Ensembl canonical isoform discovery now."
            echo "[WARN] For rare or poorly annotated species, this step can be very slow because it may need web API fallbacks."

            ANNOTATION_GTF_ARGS=()
            if [[ -n "$ANNOTATION_GTF_URL" ]]; then
                if [[ ! -f "$ANNOTATION_GTF" ]]; then
                    echo "Downloading local annotation GTF with APPRIS tags: $ANNOTATION_GTF_URL"
                    wget -c "$ANNOTATION_GTF_URL" -O "$ANNOTATION_GTF"
                else
                    echo "$ANNOTATION_GTF already exists; skipping GTF download."
                fi
                ANNOTATION_GTF_ARGS=(--annotation-gtf "$ANNOTATION_GTF")
            fi

            if [[ ! -f "$CANONICAL_FASTA" || ! -f "$CANONICAL_MAPPING" || ! -f "$CANONICAL_SELECTION_TABLE" ]]; then
                "$PYTHON_BIN" -u "$CANONICAL_METHOD_SCRIPT_PATH/select_canonical_isoforms.py" \
                    --fasta-path "$FASTA_PATH" \
                    --name "$NAME" \
                    --appris-species "$APPRIS_SPECIES" \
                    --uniprot-organism "$UNIPROT_ORGANISM" \
                    --save-fasta "$CANONICAL_FASTA" \
                    --save-mapping "$CANONICAL_MAPPING" \
                    --save-selection-table "$CANONICAL_SELECTION_TABLE" \
                    "${ANNOTATION_GTF_ARGS[@]}" \
                    --cache "$CANONICAL_API_CACHE" \
                    --sleep 1 \
                    --timeout 45 \
                    --appris-batch-size "$APPRIS_BATCH_SIZE" \
                    --final-fallback ensembl_canonical_then_longest
            else
                echo "$CANONICAL_FASTA, $CANONICAL_MAPPING, and $CANONICAL_SELECTION_TABLE already exist; skipping canonical selection."
            fi
        fi

        echo "Cleaning canonical FASTA..."
        if [[ ! -f "$CANONICAL_CLEAN_FASTA" ]]; then
            "$PYTHON_BIN" -u "$COMMON_SCRIPT_PATH/clean_fasta.py" \
                --data_path="$CANONICAL_FASTA" \
                --save_path="$CANONICAL_CLEAN_FASTA"
        else
            echo "$CANONICAL_CLEAN_FASTA already exists; skipping clean_fasta.py."
        fi

        echo "Extracting canonical protein embeddings..."
        run_esm_if_needed "$CANONICAL_CLEAN_FASTA" "$CANONICAL_ESM_DIR" "canonical"

        echo "Converting canonical protein embeddings to gene embeddings..."
        "$PYTHON_BIN" -u "$COMMON_SCRIPT_PATH/convert_protein_embeddings_to_gene_embeddings.py" \
            --embedding-dir "$CANONICAL_ESM_DIR" \
            --gene-symbol-to-protein-ids-path "$CANONICAL_MAPPING_FOR_CONVERSION" \
            --method canonical_isoform \
            --embedding-model ESM2 \
            --save-path "$GENE_EMBEDDING_PATH" \
            --save-report "$CONVERSION_REPORT"
    fi

    echo "Gene embedding output: $(realpath "$GENE_EMBEDDING_PATH")"
    echo "Conversion report: $(realpath "$CONVERSION_REPORT")"
    echo "Log file: $log_file"
    echo "End time: $(date)"
} >> "$log_file" 2>&1
