#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified converter from protein-level ESM embeddings to gene-level embeddings.

Supported gene embedding strategies:
1. aggregate: keep the original SATURN/Matching behavior, i.e. average all
   available isoform embeddings for each gene.
2. max_pooling: use all available isoform embeddings for each gene and take an
   element-wise maximum across isoforms.
3. canonical_isoform: use one selected canonical protein ID per gene.

Input mapping format:
- JSON: {"GENE": ["PROTEIN1", "PROTEIN2"]} or {"GENE": "PROTEIN1"}
- CSV/TSV: must contain a gene column and a protein column. Supported names:
  gene_symbol/gene/gene_name and protein_id/selected_id/canonical_protein_id.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import torch
from tqdm import tqdm

LAST_LAYER_BY_MODEL = {
    "ESM1b": 33,
    "MSA1b": 12,
    "ESM2": 48,
}

GENE_COLUMNS = ("gene_symbol", "gene", "gene_name", "Gene", "GeneSymbol")
PROTEIN_COLUMNS = (
    "protein_id",
    "protein_ids",
    "selected_id",
    "canonical_protein_id",
    "canonical_isoform",
    "ProteinID",
)


def normalize_protein_ids(value: object) -> List[str]:
    """Normalize JSON/CSV protein-id values into a clean list of strings."""
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        raw_values = list(value)
    else:
        text = str(value).strip()
        if not text or text.lower() in {"nan", "none", "null"}:
            return []
        # Allow comma/semicolon-separated values in CSV/TSV cells.
        raw_values = []
        for chunk in text.replace(";", ",").split(","):
            chunk = chunk.strip()
            if chunk:
                raw_values.append(chunk)
    return [str(item).strip() for item in raw_values if str(item).strip()]


def load_mapping(mapping_path: Path) -> Dict[str, List[str]]:
    """Load gene -> protein IDs mapping from JSON/CSV/TSV."""
    suffixes = "".join(mapping_path.suffixes).lower()
    if suffixes.endswith(".json"):
        with mapping_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        return {str(gene): normalize_protein_ids(ids) for gene, ids in payload.items()}

    delimiter = "\t" if suffixes.endswith(".tsv") or suffixes.endswith(".txt") else ","
    with mapping_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Mapping table has no header: {mapping_path}")
        gene_col = next((col for col in GENE_COLUMNS if col in reader.fieldnames), None)
        protein_col = next((col for col in PROTEIN_COLUMNS if col in reader.fieldnames), None)
        if gene_col is None or protein_col is None:
            raise ValueError(
                "CSV/TSV canonical mapping must contain a gene column "
                f"{GENE_COLUMNS} and a protein column {PROTEIN_COLUMNS}. "
                f"Observed columns: {reader.fieldnames}"
            )
        mapping: Dict[str, List[str]] = {}
        for row in reader:
            gene = str(row.get(gene_col, "")).strip()
            protein_ids = normalize_protein_ids(row.get(protein_col, ""))
            if gene and protein_ids:
                mapping[gene] = protein_ids
        return mapping


def strip_version(identifier: str) -> str:
    """Remove common Ensembl version suffix, e.g. ENSP000... .4 -> ENSP000..."""
    return str(identifier).split(".")[0]


def build_embedding_index(embedding_dir: Path) -> Dict[str, Path]:
    """Index .pt files by both full stem and version-stripped stem."""
    index: Dict[str, Path] = {}
    for path in embedding_dir.glob("*.pt"):
        index.setdefault(path.stem, path)
        index.setdefault(strip_version(path.stem), path)
    return index


def resolve_embedding_path(protein_id: str, embedding_dir: Path, embedding_index: Mapping[str, Path]) -> Optional[Path]:
    """Find one protein embedding file robustly."""
    protein_id = str(protein_id).strip()
    if not protein_id:
        return None
    exact = embedding_dir / f"{protein_id}.pt"
    if exact.exists():
        return exact
    base = strip_version(protein_id)
    base_path = embedding_dir / f"{base}.pt"
    if base_path.exists():
        return base_path
    return embedding_index.get(protein_id) or embedding_index.get(base)


def load_embedding(path: Path, last_layer: int) -> torch.Tensor:
    """Load ESM mean representation from one .pt file."""
    obj = torch.load(path, map_location="cpu")
    mean_representations = obj.get("mean_representations")
    if not isinstance(mean_representations, dict):
        raise KeyError(f"No mean_representations dict found in {path}")
    if last_layer in mean_representations:
        embedding = mean_representations[last_layer]
    else:
        available_layers = sorted(mean_representations)
        if not available_layers:
            raise KeyError(f"No representation layers found in {path}")
        fallback_layer = available_layers[-1]
        print(f"[WARN] Layer {last_layer} not found in {path.name}; using layer {fallback_layer} instead.")
        embedding = mean_representations[fallback_layer]
    return embedding.detach().cpu().float()


def combine_embeddings(embeddings: Sequence[torch.Tensor], method: str) -> torch.Tensor:
    """Combine isoform embeddings into one gene embedding."""
    stacked = torch.stack(list(embeddings), dim=0)
    if method == "aggregate":
        # Original aggregate strategy: mean over all available isoform embeddings.
        return torch.mean(stacked, dim=0)
    if method == "max_pooling":
        # Reviewer-requested strategy: element-wise maximum over isoform embeddings.
        return torch.max(stacked, dim=0).values
    if method == "canonical_isoform":
        # Canonical mode passes one selected protein ID per gene.
        return stacked[0]
    raise ValueError(f"Unsupported method: {method}")


def iter_selected_protein_ids(method: str, protein_ids: List[str]) -> List[str]:
    """Decide which protein IDs should be loaded for one gene."""
    if method == "canonical_isoform":
        return protein_ids[:1]
    return protein_ids


def write_report(report_path: Path, rows: Iterable[dict]) -> None:
    """Write a TSV audit report."""
    rows = list(rows)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "gene_symbol",
        "method",
        "n_protein_ids_in_mapping",
        "n_selected_protein_ids",
        "n_embeddings_found",
        "selected_protein_ids",
        "loaded_protein_ids",
        "missing_protein_ids",
        "status",
    ]
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def convert(args: argparse.Namespace) -> None:
    """Main conversion routine."""
    gene_to_protein_ids = load_mapping(args.gene_symbol_to_protein_ids_path)
    last_layer = LAST_LAYER_BY_MODEL[args.embedding_model]
    embedding_index = build_embedding_index(args.embedding_dir)

    gene_to_embedding: Dict[str, torch.Tensor] = {}
    report_rows: List[dict] = []

    print(f"[INFO] Method: {args.method}")
    print(f"[INFO] Genes in mapping: {len(gene_to_protein_ids):,}")
    print(f"[INFO] Indexed protein embeddings: {len(set(embedding_index.values())):,}")

    for gene_symbol, protein_ids in tqdm(sorted(gene_to_protein_ids.items()), desc="Converting gene embeddings"):
        protein_ids = normalize_protein_ids(protein_ids)
        selected_ids = iter_selected_protein_ids(args.method, protein_ids)
        embeddings: List[torch.Tensor] = []
        loaded_ids: List[str] = []
        missing_ids: List[str] = []

        for protein_id in selected_ids:
            embedding_path = resolve_embedding_path(protein_id, args.embedding_dir, embedding_index)
            if embedding_path is None:
                missing_ids.append(protein_id)
                continue
            try:
                embeddings.append(load_embedding(embedding_path, last_layer))
                loaded_ids.append(protein_id)
            except Exception as exc:
                print(f"[WARN] Failed to load embedding for {gene_symbol}/{protein_id}: {exc}")
                missing_ids.append(protein_id)

        if embeddings:
            gene_to_embedding[gene_symbol] = combine_embeddings(embeddings, args.method)
            status = "loaded"
        else:
            status = "missing_all_embeddings"

        report_rows.append(
            {
                "gene_symbol": gene_symbol,
                "method": args.method,
                "n_protein_ids_in_mapping": len(protein_ids),
                "n_selected_protein_ids": len(selected_ids),
                "n_embeddings_found": len(embeddings),
                "selected_protein_ids": ";".join(selected_ids),
                "loaded_protein_ids": ";".join(loaded_ids),
                "missing_protein_ids": ";".join(missing_ids),
                "status": status,
            }
        )

    args.save_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(gene_to_embedding, args.save_path)
    print(f"[INFO] Saved gene embeddings: {args.save_path}")
    print(f"[INFO] Genes with embeddings: {len(gene_to_embedding):,} / {len(gene_to_protein_ids):,}")

    if args.save_report:
        write_report(args.save_report, report_rows)
        print(f"[INFO] Saved conversion report: {args.save_report}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Unified protein-to-gene embedding converter.")
    parser.add_argument("--embedding-dir", required=True, type=Path, help="Directory containing protein-level ESM .pt files.")
    parser.add_argument("--gene-symbol-to-protein-ids-path", required=True, type=Path, help="Gene -> protein IDs mapping file.")
    parser.add_argument("--method", choices=["aggregate", "max_pooling", "canonical_isoform"], default="aggregate", help="Gene embedding strategy.")
    parser.add_argument("--embedding-model", choices=sorted(LAST_LAYER_BY_MODEL), default="ESM2", help="Embedding model used by ESM extraction.")
    parser.add_argument("--save-path", required=True, type=Path, help="Output .pt file: gene_symbol -> embedding tensor.")
    parser.add_argument("--save-report", default=None, type=Path, help="Optional TSV report path.")
    return parser.parse_args()


if __name__ == "__main__":
    convert(parse_args())
