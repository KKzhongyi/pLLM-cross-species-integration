#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a canonical-isoform FASTA from an original Ensembl pep.all.fa and a gene -> protein mapping."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

GENE_COLUMNS = ("gene_symbol", "gene", "gene_name", "Gene", "GeneSymbol")
PROTEIN_COLUMNS = (
    "protein_id",
    "protein_ids",
    "selected_id",
    "canonical_protein_id",
    "canonical_isoform",
    "ProteinID",
)


def strip_version(identifier: str) -> str:
    return str(identifier).split(".")[0]


def normalize_protein_ids(value: object) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        values = list(value)
    else:
        text = str(value).strip()
        if not text or text.lower() in {"nan", "none", "null"}:
            return []
        values = [x.strip() for x in text.replace(";", ",").split(",") if x.strip()]
    return [str(x).strip() for x in values if str(x).strip()]


def load_mapping(mapping_path: Path) -> Dict[str, List[str]]:
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
                "Canonical mapping table must contain a gene column "
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


def parse_fasta(path: Path) -> Dict[str, Tuple[str, str]]:
    """Return protein_id -> (sequence, original_description)."""
    records: Dict[str, Tuple[str, str]] = {}
    current_id = ""
    current_desc = ""
    current_seq: List[str] = []

    def flush() -> None:
        if not current_id:
            return
        seq = "".join(current_seq)
        records[current_id] = (seq, current_desc)
        records.setdefault(strip_version(current_id), (seq, current_desc))

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                current_desc = line[1:].strip()
                current_id = current_desc.split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        flush()
    return records


def wrap_sequence(sequence: str, width: int = 80) -> Iterable[str]:
    for i in range(0, len(sequence), width):
        yield sequence[i : i + width]


def write_report(path: Path, rows: Iterable[dict]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["gene_symbol", "canonical_protein_id", "status", "note"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build canonical FASTA from a gene -> canonical protein mapping.")
    parser.add_argument("--fasta-path", required=True, type=Path, help="Original Ensembl pep.all.fa.")
    parser.add_argument("--canonical-map", required=True, type=Path, help="JSON/CSV/TSV mapping gene -> canonical protein ID.")
    parser.add_argument("--save-fasta", required=True, type=Path, help="Output canonical FASTA.")
    parser.add_argument("--save-report", default=None, type=Path, help="Optional TSV report.")
    args = parser.parse_args()

    gene_to_ids = load_mapping(args.canonical_map)
    fasta_records = parse_fasta(args.fasta_path)
    args.save_fasta.parent.mkdir(parents=True, exist_ok=True)

    rows: List[dict] = []
    written_ids = set()
    n_written = 0

    with args.save_fasta.open("w", encoding="utf-8") as handle:
        for gene_symbol, protein_ids in sorted(gene_to_ids.items()):
            canonical_id = protein_ids[0] if protein_ids else ""
            record = fasta_records.get(canonical_id) or fasta_records.get(strip_version(canonical_id))
            if not canonical_id:
                rows.append({"gene_symbol": gene_symbol, "canonical_protein_id": "", "status": "missing_id", "note": "empty canonical protein ID"})
                continue
            if record is None:
                rows.append({
                    "gene_symbol": gene_symbol,
                    "canonical_protein_id": canonical_id,
                    "status": "missing_in_fasta",
                    "note": "protein ID was not found in the original Ensembl FASTA; if this is a UniProt-only sequence, rerun without --canonical-isoform-file so the canonical selector can fetch the sequence",
                })
                continue
            sequence, _description = record
            if "*" in sequence:
                rows.append({"gene_symbol": gene_symbol, "canonical_protein_id": canonical_id, "status": "skipped_stop_codon", "note": "sequence contains stop codon"})
                continue
            if canonical_id not in written_ids:
                handle.write(f">{canonical_id} gene_symbol:{gene_symbol} source:user_canonical_mapping\n")
                for chunk in wrap_sequence(sequence):
                    handle.write(chunk + "\n")
                written_ids.add(canonical_id)
                n_written += 1
            rows.append({"gene_symbol": gene_symbol, "canonical_protein_id": canonical_id, "status": "written", "note": ""})

    print(f"[INFO] Genes in canonical mapping: {len(gene_to_ids):,}")
    print(f"[INFO] Unique canonical FASTA records written: {n_written:,}")
    print(f"[INFO] Saved canonical FASTA: {args.save_fasta}")

    if args.save_report:
        write_report(args.save_report, rows)
        print(f"[INFO] Saved canonical FASTA build report: {args.save_report}")

    missing = sum(1 for row in rows if row["status"] != "written")
    if missing:
        print(f"[WARN] {missing:,} gene(s) could not be written from the provided mapping. Check the report for details.")


if __name__ == "__main__":
    main()
