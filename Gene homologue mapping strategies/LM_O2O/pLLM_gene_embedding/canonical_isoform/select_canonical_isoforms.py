#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Select one canonical protein isoform per gene for gene embedding generation.

Priority:
1. APPRIS PRINCIPAL isoform from APPRIS WebServices exporter API.
2. UniProt default/canonical sequence queried by gene symbol + organism scientific name.
3. Optional last-resort local fallback: Ensembl canonical transcript via Ensembl REST, then longest local protein.

Outputs:
- canonical FASTA: one selected protein sequence per gene.
- gene_symbol_to_protein_ID JSON: one selected sequence ID per gene, compatible with the SATURN embedding converter style.
- selection TSV: audit table describing how each gene was selected.
- API cache JSON: resumable APPRIS/UniProt/Ensembl requests.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry


APPRIS_EXPORTER_URL = "https://apprisws.bioinfo.cnio.es/rest/exporter/id/{species}/{gene_id}"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id/{gene_id}"

GENE_SYMBOL_TOKEN = "gene_symbol:"
GENE_ID_TOKEN = "gene:"
TRANSCRIPT_ID_TOKEN = "transcript:"

NAME_TO_APPRIS_SPECIES = {
    "Homo_sapiens.GRCh38.pep.all": "homo_sapiens",
    "Mus_musculus.GRCm39.pep.all": "mus_musculus",
    "Danio_rerio.GRCz11.pep.all": "danio_rerio",
    "Sus_scrofa.Sscrofa11.1.pep.all": "sus_scrofa",
    "Gorilla_gorilla.gorGor4.pep.all": "gorilla_gorilla",
    "Pan_troglodytes.Pan_tro_3.0.pep.all": "pan_troglodytes",
    "Macaca_mulatta.Mmul_10.pep.all": "macaca_mulatta",
    "Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all": "macaca_fascicularis",
    "Callithrix_jacchus.mCalJac1.pat.X.pep.all": "callithrix_jacchus",
    "Microcebus_murinus.Mmur_3.0.pep.all": "microcebus_murinus",
    "Xenopus_tropicalis.UCB_Xtro_10.0.pep.all": "xenopus_tropicalis",
}

NAME_TO_UNIPROT_ORGANISM = {
    "Homo_sapiens.GRCh38.pep.all": "Homo sapiens",
    "Mus_musculus.GRCm39.pep.all": "Mus musculus",
    "Danio_rerio.GRCz11.pep.all": "Danio rerio",
    "Sus_scrofa.Sscrofa11.1.pep.all": "Sus scrofa",
    "Gorilla_gorilla.gorGor4.pep.all": "Gorilla gorilla",
    "Pan_troglodytes.Pan_tro_3.0.pep.all": "Pan troglodytes",
    "Macaca_mulatta.Mmul_10.pep.all": "Macaca mulatta",
    "Macaca_fascicularis.Macaca_fascicularis_6.0.pep.all": "Macaca fascicularis",
    "Callithrix_jacchus.mCalJac1.pat.X.pep.all": "Callithrix jacchus",
    "Microcebus_murinus.Mmur_3.0.pep.all": "Microcebus murinus",
    "Xenopus_tropicalis.UCB_Xtro_10.0.pep.all": "Xenopus tropicalis",
}


@dataclass
class ProteinRecord:
    protein_id: str
    protein_id_base: str
    transcript_id: str
    transcript_id_base: str
    gene_id: str
    gene_id_base: str
    gene_symbol: str
    length: int
    sequence: str
    description: str


@dataclass
class Selection:
    gene_symbol: str
    gene_id: str
    selected_id: str
    selected_transcript_id: str
    selected_source: str
    selected_reason: str
    selected_length: int
    n_local_isoforms: int
    appris_reliability: str = ""
    uniprot_accession: str = ""
    warning: str = ""


def strip_version(identifier: str) -> str:
    """Remove Ensembl version suffix, e.g. ENST000... .3 -> ENST000..."""
    return str(identifier).split(".")[0]


def safe_get_token(description_tokens: List[str], prefix: str) -> str:
    """Extract a token such as gene_symbol:TP53 from a FASTA description."""
    matches = [token[len(prefix):] for token in description_tokens if token.startswith(prefix)]
    if not matches:
        return ""
    return matches[0]


def sanitize_fasta_id(text: str) -> str:
    """Keep FASTA IDs compatible with ESM output file names."""
    text = str(text).strip()
    text = re.sub(r"[^A-Za-z0-9_.:-]+", "_", text)
    return text.strip("_")


def build_session(retries: int = 3) -> requests.Session:
    """Build a robust requests Session with retry support."""
    session = requests.Session()
    retry = Retry(
        total=retries,
        connect=retries,
        read=retries,
        status=retries,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({"User-Agent": "canonical-isoform-embedding/1.0"})
    return session


def load_cache(cache_path: Path) -> Dict[str, object]:
    """Load resumable API cache."""
    if not cache_path.exists():
        return {"appris": {}, "uniprot_search": {}, "uniprot_fasta": {}, "ensembl": {}}
    try:
        with cache_path.open("r", encoding="utf-8") as handle:
            cache = json.load(handle)
    except Exception as exc:
        print(f"[WARN] Could not read cache {cache_path}: {exc}")
        cache = {}
    for key in ["appris", "uniprot_search", "uniprot_fasta", "ensembl"]:
        cache.setdefault(key, {})
    return cache


def save_cache(cache: Dict[str, object], cache_path: Path) -> None:
    """Save API cache atomically."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_path.with_suffix(cache_path.suffix + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(cache, handle, ensure_ascii=False, indent=2, sort_keys=True)
    tmp_path.replace(cache_path)


def parse_ensembl_pep_fasta(fasta_path: Path) -> Dict[str, List[ProteinRecord]]:
    """Parse Ensembl peptide FASTA and group protein isoforms by gene symbol."""
    gene_to_records: Dict[str, List[ProteinRecord]] = {}
    skipped_no_symbol = 0

    for seq in SeqIO.parse(str(fasta_path), "fasta"):
        tokens = seq.description.split()
        gene_symbol = safe_get_token(tokens, GENE_SYMBOL_TOKEN)
        gene_id = safe_get_token(tokens, GENE_ID_TOKEN)
        transcript_id = safe_get_token(tokens, TRANSCRIPT_ID_TOKEN)

        if not gene_symbol:
            skipped_no_symbol += 1
            continue

        record = ProteinRecord(
            protein_id=seq.id,
            protein_id_base=strip_version(seq.id),
            transcript_id=transcript_id,
            transcript_id_base=strip_version(transcript_id),
            gene_id=gene_id,
            gene_id_base=strip_version(gene_id),
            gene_symbol=gene_symbol,
            length=len(seq.seq),
            sequence=str(seq.seq),
            description=seq.description,
        )
        gene_to_records.setdefault(gene_symbol, []).append(record)

    print(f"[INFO] Parsed genes with symbols: {len(gene_to_records):,}")
    print(f"[INFO] FASTA records without gene_symbol skipped: {skipped_no_symbol:,}")
    return gene_to_records


def open_text_maybe_gzip(path: Path):
    """Open plain text or gzip-compressed files transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt", encoding="utf-8")


def parse_gtf_attributes(attributes: str) -> Dict[str, List[str]]:
    """Parse GTF attributes into key -> list(values), preserving repeated tag fields."""
    parsed: Dict[str, List[str]] = {}
    for match in re.finditer(r'(\S+) "([^"]*)";', attributes):
        key, value = match.group(1), match.group(2)
        parsed.setdefault(key, []).append(value)
    return parsed


def load_gtf_appris_principal_map(annotation_gtf: Optional[Path]) -> Tuple[Dict[str, dict], Dict[str, dict]]:
    """
    Load APPRIS principal isoforms from a local GENCODE/Ensembl GTF.

    GENCODE transcript records often contain tags such as appris_principal_1,
    appris_principal_2, etc. Parsing this once is much faster and more robust
    than querying APPRIS WebService for every gene.
    """
    if annotation_gtf is None:
        return {}, {}
    if not annotation_gtf.exists():
        print(f"[WARN] annotation GTF not found, ignoring: {annotation_gtf}")
        return {}, {}

    by_gene_id: Dict[str, dict] = {}
    by_gene_symbol: Dict[str, dict] = {}
    n_transcripts = 0
    n_appris = 0

    with open_text_maybe_gzip(annotation_gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            n_transcripts += 1
            attrs = parse_gtf_attributes(fields[8])
            tags = attrs.get("tag", [])
            appris_tags = [tag for tag in tags if tag.startswith("appris_principal_")]
            if not appris_tags:
                continue

            gene_id = attrs.get("gene_id", [""])[0]
            transcript_id = attrs.get("transcript_id", [""])[0]
            protein_id = attrs.get("protein_id", [""])[0]
            gene_name = attrs.get("gene_name", [""])[0]
            if not transcript_id:
                continue

            def tag_rank(tag: str) -> int:
                try:
                    return int(tag.rsplit("_", 1)[1])
                except ValueError:
                    return 99

            best_tag = sorted(appris_tags, key=tag_rank)[0]
            rank = tag_rank(best_tag)
            candidate = {
                "gene_id": gene_id,
                "gene_id_base": strip_version(gene_id),
                "gene_symbol": gene_name,
                "transcript_id": transcript_id,
                "transcript_id_base": strip_version(transcript_id),
                "protein_id": protein_id,
                "protein_id_base": strip_version(protein_id),
                "appris_tag": best_tag,
                "rank": rank,
            }
            n_appris += 1

            for key, target in [(candidate["gene_id_base"], by_gene_id), (gene_name, by_gene_symbol)]:
                if not key:
                    continue
                previous = target.get(key)
                if previous is None or (rank, transcript_id) < (previous["rank"], previous["transcript_id"]):
                    target[key] = candidate

    print(
        f"[INFO] Loaded local GTF APPRIS tags from {annotation_gtf}: "
        f"{len(by_gene_id):,} genes with APPRIS principal tags "
        f"({n_appris:,} transcript records among {n_transcripts:,} transcripts)."
    )
    return by_gene_id, by_gene_symbol


def find_record_by_gtf_candidate(records: List[ProteinRecord], candidate: dict) -> Optional[ProteinRecord]:
    """Match a GTF APPRIS candidate to local FASTA records by protein or transcript ID."""
    if not candidate:
        return None
    protein_base = candidate.get("protein_id_base", "")
    transcript_base = candidate.get("transcript_id_base", "")
    for record in records:
        if protein_base and record.protein_id_base == protein_base:
            return record
    return find_record_by_transcript(records, transcript_base)


def chunked(items: List[str], size: int) -> Iterable[List[str]]:
    """Yield fixed-size chunks from a list."""
    for i in range(0, len(items), size):
        yield items[i:i + size]


def prefetch_appris_for_genes(
    session: requests.Session,
    species: str,
    gene_ids: Iterable[str],
    cache: Dict[str, object],
    timeout: int,
    sleep_seconds: float,
    batch_size: int,
) -> None:
    """
    Batch-prefetch APPRIS annotations.

    APPRIS exporter accepts multiple Ensembl gene IDs separated by semicolons in the
    path. This avoids one HTTP request per gene and makes whole-proteome processing
    much faster while still keeping a polite sleep between batch requests.
    """
    if not species or batch_size <= 1:
        return

    appris_cache = cache["appris"]
    unique_gene_ids = sorted({strip_version(gene_id) for gene_id in gene_ids if gene_id})
    missing_gene_ids = [
        gene_id for gene_id in unique_gene_ids
        if f"{species}|{gene_id}" not in appris_cache
    ]

    if not missing_gene_ids:
        print("[INFO] APPRIS cache already covers all genes; skipping batch prefetch.")
        return

    n_batches = (len(missing_gene_ids) + batch_size - 1) // batch_size
    print(f"[INFO] Batch-prefetching APPRIS annotations for {len(missing_gene_ids):,} genes in {n_batches:,} batches of up to {batch_size}.")

    for batch_index, batch in enumerate(tqdm(list(chunked(missing_gene_ids, batch_size)), desc="Prefetching APPRIS batches"), start=1):
        query_id = ";".join(batch)
        url = APPRIS_EXPORTER_URL.format(species=species, gene_id=query_id)
        params = {"methods": "appris", "format": "json", "sc": "ensembl"}

        try:
            response = session.get(url, params=params, headers={"Accept": "application/json"}, timeout=timeout)
            time.sleep(sleep_seconds)

            if response.status_code >= 400:
                print(f"[WARN] APPRIS batch HTTP {response.status_code}; batch {batch_index}/{n_batches}. Falling back to per-gene for this batch later.")
                continue

            payload = response.json()
            if not isinstance(payload, list):
                print(f"[WARN] APPRIS batch returned non-list payload; batch {batch_index}/{n_batches}.")
                continue

            grouped: Dict[str, List[dict]] = {gene_id: [] for gene_id in batch}
            for item in payload:
                gene_id = strip_version(item.get("gene_id", "")) if isinstance(item, dict) else ""
                if gene_id in grouped:
                    grouped[gene_id].append(item)

            for gene_id, records in grouped.items():
                appris_cache[f"{species}|{gene_id}"] = records

        except Exception as exc:
            print(f"[WARN] APPRIS batch prefetch failed for batch {batch_index}/{n_batches}: {exc}")
            time.sleep(sleep_seconds)


def query_appris(session: requests.Session, species: str, gene_id: str, cache: Dict[str, object], timeout: int, sleep_seconds: float) -> List[dict]:
    """Query APPRIS exporter for a gene_id, cached."""
    gene_id_base = strip_version(gene_id)
    cache_key = f"{species}|{gene_id_base}"
    appris_cache = cache["appris"]
    if cache_key in appris_cache:
        return appris_cache[cache_key]

    url = APPRIS_EXPORTER_URL.format(species=species, gene_id=gene_id_base)
    params = {"methods": "appris", "format": "json", "sc": "ensembl"}
    try:
        response = session.get(url, params=params, headers={"Accept": "application/json"}, timeout=timeout)
        time.sleep(sleep_seconds)
        if response.status_code == 404:
            print(f"[WARN] APPRIS no record: {species} {gene_id_base}")
            appris_cache[cache_key] = []
            return []
        if response.status_code >= 400:
            print(f"[WARN] APPRIS HTTP {response.status_code}: {species} {gene_id_base}")
            appris_cache[cache_key] = []
            return []
        payload = response.json()
        if isinstance(payload, list):
            appris_cache[cache_key] = payload
            return payload
        appris_cache[cache_key] = []
        return []
    except Exception as exc:
        print(f"[WARN] APPRIS query failed for {species} {gene_id_base}: {exc}")
        time.sleep(sleep_seconds)
        appris_cache[cache_key] = []
        return []


def select_appris_principal_transcript(appris_records: List[dict]) -> Tuple[str, str, str]:
    """
    Return APPRIS principal transcript_id, reliability, reason.

    We intentionally require PRINCIPAL reliability or an appris_principal tag.
    APPRIS ALTERNATIVE/Possible Principal records are treated as non-canonical and will fall back.
    """
    principal_candidates = []
    for item in appris_records:
        if item.get("type") != "principal_isoform":
            continue
        transcript_id = item.get("transcript_id", "")
        if not transcript_id:
            continue
        reliability = str(item.get("reliability", ""))
        tag = str(item.get("tag", ""))
        annotation = str(item.get("annotation", ""))
        length_aa = int(item.get("length_aa") or 0)

        rank = 999
        if reliability.startswith("PRINCIPAL:"):
            try:
                rank = int(reliability.split(":", 1)[1])
            except ValueError:
                rank = 99
        elif "appris_principal" in tag:
            rank = 50
        elif annotation == "Principal Isoform":
            rank = 75
        else:
            continue

        principal_candidates.append((rank, -length_aa, transcript_id, reliability, tag or annotation))

    if not principal_candidates:
        return "", "", "no APPRIS PRINCIPAL isoform"

    principal_candidates.sort()
    _, _, transcript_id, reliability, reason = principal_candidates[0]
    return strip_version(transcript_id), reliability, reason


def query_uniprot_accession(session: requests.Session, gene_symbol: str, organism: str, cache: Dict[str, object], timeout: int, sleep_seconds: float) -> str:
    """Find a UniProt accession for gene_symbol + organism, preferring reviewed entries."""
    gene_symbol = str(gene_symbol).strip()
    organism = str(organism).strip()
    if not gene_symbol or not organism:
        return ""

    cache_key = f"{gene_symbol}|{organism}"
    uniprot_cache = cache["uniprot_search"]
    if cache_key in uniprot_cache:
        return uniprot_cache[cache_key]

    queries = [
        f'(gene_exact:"{gene_symbol}") AND (organism_name:"{organism}") AND reviewed:true',
        f'(gene_exact:"{gene_symbol}") AND (organism_name:"{organism}")',
        f'(gene:"{gene_symbol}") AND (organism_name:"{organism}") AND reviewed:true',
        f'(gene:"{gene_symbol}") AND (organism_name:"{organism}")',
    ]

    accession = ""
    for query in queries:
        try:
            response = session.get(
                UNIPROT_SEARCH_URL,
                params={"query": query, "format": "json", "size": 5},
                headers={"Accept": "application/json"},
                timeout=timeout,
            )
            time.sleep(sleep_seconds)
            if response.status_code >= 400:
                continue
            results = response.json().get("results", []) or []
            if results:
                accession = results[0].get("primaryAccession", "")
                break
        except Exception as exc:
            print(f"[WARN] UniProt search failed for {gene_symbol}/{organism}: {exc}")
            time.sleep(sleep_seconds)

    if not accession:
        print(f"[WARN] UniProt fallback not found: {gene_symbol}/{organism}")
    uniprot_cache[cache_key] = accession
    return accession


def fetch_uniprot_sequence(session: requests.Session, accession: str, cache: Dict[str, object], timeout: int, sleep_seconds: float) -> Tuple[str, str]:
    """Fetch UniProt canonical/default FASTA sequence for an accession."""
    accession = str(accession).strip()
    if not accession:
        return "", ""

    fasta_cache = cache["uniprot_fasta"]
    if accession in fasta_cache:
        entry = fasta_cache[accession]
        return entry.get("sequence", ""), entry.get("header", "")

    url = UNIPROT_FASTA_URL.format(accession=accession)
    try:
        response = session.get(url, headers={"Accept": "text/x-fasta,text/plain,*/*"}, timeout=timeout)
        time.sleep(sleep_seconds)
        if response.status_code >= 400:
            print(f"[WARN] UniProt FASTA HTTP {response.status_code}: {accession}")
            fasta_cache[accession] = {"sequence": "", "header": ""}
            return "", ""
        lines = [line.strip() for line in response.text.splitlines() if line.strip()]
        header = lines[0][1:] if lines and lines[0].startswith(">") else accession
        sequence = "".join(line for line in lines[1:] if not line.startswith(">"))
        fasta_cache[accession] = {"sequence": sequence, "header": header}
        return sequence, header
    except Exception as exc:
        print(f"[WARN] UniProt FASTA fetch failed for {accession}: {exc}")
        time.sleep(sleep_seconds)
        fasta_cache[accession] = {"sequence": "", "header": ""}
        return "", ""


def query_ensembl_canonical_transcript(session: requests.Session, gene_id: str, cache: Dict[str, object], timeout: int, sleep_seconds: float) -> str:
    """Optional last-resort fallback to Ensembl canonical transcript."""
    gene_id_base = strip_version(gene_id)
    ensembl_cache = cache["ensembl"]
    if gene_id_base in ensembl_cache:
        return ensembl_cache[gene_id_base]

    try:
        response = session.get(
            ENSEMBL_LOOKUP_URL.format(gene_id=gene_id_base),
            params={"expand": 1},
            headers={"Content-Type": "application/json"},
            timeout=timeout,
        )
        time.sleep(sleep_seconds)
        if response.status_code >= 400:
            ensembl_cache[gene_id_base] = ""
            return ""
        payload = response.json()
        canonical = strip_version(payload.get("canonical_transcript", ""))
        ensembl_cache[gene_id_base] = canonical
        return canonical
    except Exception as exc:
        print(f"[WARN] Ensembl canonical lookup failed for {gene_id_base}: {exc}")
        time.sleep(sleep_seconds)
        ensembl_cache[gene_id_base] = ""
        return ""


def find_record_by_transcript(records: List[ProteinRecord], transcript_base: str) -> Optional[ProteinRecord]:
    """Find local Ensembl protein record by transcript ID, ignoring version suffix."""
    if not transcript_base:
        return None
    for record in records:
        if record.transcript_id_base == transcript_base:
            return record
    return None


def choose_longest_local(records: List[ProteinRecord]) -> ProteinRecord:
    """Choose the longest local protein isoform as a transparent last resort."""
    return sorted(records, key=lambda r: (-r.length, r.protein_id))[0]


def select_for_gene(
    gene_symbol: str,
    records: List[ProteinRecord],
    session: requests.Session,
    appris_species: str,
    uniprot_organism: str,
    cache: Dict[str, object],
    timeout: int,
    sleep_seconds: float,
    final_fallback: str,
    gtf_appris_by_gene_id: Optional[Dict[str, dict]] = None,
    gtf_appris_by_gene_symbol: Optional[Dict[str, dict]] = None,
) -> Tuple[Selection, SeqRecord]:
    """Select canonical isoform sequence for one gene."""
    records = sorted(records, key=lambda r: r.protein_id)
    gene_id = records[0].gene_id

    # 0. Local GTF APPRIS PRINCIPAL tag. This is much faster than WebService calls.
    gtf_appris_by_gene_id = gtf_appris_by_gene_id or {}
    gtf_appris_by_gene_symbol = gtf_appris_by_gene_symbol or {}
    gtf_candidate = gtf_appris_by_gene_id.get(strip_version(gene_id)) or gtf_appris_by_gene_symbol.get(gene_symbol)
    gtf_record = find_record_by_gtf_candidate(records, gtf_candidate)
    if gtf_record is not None:
        selected = Selection(
            gene_symbol=gene_symbol,
            gene_id=gene_id,
            selected_id=gtf_record.protein_id,
            selected_transcript_id=gtf_record.transcript_id,
            selected_source="GTF_APPRIS_PRINCIPAL",
            selected_reason=f"Local annotation GTF {gtf_candidate.get('appris_tag', 'appris_principal')}",
            selected_length=gtf_record.length,
            n_local_isoforms=len(records),
            appris_reliability=gtf_candidate.get("appris_tag", ""),
        )
        seq = SeqRecord(Seq(gtf_record.sequence), id=gtf_record.protein_id, description=f"gene_symbol:{gene_symbol} source:GTF_APPRIS_PRINCIPAL transcript:{gtf_record.transcript_id}")
        return selected, seq

    # 1. APPRIS PRINCIPAL isoform from WebService, only for genes not covered by local GTF.
    appris_records = query_appris(session, appris_species, gene_id, cache, timeout, sleep_seconds) if appris_species else []
    appris_transcript, appris_reliability, appris_reason = select_appris_principal_transcript(appris_records)
    appris_record = find_record_by_transcript(records, appris_transcript)
    if appris_record is not None:
        selected = Selection(
            gene_symbol=gene_symbol,
            gene_id=gene_id,
            selected_id=appris_record.protein_id,
            selected_transcript_id=appris_record.transcript_id,
            selected_source="APPRIS_PRINCIPAL",
            selected_reason=f"APPRIS {appris_reliability or appris_reason}",
            selected_length=appris_record.length,
            n_local_isoforms=len(records),
            appris_reliability=appris_reliability,
        )
        seq = SeqRecord(Seq(appris_record.sequence), id=appris_record.protein_id, description=f"gene_symbol:{gene_symbol} source:APPRIS_PRINCIPAL transcript:{appris_record.transcript_id}")
        return selected, seq

    # 2. UniProt default/canonical sequence fallback.
    accession = query_uniprot_accession(session, gene_symbol, uniprot_organism, cache, timeout, sleep_seconds) if uniprot_organism else ""
    sequence, header = fetch_uniprot_sequence(session, accession, cache, timeout, sleep_seconds) if accession else ("", "")
    if sequence and "*" not in sequence:
        selected_id = sanitize_fasta_id(accession)
        selected = Selection(
            gene_symbol=gene_symbol,
            gene_id=gene_id,
            selected_id=selected_id,
            selected_transcript_id="",
            selected_source="UniProt_default_sequence",
            selected_reason=f"APPRIS unavailable ({appris_reason}); UniProt default/canonical FASTA used",
            selected_length=len(sequence),
            n_local_isoforms=len(records),
            appris_reliability=appris_reliability,
            uniprot_accession=accession,
        )
        seq = SeqRecord(Seq(sequence), id=selected_id, description=f"gene_symbol:{gene_symbol} source:UniProt_default accession:{accession} {header}")
        return selected, seq

    # 3. Optional final fallback to keep coverage explicit.
    if final_fallback == "ensembl_canonical_then_longest":
        canonical_tx = query_ensembl_canonical_transcript(session, gene_id, cache, timeout, sleep_seconds)
        canonical_record = find_record_by_transcript(records, canonical_tx)
        if canonical_record is not None:
            selected = Selection(
                gene_symbol=gene_symbol,
                gene_id=gene_id,
                selected_id=canonical_record.protein_id,
                selected_transcript_id=canonical_record.transcript_id,
                selected_source="Ensembl_canonical_last_resort",
                selected_reason="APPRIS and UniProt fallback unavailable; Ensembl canonical transcript used",
                selected_length=canonical_record.length,
                n_local_isoforms=len(records),
                warning="not APPRIS/UniProt; last-resort fallback",
            )
            seq = SeqRecord(Seq(canonical_record.sequence), id=canonical_record.protein_id, description=f"gene_symbol:{gene_symbol} source:Ensembl_canonical_last_resort transcript:{canonical_record.transcript_id}")
            return selected, seq

        print(
            f"[WARN] No canonical isoform source found for gene_symbol={gene_symbol}, gene_id={gene_id}; "
            f"APPRIS failed ({appris_reason}), UniProt failed, Ensembl canonical failed. "
            "Using longest local isoform as last-resort fallback."
        )
        longest = choose_longest_local(records)
        selected = Selection(
            gene_symbol=gene_symbol,
            gene_id=gene_id,
            selected_id=longest.protein_id,
            selected_transcript_id=longest.transcript_id,
            selected_source="local_longest_last_resort",
            selected_reason="APPRIS, UniProt, and Ensembl canonical unavailable; longest local isoform used",
            selected_length=longest.length,
            n_local_isoforms=len(records),
            warning="not APPRIS/UniProt; last-resort fallback",
        )
        seq = SeqRecord(Seq(longest.sequence), id=longest.protein_id, description=f"gene_symbol:{gene_symbol} source:local_longest_last_resort transcript:{longest.transcript_id}")
        return selected, seq

    print(
        f"[WARN] No canonical isoform source found for gene_symbol={gene_symbol}, gene_id={gene_id}; "
        f"APPRIS failed ({appris_reason}), UniProt failed, and final fallback is disabled. Skipping this gene."
    )
    selected = Selection(
        gene_symbol=gene_symbol,
        gene_id=gene_id,
        selected_id="",
        selected_transcript_id="",
        selected_source="skipped",
        selected_reason="APPRIS and UniProt fallback unavailable; final fallback disabled",
        selected_length=0,
        n_local_isoforms=len(records),
        warning="skipped",
    )
    return selected, SeqRecord(Seq(""), id=f"SKIPPED_{gene_symbol}", description="")


def infer_species_from_name(name: str, explicit_appris: str, explicit_organism: str) -> Tuple[str, str]:
    """Infer APPRIS species and UniProt organism names from Ensembl FASTA prefix."""
    appris_species = explicit_appris or NAME_TO_APPRIS_SPECIES.get(name, "")
    uniprot_organism = explicit_organism or NAME_TO_UNIPROT_ORGANISM.get(name, "")
    return appris_species, uniprot_organism


def main() -> None:
    parser = argparse.ArgumentParser(description="Select APPRIS/UniProt canonical isoform sequences for gene embeddings.")
    parser.add_argument("--fasta-path", required=True, type=Path, help="Original Ensembl pep.all.fa with full headers.")
    parser.add_argument("--name", default="", help="Dataset NAME prefix, e.g. Homo_sapiens.GRCh38.pep.all. Used to infer species.")
    parser.add_argument("--appris-species", default="", help="APPRIS species, e.g. homo_sapiens. Overrides --name inference.")
    parser.add_argument("--uniprot-organism", default="", help="UniProt organism scientific name. Overrides --name inference.")
    parser.add_argument("--save-fasta", required=True, type=Path, help="Output canonical isoform FASTA.")
    parser.add_argument("--save-mapping", required=True, type=Path, help="Output gene_symbol -> [selected_id] JSON.")
    parser.add_argument("--save-selection-table", required=True, type=Path, help="Output audit TSV table.")
    parser.add_argument("--annotation-gtf", default=None, type=Path, help="Optional local GENCODE/Ensembl GTF containing tag appris_principal_*.")
    parser.add_argument("--cache", required=True, type=Path, help="API cache JSON.")
    parser.add_argument("--timeout", type=int, default=45, help="HTTP timeout seconds.")
    parser.add_argument("--sleep", type=float, default=1.0, help="Sleep seconds after API requests.")
    parser.add_argument("--retries", type=int, default=3, help="HTTP retry count.")
    parser.add_argument("--cache-every", type=int, default=100, help="Save cache after this many genes.")
    parser.add_argument("--appris-batch-size", type=int, default=50, help="Batch size for APPRIS prefetch. Set <=1 to disable batch prefetch.")
    parser.add_argument(
        "--final-fallback",
        choices=["skip", "ensembl_canonical_then_longest"],
        default="ensembl_canonical_then_longest",
        help="What to do if both APPRIS and UniProt fallback fail.",
    )
    args = parser.parse_args()

    name = args.name or args.fasta_path.name.removesuffix(".fa")
    appris_species, uniprot_organism = infer_species_from_name(name, args.appris_species, args.uniprot_organism)
    print(f"[INFO] NAME={name}")
    print(f"[INFO] APPRIS species={appris_species or 'NA'}")
    print(f"[INFO] UniProt organism={uniprot_organism or 'NA'}")

    args.save_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.save_mapping.parent.mkdir(parents=True, exist_ok=True)
    args.save_selection_table.parent.mkdir(parents=True, exist_ok=True)

    gene_to_records = parse_ensembl_pep_fasta(args.fasta_path)
    cache = load_cache(args.cache)
    session = build_session(retries=args.retries)
    gtf_appris_by_gene_id, gtf_appris_by_gene_symbol = load_gtf_appris_principal_map(args.annotation_gtf)

    if appris_species and args.appris_batch_size > 1:
        gtf_covered_gene_ids = set(gtf_appris_by_gene_id)
        prefetch_gene_ids = [
            records[0].gene_id
            for records in gene_to_records.values()
            if records and strip_version(records[0].gene_id) not in gtf_covered_gene_ids
        ]
        prefetch_appris_for_genes(
            session=session,
            species=appris_species,
            gene_ids=prefetch_gene_ids,
            cache=cache,
            timeout=args.timeout,
            sleep_seconds=args.sleep,
            batch_size=args.appris_batch_size,
        )
        save_cache(cache, args.cache)

    selections: List[Selection] = []
    seq_records: List[SeqRecord] = []
    mapping: Dict[str, List[str]] = {}
    written_ids = set()

    for i, (gene_symbol, records) in enumerate(tqdm(sorted(gene_to_records.items()), desc="Selecting canonical isoforms"), start=1):
        selection, seq_record = select_for_gene(
            gene_symbol=gene_symbol,
            records=records,
            session=session,
            appris_species=appris_species,
            uniprot_organism=uniprot_organism,
            cache=cache,
            timeout=args.timeout,
            sleep_seconds=args.sleep,
            final_fallback=args.final_fallback,
            gtf_appris_by_gene_id=gtf_appris_by_gene_id,
            gtf_appris_by_gene_symbol=gtf_appris_by_gene_symbol,
        )
        selections.append(selection)
        if selection.selected_id:
            mapping[gene_symbol] = [selection.selected_id]
            if selection.selected_id not in written_ids and len(seq_record.seq) > 0:
                seq_records.append(seq_record)
                written_ids.add(selection.selected_id)

        if i % args.cache_every == 0:
            save_cache(cache, args.cache)

    save_cache(cache, args.cache)
    SeqIO.write(seq_records, args.save_fasta, "fasta")
    with args.save_mapping.open("w", encoding="utf-8") as handle:
        json.dump(mapping, handle, ensure_ascii=False, indent=2, sort_keys=True)
    pd.DataFrame([asdict(x) for x in selections]).to_csv(args.save_selection_table, sep="\t", index=False)

    selection_df = pd.DataFrame([asdict(x) for x in selections])
    print(f"[INFO] Wrote canonical FASTA: {args.save_fasta} ({len(seq_records):,} unique sequences)")
    print(f"[INFO] Wrote mapping JSON: {args.save_mapping} ({len(mapping):,} genes)")
    print(f"[INFO] Wrote selection audit table: {args.save_selection_table}")
    print("[INFO] Selection source counts:")
    print(selection_df["selected_source"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
