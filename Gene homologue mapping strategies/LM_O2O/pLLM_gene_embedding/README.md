# pLLM Gene Embedding Workflow

This directory generates gene-level ESM2 embeddings for the LM_O2O mapping strategy. `run_embedding_complete.sh` keeps the original entrypoint name and adds reviewer-requested alternatives while keeping aggregate as the default.

## Gene Embedding Methods

- `aggregate`: mean pooling over all available protein isoform embeddings for each gene. This reproduces the previous LM_O2O behavior.
- `max_pooling`: element-wise max pooling over all available protein isoform embeddings for each gene.
- `canonical_isoform`: one selected canonical protein isoform per gene.

## Layout

```text
pLLM_gene_embedding/
|-- run_embedding_complete.sh
|-- clean_fasta.py
|-- map_gene_symbol_to_protein_ids.py
|-- convert_protein_embeddings_to_gene_embeddings.py
|-- canonical_isoform/
|   |-- select_canonical_isoforms.py
|   `-- build_canonical_fasta_from_mapping.py
```

## Usage

Run from this directory and override paths as needed:

```bash
bash run_embedding_complete.sh \
  --name Homo_sapiens.GRCh38.pep.all \
  --gene-embedding-method aggregate
```

```bash
bash run_embedding_complete.sh \
  --name Homo_sapiens.GRCh38.pep.all \
  --gene-embedding-method max_pooling
```

```bash
bash run_embedding_complete.sh \
  --name Homo_sapiens.GRCh38.pep.all \
  --gene-embedding-method canonical_isoform \
  --canonical-isoform-file /path/to/gene_to_canonical_protein.json
```

If `--canonical-isoform-file` is omitted, `canonical_isoform/select_canonical_isoforms.py` selects one protein per gene using local GTF APPRIS tags when available, then APPRIS WebService, UniProt default sequences, and a final Ensembl canonical/longest-local fallback.

Outputs are written to:

```text
pLLM_gene_embedding/Results/<NAME>/<method>/
```

The converter also writes a TSV audit report for the protein IDs loaded, missing, and used for each gene.

## Use Embeddings In LM_O2O Matching

After generating gene embeddings, run the method-aware LM_O2O wrapper from `../LM_O2O gene homologue mapping/`:

```bash
bash "../LM_O2O gene homologue mapping/run_LM_O2O_Matching_2species.sh" \
  --name-1 Homo_sapiens.GRCh38.pep.all \
  --name-2 Mus_musculus.GRCm39.pep.all \
  --gene-embedding-method max_pooling
```

The wrapper resolves the selected `.pt` files under `pLLM_gene_embedding/Results/<NAME>/<method>/` and calls the original `correlation.py` workflow.
