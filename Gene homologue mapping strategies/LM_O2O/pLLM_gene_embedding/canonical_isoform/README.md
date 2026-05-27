# canonical_isoform gene embedding

This method uses one selected canonical protein isoform for each gene.

If `--canonical-isoform-file` is supplied, the mapping is used directly. If not, the workflow runs `select_canonical_isoforms.py`, which tries local GTF APPRIS tags, APPRIS API, UniProt default sequence, and Ensembl/longest fallback.

Run:

```bash
bash ../run_embedding_complete.sh \
  --name Homo_sapiens.GRCh38.pep.all \
  --gene-embedding-method canonical_isoform \
  --canonical-isoform-file /path/to/gene_to_canonical_protein.json
```
