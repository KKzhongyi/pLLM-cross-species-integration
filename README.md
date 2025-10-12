# Cross-Species Integration Pipeline

This repository provides a comprehensive workflow for generating **cross-species gene homologue mappings** and performing **cross-species dataset integration**.  
The workflow includes five distinct homology mapping strategies and a complete end-to-end pipeline for constructing dataset-level homolog tables, integrating data, and evaluating integration performance.

---

## 📘 Table of Contents
1. [Overview](#overview)
2. [Gene Homologue Mapping Strategies](#gene-homologue-mapping-strategies)
   - [1. ENS_M2M](#1-ens_m2m)
   - [2. ENS_O2O](#2-ens_o2o)
   - [3. HM_O2O](#3-hm_o2o)
   - [4. LM_O2O](#4-lm_o2o)
   - [5. HL_O2O](#5-hl_o2o)
3. [Cross-Species Integration Workflow](#cross-species-integration-workflow)
   - [Step 1: Generate Dataset-Level Homologue Table](#step-1-generate-dataset-level-homologue-table)
   - [Step 2: Perform Cross-Species Integration](#step-2-perform-cross-species-integration)
   - [Step 3: Evaluate Integration Performance](#step-3-evaluate-integration-performance)
4. [Example: CAA Dataset Integration](#example-caa-dataset-integration)
5. [Directory Structure](#directory-structure)

---

## Overview

Cross-species integration enables joint analysis of multi-species transcriptomic datasets by establishing reliable **gene homologue mappings**.  
This repository compares and implements **five strategies** to derive such mappings, integrating orthology-based and protein language model–based evidence.  
The resulting mappings are then used to unify gene names across species and evaluate integration quality via dimensionality reduction and quantitative metrics.

---

## Gene Homologue Mapping Strategies

### 1. ENS_M2M

- **Description:**  
  The ENS_M2M strategy directly downloads homology mappings from **Ensembl BioMart**.
- **Procedure:**
  1. On the BioMart website, under **Dataset**, select species *A*.
  2. Under **Attributes**, choose **Homologues (Max select 6 orthologues)**.
  3. In the “GENE” section, select:
     - `Gene name`
  4. In the **ORTHOLOGUES (A~Z)** section, select target species *B*, checking:
     - `B gene name`
     - `B homology type`
     - `%id. target B gene identical to query gene`
     - `%id. query gene identical to target B gene`
     - `B orthology confidence [0 low, 1 high]`
- **Output Directory:**  
  ```
  /Gene homologue mapping strategies/ENS_M2M/
  ```

---

### 2. ENS_O2O

- **Description:**  
  Derived from **ENS_M2M**, this strategy retains only **one-to-one** orthologues.
- **Procedure:**
  - Filter ENS_M2M results where `"homology type" == "ortholog_one2one"`.
  - Implemented via:  
    ```
    /Gene homologue mapping strategies/ENS_O2O/Generating ENS_O2O.ipynb
    ```
- **Output Directory:**  
  ```
  /Gene homologue mapping strategies/ENS_O2O/o2oResults/
  ```

---

### 3. HM_O2O

- **Description:**  
  A deterministic one-to-one mapping method that prioritizes **Ensembl orthology confidence** and **sequence identity**.
- **Procedure:**
  1. Extract the following Ensembl attributes:
     - `%id. target gene identical to query gene`
     - `%id. query gene identical to target gene`
     - `orthology confidence (0=low, 1=high)`
  2. Compute `identical_scores` as the mean of the two % identity attributes.
  3. Apply **global greedy selection**:
     - Priority 1: pairs with `orthology confidence = 1`
     - Priority 2: within same confidence, sort by `identical_scores (desc)`
     - Iteratively accept top pairs, removing all conflicts until resolved.
- **Output Directory:**  
  ```
  /Gene homologue mapping strategies/HM_O2O/o2oResults/
  ```

---

### 4. LM_O2O

- **Description:**  
  Uses **protein language models (ESM2)** to infer cross-species gene mappings via sequence embeddings.
- **Procedure:**
  1. Download proteome FASTA files from **Ensembl FTP**.
  2. Embed all protein sequences using **ESM2 (esm2_t48_15B_UR50D)** → 5120-dimensional embeddings.
  3. Map protein IDs to gene symbols; average protein embeddings to obtain **gene-level embeddings**.  
     See:  
     ```
     /Gene homologue mapping strategies/LM_O2O/pLLM_gene_embedding/run_embedding_complete.sh
     ```
  4. Compute **cross-species correlation matrix** between gene embeddings.
  5. Identify **double best-hit (DBH)** pairs (highest mutual correlation).
  6. Apply greedy one-to-one selection by descending correlation.
- **Key Scripts:**
  ```
  /Gene homologue mapping strategies/LM_O2O/LM_O2O gene homologue mapping/run_LM_O2O_Matching_2species.sh
  /Gene homologue mapping strategies/LM_O2O/LM_O2O gene homologue mapping/generatehomologue_forallspecies_LM_O2O.ipynb
  ```
- **Output Directory:**  
  ```
  /Gene homologue mapping strategies/LM_O2O/o2oResults/
  ```

---

### 5. HL_O2O

- **Description:**  
  Integrates **homology-based (HM_O2O)** and **language model–based (LM_O2O)** evidence for a unified one-to-one mapping.
- **Procedure:**
  1. Merge pairs from HM_O2O and LM_O2O.
  2. Collect attributes:
     - (i) `orthology confidence` (0/1)
     - (ii) `identical_scores` (Ensembl % identity mean)
     - (iii) `correlation` (from LM_O2O)
  3. Handle missing values:
     - Missing correlation → median of species pair
     - Missing confidence → 0
     - Missing identical_scores → median of candidate set
  4. Apply **min–max normalization** on continuous features.
  5. Compute combined score = normalized(identity) + normalized(correlation).
  6. Apply **global greedy one-to-one selection**:
     - Confidence = 1 prioritized over 0
     - Within each tier, sort by combined score (desc)
- **Script:**  
  ```
  /Gene homologue mapping strategies/HL_O2O/generating_HL_O2O_homologue_mapping.py
  ```
- **Output Directory:**  
  ```
  /Gene homologue mapping strategies/HL_O2O/o2oResults/
  ```

---

## Cross-Species Integration Workflow

The integration pipeline consists of **three major stages**:

### Step 1: Generate Dataset-Level Homologue Table

- Combine all pairwise species mappings within a dataset.
- Use strategy-specific notebooks under:
  ```
  /Cross_species integration workflow/Generate Homologue Table/(Gene homologue mapping strategy)/extract_homologue.ipynb
  ```
- Example output:
  ```
  /Cross_species integration workflow/Generate Homologue Table/HL_O2O/crab_eating_macaque_rhesus_macaque_human_mice_pig.csv
  ```

---

### Step 2: Perform Cross-Species Integration

- Use the generated dataset-level Homologue Table as input.
- Run integration using:
  ```
  /Cross_species integration workflow/Scripts/run_UMAP.sh
  /Cross_species integration workflow/Scripts/extract.py
  ```
- Example output:
  ```
  crab_eating_macaque_rhesus_macaque_human_mice_pig.h5
  ```

---

### Step 3: Evaluate Integration Performance

- Run UMAP and compute integration metrics:
  ```
  /Cross_species integration workflow/Scripts/UMAPgeneration.R
  /Cross_species integration workflow/Scripts/ASW.R
  ```
- Example outputs:
  ```
  crab_eating_macaque_rhesus_macaque_human_mice_pig.qs
  crab_eating_macaque_rhesus_macaque_human_mice_pigprocessed_obj_CCAIntegration.qs
  cell_type_cluster_metrics_with_CCAIntegration.txt
  ```

---

## Example: CAA Dataset Integration

Using the **Caa (Cell Atlas of Aqueous Humor)** dataset as an example:

1. Generate `Homologue Table` →  
   `/Generate Homologue Table/HL_O2O/crab_eating_macaque_rhesus_macaque_human_mice_pig.csv`
2. Integrate dataset →  
   `/Scripts/extract.py` produces `.h5` integrated file
3. Evaluate and visualize →  
   `UMAPgeneration.R` and `ASW.R` produce `.qs` files and integration metrics

This completes the **cross-species integration and evaluation** process.

---

## Directory Structure

```
├── Gene homologue mapping strategies/
│   ├── ENS_M2M/
│   ├── ENS_O2O/
│   ├── HM_O2O/
│   ├── LM_O2O/
│   └── HL_O2O/
│
└── Cross_species integration workflow/
    ├── Generate Homologue Table/
    └── Scripts/
```

---

## Citation

If you use this workflow or any derived gene homologue mapping in your work, please cite the corresponding reference or repository link.

---

© 2025 Cross-Species Integration Project
