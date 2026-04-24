# Mouse Gene Knowledge Graph

A knowledge graph for mouse (*Mus musculus*) gene biology, built with [BioCypher](https://biocypher.org/).  
Designed to support gene–phenotype link prediction using graph-based machine learning.

Mirrors the architecture of [@TonyGoncharov/phenotype_prediction_project](https://github.com/TonyGoncharov/phenotype_prediction_project) (human graph), extended with expression and protein interaction layers.

---

## Graph structure

| Node | ID | Count |
|---|---|---|
| `MouseGene` | MGI ID | 47 757 |
| `Tissue` | UBERON ID | 90 |
| `Protein` | BioGRID ID | 18 904 |

| Edge | Source | Count |
|---|---|---|
| `EXPRESSED_IN` | GXD / MGI RNA-seq | ~2 028 736 |
| `INTERACTS_WITH` | BioGRID v5.0 | 95 022 |

---

## Data sources

| Layer | Source | Notes |
|---|---|---|
| Expression | [GXD / MGI](https://www.informatics.jax.org/downloads/reports/gxdrnaseq/) | 96 RNA-seq experiments, median TPM per gene × tissue |
| Tissue ontology | [UBERON](https://obofoundry.org/ontology/uberon.html) | mapped from GXD anatomical terms |
| PPI | [BioGRID](https://thebiogrid.org/) v5.0.256, TAB3 | physical interactions only |

---

## Expression data pipeline

Raw RNA-seq data was taken from the GXD RNA-seq download folder (`/gxdrnaseq/`), which contains experiments curated and pre-processed by MGI from EMBL-EBI Expression Atlas. GXD applies controlled biological source metadata and quantile normalisation across biological replicates within each experiment, producing `avg_qnTPM` values per gene per sample set.

For this graph we used `avg_TPM` (column 17) — the average TPM across technical replicates within a sample — as the base value, and aggregated across experiments as follows:

1. **Download** — all 96 `.rpt.gz` files from the GXD RNA-seq folder (`scripts/gxd_download_aggregate.py`)
2. **Filter** — keep only rows where `Detected = Yes` and no mutant allele is present (wild-type only)
3. **Aggregate** — group by `(MGI Gene ID, Anatomical Structure)` across all experiments, compute **median `avg_TPM`**
4. **Level assignment** — assign expression level label based on median TPM:

   | Level | TPM range |
   |---|---|
   | Absent | 0 – 0.5 |
   | Low | 0.5 – 10 |
   | Medium | 10 – 100 |
   | High | > 100 |

5. **Tissue mapping** — GXD anatomical structure names mapped to UBERON IDs via OLS API + manual curation for 91 tissues (`scripts/make_uberon_mapping.py`)

Result: **2 043 209** unique gene × tissue pairs from **47 757** genes across **91** tissues.

---

## Project structure

```
├── config/
│   ├── biocypher_config.yaml
│   └── schema_config_mouse.yaml
├── data/                          # input files (not in git)
│   ├── gene_expression_summary.tsv
│   ├── gxd_uberon_mapping.tsv
│   └── biogrid_mouse_ppi.tsv
├── src/
│   ├── build_graph.py
│   └── adapters/
│       ├── base.py
│       ├── expression_adapter.py
│       └── ppi_adapter.py
└── scripts/
    ├── gxd_download_aggregate.py  # download GXD files → gene_expression_summary.tsv
    └── make_uberon_mapping.py     # GXD tissue names → UBERON IDs
```

---

## Quickstart

```bash
# 1. install dependencies
pip install biocypher pandas requests beautifulsoup4

# 2. download and aggregate GXD expression data (~96 experiments)
python scripts/gxd_download_aggregate.py
# → data/gene_expression_summary.tsv

# 3. map tissue names to UBERON IDs
python scripts/make_uberon_mapping.py
# → data/gxd_uberon_mapping.tsv

# 4. download BioGRID mouse PPI
#    https://downloads.thebiogrid.org/BioGRID/Latest-Release/
#    unzip BIOGRID-ORGANISM-LATEST.tab3.zip → rename to data/biogrid_mouse_ppi.tsv

# 5. build graph
python src/build_graph.py
# → biocypher_out/mouse/

# 6. import into Neo4j
bash biocypher_out/mouse/neo4j-admin-import-call.sh
```

---

## Planned

- Phenotype layer: `MouseGene` → `MpTopTerm` (MGI genotype–phenotype annotations)
- GO layer: `MouseGene` → `GoTerm`
- Gene–phenotype link prediction (PyKEEN / GNN)
