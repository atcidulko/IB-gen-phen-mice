# Mouse Gene Knowledge Graph

Pipeline for building a mouse gene–phenotype–expression–PPI knowledge graph
using [BioCypher](https://biocypher.org/), mirroring the architecture of
[TonyGoncharov/phenotype_prediction_project](https://github.com/TonyGoncharov/phenotype_prediction_project)
but extended with expression and PPI layers, and using MGI IDs as primary
gene identifiers.

## Graph structure

```
[MouseGene]  ──(has_mp_top_term)──▶  [MpTopTerm]   # phenotype
[MouseGene]  ──(has_go_term)───────▶  [GoTerm]      # function
[MouseGene]  ──(expressed_in)──────▶  [Tissue]      # expression (GXD)
[Protein]    ──(interacts_with)────▶  [Protein]     # PPI (BioGRID)
```

Tissue nodes carry UBERON IDs — enabling cross-species expression
comparison with a human graph sharing the same tissue node IDs.

MouseGene nodes use `MGI:xxxxxxx` as primary ID (stable across symbol
renames); gene symbol is stored as a property.

## Data sources

| File | Source | Description |
|---|---|---|
| `MGI_PhenoGenoMP.rpt` | [MGI](https://www.informatics.jax.org/) | Mouse genotype → MP annotations |
| `mp.obo` | [OBO Foundry](https://purl.obolibrary.org/obo/mp.obo) | Mammalian Phenotype Ontology |
| `go-basic.obo` | [Gene Ontology](https://purl.obolibrary.org/obo/go/go-basic.obo) | Gene Ontology |
| `gene2go.gz` | [NCBI](https://ftp.ncbi.nlm.nih.gov/gene/DATA/) | Gene → GO term mappings |
| `gene_info.gz` | [NCBI](https://ftp.ncbi.nlm.nih.gov/gene/DATA/) | Gene symbols and MGI cross-refs |
| `gene_expression_summary.tsv` | [GXD/MGI](https://www.informatics.jax.org/downloads/reports/gxdrnaseq/) | Aggregated RNA-seq (produced by `scripts/gxd_download_aggregate.py`) |
| `gxd_uberon_mapping.tsv` | manual + OLS | GXD tissue name → UBERON ID |
| `biogrid_mouse_ppi.tsv` | [BioGRID](https://thebiogrid.org/download.php) | Mouse PPI (TAB3 format) |

## Project structure

```
├── config/
│   └── schema_config_mouse.yaml   # BioCypher schema (nodes + edges)
├── src/
│   ├── adapters/
│   │   ├── base.py                # BaseAdapter — shared I/O + dedup
│   │   ├── expression_adapter.py  # Layer 3: GXD expression
│   │   └── ppi_adapter.py         # Layer 4: BioGRID PPI
│   └── layers/
│       ├── expression_export.py   # raw GXD → gene_expression_summary.tsv
│       └── ppi_export.py          # BioGRID download + filter → TSV
├── scripts/
│   ├── gxd_download_aggregate.py  # download all GXD files, compute median TPM
│   └── make_uberon_mapping.py     # map GXD tissue names → UBERON IDs
├── data/                          # input files (not in git)
└── biocypher_out/
    └── mouse/                     # Neo4j-ready CSVs + import script
```

## Quickstart

```bash
# 1. install
pip install biocypher pandas requests beautifulsoup4

# 2. download GXD expression data
python scripts/gxd_download_aggregate.py
# → produces data/gene_expression_summary.tsv

# 3. download BioGRID mouse PPI
#    go to https://thebiogrid.org/download.php
#    download BIOGRID-ORGANISM-Mus_musculus-*.tab3.zip
#    unzip and copy to data/biogrid_mouse_ppi.tsv

# 4. build UBERON mapping for 91 GXD tissues
python scripts/make_uberon_mapping.py
# → produces data/gxd_uberon_mapping.tsv

# 5. build graph
python run.py
```

## Adding a new layer

1. **Adapter** — `src/adapters/<name>_adapter.py`, inherit `BaseAdapter`, set `layer_name`
2. **Export** — `src/layers/<name>_export.py`, raw data → TSV in `data/`
3. **Registration** — add to `SPECIES_LAYERS` in `src/build_graph.py`
4. **Schema** — add nodes/edges to `config/schema_config_mouse.yaml`
