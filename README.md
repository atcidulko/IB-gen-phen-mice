# IB-gen-phen-mice

Knowledge graph linking mouse genotype, phenotype, gene expression (quantitative TPM), and protein-protein interactions. Built with [BioCypher](https://biocypher.org/).

## Data sources

| Data | Source | What it adds |
|------|--------|-------------|
| Genotype → Phenotype | [MGI PhenoGenoMP](https://www.informatics.jax.org/) | Gene–MP term edges with allelic composition, background, PMID |
| Gene Expression (TPM) | [Expression Atlas](https://www.ebi.ac.uk/gxa/) (EMBL-EBI) | Numeric TPM per gene per tissue from baseline RNA-seq |
| Protein–Protein Interactions | [STRING v12](https://string-db.org/) (mouse, taxid 10090) | Combined confidence scores for physical/functional PPI |

## Graph schema

**Nodes:** MGI Gene, Ensembl Gene, MP Phenotype, Tissue, Protein

**Edges:**
- `gene_has_phenotype` — MGI Gene → MP Phenotype (with allelic composition, PMID)
- `has_expression` — Ensembl Gene → Tissue (with TPM value, experiment ID)
- `protein_interacts_with_protein` — Protein ↔ Protein (with STRING combined score)

## Setup

```bash
# install dependencies
pip install biocypher

# download external data
python scripts/download_string.py              # STRING PPI (~400 MB)
python scripts/download_expression_atlas.py    # Expression Atlas TPM

# build the knowledge graph
python create_knowledge_graph.py
```

## Project structure

```
├── config/
│   ├── biocypher_config.yaml       # BioCypher settings
│   └── schema_config.yaml          # Graph schema definition
├── data/
│   ├── MGI_PhenoGenoMP.rpt         # MGI genotype-phenotype
│   ├── E-ERAD-169.rpt              # Expression Atlas data
│   ├── E-GEOD-63813.rpt
│   ├── E-GEOD-65775.rpt
│   └── 10090.protein.links.*.gz    # STRING PPI (downloaded)
├── template_package/
│   └── adapters/
│       ├── mgi_phenogenomp_adapter.py
│       ├── expression_atlas_adapter.py
│       └── string_ppi_adapter.py
├── scripts/
│   ├── download_string.py
│   └── download_expression_atlas.py
└── create_knowledge_graph.py
```
