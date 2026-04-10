from __future__ import annotations

from pathlib import Path
from biocypher import BioCypher

from template_package.adapters.mgi_phenogenomp_adapter import MgiGeneMpAdapter
from template_package.adapters.expression_atlas_adapter import ExpressionAtlasAdapter
from template_package.adapters.string_ppi_adapter import StringPpiAdapter


def main() -> None:
    bc = BioCypher()

    # -- 1. Genotype-Phenotype (MGI) -----------------------------------
    pheno_adapter = MgiGeneMpAdapter(
        path=Path("data/MGI_PhenoGenoMP.rpt"),
    )
    bc.write_nodes(pheno_adapter.iter_nodes())
    bc.write_edges(pheno_adapter.iter_edges())

    # -- 2. Gene Expression (Expression Atlas, numeric TPM) ------------
    tpm_path = Path("data/E-GEOD-74747-tpms.tsv")
    sdrf_path = Path("data/E-GEOD-74747.condensed-sdrf.tsv")

    if tpm_path.exists() and sdrf_path.exists():
        expr_adapter = ExpressionAtlasAdapter(
            tpm_path=tpm_path,
            sdrf_path=sdrf_path,
            experiment_id="E-GEOD-74747",
            tpm_cutoff=0.5,
        )
        bc.write_nodes(expr_adapter.iter_nodes())
        bc.write_edges(expr_adapter.iter_edges())
    else:
        print("Warning: Expression Atlas files not found in data/.")
        print("  Expected: data/E-GEOD-74747-tpms.tsv")
        print("           data/E-GEOD-74747.condensed-sdrf.tsv")

    # -- 3. Protein-Protein Interactions (STRING) ----------------------
    string_path = Path("data/10090.protein.links.v12.0.txt.gz")
    if not string_path.exists():
        string_path = Path("data/10090.protein.links.v12.0.txt")

    if string_path.exists():
        ppi_adapter = StringPpiAdapter(
            path=string_path,
            score_threshold=400,
        )
        bc.write_nodes(ppi_adapter.iter_nodes())
        bc.write_edges(ppi_adapter.iter_edges())
    else:
        print("Warning: STRING data not found. Run:")
        print("  python scripts/download_string.py")

    # -- finalise ------------------------------------------------------
    bc.write_import_call()
    bc.summary()


if __name__ == "__main__":
    main()
