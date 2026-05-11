"""src/build_graph.py

Собирает граф знаний по мыши через BioCypher.

Слои:
  1. phenotype  — MouseGene → MpTopTerm        (MGI_PhenoGenoMP.rpt + mp.obo)
  2. expression — MouseGene → Tissue           (gene_expression_summary.tsv)
  3. encodes    — MouseGene → Protein          (MGI_Gene_Model_Coord.rpt + biogrid)
  4. ppi        — Protein ↔ Protein            (biogrid_mouse_ppi.tsv)

Запуск:
  python src/build_graph.py
  python src/build_graph.py --skip-layers ppi
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from biocypher import BioCypher

from adapters.phenotype_adapter import MousePhenotypeAdapter
from adapters.expression_adapter import MouseExpressionAdapter
from adapters.encodes_adapter import MouseEncodesAdapter
from adapters.ppi_adapter import MousePPIAdapter
from adapters.go_adapter import MouseGoAdapter

LAYERS: dict[str, type] = {
    "phenotype":  MousePhenotypeAdapter,
    "expression": MouseExpressionAdapter,
    "encodes":    MouseEncodesAdapter,
    "ppi":        MousePPIAdapter,
    "go":         MouseGoAdapter,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build mouse knowledge graph")
    p.add_argument("--data-dir", default="data")
    p.add_argument("--graph-out", default="biocypher_out/mouse")
    p.add_argument("--biocypher-config", default="config/biocypher_config.yaml")
    p.add_argument("--skip-layers", default="")
    p.add_argument("--clean-graph", action="store_true")
    return p.parse_args()


def main() -> None:
    args   = parse_args()
    data_dir  = Path(args.data_dir)
    graph_out = Path(args.graph_out)
    skip      = {s.strip() for s in args.skip_layers.split(",") if s.strip()}

    if args.clean_graph and graph_out.exists():
        shutil.rmtree(graph_out)
        print(f"Удалена старая директория: {graph_out}")

    graph_out.mkdir(parents=True, exist_ok=True)

    bc = BioCypher(
        biocypher_config_path=args.biocypher_config
        if Path(args.biocypher_config).exists() else None,
        output_directory=str(graph_out),
    )

    for layer_name, adapter_cls in LAYERS.items():
        if layer_name in skip:
            print(f"[{layer_name}] пропускаем (--skip-layers)")
            continue

        print(f"\n[{layer_name}] запускаем {adapter_cls.__name__}...")
        adapter = adapter_cls(data_dir=data_dir)
        nodes = list(adapter.get_nodes())
        if nodes:
            bc.write_nodes(iter(nodes))
        bc.write_edges(adapter.get_edges())
        print(f"[{layer_name}] готово")

    bc.write_import_call()
    print(f"\nГотово. Файлы для Neo4j: {graph_out}/")
    print(f"Импорт: bash {graph_out}/neo4j-admin-import-call.sh")


if __name__ == "__main__":
    main()
