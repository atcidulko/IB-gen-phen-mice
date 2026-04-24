"""src/build_graph.py

Собирает граф знаний по мыши через BioCypher.

Слои:
  1. phenotype   — MouseGene → MpTopTerm        (MGI_PhenoGenoMP.rpt)
  2. go          — MouseGene → GoTerm           (gene2go.gz + gene_info.gz)
  3. expression  — MouseGene → Tissue           (gene_expression_summary.tsv)
  4. ppi         — Protein ↔ Protein            (biogrid_mouse_ppi.tsv)

Запуск:
  python src/build_graph.py
  python src/build_graph.py --skip-layers ppi
  python src/build_graph.py --skip-layers go,ppi
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from biocypher import BioCypher

from adapters.expression_adapter import MouseExpressionAdapter
from adapters.ppi_adapter import MousePPIAdapter

# ── регистрация всех слоёв ────────────────────────────────────────────────────
# Добавляйте новые адаптеры сюда — build_graph подхватит автоматически.

LAYERS: dict[str, type] = {
    "expression": MouseExpressionAdapter,
    "ppi":        MousePPIAdapter,
}

# Слои которые ещё не реализованы — раскомментируйте когда будут готовы:
# from adapters.phenotype_adapter import MousePhenotypeAdapter
# from adapters.go_adapter import MouseGoAdapter
# LAYERS["phenotype"] = MousePhenotypeAdapter
# LAYERS["go"]        = MouseGoAdapter


# ── аргументы CLI ─────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build mouse knowledge graph")
    p.add_argument(
        "--data-dir",
        default="data",
        help="Directory with input TSV files (default: data/)",
    )
    p.add_argument(
        "--graph-out",
        default="biocypher_out/mouse",
        help="Output directory for Neo4j import files (default: biocypher_out/mouse/)",
    )
    p.add_argument(
        "--biocypher-config",
        default="config/biocypher_config.yaml",
        help="Path to biocypher_config.yaml",
    )
    p.add_argument(
        "--skip-layers",
        default="",
        help="Comma-separated list of layer names to skip, e.g. --skip-layers ppi,go",
    )
    p.add_argument(
        "--clean-graph",
        action="store_true",
        help="Delete and recreate output directory before building",
    )
    return p.parse_args()


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    data_dir  = Path(args.data_dir)
    graph_out = Path(args.graph_out)
    skip      = {s.strip() for s in args.skip_layers.split(",") if s.strip()}

    # опционально — чистим старый вывод
    if args.clean_graph and graph_out.exists():
        shutil.rmtree(graph_out)
        print(f"Удалена старая директория: {graph_out}")

    graph_out.mkdir(parents=True, exist_ok=True)

    # инициализируем BioCypher
    bc = BioCypher(
        biocypher_config_path=args.biocypher_config
        if Path(args.biocypher_config).exists()
        else None,
        output_directory=str(graph_out),
    )

    # запускаем каждый слой
    for layer_name, adapter_cls in LAYERS.items():
        if layer_name in skip:
            print(f"[{layer_name}] пропускаем (--skip-layers)")
            continue

        print(f"\n[{layer_name}] запускаем {adapter_cls.__name__}...")
        adapter = adapter_cls(data_dir=data_dir)

        bc.write_nodes(adapter.get_nodes())
        bc.write_edges(adapter.get_edges())
        print(f"[{layer_name}] готово")

    # финализируем — генерируем скрипт импорта для Neo4j
    bc.write_import_call()
    bc.summary()

    print(f"\nГотово. Файлы для Neo4j: {graph_out}/")
    print(f"Импорт: bash {graph_out}/neo4j-admin-import-call.sh")


if __name__ == "__main__":
    main()
