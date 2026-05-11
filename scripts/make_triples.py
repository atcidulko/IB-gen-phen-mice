"""
scripts/make_triples.py

Конвертирует CSV файлы из BioCypher в единый файл троек для PyKEEN.

Запуск:
    python scripts/make_triples.py

Результат:
    data/triples.tsv
"""

import pandas as pd
from pathlib import Path

GRAPH_DIR  = Path("biocypher_out/mouse")
OUTPUT     = Path("data/triples.tsv")


def load_edges() -> pd.DataFrame:
    parts = []

    # 1. MouseGene → MpTopTerm
    for path in sorted(GRAPH_DIR.glob("MouseGeneHasMpTopTerm-part*.csv")):
        df = pd.read_csv(path, sep=";", dtype=str, header=None)
        parts.append(pd.DataFrame({
            "head":     df[0],
            "relation": "has_mp_top_term",
            "tail":     df[2],
        }))
        print(f"  {path.name}: {len(df):,} рёбер фенотипа")

    # 2. MouseGene → Tissue
    for path in sorted(GRAPH_DIR.glob("MouseGeneExpressedInTissue-part*.csv")):
        df = pd.read_csv(path, sep=";", dtype=str, header=None)
        parts.append(pd.DataFrame({
            "head":     df[0],
            "relation": "expressed_in",
            "tail":     df[5],
        }))
        print(f"  {path.name}: {len(df):,} рёбер экспрессии")

    # 3. Protein ↔ Protein
    for path in sorted(GRAPH_DIR.glob("ProteinInteractsWith-part*.csv")):
        df = pd.read_csv(path, sep=";", dtype=str, header=None)
        parts.append(pd.DataFrame({
            "head":     df[0],
            "relation": "interacts_with",
            "tail":     df[6],
        }))
        print(f"  {path.name}: {len(df):,} PPI рёбер")

    # 4. MouseGene → GoTerm
    for path in sorted(GRAPH_DIR.glob("MouseGeneHasGoTerm-part*.csv")):
        df = pd.read_csv(path, sep=";", dtype=str, header=None)
        parts.append(pd.DataFrame({
            "head":     df[0],
            "relation": "has_go_term",
            "tail":     df[4],
        }))
        print(f"  {path.name}: {len(df):,} рёбер GO")

    # 5. MouseGene → Protein (encodes)
    for path in sorted(GRAPH_DIR.glob("MouseGeneEncodesProtein-part*.csv")):
        df = pd.read_csv(path, sep=";", dtype=str, header=None)
        parts.append(pd.DataFrame({
            "head":     df[0],
            "relation": "encodes",
            "tail":     df[2],
        }))
        print(f"  {path.name}: {len(df):,} рёбер encodes")

    result = pd.concat(parts, ignore_index=True)
    result = result.dropna()
    result = result[result["head"].str.strip() != ""]
    result = result[result["tail"].str.strip() != ""]

    print(f"\nВсего троек:          {len(result):,}")
    print(f"Уникальных сущностей: {pd.concat([result['head'], result['tail']]).nunique():,}")
    print(f"Типов отношений:      {result['relation'].nunique()}")
    print(f"Распределение:")
    print(result["relation"].value_counts().to_string())
    return result


def main():
    print("── Конвертируем CSV → triples.tsv ──\n")
    df = load_edges()
    OUTPUT.parent.mkdir(exist_ok=True)
    df.to_csv(OUTPUT, sep="\t", index=False, header=False)
    print(f"\nСохранено: {OUTPUT}  ({OUTPUT.stat().st_size / 1_000_000:.1f} МБ)")


if __name__ == "__main__":
    main()
