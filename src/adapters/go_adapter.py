"""src/adapters/go_adapter.py

Layer: Gene Ontology (MouseGene -> GoTerm).

Input files:
  data/gene2go.gz    — Entrez GeneID -> GO ID (NCBI)
  data/gene_info.gz  — Entrez GeneID -> MGI ID + symbol

Логика:
  1. gene_info.gz: фильтруем taxid=10090, извлекаем GeneID -> MGI ID
  2. gene2go.gz: фильтруем taxid=10090, берём GeneID -> GO ID + term + aspect
  3. Соединяем: MGI ID -> GO ID
  4. Выдаём GoTerm узлы и MouseGene-[:has_go_term]->GoTerm рёбра
"""

from __future__ import annotations
import gzip
import itertools
from typing import Generator
import pandas as pd
from .base import BaseAdapter, EdgeTuple, NodeTuple

MOUSE_TAXID = "10090"

# Evidence codes для фильтрации — исключаем IEA (автоматические)
# None = брать все включая IEA
EXCLUDE_EVIDENCE = {"IEA"}   # поменяйте на None чтобы брать все


class MouseGoAdapter(BaseAdapter):
    """Yields GoTerm nodes and MouseGene-[:has_go_term]->GoTerm edges."""

    layer_name = "go"

    def _load_entrez_to_mgi(self) -> dict[str, str]:
        """gene_info.gz: GeneID -> MGI ID для мыши."""
        path = self.data_dir / "gene_info.gz"
        if not path.exists():
            raise FileNotFoundError(
                "gene_info.gz not found in data/\n"
                "Download: curl -O https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
            )

        mapping = {}
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                tax_id, gene_id, dbxrefs = parts[0], parts[1], parts[5]
                if tax_id != MOUSE_TAXID:
                    continue
                # dbxrefs формат: MGI:MGI:87854|Ensembl:...|...
                for xref in dbxrefs.split("|"):
                    if xref.startswith("MGI:MGI:"):
                        mgi_id = xref[4:]  # убираем первый "MGI:" → "MGI:87854"
                        mapping[gene_id] = mgi_id
                        break

        print(f"  Entrez→MGI маппинг: {len(mapping):,} записей")
        return mapping

    def _load_go_data(self, entrez_to_mgi: dict[str, str]):
        """gene2go.gz: возвращает DataFrame с колонками mgi_id, go_id, go_term, aspect, evidence."""
        path = self.data_dir / "gene2go.gz"
        if not path.exists():
            raise FileNotFoundError(
                "gene2go.gz not found in data/\n"
                "Download: curl -O https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
            )

        rows = []
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                tax_id, gene_id, go_id, evidence, _, go_term, _, category = parts[:8]
                if tax_id != MOUSE_TAXID:
                    continue
                if EXCLUDE_EVIDENCE and evidence in EXCLUDE_EVIDENCE:
                    continue
                mgi_id = entrez_to_mgi.get(gene_id)
                if not mgi_id:
                    continue
                rows.append({
                    "mgi_id":   mgi_id,
                    "go_id":    go_id,
                    "go_term":  go_term,
                    "aspect":   category,      # Function / Process / Component
                    "evidence": evidence,
                })

        df = pd.DataFrame(rows)
        print(f"  GO аннотаций (мышь, без IEA): {len(df):,}")
        print(f"  Уникальных генов:  {df['mgi_id'].nunique():,}")
        print(f"  Уникальных GO ID:  {df['go_id'].nunique():,}")
        return df

    def _load(self):
        if hasattr(self, "_cache"):
            return self._cache
        entrez_to_mgi = self._load_entrez_to_mgi()
        df = self._load_go_data(entrez_to_mgi)
        self._cache = df
        return df

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        df = self._load()

        # GoTerm узлы
        go_nodes = {}
        for _, row in df.drop_duplicates(subset=["go_id"]).iterrows():
            go_id = row["go_id"]
            if go_id not in go_nodes:
                go_nodes[go_id] = (
                    go_id,
                    "go term",
                    {
                        "go_id":     go_id,
                        "name":      row["go_term"],
                        "namespace": row["aspect"],
                    },
                )

        yield from go_nodes.values()

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        df = self._load()

        # группируем evidence codes по (mgi_id, go_id)
        grouped = (
            df.groupby(["mgi_id", "go_id", "go_term", "aspect"])["evidence"]
            .apply(lambda x: "|".join(sorted(set(x))))
            .reset_index()
        )

        seen = set()
        for _, row in grouped.iterrows():
            pair = (row["mgi_id"], row["go_id"])
            if pair in seen:
                continue
            seen.add(pair)
            yield (
                row["mgi_id"],
                row["go_id"],
                "mouse gene has go term",
                {
                    "evidence_codes": row["evidence"].split("|"),
                    "aspect":         row["aspect"],
                },
            )
