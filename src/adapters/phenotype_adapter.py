"""src/adapters/phenotype_adapter.py

Layer 1: Gene–Phenotype associations (MGI).

Input files:
  MGI_PhenoGenoMP.rpt   — MGI genotype → MP annotations
  mp.obo                — Mammalian Phenotype Ontology

Columns in MGI_PhenoGenoMP.rpt (tab-separated, no header):
  0  Allelic Composition
  1  Allele Symbol
  2  Genetic Background
  3  MP ID              ← MP term for this genotype
  4  PubMed ID
  5  MGI Gene ID        ← our gene node ID
  6  MGI Genotype ID

Logic:
  For each (MGI Gene ID, MP ID) pair:
    1. Walk up the MP ontology to find the top-level term (depth 1)
    2. Emit MpTopTerm node + MouseGene -[has_mp_top_term]-> MpTopTerm edge

  Top-level MP terms are direct children of the root MP:0000001.
"""

from __future__ import annotations

import itertools
from typing import Generator

from .base import BaseAdapter, EdgeTuple, NodeTuple

# ── OBO parser (minimal, no external deps) ────────────────────────────────────

def parse_obo(path) -> dict[str, dict]:
    """
    Parse mp.obo and return a dict:
      mp_id -> {name, is_a: [parent_ids]}
    """
    terms = {}
    current = None

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                current = {"is_a": [], "name": "", "obsolete": False}
            elif line == "" and current is not None:
                if "id" in current and not current["obsolete"]:
                    terms[current["id"]] = current
                current = None
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:].strip()
                elif line.startswith("name: "):
                    current["name"] = line[6:].strip()
                elif line.startswith("is_a: "):
                    parent = line[6:].split("!")[0].strip()
                    current["is_a"].append(parent)
                elif line.startswith("is_obsolete: true"):
                    current["obsolete"] = True

    # flush last term
    if current and "id" in current and not current.get("obsolete"):
        terms[current["id"]] = current

    return terms


def build_top_level_map(terms: dict) -> dict[str, str]:
    """
    For each MP term, find its top-level ancestor
    (direct child of MP:0000001, the root).

    Returns: mp_id -> top_level_mp_id
    """
    ROOT = "MP:0000001"

    # find all direct children of root — these are top-level terms
    top_level = {
        mp_id
        for mp_id, data in terms.items()
        if ROOT in data.get("is_a", [])
    }

    cache: dict[str, str | None] = {}

    def find_top(mp_id: str, visited: set) -> str | None:
        if mp_id in cache:
            return cache[mp_id]
        if mp_id == ROOT:
            cache[mp_id] = None
            return None
        if mp_id in top_level:
            cache[mp_id] = mp_id
            return mp_id
        if mp_id in visited:
            cache[mp_id] = None
            return None
        visited.add(mp_id)
        if mp_id not in terms:
            cache[mp_id] = None
            return None
        for parent in terms[mp_id].get("is_a", []):
            result = find_top(parent, visited)
            if result:
                cache[mp_id] = result
                return result
        cache[mp_id] = None
        return None

    return {
        mp_id: find_top(mp_id, set())
        for mp_id in terms
    }


# ── Adapter ───────────────────────────────────────────────────────────────────

class MousePhenotypeAdapter(BaseAdapter):
    """Yields MpTopTerm nodes and MouseGene-[:has_mp_top_term]->MpTopTerm edges."""

    layer_name = "phenotype"

    def __init__(self, data_dir):
        super().__init__(data_dir)
        self._cache = None  # cache parsed data

    def _load(self):
        """Load and parse both input files. Cached after first call."""
        if self._cache is not None:
            return self._cache

        import pandas as pd

        # parse ontology
        obo_path = self.data_dir / "mp.obo"
        if not obo_path.exists():
            raise FileNotFoundError(
                "mp.obo not found in data/\n"
                "Download: curl -O https://purl.obolibrary.org/obo/mp.obo"
            )
        terms    = parse_obo(obo_path)
        top_map  = build_top_level_map(terms)

        # parse MGI file
        rpt_path = self.data_dir / "MGI_PhenoGenoMP.rpt"
        if not rpt_path.exists():
            raise FileNotFoundError(
                "MGI_PhenoGenoMP.rpt not found in data/\n"
                "Download: curl -O https://www.informatics.jax.org/downloads/reports/MGI_PhenoGenoMP.rpt"
            )

        df = pd.read_csv(
            rpt_path,
            sep="\t",
            header=None,
            dtype=str,
            names=[
                "allelic_composition", "allele_symbol", "genetic_background",
                "mp_id", "pubmed_id", "mgi_gene_id", "mgi_genotype_id",
            ],
        )

        # filter out rows with missing gene or MP ID
        df = df.dropna(subset=["mgi_gene_id", "mp_id"])
        df = df[df["mp_id"].str.startswith("MP:")]

        # split rows where mgi_gene_id contains multiple IDs (e.g. MGI:123|MGI:456)
        rows = []
        for _, row in df.iterrows():
            mgi_ids = [x.strip() for x in str(row["mgi_gene_id"]).split("|")
                       if x.strip().startswith("MGI:")]
            for mgi_id in mgi_ids:
                new_row = row.copy()
                new_row["mgi_gene_id"] = mgi_id
                rows.append(new_row)

        df = pd.DataFrame(rows).reset_index(drop=True)

        self._cache = (df, top_map, terms)
        return self._cache

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        df, top_map, terms = self._load()

        gene_nodes = {}
        mp_nodes = {}

        for _, row in df.iterrows():
            mgi_id = row["mgi_gene_id"]
            if mgi_id not in gene_nodes:
                gene_nodes[mgi_id] = (
                    mgi_id,
                    "mouse gene",
                    {
                        "mgi_id":     mgi_id,
                        "symbol":     "",
                        "ensembl_id": "",
                        "organism":   "Mus musculus",
                    },
                )
            top_id = top_map.get(row["mp_id"])
            if top_id and top_id not in mp_nodes:
                name = terms.get(top_id, {}).get("name", "")
                mp_nodes[top_id] = (
                    top_id,
                    "mp top term",
                    {"mp_id": top_id, "name": name},
                )

        yield from gene_nodes.values()
        yield from mp_nodes.values()

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        df, top_map, _ = self._load()

        seen: set[tuple] = set()

        for _, row in df.iterrows():
            mgi_id = row["mgi_gene_id"]
            top_id = top_map.get(row["mp_id"])
            if not top_id:
                continue
            pair = (mgi_id, top_id)
            if pair in seen:
                continue
            seen.add(pair)
            yield (
                mgi_id,
                top_id,
                "mouse gene has mp top term",
                {},
            )
