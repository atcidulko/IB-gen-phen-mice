"""src/adapters/expression_adapter.py

Layer 3: Gene expression (GXD / MGI RNA-seq).

Input file (produced by gxd_download_aggregate.py):
  gene_expression_summary.tsv

  Columns used:
    mgi_gene_id          — MGI:xxxxxxx  (node ID for MouseGene)
    ensembl_id           — ENSMUSG...
    gene_symbol          — e.g. Trp53
    anatomical_structure — GXD tissue name
    median_avg_tpm       — aggregated expression value
    tpm_level            — Absent / Low / Medium / High
    n_experiments        — number of GXD experiments

UBERON mapping file:
  gxd_uberon_mapping.tsv

  Two columns (tab-separated, with header):
    gxd_name    — matches anatomical_structure in expression file
    uberon_id   — e.g. UBERON:0000473  (empty string if no mapping found)

  Generate this file with scripts/make_uberon_mapping.py
  or fill manually for the 91 GXD tissues.
"""

from __future__ import annotations

import itertools
from typing import Generator

from .base import BaseAdapter, EdgeTuple, NodeTuple


class MouseExpressionAdapter(BaseAdapter):
    """Yields Tissue nodes and MouseGene-[:expressed_in]->Tissue edges."""

    layer_name = "expression"

    # GXD tissues without a clean UBERON mapping get this fallback
    UBERON_FALLBACK = ""

    def _load_uberon_map(self) -> dict[str, str]:
        """Load gxd_name → uberon_id mapping. Returns empty dict if file missing."""
        path = self.data_dir / "gxd_uberon_mapping.tsv"
        if not path.exists():
            return {}
        df = self._read("gxd_uberon_mapping.tsv")
        if df.empty or "gxd_name" not in df.columns:
            return {}
        df = df.fillna("")
        return dict(zip(df["gxd_name"], df["uberon_id"]))

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        df = self._read("gene_expression_summary.tsv")
        if df.empty:
            return

        uberon_map = self._load_uberon_map()

        gene_iter = self._gene_nodes(df)
        tissue_iter = self._tissue_nodes(df, uberon_map)

        yield from self._unique_nodes(itertools.chain(gene_iter, tissue_iter))

    def _gene_nodes(self, df):
        seen = set()
        for _, row in df.iterrows():
            mgi_id = row.get("mgi_gene_id", "")
            if not mgi_id or mgi_id in seen:
                continue
            seen.add(mgi_id)
            yield (
                mgi_id,
                "mouse gene",
                {
                    "mgi_id": mgi_id,
                    "symbol": row.get("gene_symbol", ""),
                    "ensembl_id": row.get("ensembl_id", ""),
                    "organism": "Mus musculus",
                },
            )

    def _tissue_nodes(self, df, uberon_map):
        seen = set()
        for _, row in df.iterrows():
            gxd_name = row.get("anatomical_structure", "")
            if not gxd_name:
                continue
            uberon_id = uberon_map.get(gxd_name, self.UBERON_FALLBACK)
            # use UBERON ID as node ID when available, else GXD name
            node_id = uberon_id if uberon_id else f"GXD:{gxd_name}"
            if node_id in seen:
                continue
            seen.add(node_id)
            yield (
                node_id,
                "tissue",
                {
                    "uberon_id": uberon_id,
                    "name": gxd_name,
                    "gxd_name": gxd_name,
                },
            )

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        df = self._read("gene_expression_summary.tsv")
        if df.empty:
            return

        uberon_map = self._load_uberon_map()

        for _, row in df.iterrows():
            mgi_id   = row.get("mgi_gene_id", "")
            gxd_name = row.get("anatomical_structure", "")
            if not mgi_id or not gxd_name:
                continue

            uberon_id = uberon_map.get(gxd_name, self.UBERON_FALLBACK)
            target_id = uberon_id if uberon_id else f"GXD:{gxd_name}"

            tpm = row.get("median_avg_tpm", "")
            try:
                tpm_float = float(tpm)
            except (ValueError, TypeError):
                continue

            yield (
                mgi_id,
                target_id,
                "mouse gene expressed in",
                {
                    "median_avg_tpm": round(tpm_float, 4),
                    "tpm_level":      row.get("tpm_level", ""),
                    "n_experiments":  int(row.get("n_experiments", 0)),
                },
            )
