"""src/adapters/ppi_adapter.py

Layer 4: Protein–protein interactions (BioGRID, mouse).

Input file:
  biogrid_mouse_ppi.tsv

  Download from https://thebiogrid.org/download.php
  Choose: "BIOGRID-ORGANISM-Mus_musculus-*.tab3.zip"
  Unzip and rename to biogrid_mouse_ppi.tsv

  Columns used (BioGRID TAB3 format):
    Official Symbol Interactor A   — gene symbol of protein A
    Official Symbol Interactor B   — gene symbol of protein B
    BioGRID ID for Interactor A    — BioGRID internal ID
    BioGRID ID for Interactor B
    Experimental System            — detection method
    Pubmed ID                      — supporting publication
    Score                          — confidence (empty in BioGRID, filled from STRING if merged)
    Throughput                     — Low Throughput / High Throughput

  Optional companion file (for MGI ID lookup):
    mgi_symbol_to_id.tsv          — two columns: symbol, mgi_id
    Generate with scripts/make_mgi_symbol_map.py
    or download MGI_EntrezGene.rpt from MGI and parse.

Note on Protein nodes:
  BioGRID identifies proteins by gene symbol, not UniProt.
  We use "BIOGRID:<BioGRID_ID>" as the protein node ID,
  and store the gene symbol + MGI ID as properties.
  UniProt IDs can be added later via a separate mapping step.
"""

from __future__ import annotations

import itertools
from typing import Generator

from .base import BaseAdapter, EdgeTuple, NodeTuple

# BioGRID TAB3 column names we actually need
_COL_SYM_A      = "Official Symbol Interactor A"
_COL_SYM_B      = "Official Symbol Interactor B"
_COL_BGID_A     = "BioGRID ID Interactor A"       # без "for"
_COL_BGID_B     = "BioGRID ID Interactor B"       # без "for"
_COL_METHOD     = "Experimental System"
_COL_PUBMED     = "Publication Source"             # формат PUBMED:12345678
_COL_SCORE      = "Score"
_COL_THROUGHPUT = "Throughput"
_COL_UNIPROT_A  = "SWISS-PROT Accessions Interactor A"
_COL_UNIPROT_B  = "SWISS-PROT Accessions Interactor B"

# Only keep physical / direct interaction types
PHYSICAL_METHODS = {
    "Two-hybrid",
    "Co-crystal Structure",
    "Far Western",
    "Co-purification",
    "Affinity Capture-MS",
    "Affinity Capture-Western",
    "Co-fractionation",
    "Biochemical Activity",
    "Proximity Label-MS",
    "FRET",
    "PCA",
    "Co-localization",
    "Reconstituted Complex",
    "Protein-peptide",
}


class MousePPIAdapter(BaseAdapter):
    """Yields Protein nodes and Protein-[:interacts_with]->Protein edges."""

    layer_name = "ppi"

    def _load_symbol_to_mgi(self) -> dict[str, str]:
        """Load gene symbol → MGI ID mapping. Returns empty dict if file missing."""
        path = self.data_dir / "mgi_symbol_to_id.tsv"
        if not path.exists():
            return {}
        df = self._read("mgi_symbol_to_id.tsv")
        if df.empty:
            return {}
        return dict(zip(df["symbol"], df["mgi_id"]))

    def _load_ppi(self):
        """Load and filter BioGRID TAB3 file."""
        df = self._read("biogrid_mouse_ppi.tsv")
        if df.empty:
            return df

        required = [_COL_SYM_A, _COL_SYM_B, _COL_BGID_A, _COL_BGID_B, _COL_METHOD]
        missing  = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(
                f"biogrid_mouse_ppi.tsv is missing columns: {missing}\n"
                "Make sure you downloaded the TAB3 format from BioGRID."
            )

        # keep only physical interactions
        df = df[df[_COL_METHOD].isin(PHYSICAL_METHODS)].copy()

        # drop self-interactions
        df = df[df[_COL_SYM_A] != df[_COL_SYM_B]]

        # drop rows with missing IDs
        df = df.dropna(subset=[_COL_BGID_A, _COL_BGID_B])

        return df

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        df = self._load_ppi()
        if df.empty:
            return

        sym_to_mgi = self._load_symbol_to_mgi()

        def _protein_nodes(syms_ids_uniprots):
            seen = set()
            for sym, bg_id, uniprot in syms_ids_uniprots:
                node_id = f"BIOGRID:{bg_id}"
                if node_id in seen:
                    continue
                seen.add(node_id)
                # UniProt может быть несколько через |, берём первый
                uniprot_clean = str(uniprot).split("|")[0].strip() if uniprot and str(uniprot) != "nan" else ""
                yield (
                    node_id,
                    "protein",
                    {
                        "uniprot_id": uniprot_clean,
                        "mgi_id":     sym_to_mgi.get(sym, ""),
                        "symbol":     sym,
                    },
                )

        pairs_a = zip(df[_COL_SYM_A], df[_COL_BGID_A], df.get(_COL_UNIPROT_A, [""] * len(df)))
        pairs_b = zip(df[_COL_SYM_B], df[_COL_BGID_B], df.get(_COL_UNIPROT_B, [""] * len(df)))

        yield from self._unique_nodes(
            itertools.chain(_protein_nodes(pairs_a), _protein_nodes(pairs_b))
        )

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        df = self._load_ppi()
        if df.empty:
            return

        # deduplicate edges: (id_a, id_b) regardless of order
        seen: set[frozenset] = set()

        for _, row in df.iterrows():
            id_a = f"BIOGRID:{row[_COL_BGID_A]}"
            id_b = f"BIOGRID:{row[_COL_BGID_B]}"
            pair = frozenset([id_a, id_b])
            if pair in seen:
                continue
            seen.add(pair)

            score = row.get(_COL_SCORE, "")
            try:
                score_float = round(float(score), 4)
            except (ValueError, TypeError):
                score_float = 0.0

            pubmed = row.get(_COL_PUBMED, "")
            # формат: "PUBMED:12345678" или "-"
            if pubmed and str(pubmed) not in ("-", "nan"):
                pubmed_list = [p.strip().replace("PUBMED:", "") for p in str(pubmed).split("|") if "PUBMED:" in p]
            else:
                pubmed_list = []

            yield (
                id_a,
                id_b,
                "protein interacts with",
                {
                    "confidence_score":  score_float,
                    "detection_method":  row.get(_COL_METHOD, ""),
                    "source_db":         "BioGRID",
                    "pubmed_ids":        pubmed_list,
                },
            )
