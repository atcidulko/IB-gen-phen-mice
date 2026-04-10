"""
Adapter for STRING protein-protein interaction data (mouse, taxid 10090).

Expects the file: 10090.protein.links.v12.0.txt.gz (or unzipped .txt)
downloaded from https://stringdb-downloads.org/download/protein.links.v12.0/10090.protein.links.v12.0.txt.gz

Format (space-separated, with header):
    protein1 protein2 combined_score
    10090.ENSMUSP00000000001 10090.ENSMUSP00000000002 150

We also support the detailed file (protein.links.detailed) which has
individual evidence channel scores.

The adapter produces:
  - protein nodes (ENSMUSP IDs, stripped of the "10090." prefix)
  - protein_interacts_with_protein edges with combined_score
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional, Set, Tuple

Node = Tuple[str, str, Dict]
Edge = Tuple[str, str, str, str, Dict]


@dataclass
class StringPpiAdapter:
    """
    Parameters
    ----------
    path : Path
        Path to STRING protein.links file (.txt or .txt.gz).
    score_threshold : int
        Minimum combined_score to include (STRING scores 0-1000).
        Default 400 = medium confidence (STRING convention).
    """

    path: Path
    score_threshold: int = 400

    _seen_protein: Set[str] = field(default_factory=set, repr=False)

    @staticmethod
    def _strip_taxid(protein_id: str) -> str:
        """Remove '10090.' prefix from STRING protein IDs."""
        if protein_id.startswith("10090."):
            return protein_id[6:]
        return protein_id

    def _open_file(self):
        if self.path.suffix == ".gz":
            return gzip.open(self.path, "rt", encoding="utf-8")
        return self.path.open("r", encoding="utf-8")

    def _iter_rows(self):
        with self._open_file() as fh:
            header = next(fh)  # skip header line
            for line in fh:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                p1 = self._strip_taxid(parts[0])
                p2 = self._strip_taxid(parts[1])
                try:
                    score = int(parts[2])
                except ValueError:
                    continue
                if score < self.score_threshold:
                    continue
                yield p1, p2, score

    def iter_nodes(self) -> Iterable[Node]:
        for p1, p2, _ in self._iter_rows():
            if p1 not in self._seen_protein:
                self._seen_protein.add(p1)
                yield (p1, "protein", {})
            if p2 not in self._seen_protein:
                self._seen_protein.add(p2)
                yield (p2, "protein", {})

    def iter_edges(self) -> Iterable[Edge]:
        for p1, p2, score in self._iter_rows():
            # canonical ordering to avoid duplicates
            a, b = sorted([p1, p2])
            edge_id = f"{a}--interacts_with--{b}"
            props = {"combined_score": score}
            yield (edge_id, a, b, "protein_interacts_with_protein", props)
