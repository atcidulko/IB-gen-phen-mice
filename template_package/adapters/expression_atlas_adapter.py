"""
Adapter for Expression Atlas baseline RNA-Seq data (mouse).

Handles the Expression Atlas download format where:
  - TPM file has columns: GeneID, Gene Name, g1, g2, ..., gN
  - Each cell contains comma-separated replicate values (e.g. "54,54,54,54,54")
  - Condensed SDRF file maps run IDs to tissue names via "factor organism part"
  - Column order g1..gN corresponds to sorted run IDs

The adapter produces:
  - ensembl_gene nodes (Ensembl mouse gene IDs)
  - tissue nodes (organ/tissue names)
  - has_expression edges: gene -> tissue with mean TPM value
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

Node = Tuple[str, str, Dict]
Edge = Tuple[str, str, str, str, Dict]


@dataclass
class ExpressionAtlasAdapter:
    """
    Parameters
    ----------
    tpm_path : Path
        Path to the *-tpms.tsv file from Expression Atlas.
    sdrf_path : Path
        Path to the *.condensed-sdrf.tsv file (maps runs -> tissues).
    experiment_id : str
        Experiment accession (e.g. "E-GEOD-74747"). Used for provenance.
    tpm_cutoff : float
        Minimum mean TPM to create an expression edge (default 0.5).
    """

    tpm_path: Path
    sdrf_path: Path
    experiment_id: str = ""
    tpm_cutoff: float = 0.5

    _seen_gene: Set[str] = field(default_factory=set, repr=False)
    _seen_tissue: Set[str] = field(default_factory=set, repr=False)

    # ------------------------------------------------------------------
    # Parse SDRF to get ordered list of tissue names
    # ------------------------------------------------------------------
    def _parse_tissue_map(self) -> Dict[str, str]:
        """
        Parse condensed-sdrf.tsv to build {run_id: tissue_name}.
        We look for lines with 'factor' and 'organism part'.
        """
        run_to_tissue: Dict[str, str] = {}
        with self.sdrf_path.open("r", encoding="utf-8") as fh:
            for line in fh:
                cols = line.strip().split("\t")
                if len(cols) < 6:
                    continue
                # cols: experiment, ??, run_id, type, factor_name, value, [uri]
                run_id = cols[2].strip()
                entry_type = cols[3].strip()
                factor_name = cols[4].strip()
                value = cols[5].strip()

                if entry_type == "factor" and factor_name == "organism part":
                    run_to_tissue[run_id] = value

        return run_to_tissue

    def _get_ordered_tissues(self) -> List[str]:
        """
        Return tissue names in the order matching g1, g2, ..., gN columns.
        Expression Atlas orders columns by sorted run ID.
        """
        run_to_tissue = self._parse_tissue_map()
        sorted_runs = sorted(run_to_tissue.keys())
        return [run_to_tissue[r] for r in sorted_runs]

    # ------------------------------------------------------------------
    # Parse TPM values
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_cell(cell: str) -> Optional[float]:
        """
        Parse a cell like "54,54,54,54,54" -> mean TPM (54.0).
        All values are typically identical (already averaged by Atlas),
        but we take the mean to be safe.
        """
        cell = cell.strip()
        if not cell or cell in ("", "NA", "N/A", "-"):
            return None
        try:
            values = [float(v) for v in cell.split(",")]
            return sum(values) / len(values)
        except ValueError:
            return None

    def _iter_records(self):
        """
        Yield dicts: {gene_id, gene_name, tissue, tpm}
        """
        tissues = self._get_ordered_tissues()
        if not tissues:
            print(f"Warning: No tissues found in {self.sdrf_path}")
            return

        with self.tpm_path.open("r", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)  # GeneID, Gene Name, g1, g2, ...

            n_data_cols = len(header) - 2
            if n_data_cols != len(tissues):
                print(
                    f"Warning: Column mismatch: {n_data_cols} data columns "
                    f"vs {len(tissues)} tissues in SDRF. Using min."
                )

            n = min(n_data_cols, len(tissues))

            for row in reader:
                if len(row) < 3:
                    continue

                gene_id = row[0].strip()
                gene_name = row[1].strip()

                if not gene_id or not gene_id.startswith("ENSMUSG"):
                    continue

                for i in range(n):
                    tpm = self._parse_cell(row[i + 2])
                    if tpm is None:
                        continue

                    yield {
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "tissue": tissues[i],
                        "tpm": round(tpm, 2),
                    }

    # ------------------------------------------------------------------
    # BioCypher API
    # ------------------------------------------------------------------
    def iter_nodes(self) -> Iterable[Node]:
        """Yield ensembl_gene and tissue nodes (deduplicated)."""
        for rec in self._iter_records():
            gene_id = rec["gene_id"]
            if gene_id not in self._seen_gene:
                self._seen_gene.add(gene_id)
                props = {}
                if rec["gene_name"]:
                    props["symbol"] = rec["gene_name"]
                yield (gene_id, "ensembl_gene", props)

            tissue = rec["tissue"]
            if tissue not in self._seen_tissue:
                self._seen_tissue.add(tissue)
                yield (tissue, "tissue", {"name": tissue})

    def iter_edges(self) -> Iterable[Edge]:
        """Yield has_expression edges: gene -> tissue with TPM."""
        for rec in self._iter_records():
            if rec["tpm"] < self.tpm_cutoff:
                continue

            gene_id = rec["gene_id"]
            tissue = rec["tissue"]
            exp_id = self.experiment_id or self.tpm_path.stem

            edge_id = f"{gene_id}--has_expression--{tissue}--{exp_id}"
            props = {
                "tpm": rec["tpm"],
                "experiment_id": exp_id,
            }

            yield (edge_id, gene_id, tissue, "has_expression", props)
