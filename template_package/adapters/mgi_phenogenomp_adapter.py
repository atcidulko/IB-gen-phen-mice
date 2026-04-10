from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
from typing import Dict, Iterable, Iterator, Set, Tuple, Optional, List

Node = Tuple[str, str, Dict]
Edge = Tuple[str, str, str, str, Dict]


@dataclass
class MgiGeneMpAdapter:
    """
    Expects tab-delimited columns:
      0 allelic_composition
      1 allele_symbol
      2 background
      3 mp_id
      4 pmid
      5 gene_mgi_id (can be "MGI:xxx" OR "MGI:xxx|MGI:yyy|...")
      6 genotype_mgi_id
    """
    path: Path

    def iter_rows(self) -> Iterator[list[str]]:
        with self.path.open("r", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            for cols in reader:
                if not cols:
                    continue
                if cols[0].startswith("#"):
                    continue
                cols = [c.strip() for c in cols]
                if not cols[0]:
                    continue
                yield cols

    @staticmethod
    def _norm_pmid(pmid_raw: str) -> Optional[str]:
        pmid_raw = (pmid_raw or "").strip()
        if not pmid_raw:
            return None
        if pmid_raw.startswith("PMID:"):
            return pmid_raw
        if pmid_raw.isdigit():
            return f"PMID:{pmid_raw}"
        # если там сложные значения (несколько PMID и т.п.) — оставляем как есть
        return pmid_raw

    @staticmethod
    def _split_gene_ids(gene_id_raw: str) -> List[str]:
        """
        Handle cases like:
          "MGI:12345"
          "MGI:12345|MGI:67890"
        Return only valid-looking MGI CURIEs.
        """
        gene_id_raw = (gene_id_raw or "").strip()
        if not gene_id_raw:
            return []
        parts = [p.strip() for p in gene_id_raw.split("|")]
        return [p for p in parts if p.startswith("MGI:")]

    def iter_nodes(self) -> Iterable[Node]:
        seen_gene: Set[str] = set()
        seen_mp: Set[str] = set()

        for cols in self.iter_rows():
            if len(cols) < 7:
                continue

            mp_id = cols[3].strip()
            gene_ids = self._split_gene_ids(cols[5])

            if not mp_id.startswith("MP:"):
                continue
            if not gene_ids:
                continue

            # gene nodes (one per gene_id)
            for gene_id in gene_ids:
                if gene_id not in seen_gene:
                    seen_gene.add(gene_id)
                    yield (gene_id, "gene", {})  # symbol можно добавить, если сделаешь mapping MGI->symbol

            # phenotype (MP) node
            if mp_id not in seen_mp:
                seen_mp.add(mp_id)
                # Важно: node_type должен совпадать с input_label в schema_config.yaml
                yield (mp_id, "phenotype", {})  # name можно добавить, если сделаешь mapping MP->name (из mp.obo)

    def iter_edges(self) -> Iterable[Edge]:
        for cols in self.iter_rows():
            if len(cols) < 7:
                continue

            allelic_comp = cols[0].strip().replace("'", "")
            allele_symbol = cols[1].strip().replace("'", "")
            background = cols[2].strip().replace("'", "")
            mp_id = cols[3].strip()
            pmid = self._norm_pmid(cols[4])
            gene_ids = self._split_gene_ids(cols[5])
            genotype_id = cols[6].strip()

            if not mp_id.startswith("MP:"):
                continue
            if not genotype_id.startswith("MGI:"):
                continue
            if not gene_ids:
                continue

            # props общие для всех gene->MP рёбер, порожденных одной строкой
            props: Dict[str, str] = {}
            if pmid:
                props["pmid"] = pmid
            props["genotype_id"] = genotype_id  # всегда полезно для provenance

            if background:
                props["background"] = background
            if allelic_comp:
                props["allelic_composition"] = allelic_comp
            if allele_symbol:
                props["allele_symbol"] = allele_symbol

            # Отдельное ребро для каждого гена
            for gene_id in gene_ids:
                edge_id = f"{gene_id}--gene_has_phenotype--{mp_id}--{genotype_id}"
                if pmid:
                    edge_id += f"--{pmid}"

                yield (edge_id, gene_id, mp_id, "gene_has_phenotype", props)