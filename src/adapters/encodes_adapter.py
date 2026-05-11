"""src/adapters/encodes_adapter.py

Layer: Gene -> Protein (encodes edges).
"""

from __future__ import annotations
from typing import Generator
import pandas as pd
from .base import BaseAdapter, EdgeTuple, NodeTuple


class MouseEncodesAdapter(BaseAdapter):
    layer_name = "encodes"

    def _load_symbol_to_mgi(self) -> dict[str, str]:
        path = self.data_dir / "MGI_Gene_Model_Coord.rpt"
        if not path.exists():
            raise FileNotFoundError("MGI_Gene_Model_Coord.rpt not found in data/")
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, header=0, index_col=False)
        df = df[["1. MGI accession id", "3. marker symbol"]].dropna()
        mapping = dict(zip(df["3. marker symbol"].str.strip().str.upper(),
                           df["1. MGI accession id"].str.strip()))
        print(f"  Маппинг символ->MGI: {len(mapping):,} записей")
        return mapping

    def _load_symbol_to_biogrid(self) -> dict[str, str]:
        path = self.data_dir / "biogrid_mouse_ppi.tsv"
        if not path.exists():
            raise FileNotFoundError("biogrid_mouse_ppi.tsv not found in data/")
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False,
                         usecols=["Official Symbol Interactor A", "BioGRID ID Interactor A",
                                  "Official Symbol Interactor B", "BioGRID ID Interactor B"])
        pairs_a = df[["Official Symbol Interactor A", "BioGRID ID Interactor A"]].rename(
            columns={"Official Symbol Interactor A": "symbol", "BioGRID ID Interactor A": "biogrid_id"})
        pairs_b = df[["Official Symbol Interactor B", "BioGRID ID Interactor B"]].rename(
            columns={"Official Symbol Interactor B": "symbol", "BioGRID ID Interactor B": "biogrid_id"})
        combined = pd.concat([pairs_a, pairs_b]).dropna().drop_duplicates(subset=["symbol"])
        mapping = dict(zip(combined["symbol"].str.strip().str.upper(),
                           combined["biogrid_id"].str.strip()))
        print(f"  Маппинг символ->BioGRID: {len(mapping):,} записей")
        return mapping

    def get_nodes(self) -> Generator[NodeTuple, None, None]:
        yield from ()

    def get_edges(self) -> Generator[EdgeTuple, None, None]:
        print("  Загружаем маппинги для encodes рёбер...")
        sym_to_mgi = self._load_symbol_to_mgi()
        sym_to_biogrid = self._load_symbol_to_biogrid()
        edges = []
        for symbol, biogrid_id in sym_to_biogrid.items():
            mgi_id = sym_to_mgi.get(symbol)
            if not mgi_id:
                continue
            edges.append((mgi_id, f"BIOGRID:{biogrid_id}", "encodes", {}))
        print(f"  Рёбер encodes (ген -> белок): {len(edges):,}")
        yield from edges
