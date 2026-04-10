#!/usr/bin/env python3
"""
Download STRING protein-protein interaction data for mouse (taxid 10090).

Usage:
    python scripts/download_string.py

Downloads ~400 MB compressed file to data/ directory.
"""

import urllib.request
import os
from pathlib import Path

# STRING v12.0 — mouse
URL = "https://stringdb-downloads.org/download/protein.links.v12.0/10090.protein.links.v12.0.txt.gz"
DEST = Path("data/10090.protein.links.v12.0.txt.gz")


def main():
    DEST.parent.mkdir(parents=True, exist_ok=True)

    if DEST.exists():
        size_mb = DEST.stat().st_size / 1e6
        print(f"✓ File already exists: {DEST} ({size_mb:.1f} MB)")
        return

    print(f"Downloading STRING mouse PPI data...")
    print(f"  URL:  {URL}")
    print(f"  Dest: {DEST}")
    print(f"  This may take a few minutes (~400 MB)...")

    urllib.request.urlretrieve(URL, DEST)

    size_mb = DEST.stat().st_size / 1e6
    print(f"✓ Downloaded: {DEST} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
