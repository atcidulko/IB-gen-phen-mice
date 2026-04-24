"""
scripts/download_biogrid.py

Скачивает актуальные данные BioGRID по мыши (TAB3 формат)
и сохраняет в data/biogrid_mouse_ppi.tsv

Запуск:
    python scripts/download_biogrid.py

Что делает:
    1. Запрашивает BioGRID API — узнаёт последнюю версию базы
    2. Скачивает BIOGRID-ORGANISM-Mus_musculus-*.tab3.zip
    3. Распаковывает и сохраняет как data/biogrid_mouse_ppi.tsv
    4. Печатает статистику: сколько строк, какие методы детекции
"""

import io
import zipfile
from pathlib import Path

import requests

# ── настройки ────────────────────────────────────────────────────────────────

DATA_DIR    = Path("data")
OUTPUT_FILE = DATA_DIR / "biogrid_mouse_ppi.tsv"

# BioGRID REST API — узнаём актуальную версию
VERSION_URL = "https://webservice.thebiogrid.org/version"
DOWNLOAD_URL = (
    "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/"
    "BIOGRID-ORGANISM-Mus_musculus-LATEST.tab3.zip"
)


# ── скачивание ────────────────────────────────────────────────────────────────

def get_latest_version() -> str:
    try:
        r = requests.get(VERSION_URL, timeout=10)
        r.raise_for_status()
        return r.text.strip()
    except Exception:
        return "unknown"


def download_and_unzip() -> None:
    DATA_DIR.mkdir(exist_ok=True)

    version = get_latest_version()
    print(f"BioGRID версия: {version}")
    print(f"Скачиваем: {DOWNLOAD_URL}")
    print("Это может занять минуту (~50 МБ)...\n")

    r = requests.get(DOWNLOAD_URL, timeout=120, stream=True)
    r.raise_for_status()

    # читаем zip из памяти
    content = b""
    total = 0
    for chunk in r.iter_content(chunk_size=1024 * 256):
        content += chunk
        total += len(chunk)
        print(f"\r  скачано: {total / 1_000_000:.1f} МБ", end="", flush=True)
    print()

    # распаковываем
    with zipfile.ZipFile(io.BytesIO(content)) as z:
        # внутри zip один .txt файл
        names = [n for n in z.namelist() if n.endswith(".txt")]
        if not names:
            raise RuntimeError(f"В архиве нет .txt файлов: {z.namelist()}")

        tsv_name = names[0]
        print(f"Распаковываем: {tsv_name}")
        data = z.read(tsv_name)

    OUTPUT_FILE.write_bytes(data)
    print(f"Сохранено: {OUTPUT_FILE}  ({OUTPUT_FILE.stat().st_size / 1_000_000:.1f} МБ)")


# ── статистика ────────────────────────────────────────────────────────────────

def print_stats() -> None:
    import pandas as pd

    print("\nЧитаем файл для статистики...")
    df = pd.read_csv(OUTPUT_FILE, sep="\t", dtype=str, low_memory=False)

    print(f"\n── Статистика BioGRID (мышь) ──")
    print(f"  Всего взаимодействий:  {len(df):,}")
    print(f"  Уникальных генов A:    {df['Official Symbol Interactor A'].nunique():,}")

    if "Experimental System" in df.columns:
        print(f"\n  Топ-10 методов детекции:")
        top = df["Experimental System"].value_counts().head(10)
        for method, count in top.items():
            print(f"    {method:45s} {count:>8,}")

    if "Throughput" in df.columns:
        print(f"\n  Throughput:")
        for t, c in df["Throughput"].value_counts().items():
            print(f"    {t}: {c:,}")


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    if OUTPUT_FILE.exists():
        size_mb = OUTPUT_FILE.stat().st_size / 1_000_000
        ans = input(
            f"Файл уже существует ({size_mb:.1f} МБ): {OUTPUT_FILE}\n"
            "Перескачать? [y/N]: "
        ).strip().lower()
        if ans != "y":
            print("Пропускаем скачивание.")
            print_stats()
            return

    download_and_unzip()
    print_stats()


if __name__ == "__main__":
    main()
