"""
Диагностика: смотрим реальные заголовки колонок в файлах GXD.
Запустите ЭТОТ скрипт первым, посмотрите вывод,
потом скажите мне что там написано — поправим основной скрипт.
"""

import gzip
import io
from pathlib import Path

import pandas as pd

# ищем любой .rpt.gz в папке gxd_raw
raw_dir = Path("gxd_raw")
files = sorted(raw_dir.glob("*.rpt.gz"))

if not files:
    print("Папка gxd_raw пуста или не найдена.")
    print("Убедитесь что скрипт запущен из той же директории что и gxd_raw/")
else:
    path = files[0]
    print(f"Читаем файл: {path.name}\n")

    with gzip.open(path, "rb") as f:
        raw = f.read()

    # пробуем прочитать первые 3 строки как есть
    df = pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str, nrows=3)

    print("── Реальные названия колонок ──")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i:2d}  {repr(col)}")

    print(f"\n── Первая строка данных ──")
    if len(df) > 0:
        for col in df.columns:
            print(f"  {col!r:45s} → {df[col].iloc[0]!r}")
