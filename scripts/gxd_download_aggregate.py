"""
GXD RNA-seq: скачать все файлы и агрегировать медианный TPM по гену × ткани.

Источник: https://www.informatics.jax.org/downloads/reports/gxdrnaseq/

Колонки, которые используем:
  1  MGI Gene ID
  2  Ensembl ID
  3  Gene Symbol
  6  Anatomical Structure   <- ткань
  7  Theiler Stage          <- для фильтрации стадий
  10 Strain
  11 Mutant Allele Pair(s)  <- фильтруем: только wild-type (пустая строна)
  16 Detected               <- фильтруем: только "Yes"
  17 avg_TPM                <- числовое значение, с которым работаем

Результат: gene_expression_summary.tsv
  mgi_gene_id | ensembl_id | gene_symbol | anatomical_structure |
  n_experiments | median_avg_tpm | mean_avg_tpm | tpm_level
"""

import gzip
import io
import logging
import time
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── настройки ────────────────────────────────────────────────────────────────

BASE_URL    = "https://www.informatics.jax.org/downloads/reports/gxdrnaseq/"
OUT_DIR     = Path("gxd_raw")          # куда складываем скачанные файлы
OUT_SUMMARY = Path("gene_expression_summary.tsv")

# Фильтры (None = не фильтровать)
ONLY_WILD_TYPE   = True   # исключить строки с мутантными аллелями
ONLY_DETECTED    = True   # оставить только Detected == "Yes"
THEILER_STAGES   = None   # например {26, 27, 28} для взрослой мыши; None = все

DELAY_BETWEEN_REQUESTS = 0.3   # секунды между запросами (вежливость к серверу)

# ── шаг 1: получить список файлов ────────────────────────────────────────────

def list_rpt_files(base_url: str) -> list[str]:
    """Парсим HTML-листинг директории, возвращаем URL всех .rpt.gz файлов."""
    log.info("Запрашиваем листинг директории: %s", base_url)
    r = requests.get(base_url, timeout=30)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, "html.parser")
    urls = []
    for a in soup.find_all("a", href=True):
        href = a["href"]
        if href.endswith(".rpt.gz"):
            urls.append(base_url.rstrip("/") + "/" + href.lstrip("/"))
    log.info("Найдено файлов: %d", len(urls))
    return urls


# ── шаг 2: скачать файлы ─────────────────────────────────────────────────────

def download_files(urls: list[str], out_dir: Path) -> list[Path]:
    """Скачивает файлы в out_dir, пропускает уже существующие."""
    out_dir.mkdir(exist_ok=True)
    paths = []
    for i, url in enumerate(urls, 1):
        fname = url.split("/")[-1]
        dest  = out_dir / fname
        if dest.exists():
            log.info("[%d/%d] Уже есть: %s", i, len(urls), fname)
        else:
            log.info("[%d/%d] Скачиваем: %s", i, len(urls), fname)
            r = requests.get(url, timeout=60, stream=True)
            r.raise_for_status()
            dest.write_bytes(r.content)
            time.sleep(DELAY_BETWEEN_REQUESTS)
        paths.append(dest)
    return paths


# ── шаг 3: читать один файл ──────────────────────────────────────────────────

# Файлы разделены pipe (|), заголовок присутствует в каждом файле.
USECOLS = [
    "MGI Gene ID",
    "Ensembl ID",
    "Gene Symbol",
    "Anatomical Structure",
    "Theiler Stage",
    "Mutant Allele Pair(s)",
    "Detected",
    "avg_TPM",
]

def read_rpt(path: Path) -> pd.DataFrame | None:
    """Читает один .rpt.gz файл, возвращает DataFrame с нужными колонками."""
    try:
        with gzip.open(path, "rb") as f:
            raw = f.read()
        df = pd.read_csv(
            io.BytesIO(raw),
            sep="|",
            dtype=str,
            low_memory=False,
        )
        # нормализуем имена колонок (убираем лишние пробелы)
        df.columns = df.columns.str.strip()

        # проверяем наличие нужных колонок
        missing = [c for c in USECOLS if c not in df.columns]
        if missing:
            log.warning("Пропускаем %s — нет колонок: %s", path.name, missing)
            return None

        df = df[USECOLS].copy()
        df["avg_TPM"] = pd.to_numeric(df["avg_TPM"], errors="coerce")
        df["source_file"] = path.name
        return df

    except Exception as e:
        log.warning("Ошибка при чтении %s: %s", path.name, e)
        return None


# ── шаг 4: объединить и отфильтровать ────────────────────────────────────────

def load_all(paths: list[Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        df = read_rpt(path)
        if df is not None:
            frames.append(df)
    if not frames:
        raise RuntimeError("Ни один файл не прочитан — проверьте источник данных.")
    combined = pd.concat(frames, ignore_index=True)
    log.info("Всего строк после объединения: %d", len(combined))
    return combined


def apply_filters(df: pd.DataFrame) -> pd.DataFrame:
    before = len(df)

    if ONLY_DETECTED:
        df = df[df["Detected"].str.strip().str.upper() == "YES"]
        log.info("После фильтра Detected=Yes: %d (убрано %d)", len(df), before - len(df))
        before = len(df)

    if ONLY_WILD_TYPE:
        # Мутантные аллели: пустая строка или NaN означает wild-type
        wt_mask = df["Mutant Allele Pair(s)"].isna() | (df["Mutant Allele Pair(s)"].str.strip() == "")
        df = df[wt_mask]
        log.info("После фильтра wild-type: %d (убрано %d)", len(df), before - len(df))
        before = len(df)

    if THEILER_STAGES is not None:
        df["_ts_num"] = pd.to_numeric(df["Theiler Stage"].str.extract(r"(\d+)")[0], errors="coerce")
        df = df[df["_ts_num"].isin(THEILER_STAGES)]
        df = df.drop(columns=["_ts_num"])
        log.info("После фильтра Theiler Stages %s: %d (убрано %d)",
                 THEILER_STAGES, len(df), before - len(df))

    # убираем строки с NaN в avg_TPM
    df = df.dropna(subset=["avg_TPM"])
    log.info("Строк после всех фильтров: %d", len(df))
    return df


# ── шаг 5: агрегация ─────────────────────────────────────────────────────────

TPM_BINS   = [0, 0.5, 10, 100, float("inf")]
TPM_LABELS = ["Absent", "Low", "Medium", "High"]


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Группируем по (MGI Gene ID, Ensembl ID, Gene Symbol, Anatomical Structure).
    Считаем медиану avg_TPM, среднее и число экспериментов.
    """
    group_cols = ["MGI Gene ID", "Ensembl ID", "Gene Symbol", "Anatomical Structure"]

    agg = (
        df.groupby(group_cols, sort=False)["avg_TPM"]
        .agg(
            n_experiments="count",
            median_avg_tpm="median",
            mean_avg_tpm="mean",
        )
        .reset_index()
    )

    # категориальный уровень экспрессии на основе медианы
    agg["tpm_level"] = pd.cut(
        agg["median_avg_tpm"],
        bins=TPM_BINS,
        labels=TPM_LABELS,
        right=True,
    ).astype(str)

    # округляем числа
    agg["median_avg_tpm"] = agg["median_avg_tpm"].round(4)
    agg["mean_avg_tpm"]   = agg["mean_avg_tpm"].round(4)

    # переименовываем для графа
    agg = agg.rename(columns={
        "MGI Gene ID":          "mgi_gene_id",
        "Ensembl ID":           "ensembl_id",
        "Gene Symbol":          "gene_symbol",
        "Anatomical Structure": "anatomical_structure",
    })

    log.info("Уникальных пар ген × ткань: %d", len(agg))
    log.info("Уникальных генов: %d",   agg["mgi_gene_id"].nunique())
    log.info("Уникальных тканей: %d",  agg["anatomical_structure"].nunique())
    return agg


# ── шаг 6: сохранить ─────────────────────────────────────────────────────────

def save(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False)
    log.info("Сохранено: %s  (%d строк)", path, len(df))


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    # 1. листинг
    urls = list_rpt_files(BASE_URL)

    # 2. скачать
    paths = download_files(urls, OUT_DIR)

    # 3–4. читать + фильтровать
    combined = load_all(paths)
    filtered = apply_filters(combined)

    # 5. агрегация
    summary = aggregate(filtered)

    # 6. сохранить
    save(summary, OUT_SUMMARY)

    # краткая статистика
    print("\n── Распределение по уровням TPM ──")
    print(summary["tpm_level"].value_counts().to_string())
    print("\n── Топ-10 тканей по числу генов ──")
    print(
        summary.groupby("anatomical_structure")["mgi_gene_id"]
        .nunique()
        .sort_values(ascending=False)
        .head(10)
        .to_string()
    )


if __name__ == "__main__":
    main()
