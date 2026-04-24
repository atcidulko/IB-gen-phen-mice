"""
scripts/make_uberon_mapping.py

Автоматически маппирует названия тканей из GXD на UBERON ID
через OLS API (Ontology Lookup Service, EBI).

Запуск:
    python scripts/make_uberon_mapping.py

Входные данные:
    data/gene_expression_summary.tsv  — берёт уникальные значения
                                        из колонки anatomical_structure

Результат:
    data/gxd_uberon_mapping.tsv       — два столбца: gxd_name, uberon_id

Что делает скрипт:
    1. Берёт все уникальные названия тканей из вашего файла экспрессии
    2. Для каждой ткани спрашивает OLS API: "есть ли в UBERON термин
       с таким названием?"
    3. Если находит точное совпадение — записывает UBERON ID
    4. Если не находит — оставляет пустым и печатает список для
       ручной проверки

После запуска:
    Откройте data/gxd_uberon_mapping.tsv в Excel или любом редакторе.
    Найдите строки где uberon_id пустой.
    Для каждой зайдите на https://www.ebi.ac.uk/ols4/ontologies/uberon
    и найдите подходящий ID вручную.
"""

import time
from pathlib import Path

import pandas as pd
import requests

# ── настройки ────────────────────────────────────────────────────────────────

DATA_DIR        = Path("data")
EXPRESSION_FILE = DATA_DIR / "gene_expression_summary.tsv"
OUTPUT_FILE     = DATA_DIR / "gxd_uberon_mapping.tsv"

OLS_API         = "https://www.ebi.ac.uk/ols4/api/search"
DELAY           = 0.3   # секунды между запросами

# ── известные маппинги которые OLS находит неточно ───────────────────────────
# Добавляйте сюда если знаете правильный ID заранее

MANUAL_OVERRIDES: dict[str, str] = {
    # ── общие ткани ───────────────────────────────────────────────── #
    "embryo":               "UBERON:0000922",
    "metanephros":          "UBERON:0000081",
    "cerebral hemisphere":  "UBERON:0001869",
    "vomeronasal organ":    "UBERON:0002255",
    "hippocampus":          "UBERON:0001954",
    "cerebellum":           "UBERON:0002037",
    "olfactory lobe":       "UBERON:0005902",
    "medulla oblongata":    "UBERON:0001896",
    "eye":                  "UBERON:0000970",
    "heart":                "UBERON:0000948",
    "liver":                "UBERON:0002107",
    "lung":                 "UBERON:0002048",
    "kidney":               "UBERON:0002113",
    "spleen":               "UBERON:0002106",
    "testis":               "UBERON:0000473",
    "ovary":                "UBERON:0000992",
    "brain":                "UBERON:0000955",
    "skin":                 "UBERON:0002097",
    "bone marrow":          "UBERON:0002371",
    "thymus":               "UBERON:0002370",
    "adipose tissue":       "UBERON:0001013",
    "muscle organ":         "UBERON:0001630",
    "pancreas":             "UBERON:0001264",
    "stomach":              "UBERON:0000945",
    "intestine":            "UBERON:0000160",
    "colon":                "UBERON:0001155",
    "blood":                "UBERON:0000178",
    "placenta":             "UBERON:0001987",
    "uterus":               "UBERON:0000995",
    "prostate gland":       "UBERON:0002367",
    "thyroid gland":        "UBERON:0002046",
    "adrenal gland":        "UBERON:0002369",
    "spinal cord":          "UBERON:0002240",
    "mammary gland":        "UBERON:0001911",

    # ── ранние эмбриональные стадии ───────────────────────────────── #
    "1-cell stage conceptus":  "UBERON:0000106",  # zygote
    "2-cell stage conceptus":  "UBERON:0000107",  # 2-cell stage embryo

    # ── жаберные дуги (краниофациальное развитие) ─────────────────── #
    "1st branchial arch mandibular component": "UBERON:0009752",
    "1st branchial arch maxillary component":  "UBERON:0004079",

    # ── мозг ──────────────────────────────────────────────────────── #
    "arcuate nucleus":                "UBERON:0001581",
    "dentate gyrus":                  "UBERON:0001885",
    "thalamus":                       "UBERON:0001897",
    "telencephalon subventricular zone": "UBERON:0004922",
    "future brain floor plate":       "UBERON:0003850",

    # ── жировая ткань ─────────────────────────────────────────────── #
    "brown fat":   "UBERON:0001348",  # brown adipose tissue
    "white fat":   "UBERON:0001347",  # white adipose tissue

    # ── кишечник ──────────────────────────────────────────────────── #
    "cecum":                          "UBERON:0001153",
    "large intestine crypt of lieberkuhn": "UBERON:0013475",
    "islets of Langerhans":           "UBERON:0000006",

    # ── мышцы ─────────────────────────────────────────────────────── #
    "gastrocnemius muscle":     "UBERON:0001370",
    "quadriceps femoris muscle": "UBERON:0001377",
    "soleus":                   "UBERON:0001378",
    "upper leg muscle":         "UBERON:0004262",

    # ── кость ─────────────────────────────────────────────────────── #
    "femur diaphysis":   "UBERON:0006800",
    "femur metaphysis":  "UBERON:0006793",

    # ── сердце ────────────────────────────────────────────────────── #
    "heart ventricle":   "UBERON:0002082",

    # ── репродуктивные клетки ─────────────────────────────────────── #
    "oocyte":           "UBERON:0000023",
    "secondary oocyte": "UBERON:0002174",
    "spermatid":        "UBERON:0004235",
    "spermatocyte":     "UBERON:0004230",

    # ── лёгкое ────────────────────────────────────────────────────── #
    "right lung middle lobe": "UBERON:0002174",  # right lung lobe

    # ── ухо ───────────────────────────────────────────────────────── #
    "organ of Corti":   "UBERON:0002227",

    # ── эмбриональные структуры ───────────────────────────────────── #
    "head surface ectoderm": "UBERON:0016887",
    "latero-nasal process":  "UBERON:0006234",
    "medial-nasal process":  "UBERON:0006235",
    "nasal pit":             "UBERON:0006236",
}


# ── OLS lookup ────────────────────────────────────────────────────────────────

def lookup_uberon(tissue_name: str) -> str:
    """
    Запрашивает OLS API и возвращает UBERON ID для точного совпадения.
    Возвращает пустую строку если ничего не найдено.
    """
    params = {
        "q":          tissue_name,
        "ontology":   "uberon",
        "exact":      "true",
        "fieldList":  "iri,label,obo_id",
        "rows":       5,
    }
    try:
        r = requests.get(OLS_API, params=params, timeout=10)
        r.raise_for_status()
        docs = r.json().get("response", {}).get("docs", [])
        for doc in docs:
            label   = doc.get("label", "").lower()
            obo_id  = doc.get("obo_id", "")
            if label == tissue_name.lower() and obo_id.startswith("UBERON:"):
                return obo_id
    except Exception as e:
        print(f"  OLS error for '{tissue_name}': {e}")
    return ""


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    # 1. загружаем уникальные ткани
    if not EXPRESSION_FILE.exists():
        print(f"Файл не найден: {EXPRESSION_FILE}")
        print("Сначала запустите scripts/gxd_download_aggregate.py")
        return

    df_expr = pd.read_csv(EXPRESSION_FILE, sep="\t", dtype=str,
                          usecols=["anatomical_structure"])
    tissues = sorted(df_expr["anatomical_structure"].dropna().unique())
    print(f"Уникальных тканей в файле экспрессии: {len(tissues)}\n")

    # 2. маппинг
    rows = []
    not_found = []

    for i, name in enumerate(tissues, 1):
        # сначала проверяем ручные переопределения
        if name.lower() in {k.lower(): v for k, v in MANUAL_OVERRIDES.items()}:
            uberon_id = MANUAL_OVERRIDES.get(name, MANUAL_OVERRIDES.get(name.lower(), ""))
            status = "manual"
        else:
            uberon_id = lookup_uberon(name)
            status = "ols" if uberon_id else "NOT FOUND"
            time.sleep(DELAY)

        print(f"[{i:3d}/{len(tissues)}] {name:40s} → {uberon_id or '—':20s} ({status})")

        if not uberon_id:
            not_found.append(name)

        rows.append({"gxd_name": name, "uberon_id": uberon_id})

    # 3. сохраняем
    DATA_DIR.mkdir(exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\nСохранено: {OUTPUT_FILE}")

    # 4. отчёт о не найденных
    if not_found:
        print(f"\n── {len(not_found)} тканей без UBERON ID (нужна ручная проверка) ──")
        for name in not_found:
            print(f"  {name}")
        print("\nОткройте data/gxd_uberon_mapping.tsv и заполните uberon_id вручную.")
        print("Ищите на: https://www.ebi.ac.uk/ols4/ontologies/uberon")
    else:
        print("\nВсе ткани успешно замаплены!")


if __name__ == "__main__":
    main()
