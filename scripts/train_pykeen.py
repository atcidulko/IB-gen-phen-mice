"""
scripts/train_pykeen.py

Обучает модель предсказания связей через PyKEEN.
Читает готовый файл троек data/triples.tsv (сделайте его через make_triples.py).

Запуск на сервере:
    python scripts/train_pykeen.py

Результат:
    models/pykeen/   — обученная модель
    results/         — метрики и предсказания ген → фенотип
"""

from pathlib import Path
import pandas as pd
import torch
from pykeen.pipeline import pipeline
from pykeen.triples import TriplesFactory

# ── настройки ────────────────────────────────────────────────────────────────

TRIPLES_FILE  = Path("data/triples.tsv")
OUT_DIR       = Path("models/pykeen")
RESULTS_DIR   = Path("results")

MODEL         = "DistMult"    # TransE / RotatE / DistMult / ComplEx
EPOCHS        = 30
EMBEDDING_DIM = 128
BATCH_SIZE    = 1024        # увеличьте если есть GPU и много RAM

# автоматически выбираем GPU если есть
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def main():
    print(f"Устройство: {DEVICE}")
    if DEVICE == "cuda":
        print(f"GPU: {torch.cuda.get_device_name(0)}")

    # ── загрузка троек ────────────────────────────────────────────────────────
    if not TRIPLES_FILE.exists():
        raise FileNotFoundError(
            f"{TRIPLES_FILE} не найден.\n"
            "Сначала запустите: python scripts/make_triples.py"
        )

    print(f"\nЗагружаем тройки из {TRIPLES_FILE}...")
    triples_df = pd.read_csv(TRIPLES_FILE, sep="\t", header=None,
                             names=["head", "relation", "tail"], dtype=str)
    triples_df = triples_df.dropna()

    print(f"Всего троек:          {len(triples_df):,}")
    print(f"Уникальных сущностей: {pd.concat([triples_df['head'], triples_df['tail']]).nunique():,}")
    print(f"Типов отношений:      {triples_df['relation'].nunique()}")
    print(triples_df["relation"].value_counts().to_string())

    # ── создаём TriplesFactory ────────────────────────────────────────────────
    tf = TriplesFactory.from_labeled_triples(
        triples=triples_df[["head", "relation", "tail"]].values,
        create_inverse_triples=True,  # обратные рёбра улучшают качество
    )

    training, testing, validation = tf.split([0.8, 0.1, 0.1], random_state=42)

    print(f"\nТренировочных троек: {training.num_triples:,}")
    print(f"Тестовых:            {testing.num_triples:,}")
    print(f"Валидационных:       {validation.num_triples:,}")

    # ── обучение ─────────────────────────────────────────────────────────────
    print(f"\nЗапускаем {MODEL} на {DEVICE}...")

    result = pipeline(
        training=training,
        testing=testing,
        validation=validation,
        model=MODEL,
        model_kwargs=dict(embedding_dim=EMBEDDING_DIM),
        training_kwargs=dict(
            num_epochs=EPOCHS,
            batch_size=BATCH_SIZE,
        ),
        optimizer="Adam",
        optimizer_kwargs=dict(lr=0.001),
        negative_sampler="basic",
        negative_sampler_kwargs=dict(num_negs_per_pos=5),
        evaluator_kwargs=dict(filtered=True),
        random_seed=42,
        device=DEVICE,
    )

    # ── сохранение ───────────────────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    result.save_to_directory(OUT_DIR)
    print(f"\nМодель сохранена: {OUT_DIR}/")

    metrics = result.metric_results.to_df()
    metrics.to_csv(RESULTS_DIR / "metrics.csv", index=False)
    print("\n── Метрики ──")
    print(metrics.to_string(index=False))

    # ── предсказание ген → фенотип ────────────────────────────────────────────
    print("\n── Предсказываем фенотипы ──")

    pheno_df = triples_df[triples_df["relation"] == "has_mp_top_term"]
    genes    = pheno_df["head"].unique()
    mp_terms = pheno_df["tail"].unique()

    entity_to_id   = result.training.entity_to_id
    relation_to_id = result.training.relation_to_id
    model          = result.model
    model.eval()

    rel_id = relation_to_id.get("has_mp_top_term")

    rows = []
    with torch.no_grad():
        for gene in genes:
            if gene not in entity_to_id:
                continue
            gene_id = entity_to_id[gene]
            for mp in mp_terms:
                if mp not in entity_to_id:
                    continue
                mp_id = entity_to_id[mp]
                score = model.score_hrt(
                    torch.tensor([[gene_id, rel_id, mp_id]], device=DEVICE)
                ).item()
                rows.append({"gene": gene, "mp_term": mp, "score": score})

    pred_df = pd.DataFrame(rows)
    pred_df = pred_df.sort_values(["gene", "score"], ascending=[True, False])
    top_preds = pred_df.groupby("gene").head(5)
    top_preds.to_csv(RESULTS_DIR / "phenotype_predictions.tsv", sep="\t", index=False)
    print(f"Сохранено: {RESULTS_DIR}/phenotype_predictions.tsv")
    print(top_preds.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
