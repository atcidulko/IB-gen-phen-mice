"""
scripts/train_pykeen_neuro.py

Обучает модель предсказания связей для одного фенотипа:
MP:0003631 — nervous system phenotype (5762 гена)

Стратегия:
  1. Берём гены с известной связью с MP:0003631
  2. Добавляем их экспрессию и PPI как контекст
  3. Обучаем TransE предсказывать ген → MP:0003631

Запуск:
    python scripts/train_pykeen_neuro.py

Результат:
    models/neuro/   — обученная модель
    results/neuro/  — метрики и предсказания
"""

from pathlib import Path
import pandas as pd
import torch
from pykeen.pipeline import pipeline
from pykeen.triples import TriplesFactory

# ── настройки ────────────────────────────────────────────────────────────────

TRIPLES_FILE  = Path("data/triples.tsv")
TARGET_MP     = "MP:0003631"   # nervous system phenotype
OUT_DIR       = Path("models/neuro")
RESULTS_DIR   = Path("results/neuro")

MODEL         = "DistMult"
EPOCHS        = 300
EMBEDDING_DIM = 128
BATCH_SIZE    = 512
DEVICE        = "cuda" if torch.cuda.is_available() else "cpu"


def main():
    print(f"Устройство: {DEVICE}")
    print(f"Целевой фенотип: {TARGET_MP} (nervous system phenotype)\n")

    # ── загрузка всех троек ───────────────────────────────────────────────────
    print("Загружаем triples.tsv...")
    all_df = pd.read_csv(TRIPLES_FILE, sep="\t", header=None,
                         names=["head", "relation", "tail"], dtype=str)
    all_df = all_df.dropna()

    # ── шаг 1: берём гены связанные с TARGET_MP, делим на train/holdout ────────
    import numpy as np
    pheno_df = all_df[all_df["relation"] == "has_mp_top_term"]
    target_pheno = pheno_df[pheno_df["tail"] == TARGET_MP]
    target_genes = list(target_pheno["head"].unique())

    rng = np.random.default_rng(42)
    rng.shuffle(target_genes)
    split = int(len(target_genes) * 0.8)
    train_genes  = set(target_genes[:split])
    holdout_genes = set(target_genes[split:])

    print(f"Генов с {TARGET_MP}:    {len(target_genes):,}")
    print(f"  → train:             {len(train_genes):,}")
    print(f"  → holdout (predict): {len(holdout_genes):,}")

    # ── шаг 2: фильтруем все рёбра — только train гены в фенотипе ────────────
    # фенотипные рёбра: только для train генов (holdout скрываем)
    pheno_triples = target_pheno[target_pheno["head"].isin(train_genes)][["head", "relation", "tail"]].copy()

    # экспрессия: для ВСЕХ генов (train + holdout) — это контекст
    expr_df = all_df[all_df["relation"] == "expressed_in"]
    expr_triples = expr_df[expr_df["head"].isin(set(target_genes))][["head", "relation", "tail"]]

    # PPI: взаимодействия где хотя бы один партнёр среди всех генов
    ppi_df = all_df[all_df["relation"] == "interacts_with"]
    ppi_triples = ppi_df[
        ppi_df["head"].isin(set(target_genes)) | ppi_df["tail"].isin(set(target_genes))
    ][["head", "relation", "tail"]]

    # объединяем
    triples_df = pd.concat([pheno_triples, expr_triples, ppi_triples], ignore_index=True)
    triples_df = triples_df.drop_duplicates()

    print(f"\nИтоговый граф для обучения:")
    print(f"  Фенотип рёбра:    {len(pheno_triples):,}")
    print(f"  Экспрессия рёбра: {len(expr_triples):,}")
    print(f"  PPI рёбра:        {len(ppi_triples):,}")
    print(f"  Всего троек:      {len(triples_df):,}")
    print(f"  Уникальных сущностей: {pd.concat([triples_df['head'], triples_df['tail']]).nunique():,}")

    # ── шаг 3: обучение ───────────────────────────────────────────────────────
    tf = TriplesFactory.from_labeled_triples(
        triples=triples_df[["head", "relation", "tail"]].values,
        create_inverse_triples=True,
    )

    training, testing, validation = tf.split([0.8, 0.1, 0.1], random_state=42)

    print(f"\nТренировочных: {training.num_triples:,}")
    print(f"Тестовых:      {testing.num_triples:,}")
    print(f"Валидационных: {validation.num_triples:,}")
    print(f"\nЗапускаем {MODEL}...")

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

    # ── сохранение ────────────────────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    result.save_to_directory(OUT_DIR)

    metrics = result.metric_results.to_df()
    metrics.to_csv(RESULTS_DIR / "metrics.csv", index=False)
    print("\n── Метрики ──")
    print(metrics.to_string(index=False))

    # ── предсказание: holdout гены ───────────────────────────────────────────
    print(f"\n── Предсказываем для holdout генов ({len(holdout_genes):,}) ──")

    genes_no_pheno = holdout_genes

    entity_to_id   = result.training.entity_to_id
    relation_to_id = result.training.relation_to_id
    model          = result.model
    model.eval()

    rel_id = relation_to_id.get("has_mp_top_term")
    mp_eid = entity_to_id.get(TARGET_MP)

    if rel_id is None or mp_eid is None:
        print("Не удалось найти relation/entity ID — проверьте обучение")
        return

    rows = []
    with torch.no_grad():
        for gene in genes_no_pheno:
            if gene not in entity_to_id:
                continue
            gene_eid = entity_to_id[gene]
            score = model.score_hrt(
                torch.tensor([[gene_eid, rel_id, mp_eid]], device=DEVICE)
            ).item()
            rows.append({"gene": gene, "mp_term": TARGET_MP, "score": score})

    if not rows:
        print("Нет генов для предсказания — все гены уже имеют этот фенотип в графе обучения.")
        return

    pred_df = pd.DataFrame(rows).sort_values("score", ascending=False)
    pred_df.to_csv(RESULTS_DIR / "new_gene_predictions.tsv", sep="\t", index=False)

    print(f"\nТоп-20 новых генов предсказанных для {TARGET_MP}:")
    print(pred_df.head(20).to_string(index=False))
    print(f"\nСохранено: {RESULTS_DIR}/new_gene_predictions.tsv")


if __name__ == "__main__":
    main()
