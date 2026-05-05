import pandas as pd
import numpy as np
from pathlib import Path

CSV_PATH = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"

FITNESS_COL = "ddelta"
SD_FLOOR = 0.1
OUTPUT_PREFIX = "weighted_optimal_path_ddelta"

START = "WWWWWWWW"
END = "OOOOOOOO"

def hamming_distance_wo(a: str, b: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def step_direction_by_load(a_row, b_row):
    if a_row["mut_count"] < b_row["mut_count"]:
        return a_row, b_row
    elif b_row["mut_count"] < a_row["mut_count"]:
        return b_row, a_row
    else:
        return None, None

def reconstruct_path(pred, start, end):
    if end not in pred and end != start:
        return None
    path = [end]
    cur = end
    while cur != start:
        cur = pred[cur]
        path.append(cur)
    path.reverse()
    return path

def path_to_edge_table(path, edges_lookup_df):
    if path is None or len(path) < 2:
        return pd.DataFrame()
    recs = []
    for a, b in zip(path[:-1], path[1:]):
        hit = edges_lookup_df[(edges_lookup_df["source"] == a) & (edges_lookup_df["target"] == b)]
        if not hit.empty:
            recs.append(hit.iloc[0].to_dict())
    return pd.DataFrame(recs)

df = pd.read_csv(CSV_PATH)

required = {"library_id", "wo", "mut_count", FITNESS_COL}
missing = required - set(df.columns)
if missing:
    raise KeyError(f"Missing required columns: {missing}")

df = df[["library_id", "wo", "mut_count", FITNESS_COL]].copy()
df = df.replace([np.inf, -np.inf], np.nan)
df[FITNESS_COL] = pd.to_numeric(df[FITNESS_COL], errors="coerce")
df["mut_count"] = pd.to_numeric(df["mut_count"], errors="coerce")
df = df.dropna(subset=["library_id", "wo", "mut_count", FITNESS_COL]).copy()

avg_df = (
    df.groupby(["wo", "mut_count"], as_index=False)
      .agg(mean_fitness=(FITNESS_COL, "mean"))
)

if START not in set(avg_df["wo"]) or END not in set(avg_df["wo"]):
    raise RuntimeError("Start and/or end genotype missing from averaged dataframe.")

lib_lookup = {}
for lib in sorted(df["library_id"].unique()):
    lib_lookup[lib] = df[df["library_id"] == lib].set_index("wo")

genotypes = sorted(avg_df["wo"].unique())
edge_records = []

for i, g1 in enumerate(genotypes):
    for g2 in genotypes[i + 1:]:
        if hamming_distance_wo(g1, g2) != 1:
            continue

        row1 = avg_df.loc[avg_df["wo"] == g1].iloc[0]
        row2 = avg_df.loc[avg_df["wo"] == g2].iloc[0]
        low_row, high_row = step_direction_by_load(row1, row2)

        if low_row is None or high_row is None:
            continue

        low_g, high_g = low_row["wo"], high_row["wo"]

        step_effects = []
        for _, lookup in lib_lookup.items():
            if low_g in lookup.index and high_g in lookup.index:
                low_val = lookup.loc[low_g, FITNESS_COL]
                high_val = lookup.loc[high_g, FITNESS_COL]
                step_effects.append(high_val - low_val)

        if len(step_effects) == 0:
            continue

        step_effects = np.asarray(step_effects, dtype=float)
        mean_step = float(step_effects.mean())
        sd_step = float(step_effects.std(ddof=1)) if len(step_effects) > 1 else 0.0
        raw_repro = abs(mean_step) / (sd_step + SD_FLOOR)
        log_repro = np.log10(raw_repro + 1.0)
        weighted_score = mean_step * log_repro

        edge_records.append({
            "source": low_g,
            "target": high_g,
            "mut_source": int(low_row["mut_count"]),
            "mut_target": int(high_row["mut_count"]),
            "mean_step": mean_step,
            "sd_step": sd_step,
            "raw_reproducibility": raw_repro,
            "log_reproducibility": log_repro,
            "weighted_score": weighted_score,
            "source_mean_fitness": float(low_row["mean_fitness"]),
            "target_mean_fitness": float(high_row["mean_fitness"]),
        })

edges_df = pd.DataFrame(edge_records)

if edges_df.empty:
    raise RuntimeError("No valid directed one-step edges were constructed.")

order_nodes = avg_df.sort_values(["mut_count", "wo"])["wo"].tolist()

adj = {}
for _, row in edges_df.iterrows():
    adj.setdefault(row["source"], []).append(row.to_dict())

best_score = {node: -np.inf for node in order_nodes}
best_score[START] = 0.0
pred = {}

for node in order_nodes:
    if node not in adj or best_score[node] == -np.inf:
        continue
    for e in adj[node]:
        cand = best_score[node] + e["weighted_score"]
        if cand > best_score[e["target"]]:
            best_score[e["target"]] = cand
            pred[e["target"]] = node

opt_path = reconstruct_path(pred, START, END)
opt_score = best_score[END] if best_score[END] != -np.inf else None
opt_edges = path_to_edge_table(opt_path, edges_df)

if opt_edges is not None and not opt_edges.empty:
    n_positive = int((opt_edges["mean_step"] > 0).sum())
    n_negative = int((opt_edges["mean_step"] < 0).sum())
    n_nonnegative = int((opt_edges["mean_step"] >= 0).sum())
    min_step = float(opt_edges["mean_step"].min())
    mean_step_path = float(opt_edges["mean_step"].mean())
    mean_repro_path = float(opt_edges["log_reproducibility"].mean())
else:
    n_positive = n_negative = n_nonnegative = 0
    min_step = np.nan
    mean_step_path = np.nan
    mean_repro_path = np.nan

summary = {
    "start": START,
    "end": END,
    "path_exists": opt_path is not None,
    "path_length": len(opt_path) - 1 if opt_path is not None else np.nan,
    "cumulative_weighted_score": opt_score if opt_score is not None else np.nan,
    "n_positive_steps": n_positive,
    "n_nonnegative_steps": n_nonnegative,
    "n_negative_steps": n_negative,
    "minimum_step_effect_on_path": min_step,
    "mean_step_effect_on_path": mean_step_path,
    "mean_log_reproducibility_on_path": mean_repro_path,
    "sd_floor": SD_FLOOR,
    "score_definition": "mean_step * log10(|mean_step|/(SD + floor) + 1)",
}
summary_df = pd.DataFrame([summary])

outdir = Path(".")
all_edges_path = outdir / f"{OUTPUT_PREFIX}_all_directed_edges.csv"
path_edges_path = outdir / f"{OUTPUT_PREFIX}_optimal_path_edges.csv"
path_nodes_path = outdir / f"{OUTPUT_PREFIX}_optimal_path_nodes.txt"
summary_path = outdir / f"{OUTPUT_PREFIX}_summary.csv"

edges_df.to_csv(all_edges_path, index=False)
opt_edges.to_csv(path_edges_path, index=False)
summary_df.to_csv(summary_path, index=False)

with open(path_nodes_path, "w") as f:
    if opt_path is None:
        f.write("No Wuhan-to-Omicron path identified.\n")
    else:
        f.write(" -> ".join(opt_path) + "\n")

print("=" * 96)
print("WEIGHTED OPTIMAL PATH ANALYSIS: DDELTA NETWORK")
print("=" * 96)
print(f"Edge score = mean_step * log10(|mean_step| / (SD + {SD_FLOOR}) + 1)")
print()

if opt_path is None:
    print("No Wuhan-to-Omicron path identified.")
else:
    print(f"Path length: {len(opt_path)-1}")
    print(f"Cumulative weighted score: {opt_score:.6f}")
    print(f"Positive steps: {n_positive}")
    print(f"Nonnegative steps: {n_nonnegative}")
    print(f"Negative steps: {n_negative}")
    print(f"Minimum step effect on path: {min_step:.6f}")
    print(f"Mean step effect on path: {mean_step_path:.6f}")
    print(f"Mean log reproducibility on path: {mean_repro_path:.6f}")
    print()
    print("Optimal path:")
    print(" -> ".join(opt_path))
    print()
    print("Path edge summary:")
    print(opt_edges[[
        "source", "target", "mut_source", "mut_target",
        "mean_step", "sd_step", "log_reproducibility", "weighted_score"
    ]].to_string(index=False))

print()
print("=" * 96)
print("FILES WRITTEN")
print("=" * 96)
for p in [all_edges_path, path_edges_path, path_nodes_path, summary_path]:
    print(p.resolve())
