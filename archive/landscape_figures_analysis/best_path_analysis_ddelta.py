import pandas as pd
import numpy as np
from pathlib import Path

CSV_PATH = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"

FITNESS_COL = "ddelta"
SD_FLOOR = 0.1
REPRO_THRESH = 0.30
OUTPUT_PREFIX = "best_path_analysis_ddelta"

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

if "WWWWWWWW" not in set(avg_df["wo"]) or "OOOOOOOO" not in set(avg_df["wo"]):
    raise RuntimeError("Wuhan and/or Omicron BA.1 genotype missing from averaged dataframe.")

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
        uphill_component = max(mean_step, 0.0)
        reliable_uphill_score = uphill_component * log_repro

        edge_records.append({
            "source": low_g,
            "target": high_g,
            "mut_source": int(low_row["mut_count"]),
            "mut_target": int(high_row["mut_count"]),
            "mean_step": mean_step,
            "sd_step": sd_step,
            "raw_reproducibility": raw_repro,
            "log_reproducibility": log_repro,
            "reliable_uphill_score": reliable_uphill_score,
            "source_mean_fitness": float(low_row["mean_fitness"]),
            "target_mean_fitness": float(high_row["mean_fitness"]),
        })

edges_df = pd.DataFrame(edge_records)

if edges_df.empty:
    raise RuntimeError("No valid directed one-step edges were constructed.")

filtered_edges = edges_df[
    (edges_df["mean_step"] > 0) &
    (edges_df["log_reproducibility"] >= REPRO_THRESH) &
    (edges_df["reliable_uphill_score"] > 0)
].copy()

adj = {}
for _, row in filtered_edges.iterrows():
    adj.setdefault(row["source"], []).append(row.to_dict())

start = "WWWWWWWW"
end = "OOOOOOOO"
order_nodes = avg_df.sort_values(["mut_count", "wo"])["wo"].tolist()

best_sum = {node: -np.inf for node in order_nodes}
best_sum[start] = 0.0
pred_sum = {}

for node in order_nodes:
    if node not in adj or best_sum[node] == -np.inf:
        continue
    for e in adj[node]:
        cand = best_sum[node] + e["reliable_uphill_score"]
        if cand > best_sum[e["target"]]:
            best_sum[e["target"]] = cand
            pred_sum[e["target"]] = node

sum_path = reconstruct_path(pred_sum, start, end)
sum_score = best_sum[end] if best_sum[end] != -np.inf else None

best_bottleneck = {node: -np.inf for node in order_nodes}
best_bottleneck[start] = np.inf
pred_bottle = {}

for node in order_nodes:
    if node not in adj or best_bottleneck[node] == -np.inf:
        continue
    for e in adj[node]:
        cand = min(best_bottleneck[node], e["reliable_uphill_score"])
        if cand > best_bottleneck[e["target"]]:
            best_bottleneck[e["target"]] = cand
            pred_bottle[e["target"]] = node

bottle_path = reconstruct_path(pred_bottle, start, end)
bottle_score = best_bottleneck[end] if best_bottleneck[end] != -np.inf else None

def path_to_edge_table(path, edges_lookup_df):
    if path is None or len(path) < 2:
        return pd.DataFrame()
    recs = []
    for a, b in zip(path[:-1], path[1:]):
        hit = edges_lookup_df[(edges_lookup_df["source"] == a) & (edges_lookup_df["target"] == b)]
        if not hit.empty:
            recs.append(hit.iloc[0].to_dict())
    return pd.DataFrame(recs)

sum_path_edges = path_to_edge_table(sum_path, filtered_edges)
bottle_path_edges = path_to_edge_table(bottle_path, filtered_edges)

summary = {
    "n_total_directed_one_step_edges": len(edges_df),
    "n_filtered_uphill_repro_edges": len(filtered_edges),
    "fraction_filtered": len(filtered_edges) / len(edges_df),
    "reproducibility_threshold": REPRO_THRESH,
    "sd_floor": SD_FLOOR,
    "sum_path_exists": sum_path is not None,
    "sum_path_score": sum_score if sum_score is not None else np.nan,
    "sum_path_length": len(sum_path) - 1 if sum_path is not None else np.nan,
    "bottleneck_path_exists": bottle_path is not None,
    "bottleneck_path_score": bottle_score if bottle_score is not None else np.nan,
    "bottleneck_path_length": len(bottle_path) - 1 if bottle_path is not None else np.nan,
}
summary_df = pd.DataFrame([summary])

outdir = Path(".")
edges_all_path = outdir / f"{OUTPUT_PREFIX}_all_directed_edges.csv"
edges_filtered_path = outdir / f"{OUTPUT_PREFIX}_filtered_edges.csv"
summary_path = outdir / f"{OUTPUT_PREFIX}_summary.csv"
sum_edges_path = outdir / f"{OUTPUT_PREFIX}_max_cumulative_path_edges.csv"
bottle_edges_path = outdir / f"{OUTPUT_PREFIX}_max_bottleneck_path_edges.csv"
sum_nodes_path = outdir / f"{OUTPUT_PREFIX}_max_cumulative_path_nodes.txt"
bottle_nodes_path = outdir / f"{OUTPUT_PREFIX}_max_bottleneck_path_nodes.txt"

edges_df.to_csv(edges_all_path, index=False)
filtered_edges.to_csv(edges_filtered_path, index=False)
summary_df.to_csv(summary_path, index=False)
sum_path_edges.to_csv(sum_edges_path, index=False)
bottle_path_edges.to_csv(bottle_edges_path, index=False)

with open(sum_nodes_path, "w") as f:
    if sum_path is None:
        f.write("No full Wuhan-to-Omicron path found under current filter.\n")
    else:
        f.write(" -> ".join(sum_path) + "\n")

with open(bottle_nodes_path, "w") as f:
    if bottle_path is None:
        f.write("No full Wuhan-to-Omicron path found under current filter.\n")
    else:
        f.write(" -> ".join(bottle_path) + "\n")

print("=" * 96)
print("BEST PATH ANALYSIS: DDELTA NETWORK")
print("=" * 96)
print(f"Total directed one-step edges: {len(edges_df)}")
print(f"Filtered uphill + reproducible edges: {len(filtered_edges)}")
print(f"Fraction retained: {len(filtered_edges)/len(edges_df):.4f}")
print(f"Reproducibility threshold: {REPRO_THRESH}")
print(f"SD floor: {SD_FLOOR}")
print()
print("-" * 96)
print("MAXIMUM CUMULATIVE SCORE PATH")
print("-" * 96)
if sum_path is None:
    print("No full Wuhan-to-Omicron path found under current filter.")
else:
    print(f"Path length: {len(sum_path)-1}")
    print(f"Total cumulative score: {sum_score:.6f}")
    print("Path:")
    print(" -> ".join(sum_path))
print()
print("-" * 96)
print("MAXIMUM BOTTLENECK PATH")
print("-" * 96)
if bottle_path is None:
    print("No full Wuhan-to-Omicron path found under current filter.")
else:
    print(f"Path length: {len(bottle_path)-1}")
    print(f"Bottleneck score: {bottle_score:.6f}")
    print("Path:")
    print(" -> ".join(bottle_path))
print()
print("=" * 96)
print("FILES WRITTEN")
print("=" * 96)
for p in [edges_all_path, edges_filtered_path, summary_path, sum_edges_path, bottle_edges_path, sum_nodes_path, bottle_nodes_path]:
    print(p.resolve())
