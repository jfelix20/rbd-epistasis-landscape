import pandas as pd
import numpy as np
from pathlib import Path

# =========================================================
# SETTINGS
# =========================================================
CSV_PATH = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"

FITNESS_COL = "ddelta"
SD_FLOOR = 0.1
REPRO_THRESH = 0.30
OUTPUT_PREFIX = "ddelta_path_diagnostics"

START = "WWWWWWWW"
END = "OOOOOOOO"

# =========================================================
# HELPERS
# =========================================================
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

# =========================================================
# LOAD DATA
# =========================================================
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

# =========================================================
# BUILD DIRECTED ONE-STEP EDGES
# =========================================================
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
        reliable_uphill_score = max(mean_step, 0.0) * log_repro

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

# =========================================================
# BASIC ENDPOINT DIAGNOSTICS
# =========================================================
outgoing_start_all = edges_df[edges_df["source"] == START].copy()
outgoing_start_filtered = filtered_edges[filtered_edges["source"] == START].copy()

incoming_end_all = edges_df[edges_df["target"] == END].copy()
incoming_end_filtered = filtered_edges[filtered_edges["target"] == END].copy()

# =========================================================
# REACHABILITY FROM START UNDER CURRENT FILTER
# =========================================================
adj = {}
for _, row in filtered_edges.iterrows():
    adj.setdefault(row["source"], []).append(row.to_dict())

reachable = set([START])
frontier = [START]

while frontier:
    current = frontier.pop()
    for e in adj.get(current, []):
        tgt = e["target"]
        if tgt not in reachable:
            reachable.add(tgt)
            frontier.append(tgt)

reachable_df = avg_df[avg_df["wo"].isin(reachable)].copy()

# Furthest reachable by mutational load
max_reachable_load = reachable_df["mut_count"].max()
furthest_nodes = reachable_df[reachable_df["mut_count"] == max_reachable_load].copy()
furthest_nodes["hamming_to_end"] = furthest_nodes["wo"].apply(lambda x: hamming_distance_wo(x, END))
furthest_nodes = furthest_nodes.sort_values(["hamming_to_end", "mean_fitness"], ascending=[True, False])

# =========================================================
# BEST PARTIAL PATH
# Rank reachable endpoints by:
# 1) greatest mutational load reached
# 2) smallest Hamming distance to END
# 3) highest cumulative reliable-uphill score from START
# =========================================================
order_nodes = avg_df.sort_values(["mut_count", "wo"])["wo"].tolist()

best_sum = {node: -np.inf for node in order_nodes}
best_sum[START] = 0.0
pred_sum = {}

for node in order_nodes:
    if node not in adj or best_sum[node] == -np.inf:
        continue
    for e in adj[node]:
        cand = best_sum[node] + e["reliable_uphill_score"]
        if cand > best_sum[e["target"]]:
            best_sum[e["target"]] = cand
            pred_sum[e["target"]] = node

reachable_summary = reachable_df.copy()
reachable_summary["best_sum_score"] = reachable_summary["wo"].map(best_sum)
reachable_summary["hamming_to_end"] = reachable_summary["wo"].apply(lambda x: hamming_distance_wo(x, END))
reachable_summary = reachable_summary[reachable_summary["best_sum_score"] > -np.inf].copy()

best_partial_target = (
    reachable_summary.sort_values(
        ["mut_count", "hamming_to_end", "best_sum_score"],
        ascending=[False, True, False]
    )
    .iloc[0]["wo"]
)

best_partial_path = reconstruct_path(pred_sum, START, best_partial_target)
best_partial_edges = path_to_edge_table(best_partial_path, filtered_edges)

# =========================================================
# SUMMARY TABLE
# =========================================================
summary = {
    "n_total_directed_edges": len(edges_df),
    "n_filtered_uphill_repro_edges": len(filtered_edges),
    "fraction_filtered": len(filtered_edges) / len(edges_df),
    "reproducibility_threshold": REPRO_THRESH,
    "n_outgoing_from_start_all": len(outgoing_start_all),
    "n_outgoing_from_start_filtered": len(outgoing_start_filtered),
    "n_incoming_to_end_all": len(incoming_end_all),
    "n_incoming_to_end_filtered": len(incoming_end_filtered),
    "n_reachable_nodes_from_start": len(reachable),
    "max_reachable_mut_load": int(max_reachable_load),
    "best_partial_target": best_partial_target,
    "best_partial_target_mut_load": int(reachable_summary[reachable_summary["wo"] == best_partial_target]["mut_count"].iloc[0]),
    "best_partial_target_hamming_to_end": int(reachable_summary[reachable_summary["wo"] == best_partial_target]["hamming_to_end"].iloc[0]),
    "best_partial_sum_score": float(reachable_summary[reachable_summary["wo"] == best_partial_target]["best_sum_score"].iloc[0]),
    "full_path_exists": END in reachable,
}
summary_df = pd.DataFrame([summary])

# =========================================================
# WRITE OUTPUTS
# =========================================================
outdir = Path(".")
summary_path = outdir / f"{OUTPUT_PREFIX}_summary.csv"
start_all_path = outdir / f"{OUTPUT_PREFIX}_start_outgoing_all.csv"
start_filtered_path = outdir / f"{OUTPUT_PREFIX}_start_outgoing_filtered.csv"
end_all_path = outdir / f"{OUTPUT_PREFIX}_end_incoming_all.csv"
end_filtered_path = outdir / f"{OUTPUT_PREFIX}_end_incoming_filtered.csv"
reachable_path = outdir / f"{OUTPUT_PREFIX}_reachable_nodes.csv"
furthest_path = outdir / f"{OUTPUT_PREFIX}_furthest_reachable_nodes.csv"
best_partial_edges_path = outdir / f"{OUTPUT_PREFIX}_best_partial_path_edges.csv"
best_partial_nodes_path = outdir / f"{OUTPUT_PREFIX}_best_partial_path_nodes.txt"

summary_df.to_csv(summary_path, index=False)
outgoing_start_all.to_csv(start_all_path, index=False)
outgoing_start_filtered.to_csv(start_filtered_path, index=False)
incoming_end_all.to_csv(end_all_path, index=False)
incoming_end_filtered.to_csv(end_filtered_path, index=False)
reachable_summary.to_csv(reachable_path, index=False)
furthest_nodes.to_csv(furthest_path, index=False)
best_partial_edges.to_csv(best_partial_edges_path, index=False)

with open(best_partial_nodes_path, "w") as f:
    if best_partial_path is None:
        f.write("No partial path found.\n")
    else:
        f.write(" -> ".join(best_partial_path) + "\n")

# =========================================================
# PRINT REPORT
# =========================================================
print("=" * 96)
print("DDELTA PATH DIAGNOSTICS")
print("=" * 96)
print(f"Total directed one-step edges: {len(edges_df)}")
print(f"Filtered uphill + reproducible edges: {len(filtered_edges)}")
print(f"Fraction retained: {len(filtered_edges) / len(edges_df):.4f}")
print(f"Reproducibility threshold: {REPRO_THRESH}")
print()

print("-" * 96)
print("START NODE (WUHAN) DIAGNOSTICS")
print("-" * 96)
print(f"Outgoing edges from {START} (all): {len(outgoing_start_all)}")
print(f"Outgoing edges from {START} (filtered): {len(outgoing_start_filtered)}")
print()

print("-" * 96)
print("END NODE (OMICRON BA.1) DIAGNOSTICS")
print("-" * 96)
print(f"Incoming edges to {END} (all): {len(incoming_end_all)}")
print(f"Incoming edges to {END} (filtered): {len(incoming_end_filtered)}")
print()

print("-" * 96)
print("REACHABILITY")
print("-" * 96)
print(f"Reachable nodes from {START}: {len(reachable)}")
print(f"Maximum reachable mutational load: {int(max_reachable_load)}")
print(f"Full path to {END} exists under filter: {END in reachable}")
print()

print("-" * 96)
print("BEST PARTIAL TARGET")
print("-" * 96)
print(f"Best reachable target: {best_partial_target}")
print(f"Mutational load: {int(reachable_summary[reachable_summary['wo'] == best_partial_target]['mut_count'].iloc[0])}")
print(f"Hamming distance to {END}: {int(reachable_summary[reachable_summary['wo'] == best_partial_target]['hamming_to_end'].iloc[0])}")
print(f"Cumulative score: {float(reachable_summary[reachable_summary['wo'] == best_partial_target]['best_sum_score'].iloc[0]):.6f}")
print()
if best_partial_path is not None:
    print("Best partial path:")
    print(" -> ".join(best_partial_path))
print()

print("=" * 96)
print("FILES WRITTEN")
print("=" * 96)
for p in [
    summary_path, start_all_path, start_filtered_path, end_all_path, end_filtered_path,
    reachable_path, furthest_path, best_partial_edges_path, best_partial_nodes_path
]:
    print(p.resolve())
