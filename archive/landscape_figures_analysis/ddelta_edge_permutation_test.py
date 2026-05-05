import pandas as pd
import numpy as np
from pathlib import Path

# =========================================================
# SETTINGS
# =========================================================
CSV_PATH = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"

FITNESS_COL = "ddelta"
SD_FLOOR = 0.1
THRESHOLDS = [0.30, 0.50]
N_PERMUTATIONS = 1000
RANDOM_SEED = 42
OUTPUT_PREFIX = "ddelta_edge_permutation_test"

# =========================================================
# HELPERS
# =========================================================
def hamming_distance_wo(a: str, b: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def step_direction(a_row, b_row):
    # Orient edges from lower mutational load -> higher mutational load
    if a_row["mut_count"] < b_row["mut_count"]:
        return a_row, b_row
    elif b_row["mut_count"] < a_row["mut_count"]:
        return b_row, a_row
    else:
        return (a_row, b_row) if a_row["wo"] < b_row["wo"] else (b_row, a_row)

def build_edge_table(df: pd.DataFrame, fitness_col: str, sd_floor: float) -> pd.DataFrame:
    """
    Build one-step edges and compute reproducibility score:
    log10(|mean step| / (SD across libraries + sd_floor) + 1)
    """
    sub = df[["library_id", "wo", "mut_count", fitness_col]].copy()
    sub = sub.replace([np.inf, -np.inf], np.nan)
    sub[fitness_col] = pd.to_numeric(sub[fitness_col], errors="coerce")
    sub["mut_count"] = pd.to_numeric(sub["mut_count"], errors="coerce")
    sub = sub.dropna(subset=["library_id", "wo", "mut_count", fitness_col]).copy()

    avg_df = (
        sub.groupby(["wo", "mut_count"], as_index=False)
           .agg(mean_fitness=(fitness_col, "mean"))
    )

    lib_lookup = {}
    for lib in sorted(sub["library_id"].unique()):
        lib_lookup[lib] = sub[sub["library_id"] == lib].set_index("wo")

    genotypes = sorted(avg_df["wo"].unique())
    edge_records = []

    for i, g1 in enumerate(genotypes):
        for g2 in genotypes[i + 1:]:
            if hamming_distance_wo(g1, g2) != 1:
                continue

            row1 = avg_df.loc[avg_df["wo"] == g1].iloc[0]
            row2 = avg_df.loc[avg_df["wo"] == g2].iloc[0]
            low_row, high_row = step_direction(row1, row2)
            low_g, high_g = low_row["wo"], high_row["wo"]

            step_effects = []
            for _, lookup in lib_lookup.items():
                if low_g in lookup.index and high_g in lookup.index:
                    low_val = lookup.loc[low_g, fitness_col]
                    high_val = lookup.loc[high_g, fitness_col]
                    step_effects.append(high_val - low_val)

            if len(step_effects) == 0:
                continue

            step_effects = np.asarray(step_effects, dtype=float)
            mean_step = float(step_effects.mean())
            sd_step = float(step_effects.std(ddof=1)) if len(step_effects) > 1 else 0.0
            raw_repro = abs(mean_step) / (sd_step + sd_floor)
            log_repro = np.log10(raw_repro + 1.0)

            edge_records.append({
                "wo_a": low_g,
                "wo_b": high_g,
                "mut_a": float(low_row["mut_count"]),
                "mut_b": float(high_row["mut_count"]),
                "mean_step": mean_step,
                "sd_step": sd_step,
                "raw_reproducibility": raw_repro,
                "log_reproducibility": log_repro,
                "edge_mean_fitness": float((low_row["mean_fitness"] + high_row["mean_fitness"]) / 2.0),
            })

    edges = pd.DataFrame(edge_records)
    if edges.empty:
        raise RuntimeError("No valid one-step edges were constructed.")
    return edges

def summarize_edges(edges: pd.DataFrame, thresholds):
    out = {
        "n_edges_total": len(edges),
        "median_log_reproducibility": float(edges["log_reproducibility"].median()),
        "max_log_reproducibility": float(edges["log_reproducibility"].max()),
        "mean_log_reproducibility": float(edges["log_reproducibility"].mean()),
    }
    for t in thresholds:
        out[f"n_edges_ge_{t:.2f}"] = int((edges["log_reproducibility"] >= t).sum())
        out[f"frac_edges_ge_{t:.2f}"] = float((edges["log_reproducibility"] >= t).mean())
    return out

# =========================================================
# MAIN
# =========================================================
rng = np.random.default_rng(RANDOM_SEED)

df = pd.read_csv(CSV_PATH)

required = {"library_id", "wo", "mut_count", FITNESS_COL}
missing = required - set(df.columns)
if missing:
    raise KeyError(f"Missing required columns: {missing}")

# Observed data
obs_edges = build_edge_table(df, FITNESS_COL, SD_FLOOR)
obs_summary = summarize_edges(obs_edges, THRESHOLDS)

# Permutation null:
# shuffle ddelta values within each library, preserving the per-library value distribution
perm_records = []

for perm_idx in range(N_PERMUTATIONS):
    perm_df = df.copy()

    shuffled_parts = []
    for lib, sub in perm_df.groupby("library_id", sort=True):
        sub = sub.copy()
        vals = sub[FITNESS_COL].to_numpy().copy()
        rng.shuffle(vals)
        sub[FITNESS_COL] = vals
        shuffled_parts.append(sub)

    perm_df = pd.concat(shuffled_parts, ignore_index=True)

    perm_edges = build_edge_table(perm_df, FITNESS_COL, SD_FLOOR)
    perm_summary = summarize_edges(perm_edges, THRESHOLDS)
    perm_summary["permutation"] = perm_idx + 1
    perm_records.append(perm_summary)

perm_stats = pd.DataFrame(perm_records)

# Empirical p-values:
# probability null is >= observed
results = []
for key, obs_value in obs_summary.items():
    if key == "n_edges_total":
        continue
    null_vals = perm_stats[key].to_numpy()
    p_empirical = (np.sum(null_vals >= obs_value) + 1) / (len(null_vals) + 1)
    results.append({
        "statistic": key,
        "observed_value": obs_value,
        "null_mean": float(np.mean(null_vals)),
        "null_median": float(np.median(null_vals)),
        "null_sd": float(np.std(null_vals, ddof=1)),
        "null_min": float(np.min(null_vals)),
        "null_max": float(np.max(null_vals)),
        "empirical_p_ge_observed": float(p_empirical),
    })

results_df = pd.DataFrame(results)

# Write outputs
outdir = Path(".")
obs_edges_path = outdir / f"{OUTPUT_PREFIX}_observed_edges.csv"
perm_stats_path = outdir / f"{OUTPUT_PREFIX}_permutation_summary.csv"
results_path = outdir / f"{OUTPUT_PREFIX}_observed_vs_null.csv"

obs_edges.to_csv(obs_edges_path, index=False)
perm_stats.to_csv(perm_stats_path, index=False)
results_df.to_csv(results_path, index=False)

# Print concise report
print("=" * 96)
print("OBSERVED DDELTA EDGE REPRODUCIBILITY SUMMARY")
print("=" * 96)
for k, v in obs_summary.items():
    print(f"{k}: {v}")
print()

print("=" * 96)
print(f"PERMUTATION TEST RESULTS ({N_PERMUTATIONS} shuffles)")
print("=" * 96)
print(results_df.to_string(index=False))
print()

print("=" * 96)
print("FILES WRITTEN")
print("=" * 96)
print(obs_edges_path.resolve())
print(perm_stats_path.resolve())
print(results_path.resolve())
print()

print("Interpretation guide:")
print("- Very small empirical p-values mean the observed reproducibility structure exceeds chance.")
print("- If counts above threshold (e.g., >=0.30 or >=0.50) are far above null, stochasticity alone is unlikely.")
