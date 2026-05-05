import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

CSV_PATH = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"

# Metrics to compare
METRICS = {
    "PosPre": "log2fc_pos_pre",
    "DeltaDelta": "ddelta",
}

# Reproducibility settings
SD_FLOOR = 0.1
DISAGREE_THRESH = 0.30
BORDERLINE_THRESH = 0.50

# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------
def hamming_distance_wo(a: str, b: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def step_direction(a_row, b_row):
    # Orient edge from lower mutational load to higher mutational load
    if a_row["mut_count"] < b_row["mut_count"]:
        return a_row, b_row
    elif b_row["mut_count"] < a_row["mut_count"]:
        return b_row, a_row
    else:
        return (a_row, b_row) if a_row["wo"] < b_row["wo"] else (b_row, a_row)

def build_edge_table(df: pd.DataFrame, fitness_col: str) -> pd.DataFrame:
    required = {"library_id", "wo", "mut_count", fitness_col}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns for {fitness_col}: {missing}")

    sub = df[["library_id", "wo", "mut_count", fitness_col]].copy()
    sub = sub.replace([np.inf, -np.inf], np.nan)
    sub[fitness_col] = pd.to_numeric(sub[fitness_col], errors="coerce")
    sub["mut_count"] = pd.to_numeric(sub["mut_count"], errors="coerce")
    sub = sub.dropna(subset=["library_id", "wo", "mut_count", fitness_col]).copy()

    # Average node value across libraries
    avg_df = (
        sub.groupby(["wo", "mut_count"], as_index=False)
           .agg(mean_fitness=(fitness_col, "mean"))
    )

    # Per-library lookup
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

            raw_reproducibility = abs(mean_step) / (sd_step + SD_FLOOR)
            log_reproducibility = np.log10(raw_reproducibility + 1.0)

            edge_records.append({
                "wo_a": low_g,
                "wo_b": high_g,
                "mut_a": float(low_row["mut_count"]),
                "mut_b": float(high_row["mut_count"]),
                "mean_step": mean_step,
                "abs_mean_step": abs(mean_step),
                "sd_step": sd_step,
                "raw_reproducibility": raw_reproducibility,
                "log_reproducibility": log_reproducibility,
                "mean_fitness_a": float(low_row["mean_fitness"]),
                "mean_fitness_b": float(high_row["mean_fitness"]),
                "edge_mean_fitness": float((low_row["mean_fitness"] + high_row["mean_fitness"]) / 2.0),
                "n_libs": len(step_effects),
            })

    edges = pd.DataFrame(edge_records)
    if edges.empty:
        raise RuntimeError(f"No valid one-step edges constructed for {fitness_col}.")
    return edges

def summarize_metric(edges: pd.DataFrame, metric_name: str):
    reproducible = edges[edges["log_reproducibility"] >= DISAGREE_THRESH].copy()
    strongly_repro = edges[edges["log_reproducibility"] >= BORDERLINE_THRESH].copy()

    summary = {
        "metric": metric_name,
        "n_edges_total": len(edges),
        "n_edges_reproducible_0.30": len(reproducible),
        "frac_edges_reproducible_0.30": len(reproducible) / len(edges),
        "n_edges_reproducible_0.50": len(strongly_repro),
        "frac_edges_reproducible_0.50": len(strongly_repro) / len(edges),
        "median_log_repro_all": edges["log_reproducibility"].median(),
        "median_log_repro_reproducible": reproducible["log_reproducibility"].median() if len(reproducible) else np.nan,
        "median_edge_mean_fitness_reproducible": reproducible["edge_mean_fitness"].median() if len(reproducible) else np.nan,
        "median_edge_mean_fitness_all": edges["edge_mean_fitness"].median(),
    }
    return summary, reproducible, strongly_repro

# ---------------------------------------------------------
# Main
# ---------------------------------------------------------
df = pd.read_csv(CSV_PATH)

edge_tables = {}
summaries = []

for metric_name, fitness_col in METRICS.items():
    edges = build_edge_table(df, fitness_col)
    edge_tables[metric_name] = edges
    summary, reproducible, strongly_repro = summarize_metric(edges, metric_name)
    summaries.append(summary)

summary_df = pd.DataFrame(summaries)

# ---------------------------------------------------------
# Statistical comparisons between PosPre and DeltaDelta
# ---------------------------------------------------------
edges_pos = edge_tables["PosPre"]
edges_dd = edge_tables["DeltaDelta"]

# 1) Overall reproducibility distribution comparison
mw_repro = mannwhitneyu(
    edges_pos["log_reproducibility"],
    edges_dd["log_reproducibility"],
    alternative="two-sided"
)

# 2) Among reproducible edges (>= 0.30), compare edge_mean_fitness
pos_repro = edges_pos[edges_pos["log_reproducibility"] >= DISAGREE_THRESH].copy()
dd_repro = edges_dd[edges_dd["log_reproducibility"] >= DISAGREE_THRESH].copy()

if len(pos_repro) > 0 and len(dd_repro) > 0:
    mw_edge_fitness = mannwhitneyu(
        pos_repro["edge_mean_fitness"],
        dd_repro["edge_mean_fitness"],
        alternative="two-sided"
    )
else:
    mw_edge_fitness = None

# 3) Within each metric, ask whether reproducible edges are concentrated at higher fitness
# Compare edge_mean_fitness for reproducible edges vs all edges
mw_pos_internal = mannwhitneyu(
    pos_repro["edge_mean_fitness"],
    edges_pos["edge_mean_fitness"],
    alternative="greater"
) if len(pos_repro) > 0 else None

mw_dd_internal = mannwhitneyu(
    dd_repro["edge_mean_fitness"],
    edges_dd["edge_mean_fitness"],
    alternative="greater"
) if len(dd_repro) > 0 else None

# 4) Correlation between reproducibility and edge fitness
rho_pos, p_pos = spearmanr(edges_pos["log_reproducibility"], edges_pos["edge_mean_fitness"])
rho_dd, p_dd = spearmanr(edges_dd["log_reproducibility"], edges_dd["edge_mean_fitness"])

# ---------------------------------------------------------
# Print results
# ---------------------------------------------------------
print("=" * 88)
print("EDGE-LEVEL SUMMARY")
print("=" * 88)
print(summary_df.to_string(index=False))
print()

print("=" * 88)
print("MANN–WHITNEY COMPARISON OF EDGE REPRODUCIBILITY DISTRIBUTIONS")
print("=" * 88)
print(f"PosPre vs DeltaDelta: U = {mw_repro.statistic:.4f}, p = {mw_repro.pvalue:.4e}")
print()

print("=" * 88)
print("MANN–WHITNEY COMPARISON OF EDGE FITNESS AMONG REPRODUCIBLE EDGES (>= 0.30)")
print("=" * 88)
if mw_edge_fitness is not None:
    print(f"PosPre reproducible vs DeltaDelta reproducible edge_mean_fitness: "
          f"U = {mw_edge_fitness.statistic:.4f}, p = {mw_edge_fitness.pvalue:.4e}")
else:
    print("Insufficient reproducible edges for comparison.")
print()

print("=" * 88)
print("WITHIN-METRIC TEST: ARE REPRODUCIBLE EDGES BIASED TOWARD HIGHER-FITNESS REGIONS?")
print("=" * 88)
if mw_pos_internal is not None:
    print(f"PosPre: reproducible edge_mean_fitness > all edge_mean_fitness: "
          f"U = {mw_pos_internal.statistic:.4f}, p = {mw_pos_internal.pvalue:.4e}")
else:
    print("PosPre: insufficient reproducible edges.")
if mw_dd_internal is not None:
    print(f"DeltaDelta: reproducible edge_mean_fitness > all edge_mean_fitness: "
          f"U = {mw_dd_internal.statistic:.4f}, p = {mw_dd_internal.pvalue:.4e}")
else:
    print("DeltaDelta: insufficient reproducible edges.")
print()

print("=" * 88)
print("SPEARMAN CORRELATION: EDGE REPRODUCIBILITY VS EDGE FITNESS")
print("=" * 88)
print(f"PosPre: rho = {rho_pos:.4f}, p = {p_pos:.4e}")
print(f"DeltaDelta: rho = {rho_dd:.4f}, p = {p_dd:.4e}")
print()

# Optional: write summary tables
summary_out = Path("edge_reproducibility_summary.csv")
summary_df.to_csv(summary_out, index=False)

edges_pos_out = Path("edges_pospre_detailed.csv")
edges_dd_out = Path("edges_ddelta_detailed.csv")
edges_pos.to_csv(edges_pos_out, index=False)
edges_dd.to_csv(edges_dd_out, index=False)

print("=" * 88)
print("FILES WRITTEN")
print("=" * 88)
print(summary_out.resolve())
print(edges_pos_out.resolve())
print(edges_dd_out.resolve())
