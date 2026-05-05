import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# =========================================================
# SETTINGS
# =========================================================
CSV_PATH = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"

# Node metric (z-axis + node color)
FITNESS_COL = "log2fc_pos_pre"

# Metric used to order genotypes left-to-right within each mutational-load row
ORDER_COL = "log2fc_pos_pre"
ASCENDING = True

# Reproducibility settings
SD_FLOOR = 0.1
DISAGREE_THRESH = 0.30     # hide
BORDERLINE_THRESH = 0.50   # faint below this, stronger above

# Styling
NODE_CMAP = "turbo"
EDGE_CMAP = "coolwarm"     # blue=downhill, white=neutral, red=uphill
POINT_SIZE = 36
EDGE_WIDTH_MIN = 0.4
EDGE_WIDTH_MAX = 3.2

# =========================================================
# 1. LOAD AND CLEAN DATA
# =========================================================
df = pd.read_csv(CSV_PATH)

required = {"library_id", "wo", "mut_count", FITNESS_COL, ORDER_COL}
missing = required - set(df.columns)
if missing:
    raise KeyError(f"Missing required columns: {missing}")

keep_cols = ["library_id", "wo", "mut_count"]
for col in [FITNESS_COL, ORDER_COL]:
    if col not in keep_cols:
        keep_cols.append(col)

df = df[keep_cols].copy()
df = df.replace([np.inf, -np.inf], np.nan)
df[FITNESS_COL] = pd.to_numeric(df[FITNESS_COL], errors="coerce")
df[ORDER_COL] = pd.to_numeric(df[ORDER_COL], errors="coerce")
df["mut_count"] = pd.to_numeric(df["mut_count"], errors="coerce")
df = df.dropna(subset=["library_id", "wo", "mut_count", FITNESS_COL, ORDER_COL]).copy()

# =========================================================
# 2. BUILD CONSENSUS X ORDER WITHIN EACH MUTATIONAL LOAD
# =========================================================
order_df = (
    df.groupby(["wo", "mut_count"], as_index=False)[ORDER_COL]
      .mean()
      .rename(columns={ORDER_COL: "order_value"})
)

order_df = order_df.sort_values(
    by=["mut_count", "order_value", "wo"],
    ascending=[True, ASCENDING, True]
).copy()

order_df["rank_within_load"] = order_df.groupby("mut_count").cumcount()
row_sizes = order_df.groupby("mut_count")["wo"].count().to_dict()

def centered_x(row):
    n = row_sizes[row["mut_count"]]
    return row["rank_within_load"] - (n - 1) / 2.0

order_df["x"] = order_df.apply(centered_x, axis=1)
order_df["y"] = order_df["mut_count"]

# =========================================================
# 3. AVERAGE ACROSS LIBRARIES BY GENOTYPE (for nodes)
# =========================================================
avg_df = (
    df.groupby(["wo", "mut_count"], as_index=False)
      .agg(mean_fitness=(FITNESS_COL, "mean"))
)

avg_df = avg_df.merge(
    order_df[["wo", "mut_count", "x", "y"]],
    on=["wo", "mut_count"],
    how="left"
)

node_pos = dict(zip(avg_df["wo"], zip(avg_df["x"], avg_df["y"], avg_df["mean_fitness"])))

# =========================================================
# 4. EDGE HELPERS
# =========================================================
def hamming_distance_wo(a, b):
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def step_direction(a_row, b_row):
    # orient low mutational load -> high mutational load
    if a_row["mut_count"] < b_row["mut_count"]:
        return a_row, b_row
    elif b_row["mut_count"] < a_row["mut_count"]:
        return b_row, a_row
    else:
        return (a_row, b_row) if a_row["wo"] < b_row["wo"] else (b_row, a_row)

lib_lookup = {}
for lib in sorted(df["library_id"].unique()):
    sub = df[df["library_id"] == lib].copy()
    lib_lookup[lib] = sub.set_index("wo")

# =========================================================
# 5. BUILD ONE-STEP EDGES WITH:
#    - color = accessibility (signed mean step effect)
#    - alpha = reproducibility
#    - width = |mean step effect|
# =========================================================
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
        for lib, lookup in lib_lookup.items():
            if low_g in lookup.index and high_g in lookup.index:
                low_val = lookup.loc[low_g, FITNESS_COL]
                high_val = lookup.loc[high_g, FITNESS_COL]
                step_effects.append(high_val - low_val)

        if len(step_effects) == 0:
            continue

        step_effects = np.asarray(step_effects, dtype=float)
        mean_step = step_effects.mean()
        sd_step = step_effects.std(ddof=1) if len(step_effects) > 1 else 0.0

        raw_reproducibility = np.abs(mean_step) / (sd_step + SD_FLOOR)
        log_reproducibility = np.log10(raw_reproducibility + 1.0)

        edge_records.append({
            "wo_a": low_g,
            "wo_b": high_g,
            "mean_step": mean_step,
            "abs_mean_step": abs(mean_step),
            "sd_step": sd_step,
            "log_reproducibility": log_reproducibility,
        })

edges_df = pd.DataFrame(edge_records)

if edges_df.empty:
    raise RuntimeError("No valid one-step edges were constructed.")

# width scaled by absolute mean step effect
step_min = edges_df["abs_mean_step"].min()
step_max = edges_df["abs_mean_step"].max()
if step_max == step_min:
    edges_df["edge_width"] = (EDGE_WIDTH_MIN + EDGE_WIDTH_MAX) / 2.0
else:
    edges_df["edge_width"] = EDGE_WIDTH_MIN + (
        (edges_df["abs_mean_step"] - step_min) / (step_max - step_min)
    ) * (EDGE_WIDTH_MAX - EDGE_WIDTH_MIN)

# color scaled by signed mean step effect
color_lim = np.nanmax(np.abs(edges_df["mean_step"].values))
edge_norm = colors.TwoSlopeNorm(vmin=-color_lim, vcenter=0.0, vmax=color_lim)
edge_cmap = plt.colormaps[EDGE_CMAP]

# =========================================================
# 6. PLOT 3D NETWORK
# =========================================================
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection="3d")

for _, edge in edges_df.iterrows():
    val = edge["log_reproducibility"]

    # silent at strong disagreement
    if val < DISAGREE_THRESH:
        continue

    # faint at borderline, stronger above
    if val < BORDERLINE_THRESH:
        alpha = 0.10
    else:
        alpha = 0.10 + 0.90 * ((val - BORDERLINE_THRESH) / (edges_df["log_reproducibility"].max() - BORDERLINE_THRESH + 1e-12))
        alpha = min(max(alpha, 0.10), 1.0)

    x1, y1, z1 = node_pos[edge["wo_a"]]
    x2, y2, z2 = node_pos[edge["wo_b"]]
    ax.plot(
        [x1, x2], [y1, y2], [z1, z2],
        color=edge_cmap(edge_norm(edge["mean_step"])),
        alpha=alpha,
        linewidth=edge["edge_width"],
        zorder=1
    )

# nodes
sc = ax.scatter(
    avg_df["x"], avg_df["y"], avg_df["mean_fitness"],
    c=avg_df["mean_fitness"], cmap=NODE_CMAP,
    s=POINT_SIZE, edgecolors="black", linewidths=0.25,
    depthshade=False, zorder=2
)

ax.set_xlabel("Genotypes ordered within mutational load")
ax.set_ylabel("Mutational load")
ax.set_zlabel(FITNESS_COL)
ax.set_title("Average 3D Genotype Network Across Libraries\nEdge color = accessibility, alpha = reproducibility, width = |step effect|")

ax.set_xticks([])
ax.set_yticks(range(0, 9))
ax.set_ylim(8.2, -0.2)  # Wuhan -> Omicron BA.1

# label reference genotypes
labels = {
    "WWWWWWWW": "Wuhan",
    "OOOOOOOO": "Omicron BA.1",
}
for genotype, label in labels.items():
    ref = avg_df[avg_df["wo"] == genotype]
    if not ref.empty:
        row = ref.iloc[0]
        ax.scatter(
            [row["x"]], [row["y"]], [row["mean_fitness"]],
            s=82, c="black", depthshade=False, zorder=3
        )
        ax.text(
            row["x"], row["y"], row["mean_fitness"],
            f"  {label}",
            color="black", zorder=4
        )

# node colorbar
cbar_nodes = fig.colorbar(sc, ax=ax, shrink=0.58, aspect=16, pad=0.08)
cbar_nodes.set_label(FITNESS_COL)

# edge colorbar for accessibility
sm = plt.cm.ScalarMappable(norm=edge_norm, cmap=edge_cmap)
sm.set_array([])
cbar_edges = fig.colorbar(sm, ax=ax, shrink=0.58, aspect=16, pad=0.12)
cbar_edges.set_label("Mean one-step effect (uphill ↔ downhill)")

plt.tight_layout()
plt.show()
