import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

CSV_PATH = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"

FITNESS_COL = "log2fc_pos_pre"
ORDER_COL = "log2fc_pos_pre"
ASCENDING = True

EPS = 1e-6
EDGE_WIDTH_MIN = 0.3
EDGE_WIDTH_MAX = 3.0
EDGE_ALPHA = 0.85
EDGE_CMAP = "Greys"
NODE_CMAP = "turbo"
POINT_SIZE = 36

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

# FIXED LINE (no stray quote)
df = df.dropna(subset=["library_id", "wo", "mut_count", FITNESS_COL, ORDER_COL]).copy()

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

def hamming_distance_wo(a, b):
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def step_direction(a_row, b_row):
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
        reproducibility = np.abs(mean_step) / (sd_step + EPS)

        edge_records.append({
            "wo_a": low_g,
            "wo_b": high_g,
            "reproducibility": reproducibility,
        })

edges_df = pd.DataFrame(edge_records)

rep_min = edges_df["reproducibility"].min()
rep_max = edges_df["reproducibility"].max()

edges_df["edge_width"] = EDGE_WIDTH_MIN + (
    (edges_df["reproducibility"] - rep_min) / (rep_max - rep_min + 1e-12)
) * (EDGE_WIDTH_MAX - EDGE_WIDTH_MIN)

rep_norm = colors.Normalize(vmin=rep_min, vmax=rep_max)
edge_cmap = plt.colormaps[EDGE_CMAP]

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection="3d")

for _, edge in edges_df.iterrows():
    x1, y1, z1 = node_pos[edge["wo_a"]]
    x2, y2, z2 = node_pos[edge["wo_b"]]
    ax.plot([x1, x2], [y1, y2], [z1, z2],
            color=edge_cmap(rep_norm(edge["reproducibility"])),
            alpha=EDGE_ALPHA,
            linewidth=edge["edge_width"])

sc = ax.scatter(
    avg_df["x"], avg_df["y"], avg_df["mean_fitness"],
    c=avg_df["mean_fitness"], cmap=NODE_CMAP,
    s=POINT_SIZE, edgecolors="black", linewidths=0.25
)

ax.set_ylim(8.2, -0.2)

cbar_nodes = fig.colorbar(sc, ax=ax)
cbar_nodes.set_label(FITNESS_COL)

sm = plt.cm.ScalarMappable(norm=rep_norm, cmap=edge_cmap)
sm.set_array([])
cbar_edges = fig.colorbar(sm, ax=ax)
cbar_edges.set_label("Edge reproducibility")

plt.show()
