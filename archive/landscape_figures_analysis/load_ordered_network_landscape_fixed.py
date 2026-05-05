import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# SETTINGS
# =========================================================
CSV_PATH = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"

# Metric used for node coloring in each library panel
COLOR_COL = "log2fc_pos_pre"

# Metric used to order genotypes within each mutational-load row.
# Using a cross-library mean keeps the x-order consistent across panels.
ORDER_COL = "log2fc_pos_pre"

# Sort direction within each mutational-load row:
# True  -> lowest on left, highest on right
# False -> highest on left, lowest on right
ASCENDING = True

# Colormap and point style
CMAP = "turbo"   # change to "cividis" for colorblind-safe
POINT_SIZE = 40

# =========================================================
# 1. LOAD DATA
# =========================================================
df = pd.read_csv(CSV_PATH)

required = {"library_id", "wo", "mut_count", COLOR_COL, ORDER_COL}
missing = required - set(df.columns)
if missing:
    raise KeyError(f"Missing required columns: {missing}")

# Avoid duplicate columns if COLOR_COL == ORDER_COL
keep_cols = ["library_id", "wo", "mut_count"]
for col in [COLOR_COL, ORDER_COL]:
    if col not in keep_cols:
        keep_cols.append(col)

df = df[keep_cols].copy()
df = df.replace([np.inf, -np.inf], np.nan)

# Make sure plotting/ordering columns are numeric
df[COLOR_COL] = pd.to_numeric(df[COLOR_COL], errors="coerce")
df[ORDER_COL] = pd.to_numeric(df[ORDER_COL], errors="coerce")

df = df.dropna(subset=["library_id", "wo", "mut_count", COLOR_COL, ORDER_COL]).copy()

# =========================================================
# 2. BUILD A CONSENSUS ORDER WITHIN EACH MUTATIONAL LOAD
#    so all library panels use the same left-to-right layout
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

pos = dict(zip(order_df["wo"], zip(order_df["x"], order_df["y"])))

# =========================================================
# 3. EDGES: connect genotypes differing by exactly one W/O bit
# =========================================================
def hamming_distance_wo(a: str, b: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def get_edges_for_present_genotypes(df_lib):
    genotypes = sorted(set(df_lib["wo"].tolist()))
    edges = []
    for i, g1 in enumerate(genotypes):
        for g2 in genotypes[i + 1:]:
            if hamming_distance_wo(g1, g2) == 1:
                edges.append((g1, g2))
    return edges

# =========================================================
# 4. COLOR SCALE
# =========================================================
vmin = df[COLOR_COL].quantile(0.05)
vmax = df[COLOR_COL].quantile(0.95)

# =========================================================
# 5. PLOT
# =========================================================
libraries = sorted(df["library_id"].unique())

fig, axes = plt.subplots(
    1, len(libraries),
    figsize=(5 * len(libraries), 6),
    sharey=True
)

if len(libraries) == 1:
    axes = [axes]

last_scatter = None
max_half_width = max((n - 1) / 2.0 for n in row_sizes.values())

for ax, lib in zip(axes, libraries):
    df_lib = df[df["library_id"] == lib].copy()

    edges = get_edges_for_present_genotypes(df_lib)
    for wo_a, wo_b in edges:
        x1, y1 = pos[wo_a]
        x2, y2 = pos[wo_b]
        ax.plot(
            [x1, x2], [y1, y2],
            color="lightgray",
            linewidth=0.5,
            alpha=0.35,
            zorder=1
        )

    xs = df_lib["wo"].map(lambda w: pos[w][0]).values
    ys = df_lib["wo"].map(lambda w: pos[w][1]).values
    cvals = df_lib[COLOR_COL].values

    last_scatter = ax.scatter(
        xs, ys,
        c=cvals,
        cmap=CMAP,
        vmin=vmin,
        vmax=vmax,
        s=POINT_SIZE,
        edgecolors="black",
        linewidths=0.3,
        zorder=2
    )

    ax.set_title(f"Library {lib}")
    ax.set_xlabel("Genotypes ordered within mutational load")
    ax.set_xticks([])
    ax.set_ylim(-0.5, 8.5)
    ax.set_yticks(range(0, 9))
    ax.set_yticklabels(range(0, 9))
    ax.set_xlim(-max_half_width - 2, max_half_width + 2)

    if ax is axes[0]:
        ax.set_ylabel("Mutational load (Omicron-matching sites)")

cbar = fig.colorbar(
    last_scatter,
    ax=axes[-1],
    location="right",
    shrink=0.8,
    pad=0.05
)
cbar.set_label(f"{COLOR_COL}")

plt.tight_layout()
plt.show()

print("Plotted load-ordered genotype network for all libraries.")
print(f"Ordering metric: {ORDER_COL} (cross-library mean)")
print(f"Color metric:    {COLOR_COL}")
print(f"Ascending order: {ASCENDING} -> highest values on the {'right' if ASCENDING else 'left'}")
