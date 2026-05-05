import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# =========================================================
# SETTINGS
# =========================================================
CSV_PATH = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"

# Metric used for the z-axis and color
FITNESS_COL = "log2fc_pos_pre"

# Metric used to order genotypes within each mutational-load row
# (cross-library mean of this metric determines left-to-right position)
ORDER_COL = "log2fc_pos_pre"

# Sort direction within each mutational-load row:
# True  -> lowest on left, highest on right
# False -> highest on left, lowest on right
ASCENDING = True

CMAP = "turbo"

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
df = df.dropna(subset=["wo", "mut_count", FITNESS_COL, ORDER_COL]).copy()

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
# 3. AVERAGE ACROSS LIBRARIES BY GENOTYPE
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

# =========================================================
# 4. INTERPOLATE A 3D SURFACE
# =========================================================
x_min, x_max = avg_df["x"].min(), avg_df["x"].max()
y_min, y_max = avg_df["y"].min(), avg_df["y"].max()

grid_x, grid_y = np.mgrid[x_min:x_max:220j, y_min:y_max:220j]

grid_z = griddata(
    points=avg_df[["x", "y"]].values,
    values=avg_df["mean_fitness"].values,
    xi=(grid_x, grid_y),
    method="linear"
)

grid_z = np.ma.array(grid_z, mask=np.isnan(grid_z))

# =========================================================
# 5. PLOT 3D LANDSCAPE
# =========================================================
fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection="3d")

surf = ax.plot_surface(
    grid_x,
    grid_y,
    grid_z,
    cmap=CMAP,
    linewidth=0,
    antialiased=True,
    alpha=0.96
)

ax.set_xlabel("Genotypes ordered within mutational load")
ax.set_ylabel("Mutational load")
ax.set_zlabel(FITNESS_COL)
ax.set_title("Average 3D Landscape Across Libraries")

# cleaner x axis
ax.set_xticks([])
ax.set_yticks(range(0, 9))
ax.set_ylim(8.2, -0.2)  # invert to run from Wuhan (0) to Omicron (8)

# Label reference genotypes
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
            s=70, c="black", depthshade=False
        )
        ax.text(
            row["x"], row["y"], row["mean_fitness"],
            f"  {label}",
            color="black"
        )

cbar = fig.colorbar(surf, shrink=0.65, aspect=16, pad=0.08)
cbar.set_label(FITNESS_COL)

plt.tight_layout()
plt.show()

print("Plotted averaged 3D landscape with z-axis =", FITNESS_COL)
print("Mutational load axis inverted: Wuhan -> Omicron BA.1")
print(f"Ordering metric: {ORDER_COL} (cross-library mean)")
print(f"Ascending order: {ASCENDING} -> highest values on the {'right' if ASCENDING else 'left'}")
