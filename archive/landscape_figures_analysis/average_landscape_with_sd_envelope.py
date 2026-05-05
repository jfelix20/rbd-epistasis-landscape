import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ===========================
# 1. Load data
# ===========================
csv_path = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"
df = pd.read_csv(csv_path)

# ===========================
# 2. Choose fitness metric
# ===========================
FITNESS_COL = "ddelta"

if FITNESS_COL not in df.columns:
    raise KeyError(f"{FITNESS_COL} not found in the dataset.")
if "wo" not in df.columns:
    raise KeyError("'wo' column not found in the dataset.")

# ===========================
# 3. Convert genotype to binary / xy coordinates
# ===========================
def wo_to_binary8(wo_str: str) -> str:
    wo_str = str(wo_str).strip().upper()
    if len(wo_str) != 8:
        raise ValueError(f"Expected 8-character genotype string in 'wo', got: {wo_str}")
    mapping = {"W": "0", "O": "1"}
    try:
        return "".join(mapping[ch] for ch in wo_str)
    except KeyError as e:
        raise ValueError(f"Unexpected character in wo string {wo_str!r}: {e}")

def binary_to_xy(bin_str: str):
    bin_str = str(bin_str).zfill(8)
    x = int(bin_str[:4], 2)
    y = int(bin_str[4:], 2)
    return x, y

df["binary_8"] = df["wo"].apply(wo_to_binary8)
df["x"], df["y"] = zip(*df["binary_8"].apply(binary_to_xy))

# ===========================
# 4. Aggregate mean and between-library SD by genotype
# ===========================
summary_df = (
    df.groupby(["wo", "binary_8", "x", "y"], as_index=False)
      .agg(
          mean_fitness=(FITNESS_COL, "mean"),
          sd_fitness=(FITNESS_COL, "std"),
          n_libraries=("library_id", "nunique")
      )
)

summary_df["sd_fitness"] = summary_df["sd_fitness"].fillna(0.0)
summary_df["upper_fitness"] = summary_df["mean_fitness"] + summary_df["sd_fitness"]
summary_df["lower_fitness"] = summary_df["mean_fitness"] - summary_df["sd_fitness"]

# ===========================
# 5. Interpolate onto a regular grid
# ===========================
x_min, x_max = summary_df["x"].min(), summary_df["x"].max()
y_min, y_max = summary_df["y"].min(), summary_df["y"].max()

grid_x, grid_y = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]

def make_grid(value_col: str):
    grid = griddata(
        points=summary_df[["x", "y"]].values,
        values=summary_df[value_col].values,
        xi=(grid_x, grid_y),
        method="linear"
    )
    return np.ma.array(grid, mask=np.isnan(grid))

grid_mean = make_grid("mean_fitness")
grid_upper = make_grid("upper_fitness")
grid_lower = make_grid("lower_fitness")

# ===========================
# 6. Plot mean surface with thin SD envelope
# ===========================
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# Main mean surface
surf_mean = ax.plot_surface(
    grid_x,
    grid_y,
    grid_mean,
    cmap="viridis",
    linewidth=0,
    antialiased=True,
    alpha=0.95
)

# Thin "cloud" / envelope surfaces
ax.plot_surface(
    grid_x,
    grid_y,
    grid_upper,
    linewidth=0,
    antialiased=True,
    alpha=0.18
)

ax.plot_surface(
    grid_x,
    grid_y,
    grid_lower,
    linewidth=0,
    antialiased=True,
    alpha=0.18
)

ax.set_xlabel("Positions 477/478/484/493")
ax.set_ylabel("Positions 496/498/501/505")
ax.set_zlabel(f"Mean {FITNESS_COL}")
ax.set_title("Average Experimental Fitness Landscape\nwith Between-Library SD Envelope")

cbar = fig.colorbar(surf_mean, shrink=0.6, aspect=14, pad=0.08)
cbar.set_label(f"Mean {FITNESS_COL}")

plt.tight_layout()
plt.show()

print("Generated mean landscape with thin between-library SD envelope.")
