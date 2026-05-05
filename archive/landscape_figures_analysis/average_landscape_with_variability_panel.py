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

# with 3 libraries this should usually be defined, but fill any NaNs just in case
summary_df["sd_fitness"] = summary_df["sd_fitness"].fillna(0.0)

# ===========================
# 5. Interpolate onto a regular grid
# ===========================
x_min, x_max = summary_df["x"].min(), summary_df["x"].max()
y_min, y_max = summary_df["y"].min(), summary_df["y"].max()

grid_x, grid_y = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]

grid_mean = griddata(
    points=summary_df[["x", "y"]].values,
    values=summary_df["mean_fitness"].values,
    xi=(grid_x, grid_y),
    method="linear"
)

grid_sd = griddata(
    points=summary_df[["x", "y"]].values,
    values=summary_df["sd_fitness"].values,
    xi=(grid_x, grid_y),
    method="linear"
)

grid_mean = np.ma.array(grid_mean, mask=np.isnan(grid_mean))
grid_sd = np.ma.array(grid_sd, mask=np.isnan(grid_sd))

# ===========================
# 6. Plot: mean landscape + variability panel
# ===========================
fig = plt.figure(figsize=(14, 6))

# ---- Panel A: mean 3D surface ----
ax1 = fig.add_subplot(1, 2, 1, projection="3d")
surf = ax1.plot_surface(
    grid_x,
    grid_y,
    grid_mean,
    cmap="viridis",
    linewidth=0,
    antialiased=True,
    alpha=0.95
)

ax1.set_xlabel("Positions 477/478/484/493")
ax1.set_ylabel("Positions 496/498/501/505")
ax1.set_zlabel(f"Mean {FITNESS_COL}")
ax1.set_title("Average Experimental Fitness Landscape")

cbar1 = fig.colorbar(surf, ax=ax1, shrink=0.65, aspect=14, pad=0.08)
cbar1.set_label(f"Mean {FITNESS_COL}")

# ---- Panel B: between-library SD heatmap ----
ax2 = fig.add_subplot(1, 2, 2)
im = ax2.imshow(
    grid_sd.T,
    origin="lower",
    extent=(x_min, x_max, y_min, y_max),
    aspect="equal",
    cmap="magma"
)

ax2.set_xlabel("Positions 477/478/484/493")
ax2.set_ylabel("Positions 496/498/501/505")
ax2.set_title("Between-Library Variability (SD)")

# genotype vertex ticks
ax2.set_xticks(range(int(x_min), int(x_max) + 1, 3))
ax2.set_yticks(range(int(y_min), int(y_max) + 1, 3))

cbar2 = fig.colorbar(im, ax=ax2, shrink=0.85, aspect=20, pad=0.02)
cbar2.set_label(f"SD of {FITNESS_COL} across libraries")

plt.tight_layout()
plt.show()

print("Generated mean landscape with variability panel.")
