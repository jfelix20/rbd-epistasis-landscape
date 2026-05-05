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
# Change this if you want a different landscape metric
FITNESS_COL = "ddelta"

if FITNESS_COL not in df.columns:
    raise KeyError(f"{FITNESS_COL} not found in the dataset.")

if "wo" not in df.columns:
    raise KeyError("'wo' column not found in the dataset.")

# ===========================
# 3. Convert wo strings to binary_8 and x/y coordinates
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
    x = int(bin_str[:4], 2)   # positions 1–4
    y = int(bin_str[4:], 2)   # positions 5–8
    return x, y

df["binary_8"] = df["wo"].apply(wo_to_binary8)
df["x"], df["y"] = zip(*df["binary_8"].apply(binary_to_xy))

# ===========================
# 4. Average fitness across libraries by genotype
# ===========================
avg_df = (
    df.groupby(["wo", "binary_8", "x", "y"], as_index=False)
      .agg(
          mean_fitness=(FITNESS_COL, "mean"),
          n_libraries=("library_id", "nunique")
      )
)

# ===========================
# 5. Build smooth 3D surface
# ===========================
x_min, x_max = avg_df["x"].min(), avg_df["x"].max()
y_min, y_max = avg_df["y"].min(), avg_df["y"].max()

grid_x, grid_y = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]

grid_z = griddata(
    points=avg_df[["x", "y"]].values,
    values=avg_df["mean_fitness"].values,
    xi=(grid_x, grid_y),
    method="linear"
)

grid_z = np.ma.array(grid_z, mask=np.isnan(grid_z))

# ===========================
# 6. Plot averaged landscape (no dots)
# ===========================
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

surf = ax.plot_surface(
    grid_x,
    grid_y,
    grid_z,
    cmap="viridis",
    linewidth=0,
    antialiased=True,
    alpha=0.95
)

ax.set_xlabel("Positions 477/478/484/493")
ax.set_ylabel("Positions 496/498/501/505")
ax.set_zlabel(f"Mean fitness ({FITNESS_COL})")
ax.set_title("Average Experimental Fitness Landscape Across Libraries")

fig.colorbar(surf, shrink=0.55, aspect=12, label=f"Mean fitness ({FITNESS_COL})")

plt.tight_layout()
plt.show()

print("Averaged 3-library landscape plotted successfully.")
