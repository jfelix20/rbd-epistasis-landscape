import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D  # needed for 3D plotting
from scipy.interpolate import griddata   # pip install scipy if needed

# ===========================
# 1. Load data
# ===========================
excel_path = "linear_reg_human.xlsx"
df = pd.read_excel(excel_path)

# ===========================
# 2. Choose fitness metric
# ===========================
# You can change this to "delta_log2fc_pos_minus_neg" if you want that instead
FITNESS_COL = "DeltaG_pred_docking"

if FITNESS_COL not in df.columns:
    raise KeyError(f"{FITNESS_COL} not found! Check column names in the Excel file.")

# ===========================
# 3. Make sure we have binary_8 and map to x,y
# ===========================
if "binary_8" not in df.columns:
    if "binary" not in df.columns:
        raise KeyError("Neither 'binary_8' nor 'binary' column found.")
    df["binary_8"] = df["binary"].astype(str).str.zfill(8)
else:
    df["binary_8"] = df["binary_8"].astype(str).str.zfill(8)

def binary_to_xy(bin_str):
    bin_str = str(bin_str).zfill(8)
    x = int(bin_str[:4], 2)   # bits 1-4
    y = int(bin_str[4:], 2)   # bits 5-8
    return x, y

df["x"], df["y"] = zip(*df["binary_8"].apply(binary_to_xy))
df["fitness"] = df[FITNESS_COL]

# ===========================
# 4. Build grid for smooth 3D surface
# ===========================
x_min, x_max = df["x"].min(), df["x"].max()
y_min, y_max = df["y"].min(), df["y"].max()

grid_x, grid_y = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]

grid_z = griddata(
    points=df[["x", "y"]].values,
    values=df["fitness"].values,
    xi=(grid_x, grid_y),
    method="linear"
)

# Mask NaNs so they don't break the surface plot
grid_z = np.ma.array(grid_z, mask=np.isnan(grid_z))

# ===========================
# 5. 3D surface plot (landscape)
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
    alpha=0.9
)

# Label axes
ax.set_xlabel("Bits 1–4 (x)")
ax.set_ylabel("Bits 5–8 (y)")
ax.set_zlabel(f"Fitness ({FITNESS_COL})")
ax.set_title("Experimental Fitness Landscape (3D)")

# Colorbar
fig.colorbar(surf, shrink=0.5, aspect=10, label=f"Fitness ({FITNESS_COL})")

# ===========================
# 6. Mark Wuhan and Omicron
# ===========================
wuhan = df[df["binary_8"] == "00000000"]
omicron = df[df["binary_8"] == "11111111"]

if not wuhan.empty:
    ax.scatter(
        wuhan["x"], wuhan["y"], wuhan["fitness"],
        s=80, c="red", label="Wuhan"
    )

if not omicron.empty:
    ax.scatter(
        omicron["x"], omicron["y"], omicron["fitness"],
        s=80, c="cyan", label="Omicron"
    )

if (not wuhan.empty) or (not omicron.empty):
    ax.legend()

plt.tight_layout()
plt.show()

print("3D experimental landscape plotted.")