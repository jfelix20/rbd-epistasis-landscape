#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (needed for 3D plotting)


# ===========================
# 1. Load data
# ===========================
excel_path = "human_with_experimental_landscape_jons_data.xlsx"
df = pd.read_excel(excel_path)

# ===========================
# 2. Choose fitness metric
# ===========================
# Default: ACE2-specific enrichment
FITNESS_COL = "delta_log2fc_pos_minus_neg"  # log2(Pos/Pre) - log2(Neg/Pre)

if FITNESS_COL not in df.columns:
    raise KeyError(f"{FITNESS_COL} not found! Check column names in the Excel file.")

# ===========================
# 3. Make sure we have binary_8 and map to x, y
# ===========================
if "binary_8" in df.columns:
    df["binary_8"] = df["binary_8"].astype(str).str.zfill(8)
elif "binary" in df.columns:
    df["binary_8"] = df["binary"].astype(str).str.zfill(8)
else:
    raise KeyError("Neither 'binary_8' nor 'binary' column found in the dataframe.")


def binary_to_xy(bin_str: str):
    """Map 8-bit genotype to (x, y) on a 16×16 grid."""
    bin_str = str(bin_str).zfill(8)
    x = int(bin_str[:4], 2)   # bits 1–4 -> positions 477, 478, 484, 493
    y = int(bin_str[4:], 2)   # bits 5–8 -> positions 496, 498, 501, 505
    return x, y


df["x"], df["y"] = zip(*df["binary_8"].apply(binary_to_xy))
df["fitness"] = df[FITNESS_COL].astype(float)

# ===========================
# 4. Build discrete 16×16 fitness grid
# ===========================
# Grid indices 0–15 in both x and y
GRID_SIZE = 16
F_grid = np.full((GRID_SIZE, GRID_SIZE), np.nan)

# If there happen to be duplicates, average them
for (x, y), sub in df.groupby(["x", "y"]):
    if 0 <= x < GRID_SIZE and 0 <= y < GRID_SIZE:
        F_grid[x, y] = sub["fitness"].mean()

# Optional: fix color range for comparability
# Adjust these if you want a different global scale
vmin = -6.0
vmax = 4.0

# Mask NaNs for plotting
F_masked = np.ma.masked_invalid(F_grid)

# Locate Wuhan and Omicron
wuhan_bin = "00000000"
omicron_bin = "11111111"
wuhan = df[df["binary_8"] == wuhan_bin]
omicron = df[df["binary_8"] == omicron_bin]

# ===========================
# 5. 3D discrete surface + data points
# ===========================
x_vals = np.arange(GRID_SIZE)
y_vals = np.arange(GRID_SIZE)
X, Y = np.meshgrid(x_vals, y_vals, indexing="ij")

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

surf = ax.plot_surface(
    X,
    Y,
    F_masked,
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
    rstride=1,
    cstride=1,
    linewidth=0,
    antialiased=False,
    alpha=0.9,
)

# Overlay actual data points
ax.scatter(
    df["x"],
    df["y"],
    df["fitness"],
    s=15,
    c="k",
    alpha=0.7,
    label="Measured genotypes",
)

# Highlight Wuhan and Omicron if present
if not wuhan.empty:
    ax.scatter(
        wuhan["x"],
        wuhan["y"],
        wuhan["fitness"],
        s=80,
        c="red",
        label="Wuhan (00000000)",
    )

if not omicron.empty:
    ax.scatter(
        omicron["x"],
        omicron["y"],
        omicron["fitness"],
        s=80,
        c="cyan",
        label="Omicron (11111111)",
    )

# Axis labels: map bits to positions
ax.set_xlabel("Bits 1–4 (x): 477, 478, 484, 493")
ax.set_ylabel("Bits 5–8 (y): 496, 498, 501, 505")
ax.set_zlabel(f"Fitness = {FITNESS_COL}")

ax.set_title("Experimental fitness landscape (ACE2-specific enrichment)")

# Cleaner ticks: show only 0 and 15 with bit labels
ax.set_xticks([0, 15])
ax.set_xticklabels(["0000 (all Wuhan-like)", "1111 (all Omicron-like)"], rotation=20)
ax.set_yticks([0, 15])
ax.set_yticklabels(["0000", "1111"], rotation=20)

fig.colorbar(
    surf,
    shrink=0.5,
    aspect=10,
    label=f"Fitness ({FITNESS_COL})",
)

ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
plt.show()

print("3D discrete experimental landscape plotted.")

# ===========================
# 6. Top-down 2D heatmap + data points
# ===========================
fig2, ax2 = plt.subplots(figsize=(8, 6))

# imshow expects [row, col] -> [y, x], so transpose and use extent
im = ax2.imshow(
    F_masked.T,
    origin="lower",
    extent=[0, GRID_SIZE - 1, 0, GRID_SIZE - 1],
    cmap="viridis",
    vmin=vmin,
    vmax=vmax,
    aspect="equal",
)

# Overlay the actual data locations
ax2.scatter(df["x"], df["y"], s=15, c="k", alpha=0.5)

# Highlight Wuhan and Omicron in 2D
if not wuhan.empty:
    ax2.scatter(wuhan["x"], wuhan["y"], s=60, c="red", label="Wuhan (00000000)")
if not omicron.empty:
    ax2.scatter(omicron["x"], omicron["y"], s=60, c="cyan", label="Omicron (11111111)")

ax2.set_xlabel("Bits 1–4 (x): 477, 478, 484, 493")
ax2.set_ylabel("Bits 5–8 (y): 496, 498, 501, 505")
ax2.set_title("Experimental fitness landscape (top-down view)")

# Same simplified ticks as 3D
ax2.set_xticks([0, 15])
ax2.set_xticklabels(["0000 (Wuhan-like)", "1111 (Omicron-like)"])
ax2.set_yticks([0, 15])
ax2.set_yticklabels(["0000", "1111"])

cbar = fig2.colorbar(im, ax=ax2)
cbar.set_label(f"Fitness ({FITNESS_COL})")

ax2.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
plt.show()

print("2D discrete heatmap of experimental landscape plotted.")

