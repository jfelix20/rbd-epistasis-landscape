import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ===========================
# 1. Load data
# ===========================
excel_path = "human_with_experimental_landscape_jons_data.xlsx"
df = pd.read_excel(excel_path)

FITNESS_COL = "delta_log2fc_pos_minus_neg"
if FITNESS_COL not in df.columns:
    raise KeyError(f"{FITNESS_COL} not found! Check column names in the Excel file.")

# ===========================
# 2. Ensure 8-bit binary string
# ===========================
if "binary_8" in df.columns:
    df["binary_8"] = df["binary_8"].astype(str).str.zfill(8)
elif "binary" in df.columns:
    df["binary_8"] = df["binary"].astype(str).str.zfill(8)
else:
    raise KeyError("Neither 'binary_8' nor 'binary' column found.")

df["fitness"] = pd.to_numeric(df[FITNESS_COL], errors="coerce")

# ===========================
# 3. Gray-code ordering utilities
# ===========================
def gray_encode(i: int) -> int:
    return i ^ (i >> 1)

# Gray sequence for 4 bits: position 0..15 corresponds to bitpattern gray_encode(position)
gray_seq = [gray_encode(i) for i in range(16)]

# Inverse map: for a given 4-bit pattern (as int), what is its index in the Gray sequence?
inv_gray = {g: i for i, g in enumerate(gray_seq)}

def bits4_to_int(bits4: str) -> int:
    return int(bits4, 2)

def int_to_bits4(x: int) -> str:
    return format(x, "04b")

# Tick labels in Gray order (these are the actual 4-bit patterns along the axis)
x_tick_labels = [int_to_bits4(g) for g in gray_seq]
y_tick_labels = [int_to_bits4(g) for g in gray_seq]

# ===========================
# 4. Map each genotype to (x_idx, y_idx) using Gray-order indices
# ===========================
def binary8_to_gray_xy_index(bin8: str):
    bin8 = str(bin8).zfill(8)
    x_bits = bin8[:4]   # bits 1-4
    y_bits = bin8[4:]   # bits 5-8

    x_int = bits4_to_int(x_bits)
    y_int = bits4_to_int(y_bits)

    # Place the bitpattern at its position in the Gray ordering
    x_idx = inv_gray[x_int]
    y_idx = inv_gray[y_int]
    return x_idx, y_idx

df["x_idx"], df["y_idx"] = zip(*df["binary_8"].apply(binary8_to_gray_xy_index))

# ===========================
# 5. Build a 16x16 grid with NO interpolation
# ===========================
Z = np.full((16, 16), np.nan, dtype=float)

for _, row in df.iterrows():
    x = int(row["x_idx"])
    y = int(row["y_idx"])
    val = row["fitness"]
    if np.isfinite(val):
        # If duplicates exist, keep the mean (shouldn't happen in a full 256-variant table)
        if np.isnan(Z[y, x]):
            Z[y, x] = val
        else:
            Z[y, x] = np.nanmean([Z[y, x], val])

# ===========================
# 6. Plot as a 2D heatmap (recommended)
# ===========================
fig, ax = plt.subplots(figsize=(9, 8))

im = ax.imshow(Z, origin="lower", aspect="equal")  # no custom colors specified
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label(f"Fitness ({FITNESS_COL})")

ax.set_title("Experimental fitness landscape (Gray-code ordered, no interpolation)")
ax.set_xlabel("Bits 1–4 (x): 477, 478, 484, 493")
ax.set_ylabel("Bits 5–8 (y): 496, 498, 501, 505")

ax.set_xticks(range(16))
ax.set_yticks(range(16))
ax.set_xticklabels(x_tick_labels, rotation=45, ha="right")
ax.set_yticklabels(y_tick_labels)

# Mark Wuhan and Omicron
# Wuhan = 00000000, Omicron = 11111111
wuhan_x, wuhan_y = binary8_to_gray_xy_index("00000000")
omic_x, omic_y   = binary8_to_gray_xy_index("11111111")

ax.scatter([wuhan_x], [wuhan_y], s=80, marker="o", label="Wuhan (00000000)")
ax.scatter([omic_x],  [omic_y],  s=80, marker="o", label="Omicron (11111111)")

ax.legend(loc="upper right", frameon=True)
plt.tight_layout()
plt.show()

print("Done: Gray-code heatmap plotted (no interpolation).")

