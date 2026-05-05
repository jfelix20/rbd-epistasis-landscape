
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.interpolate import griddata

# ==========================================================
# Configuration
# ==========================================================
CSV_PATH = "/Users/jonathanfelix/Documents/Yeast_Lib_Final_Deduplicated/Summary_Figures/WO_landscape_ddelta_posneg_volcano_source.csv"

# Change this if you want a different landscape metric
FITNESS_COL = "ddelta"   # alternatives: "log2fc_pos_pre", "log2fc_neg_pre", "delta_pos", "delta_neg"

# Interpolation method for surface
INTERP_METHOD = "linear"   # "linear", "nearest", or "cubic"

# If True, save the figure in the same folder as the input csv
SAVE_FIGURE = False
OUTPUT_FIG = "landscape_3libraries_ddelta.png"


# ==========================================================
# Helpers
# ==========================================================
def wo_to_binary8(wo_str: str) -> str:
    """
    Convert a genotype string like 'OWWWOWOO' into an 8-bit binary string
    where O -> 1 and W -> 0.
    """
    wo_str = str(wo_str).strip().upper()
    if len(wo_str) != 8 or any(ch not in {"O", "W"} for ch in wo_str):
        raise ValueError(f"Invalid wo genotype: {wo_str}")
    return "".join("1" if ch == "O" else "0" for ch in wo_str)


def binary_to_xy(bin_str: str):
    """
    Map an 8-bit binary genotype to 2D coordinates:
    first 4 bits -> x
    last 4 bits -> y
    """
    bin_str = str(bin_str).zfill(8)
    x = int(bin_str[:4], 2)
    y = int(bin_str[4:], 2)
    return x, y


def get_binary8(df: pd.DataFrame) -> pd.Series:
    """
    Build binary_8 using the best available column.
    Priority:
      1) binary_8
      2) wo
      3) binary
    """
    if "binary_8" in df.columns:
        return df["binary_8"].astype(str).str.zfill(8)

    if "wo" in df.columns:
        return df["wo"].apply(wo_to_binary8)

    if "binary" in df.columns:
        return df["binary"].astype(str).str.zfill(8)

    raise KeyError("Could not find 'binary_8', 'wo', or 'binary' column.")


def plot_landscape_for_library(ax, subdf, fitness_col, title):
    """
    Plot one 3D landscape for a single library on the provided axis.
    """
    if subdf.empty:
        ax.set_title(f"{title}\n(no data)")
        return None

    subdf = subdf.copy()
    subdf["x"], subdf["y"] = zip(*subdf["binary_8"].apply(binary_to_xy))
    subdf["fitness"] = subdf[fitness_col]

    # Build regular grid spanning the full 4-bit x 4-bit layout (0-15)
    grid_x, grid_y = np.mgrid[0:15:100j, 0:15:100j]

    grid_z = griddata(
        points=subdf[["x", "y"]].values,
        values=subdf["fitness"].values,
        xi=(grid_x, grid_y),
        method=INTERP_METHOD
    )

    # Fall back to nearest if the chosen interpolation leaves too many gaps
    if grid_z is None or np.all(np.isnan(grid_z)):
        grid_z = griddata(
            points=subdf[["x", "y"]].values,
            values=subdf["fitness"].values,
            xi=(grid_x, grid_y),
            method="nearest"
        )

    grid_z = np.ma.array(grid_z, mask=np.isnan(grid_z))

    surf = ax.plot_surface(
        grid_x,
        grid_y,
        grid_z,
        cmap="viridis",
        linewidth=0,
        antialiased=True,
        alpha=0.9
    )

    # Overlay actual measured genotype points
    ax.scatter(
        subdf["x"],
        subdf["y"],
        subdf["fitness"],
        c=subdf["fitness"],
        cmap="viridis",
        s=18,
        edgecolor="k",
        linewidth=0.25,
        alpha=0.95
    )

    # Mark Wuhan and Omicron if present
    wuhan = subdf[subdf["binary_8"] == "00000000"]
    omicron = subdf[subdf["binary_8"] == "11111111"]

    if not wuhan.empty:
        ax.scatter(
            wuhan["x"], wuhan["y"], wuhan["fitness"],
            s=90, c="red", label="Wuhan", depthshade=False
        )

    if not omicron.empty:
        ax.scatter(
            omicron["x"], omicron["y"], omicron["fitness"],
            s=90, c="cyan", label="Omicron", depthshade=False
        )

    ax.set_title(title)
    ax.set_xlabel("477/478/484/493")
    ax.set_ylabel("496/498/501/505")
    ax.set_zlabel(fitness_col)
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 15)

    if (not wuhan.empty) or (not omicron.empty):
        ax.legend(loc="upper left", fontsize=8)

    return surf


# ==========================================================
# Main
# ==========================================================
def main():
    # 1) Load data
    df = pd.read_csv(CSV_PATH)

    # 2) Basic checks
    required_cols = ["library_id", FITNESS_COL]
    for col in required_cols:
        if col not in df.columns:
            raise KeyError(f"Required column '{col}' not found in input CSV.")

    # 3) Clean and derive binary_8
    df[FITNESS_COL] = pd.to_numeric(df[FITNESS_COL], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["library_id", FITNESS_COL]).copy()

    df["binary_8"] = get_binary8(df)

    # 4) Sort libraries
    libraries = sorted(df["library_id"].dropna().unique())
    n_libs = len(libraries)

    if n_libs == 0:
        raise ValueError("No libraries found after cleaning.")

    # 5) Plot
    fig = plt.figure(figsize=(6 * n_libs, 6))

    surfaces = []
    for i, lib in enumerate(libraries, start=1):
        ax = fig.add_subplot(1, n_libs, i, projection="3d")
        subdf = df[df["library_id"] == lib].copy()
        surf = plot_landscape_for_library(
            ax=ax,
            subdf=subdf,
            fitness_col=FITNESS_COL,
            title=f"Library {lib}"
        )
        if surf is not None:
            surfaces.append(surf)

    fig.suptitle(f"Experimental Fitness Landscapes by Library ({FITNESS_COL})", y=0.98)

    # One shared colorbar using the first valid surface
    if surfaces:
        fig.colorbar(
            surfaces[0],
            ax=fig.axes,
            shrink=0.6,
            aspect=22,
            pad=0.03,
            label=FITNESS_COL
        )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    if SAVE_FIGURE:
        import os
        output_dir = os.path.dirname(CSV_PATH)
        output_path = os.path.join(output_dir, OUTPUT_FIG)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Saved figure to: {output_path}")

    plt.show()
    print(f"3D experimental landscapes plotted for libraries: {libraries}")


if __name__ == "__main__":
    main()
