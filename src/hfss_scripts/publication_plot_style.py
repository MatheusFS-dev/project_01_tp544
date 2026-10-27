"""Reusable publication-style plotting helpers.

This module defines a compact Matplotlib visual language for scientific figures.
It is intended as a reusable style reference, not as a mandatory framework.
"""

from pathlib import Path

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np

PALETTE = {
    "background": "#ffffff",
    "grid": "#d7dbe2",
    "axis": "#20242b",
    "muted": "#9aa4b2",
    "primary": "#1f4e79",
    "secondary": "#5b6c5d",
    "highlight": "#b54708",
    "selection": "#7a1f3d",
    "reference": "#3f3f46",
    "box_fill": "#d9e6f2",
}


def apply_publication_style() -> None:
    """Apply the shared publication-style Matplotlib theme.

    Args:
        None.

    Returns:
        None.

    Raises:
        None.
    """
    plt.rcParams.update(
        {
            "figure.facecolor": PALETTE["background"],
            "axes.facecolor": PALETTE["background"],
            "axes.edgecolor": PALETTE["axis"],
            "axes.linewidth": 0.8,
            "axes.labelcolor": PALETTE["axis"],
            "axes.labelsize": 16,
            "axes.titlesize": 20,
            "axes.titleweight": "semibold",
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 16,
            "grid.color": PALETTE["grid"],
            "grid.linewidth": 0.6,
            "grid.alpha": 0.65,
            "legend.frameon": False,
            "legend.fontsize": 16,
            "legend.title_fontsize": 16,
            "xtick.color": PALETTE["axis"],
            "ytick.color": PALETTE["axis"],
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "savefig.facecolor": PALETTE["background"],
        }
    )


def create_figure(width: float = 6.4, height: float = 4.4) -> tuple[plt.Figure, plt.Axes]:
    """Create a single-panel figure with the shared style applied.

    Args:
        width: Figure width in inches.
        height: Figure height in inches.

    Returns:
        A tuple containing the created figure and axes.

    Raises:
        None.
    """
    apply_publication_style()
    figure, axes = plt.subplots(figsize=(width, height), constrained_layout=True)
    return figure, axes


def create_two_panel_figure(width: float = 10.2, height: float = 4.2) -> tuple[plt.Figure, np.ndarray]:
    """Create a two-panel figure with the shared style applied.

    Args:
        width: Figure width in inches.
        height: Figure height in inches.

    Returns:
        A tuple containing the created figure and an array of two axes.

    Raises:
        None.
    """
    apply_publication_style()
    figure, axes = plt.subplots(1, 2, figsize=(width, height), constrained_layout=True)
    return figure, axes


def style_axes(
    axes: plt.Axes,
    xlabel: str,
    ylabel: str,
    title: str | None = None,
    grid_axis: str = "both",
) -> None:
    """Apply consistent axis-level styling.

    Args:
        axes: Axes to style.
        xlabel: X-axis label.
        ylabel: Y-axis label.
        title: Optional title. Use None for in-paper panels where captions provide context.
        grid_axis: Grid direction, one of "x", "y", or "both".

    Returns:
        None.

    Raises:
        ValueError: If grid_axis is not one of "x", "y", or "both".
    """
    if grid_axis not in {"x", "y", "both"}:
        raise ValueError(f"grid_axis must be 'x', 'y', or 'both', got {grid_axis!r}.")

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    if title is not None:
        axes.set_title(title, pad=10.0)
    axes.grid(True, axis=grid_axis, linestyle="--", dashes=(2.0, 2.4))
    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)
    axes.tick_params(length=4.0, width=0.8)


def style_legend(axes: plt.Axes, location: str = "best") -> None:
    """Add a frameless legend when labeled artists exist.

    Args:
        axes: Axes containing labeled artists.
        location: Legend location accepted by Matplotlib.

    Returns:
        None.

    Raises:
        None.
    """
    handles, labels = axes.get_legend_handles_labels()
    if handles:
        axes.legend(handles, labels, loc=location, frameon=False)


def save_figure(figure: plt.Figure, output_path: str | Path, dpi: int = 300) -> None:
    """Save a figure with publication-friendly defaults.

    Args:
        figure: Figure to save.
        output_path: Destination file path.
        dpi: Dots per inch used for raster output.

    Returns:
        None.

    Raises:
        OSError: If the parent directory cannot be created or the file cannot be written.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")


def draw_reference_boxplot(axes: plt.Axes, groups: list[np.ndarray], labels: list[str]) -> None:
    """Draw a compact publication-style boxplot.

    Args:
        axes: Target axes.
        groups: Numeric arrays, one per category.
        labels: Category labels, one per group.

    Returns:
        None.

    Raises:
        ValueError: If the number of groups and labels differs.
    """
    if len(groups) != len(labels):
        raise ValueError("groups and labels must have the same length.")

    axes.boxplot(
        groups,
        labels=labels,
        patch_artist=True,
        widths=0.55,
        medianprops={"color": PALETTE["highlight"], "linewidth": 1.4},
        boxprops={"facecolor": PALETTE["box_fill"], "edgecolor": PALETTE["axis"], "linewidth": 0.9},
        whiskerprops={"color": PALETTE["axis"], "linewidth": 0.9},
        capprops={"color": PALETTE["axis"], "linewidth": 0.9},
        flierprops={
            "marker": "o",
            "markersize": 4.2,
            "markerfacecolor": PALETTE["muted"],
            "markeredgecolor": PALETTE["muted"],
            "alpha": 0.7,
        },
    )


def add_boxplot_legend(axes: plt.Axes, location: str = "upper right") -> None:
    """Add a compact visual legend for the reference boxplot style.

    Args:
        axes: Axes containing the boxplot.
        location: Legend location accepted by Matplotlib.

    Returns:
        None.

    Raises:
        None.
    """
    handles = [
        Patch(facecolor=PALETTE["box_fill"], edgecolor=PALETTE["axis"], linewidth=0.9, label="Interquartile range"),
        Line2D([0], [0], color=PALETTE["highlight"], linewidth=1.4, label="Median"),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markersize=4.2,
            markerfacecolor=PALETTE["muted"],
            markeredgecolor=PALETTE["muted"],
            alpha=0.7,
            label="Outlier",
        ),
    ]
    axes.legend(handles=handles, loc=location, frameon=False)
