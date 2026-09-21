"""Colormap utilities for scalar color legends and ``--cmap`` atom coloring."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from xyzrender.colors import PALETTES, Color, palette_color


def cmap_value_range(
    values: Iterable[float],
    *,
    cmap_range: tuple[float, float] | None,
    cmap_symm: bool,
) -> tuple[float, float]:
    """Resolve vmin/vmax for scalar colormaps (atoms, bonds, etc.)."""
    if cmap_range is not None and cmap_symm:
        msg = "--cmap-range and --cmap-symm are mutually exclusive"
        raise ValueError(msg)
    if cmap_range is not None:
        return cmap_range
    vals = list(values)
    if cmap_symm:
        vmax = max(abs(v) for v in vals)
        return -vmax, vmax
    return min(vals), max(vals)


def bond_color_hex(value: float, palette: str, vmin: float, vmax: float) -> str:
    """Map a scalar bond property to a palette hex color."""
    vrange = max(vmax - vmin, 1e-10)
    return palette_color(palette, (value - vmin) / vrange).hex


def build_palette_lut(palette: str, size: int = 256) -> np.ndarray:
    """Return an RGB LUT sampled from a named palette."""
    lut = np.zeros((size, 3), dtype=np.uint8)
    scale = max(size - 1, 1)
    for i in range(size):
        c = palette_color(palette, i / scale)
        lut[i] = (c.r, c.g, c.b)
    return lut


# ---------------------------------------------------------------------------
# Atom color list
# ---------------------------------------------------------------------------


def atom_colors(
    atom_cmap: dict[int, float],
    n: int,
    palette: str,
    vmin: float,
    vmax: float,
    unlabeled_hex: str,
) -> list[Color]:
    """Return per-atom Color list; atoms absent from atom_cmap get unlabeled_hex."""
    unlabeled = Color.from_hex(unlabeled_hex)
    vrange = max(vmax - vmin, 1e-10)
    return [
        palette_color(palette, (atom_cmap[ai] - vmin) / vrange) if ai in atom_cmap else unlabeled for ai in range(n)
    ]


# ---------------------------------------------------------------------------
# Colorbar SVG
# ---------------------------------------------------------------------------

_BAR_W = 30.0
_MARGIN = 16.0
_TICK_GAP = 16.0
_CBAR_FONT = "DejaVu Sans Mono"
_CBAR_TICK_COLOR = "#000000"


def _colorbar_tick_label(
    x: float,
    y: float,
    text: str,
    fs: float,
    *,
    anchor: str = "start",
) -> list[str]:
    """Tick label with white halo (GIF/PNG via resvg) and explicit DejaVu font."""
    attrs = (
        f'x="{x:.1f}" y="{y:.1f}" font-family="{_CBAR_FONT}, monospace" font-size="{fs:.1f}px" '
        f'font-weight="bold" text-anchor="{anchor}" dominant-baseline="central"'
    )
    sw = fs * 0.35
    return [
        f'  <text {attrs} fill="#ffffff" stroke="#ffffff" '
        f'stroke-width="{sw:.1f}" stroke-linejoin="round">{text}</text>',
        f'  <text {attrs} fill="{_CBAR_TICK_COLOR}">{text}</text>',
    ]


def colorbar_extra_width(
    vmin: float,
    vmax: float,
    fs: float,
    unit: str | None = None,
) -> int:
    """Extra SVG canvas width needed to fit the colorbar + labels."""
    fs = min(fs, 40.0)
    char_w = fs * 0.62
    mid = (vmin + vmax) / 2
    max_label_chars = max(len(f"{v:.3f}".replace("-", "\u2212")) for v in (vmin, mid, vmax))
    unit_chars = len(unit) if unit else 0
    label_chars = max(max_label_chars, unit_chars)
    return int(_MARGIN + _BAR_W + _TICK_GAP + 3 + label_chars * char_w + 10)


def colorbar_svg(
    vmin: float,
    vmax: float,
    palette: str,
    mol_canvas_w: float,
    canvas_h: float,
    font_size: float,
    label_color: str,
    unit: str | None = None,
) -> list[str]:
    """Return SVG element strings for a vertical colorbar to the right of the molecule."""
    stops = PALETTES[palette]

    bar_x = mol_canvas_w + _MARGIN
    bar_h = max(min(canvas_h * 0.80, 400.0), 60.0)
    bar_top = (canvas_h - bar_h) / 2
    bar_bot = bar_top + bar_h

    # Gradient: top = vmax, bottom = vmin.
    n = len(stops)
    grad_stops = "".join(
        f'<stop offset="{int(i / (n - 1) * 100)}%" stop-color="{c.hex}"/>' for i, c in enumerate(reversed(stops))
    )
    tick_color = _CBAR_TICK_COLOR
    elems = [
        f'  <defs><linearGradient id="_cbg" x1="0" y1="0" x2="0" y2="1">{grad_stops}</linearGradient></defs>',
        f'  <rect x="{bar_x:.1f}" y="{bar_top:.1f}" width="{_BAR_W:.1f}" height="{bar_h:.1f}" '
        f'fill="url(#_cbg)" stroke="{tick_color}" stroke-width="5"/>',
    ]

    tick_x1 = bar_x + _BAR_W
    label_x = tick_x1 + _TICK_GAP + 3
    fs = min(font_size, 40.0)

    ticks = [
        (bar_top, vmax),
        ((bar_top + bar_bot) / 2, (vmin + vmax) / 2),
        (bar_bot, vmin),
    ]

    for ty, val in ticks:
        s = f"{val:.3f}".replace("-", "\u2212")
        elems.append(
            f'  <line x1="{tick_x1:.1f}" y1="{ty:.1f}" x2="{tick_x1 + _TICK_GAP:.1f}" y2="{ty:.1f}" '
            f'stroke="{tick_color}" stroke-width="5"/>'
        )
        elems.extend(_colorbar_tick_label(label_x, ty, s, fs))

    if unit:
        unit_fs = fs * 0.85
        unit_y = min(bar_bot + unit_fs * 1.4, canvas_h - unit_fs * 0.6)
        elems.extend(
            _colorbar_tick_label(bar_x + _BAR_W / 2, unit_y, unit, unit_fs, anchor="middle"),
        )

    return elems
