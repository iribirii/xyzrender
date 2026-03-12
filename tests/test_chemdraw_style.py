from __future__ import annotations

from pathlib import Path

from xyzrender import SVGResult, load, render


def _first_svg_tag(svg: str, tag: str) -> bool:
    return f"<{tag}" in svg


def test_chemdraw_style_no_atom_circles(tmp_path: Path) -> None:
    mol = load("examples/structures/ethanol.xyz")
    out = tmp_path / "ethanol_chemdraw.svg"
    result = render(mol, config="chemdraw", output=out)
    assert isinstance(result, SVGResult)
    svg = out.read_text()
    # No atom circles should be drawn for Chemdraw style
    assert not _first_svg_tag(svg, "circle")


def test_chemdraw_style_non_carbon_labels(tmp_path: Path) -> None:
    mol = load("examples/structures/ethanol.xyz")
    out = tmp_path / "ethanol_chemdraw_labels.svg"
    render(mol, config="chemdraw", output=out)
    svg = out.read_text()
    # Expect at least one text label for O and H atoms
    assert "<text" in svg
    assert "O" in svg

