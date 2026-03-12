from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np

from xyzrender import SVGResult, ensemble
from xyzrender.api import _build_ensemble_molecule


def _write_multiframe_xyz(path: Path, frames: list[list[tuple[str, tuple[float, float, float]]]]) -> None:
    """Write a simple multi-frame XYZ file for testing."""
    lines: list[str] = []
    for frame in frames:
        lines.append(f"{len(frame)}\n")
        lines.append("test frame\n")
        for sym, (x, y, z) in frame:
            lines.append(f"{sym:<3} {x:15.8f} {y:15.8f} {z:15.8f}\n")
    path.write_text("".join(lines))


def test_build_ensemble_molecule(tmp_path: Path) -> None:
    base = [
        ("H", (0.0, 0.0, 0.0)),
        ("O", (0.0, 0.0, 1.0)),
    ]
    # Two additional frames: small displacements
    frames = [
        base,
        [("H", (0.1, 0.0, 0.0)), ("O", (0.0, 0.1, 1.0))],
        [("H", (-0.1, 0.0, 0.0)), ("O", (0.0, -0.1, 1.0))],
    ]
    xyz_path = tmp_path / "traj.xyz"
    _write_multiframe_xyz(xyz_path, frames)

    mol = _build_ensemble_molecule(xyz_path)
    g = mol.graph

    # Three conformers × 2 atoms each
    assert g.number_of_nodes() == 6
    assert g.number_of_edges() == 3  # one bond per conformer

    # All atoms carry molecule_index and no overlay flag
    molecule_indices = {data["molecule_index"] for _, data in g.nodes(data=True)}
    assert molecule_indices == {0, 1, 2}
    assert all("overlay" not in data for _, data in g.nodes(data=True))

    # All bonds share the same molecule_index as their atoms
    for i, j, d in g.edges(data=True):
        mi = g.nodes[i]["molecule_index"]
        mj = g.nodes[j]["molecule_index"]
        assert mi == mj == d["molecule_index"]
        # No overlay-specific bond colour override
        assert "bond_color_override" not in d


def test_ensemble_api_returns_svg(tmp_path: Path) -> None:
    base = [
        ("H", (0.0, 0.0, 0.0)),
        ("O", (0.0, 0.0, 1.0)),
    ]
    frames = [
        base,
        [("H", (0.1, 0.0, 0.0)), ("O", (0.0, 0.1, 1.0))],
    ]
    xyz_path = tmp_path / "traj.xyz"
    _write_multiframe_xyz(xyz_path, frames)

    result = ensemble(xyz_path, output=tmp_path / "out.svg")
    assert isinstance(result, SVGResult)
    out_path = tmp_path / "out.svg"
    assert out_path.exists()
    svg = out_path.read_text()
    assert "<svg" in svg

