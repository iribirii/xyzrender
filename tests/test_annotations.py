"""Tests for annotation / colormap file loaders."""

from pathlib import Path

import pytest

from xyzrender import load
from xyzrender.annotations import load_bond_cmap

STRUCTURES = Path(__file__).parent.parent / "examples" / "structures"


@pytest.fixture(scope="module")
def ethanol():
    return load(STRUCTURES / "ethanol.xyz")


def test_load_bond_cmap_valid(ethanol, tmp_path):
    path = tmp_path / "bonds.txt"
    path.write_text("2 3 0.5\n# comment\n1 2 1.0\n")
    result = load_bond_cmap(str(path), ethanol.graph)
    assert result[(1, 2)] == 0.5  # 0-indexed 1-2 from line "2 3"
    assert result[(0, 1)] == 1.0


def test_load_bond_cmap_canonicalizes_pair(ethanol, tmp_path):
    path = tmp_path / "bonds.txt"
    path.write_text("3 2 0.25\n")
    assert load_bond_cmap(str(path), ethanol.graph)[(1, 2)] == 0.25


def test_load_bond_cmap_missing_atom(ethanol, tmp_path):
    path = tmp_path / "bonds.txt"
    path.write_text("2 99 1.0\n")
    with pytest.raises(ValueError, match="not found"):
        load_bond_cmap(str(path), ethanol.graph)


def test_load_bond_cmap_not_bonded(ethanol, tmp_path):
    path = tmp_path / "bonds.txt"
    path.write_text("3 4 1.0\n")
    with pytest.raises(ValueError, match="not bonded"):
        load_bond_cmap(str(path), ethanol.graph)


def test_load_bond_cmap_same_atom(ethanol, tmp_path):
    path = tmp_path / "bonds.txt"
    path.write_text("2 2 1.0\n")
    with pytest.raises(ValueError, match="must differ"):
        load_bond_cmap(str(path), ethanol.graph)
