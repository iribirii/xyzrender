"""Ensemble overlay: align and merge multiple conformers into one graph.

This module is intentionally independent from :mod:`xyzrender.overlay`.
It reuses the same Kabsch-based alignment idea, but does *not* apply any
overlay-specific styling — colours are left to the normal CPK palette.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import networkx as nx


# Tiny z-offset between conformers to avoid z-fighting in SVG rendering.
_Z_NUDGE: float = -1e-3


def _node_list(graph: nx.Graph) -> list:
    return list(graph.nodes())


def _positions_from_frame(frame: dict) -> np.ndarray:
    """Return positions array for a trajectory *frame* dict."""
    return np.array(frame["positions"], dtype=float)


def _kabsch_rotation(p_centered: np.ndarray, q_centered: np.ndarray) -> np.ndarray:
    """Kabsch rotation matrix rot s.t. q_centered @ rot.T ≈ p_centered.

    Both arrays must already be mean-centred.  Ensures det(rot) = +1 (proper
    rotation, no reflection).
    """
    h = q_centered.T @ p_centered
    u, _, vt = np.linalg.svd(h)
    det = np.linalg.det(vt.T @ u.T)
    d_mat = np.diag([1.0, 1.0, det])
    return vt.T @ d_mat @ u.T


def align(
    frames: list[dict],
    *,
    reference_frame: int = 0,
) -> list[np.ndarray]:
    """Align all trajectory *frames* onto *reference_frame*.

    Parameters
    ----------
    frames:
        List of ``{"symbols": [...], "positions": [[x,y,z], ...]}`` dicts as
        returned by :func:`xyzrender.readers.load_trajectory_frames`.
    reference_frame:
        Index of the reference frame.  All other frames are RMSD-aligned
        onto this frame via the Kabsch algorithm.

    Returns
    -------
    list of np.ndarray
        One array per frame with aligned 3-D positions, in the same order as
        *frames*.  The reference frame positions are returned unchanged.
    """
    if not frames:
        msg = "ensemble.align: no frames provided"
        raise ValueError(msg)
    if not (0 <= reference_frame < len(frames)):
        msg = f"ensemble.align: reference_frame {reference_frame} out of range for {len(frames)} frames"
        raise ValueError(msg)

    ref = frames[reference_frame]
    ref_pos = _positions_from_frame(ref)
    n_atoms = ref_pos.shape[0]

    aligned: list[np.ndarray] = []
    c_ref = ref_pos.mean(axis=0)

    for idx, frame in enumerate(frames):
        pos = _positions_from_frame(frame)
        if pos.shape != ref_pos.shape:
            msg = (
                f"ensemble.align: frame {idx} has shape {pos.shape}, "
                f"expected {ref_pos.shape} from reference frame"
            )
            raise ValueError(msg)
        if idx == reference_frame:
            aligned.append(ref_pos.copy())
            continue
        c = pos.mean(axis=0)
        rot = _kabsch_rotation(ref_pos - c_ref, pos - c)
        aligned.append((pos - c) @ rot.T + c_ref)

    assert len(aligned) == len(frames)
    assert all(a.shape == (n_atoms, 3) for a in aligned)
    return aligned


def merge_graphs(reference_graph: nx.Graph, aligned_positions: list[np.ndarray]) -> nx.Graph:
    """Merge *reference_graph* with additional conformers into a single graph.

    No overlay-specific styling attributes are added; the renderer will use
    the normal CPK palette based on element symbols.

    Node attributes:
    - ``molecule_index``: conformer index (0 for reference frame).

    Edge attributes:
    - ``molecule_index``: conformer index for that bond set.
    """
    import networkx as nx

    if not aligned_positions:
        msg = "ensemble.merge_graphs: aligned_positions must contain at least one frame"
        raise ValueError(msg)

    base_nodes = _node_list(reference_graph)
    n_base = len(base_nodes)
    n_frames = len(aligned_positions)

    if aligned_positions[0].shape[0] != n_base:
        msg = (
            "ensemble.merge_graphs: position array length does not match "
            f"reference graph (got {aligned_positions[0].shape[0]}, expected {n_base})"
        )
        raise ValueError(msg)

    merged = nx.Graph()
    merged.graph.update(reference_graph.graph)

    # Reference conformer (index 0): keep original node IDs.
    pos0 = aligned_positions[0]
    for k, nid in enumerate(base_nodes):
        data = dict(reference_graph.nodes[nid])
        data["molecule_index"] = 0
        x, y, z = pos0[k]
        data["position"] = (float(x), float(y), float(z))
        merged.add_node(nid, **data)

    for i, j, d in reference_graph.edges(data=True):
        merged.add_edge(i, j, **dict(d), molecule_index=0)

    # Additional conformers: copy node/edge attributes, renumbering node IDs.
    next_id = n_base
    for conf_idx in range(1, n_frames):
        pos = aligned_positions[conf_idx]
        if pos.shape[0] != n_base:
            msg = (
                "ensemble.merge_graphs: position array length does not match "
                f"reference graph (got {pos.shape[0]}, expected {n_base})"
            )
            raise ValueError(msg)

        id_map = {old: next_id + i for i, old in enumerate(base_nodes)}

        for k, old_id in enumerate(base_nodes):
            data = dict(reference_graph.nodes[old_id])
            data["molecule_index"] = conf_idx
            x, y, z = pos[k]
            # Nudge z slightly so conformers don't z-fight in SVG rendering.
            data["position"] = (float(x), float(y), float(z) + conf_idx * _Z_NUDGE)
            merged.add_node(id_map[old_id], **data)

        for i, j, d in reference_graph.edges(data=True):
            merged.add_edge(id_map[i], id_map[j], **dict(d), molecule_index=conf_idx)

        next_id += n_base

    return merged

