"""Standalone traversal-based alignment of nematic directors (TADA)."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from scipy.spatial import KDTree

def normalize_directors(directors):
    directors = np.asarray(directors, float)
    if directors.ndim != 2 or directors.shape[1] != 3:
        raise ValueError("directors must have shape (N, 3)")
    norms = np.linalg.norm(directors, axis=1)
    if np.any(~np.isfinite(directors)) or np.any(norms <= np.finfo(float).eps):
        raise ValueError("directors must be finite and non-zero")
    return directors / norms[:, None]

def sorted_neighbors(positions, cutoff):
    tree = KDTree(positions)
    result = []
    for i, indices in enumerate(tree.query_ball_point(positions, cutoff, workers=-1)):
        ordered = sorted(indices, key=lambda j: np.linalg.norm(positions[i] - positions[j]))
        result.append(np.asarray(ordered, dtype=np.int64))
    return result

def run_tada(positions, directors, cutoff, seed=0):
    """Orient directors by growing from a selected seed through the neighbor graph."""
    positions, directors = np.asarray(positions, float), normalize_directors(directors)
    if positions.shape != directors.shape:
        raise ValueError("positions and directors must both have shape (N, 3)")
    if cutoff <= 0:
        raise ValueError("cutoff must be greater than zero")
    n = len(positions)
    if not n:
        return directors.copy(), np.empty(0, dtype=np.int32)
    if seed < 0 or seed >= n:
        raise ValueError("seed index is outside the particle range")
    neighbors = sorted_neighbors(positions, cutoff)
    aligned, order = np.zeros_like(directors), np.zeros(n, dtype=np.int32)
    next_order, center = 1, seed
    aligned[center], order[center] = directors[center], next_order
    while next_order < n:
        for candidate in neighbors[center]:
            if order[candidate] != 0:
                continue
            assigned = [j for j in neighbors[candidate] if order[j] != 0]
            if not assigned:
                continue
            closest = assigned[0]
            sign = 1.0 if np.dot(aligned[closest], directors[candidate]) >= 0 else -1.0
            aligned[candidate] = sign * directors[candidate]
            next_order += 1
            order[candidate] = next_order
        candidates = np.flatnonzero(order == min(next_order, n))
        if len(candidates):
            center = int(candidates[0])
        else:
            unassigned = np.flatnonzero(order == 0)
            if not len(unassigned):
                break
            center = int(unassigned[0])
            next_order += 1
            aligned[center], order[center] = directors[center], next_order
    return aligned, order

def main():
    parser = argparse.ArgumentParser(description="Apply TADA to an x y z nx ny nz table")
    parser.add_argument("input", type=Path)
    parser.add_argument("--cutoff", type=float, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=Path("tada_output.txt"))
    args = parser.parse_args()
    data = np.loadtxt(args.input)
    if data.ndim != 2 or data.shape[1] != 6:
        raise ValueError("input must contain x y z nx ny nz columns")
    aligned, order = run_tada(data[:, :3], data[:, 3:], args.cutoff, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output, np.column_stack((data[:, :3], aligned, order)), header="x y z nx ny nz assignment_order")
    print(f"Wrote {args.output} ({len(data)} particles)")

if __name__ == "__main__":
    main()
