"""Erdős–Rényi construction matching NetLogo `build-er-network` draw order."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.random import Generator


@dataclass(frozen=True)
class ERGraph:
    """Undirected simple graph as sorted adjacency lists."""

    n: int
    neighbors: list[list[int]]
    edges: list[tuple[int, int]]


def build_er_netlogo_order(num_nodes: int, avg_degree: float, rng: Generator) -> ERGraph:
    """
    G(n,p) with p = avg_degree / (n-1), same as NetLogo.

    Random draws follow the same nested loop as NetLogo: for each i in 0..n-1,
    for each j > i, one Bernoulli(p) trial.
    """
    n = int(num_nodes)
    if n < 2:
        return ERGraph(n, [[] for _ in range(n)], [])

    p = float(avg_degree) / float(n - 1)
    neighbors: list[list[int]] = [[] for _ in range(n)]
    edges: list[tuple[int, int]] = []

    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p:
                neighbors[i].append(j)
                neighbors[j].append(i)
                edges.append((i, j))

    for i in range(n):
        neighbors[i].sort()

    return ERGraph(n=n, neighbors=neighbors, edges=edges)
