"""
Spine (12 segments) + spinePairs (5 width pairs) from a binary larva mask.

Uses skeletonization + shortest path on skeleton graph (approximation of larvaDistanceMap
interior path), then samples perpendicular width pairs at mid-segments.
"""

from __future__ import annotations

from collections import deque
from typing import List, Tuple

import cv2
import numpy as np

try:
    from skimage.morphology import skeletonize
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "imba_tracker spine extraction requires scikit-image (skimage.morphology.skeletonize)."
    ) from e


def _neighbors8(y: int, x: int, h: int, w: int) -> List[Tuple[int, int]]:
    out = []
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy == 0 and dx == 0:
                continue
            yy, xx = y + dy, x + dx
            if 0 <= yy < h and 0 <= xx < w:
                out.append((yy, xx))
    return out


def _bfs_parents(
    start: Tuple[int, int],
    adj: dict[Tuple[int, int], list[Tuple[int, int]]],
) -> dict[Tuple[int, int], Tuple[int, int] | None]:
    parent: dict[Tuple[int, int], Tuple[int, int] | None] = {start: None}
    q = deque([start])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in parent:
                parent[v] = u
                q.append(v)
    return parent


def _skeleton_path(
    sk: np.ndarray,
) -> Tuple[np.ndarray, Tuple[int, int], Tuple[int, int]] | None:
    """Return ordered path as (N,2) float32 in (x,y) pixel coords, plus endpoints."""
    ys, xs = np.where(sk)
    if len(xs) < 3:
        return None
    pts = list(zip(xs.tolist(), ys.tolist()))
    h, w = sk.shape
    adj: dict[Tuple[int, int], list[Tuple[int, int]]] = {p: [] for p in pts}
    for x, y in pts:
        for nx, ny in _neighbors8(y, x, h, w):
            if sk[ny, nx] and (nx, ny) in adj:
                adj[(x, y)].append((nx, ny))

    def bfs_dist(start: Tuple[int, int]) -> dict[Tuple[int, int], int]:
        d: dict[Tuple[int, int], int] = {start: 0}
        q = deque([start])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v not in d:
                    d[v] = d[u] + 1
                    q.append(v)
        return d

    deg = {p: len(adj[p]) for p in pts}
    ends = [p for p in pts if deg[p] == 1]
    if len(ends) >= 2:
        best_pair = None
        best_d = -1
        for i, a in enumerate(ends):
            da = bfs_dist(a)
            for b in ends[i + 1 :]:
                if b in da and da[b] > best_d:
                    best_d = da[b]
                    best_pair = (a, b)
        if best_pair is None:
            return None
        a, b = best_pair
    else:
        best_pair = None
        best_d = -1
        samp = pts[:: max(1, len(pts) // 64)]
        for i, p1 in enumerate(samp):
            da = bfs_dist(p1)
            for p2 in samp[i + 1 :]:
                if p2 in da and da[p2] > best_d:
                    best_d = da[p2]
                    best_pair = (p1, p2)
        if best_pair is None:
            return None
        a, b = best_pair

    parent = _bfs_parents(a, adj)
    if b not in parent:
        return None
    path_pts: List[Tuple[int, int]] = []
    cur: Tuple[int, int] | None = b
    while cur is not None:
        path_pts.append(cur)
        cur = parent[cur]
    path_pts.reverse()
    arr = np.array([[float(p[0]), float(p[1])] for p in path_pts], dtype=np.float32)
    return arr, a, b


def _resample_path(path_xy: np.ndarray, n: int) -> np.ndarray:
    """path_xy (N,2) in x,y order; resample n points along polyline length."""
    if len(path_xy) < 2:
        return np.tile(path_xy[:1], (n, 1))
    seg = np.sqrt(np.sum(np.diff(path_xy, axis=0) ** 2, axis=1))
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    total = float(cum[-1])
    if total < 1e-6:
        return np.tile(path_xy[:1], (n, 1))
    targets = np.linspace(0.0, total, n, dtype=np.float64)
    out = np.zeros((n, 2), dtype=np.float32)
    j = 0
    for i, t in enumerate(targets):
        while j + 1 < len(cum) and cum[j + 1] < t:
            j += 1
        j = min(j, len(path_xy) - 2)
        t0, t1 = cum[j], cum[j + 1]
        a = (t - t0) / max(t1 - t0, 1e-9)
        out[i] = (1 - a) * path_xy[j] + a * path_xy[j + 1]
    return out


def _contour_poly(mask_u8: np.ndarray) -> np.ndarray | None:
    cts, _ = cv2.findContours(mask_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not cts:
        return None
    c = max(cts, key=cv2.contourArea)
    return c.reshape(-1, 2).astype(np.float32)


def _width_pair_from_contour(
    mid: np.ndarray,
    perp: np.ndarray,
    contour: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Two contour points with min/max projection onto perp (width segment)."""
    proj = np.dot(contour - mid.reshape(1, 2), perp.astype(np.float64))
    i_min = int(np.argmin(proj))
    i_max = int(np.argmax(proj))
    return contour[i_min].copy(), contour[i_max].copy()


def _perpendicular(spine: np.ndarray, seg_idx: int) -> np.ndarray:
    """Unit perpendicular to spine at segment mid between seg_idx and seg_idx+1."""
    if seg_idx + 1 >= len(spine):
        seg_idx = len(spine) - 2
    t = spine[seg_idx].astype(np.float64)
    t2 = spine[seg_idx + 1].astype(np.float64)
    tg = t2 - t
    n = np.array([-tg[1], tg[0]], dtype=np.float64)
    return (n / (np.linalg.norm(n) + 1e-9)).astype(np.float32)


def compute_spine_and_pairs(
    mask_bool: np.ndarray,
    head_xy: np.ndarray | None,
) -> Tuple[
    np.ndarray,
    List[Tuple[np.ndarray, np.ndarray]],
    Tuple[np.ndarray, np.ndarray],
    bool,
]:
    """
    Returns:
      spine_xy (12,2) in ROI pixel coords (x,y) with spine[0]=head
      spine_pairs: 5 pairs (left, right) as (2,) float32 in ROI coords
      (tail_tip, head_tip): nearest contour points to spine ends for legacy rev logic
      head_tail_confusion
    """
    h, w = mask_bool.shape[:2]
    m = (mask_bool.astype(np.uint8)) * 255
    sk = skeletonize(mask_bool.astype(bool)).astype(np.uint8)
    path_pack = _skeleton_path(sk)
    if path_pack is None:
        cxy = np.array(
            [np.mean(np.where(mask_bool)[1]), np.mean(np.where(mask_bool)[0])],
            dtype=np.float32,
        )
        sp = np.tile(cxy, (12, 1))
        return sp, [], (cxy, cxy), True

    path_xy, ep_a, ep_b = path_pack
    # path_xy columns are x,y
    p_a = np.array([float(ep_a[0]), float(ep_a[1])], dtype=np.float32)
    p_b = np.array([float(ep_b[0]), float(ep_b[1])], dtype=np.float32)

    contour = _contour_poly(m)
    if contour is None:
        sp = _resample_path(path_xy, 12)
        return sp, [], (p_a, p_b), True

    def nearest_cont(pt: np.ndarray) -> np.ndarray:
        d2 = np.sum((contour - pt.reshape(1, 2)) ** 2, axis=1)
        return contour[int(np.argmin(d2))].copy()

    # Orient path: head = end closer to previous centroid (ROI coords)
    if head_xy is not None:
        ha = np.sum((p_a - head_xy) ** 2)
        hb = np.sum((p_b - head_xy) ** 2)
        if ha > hb:
            path_xy = path_xy[::-1].copy()
            p_a, p_b = p_b, p_a

    spine12 = _resample_path(path_xy, 12)
    # Contour point nearest physical tail (spine end spine[11])
    tail_tip = nearest_cont(spine12[11])
    head_tip = nearest_cont(spine12[0])

    # spinePairs at mid of segments 0-1, 2-3, 4-5, 6-7, 8-9
    pairs: List[Tuple[np.ndarray, np.ndarray]] = []
    for k in (0, 2, 4, 6, 8):
        mid = 0.5 * (spine12[k] + spine12[k + 1])
        perp = _perpendicular(spine12, k)
        left, right = _width_pair_from_contour(mid, perp, contour)
        pairs.append((left, right))

    d0 = float(np.linalg.norm(spine12[0] - tail_tip))
    d1 = float(np.linalg.norm(spine12[11] - tail_tip))
    confusion = abs(d0 - d1) < 0.25 * max(d0, d1, 1.0)

    return spine12, pairs, (tail_tip, head_tip), confusion
