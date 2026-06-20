"""Greedy / Hungarian assignment of blob masks to persistent track IDs."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

try:
    from scipy.optimize import linear_sum_assignment
except ImportError:  # pragma: no cover
    linear_sum_assignment = None

from imba_tracker.types import TrackedBlob


@dataclass
class TrackState:
    track_id: int
    bbox: tuple[int, int, int, int]
    mask_crop: np.ndarray  # uint8 0/1, shape matches bbox
    centroid_xy: np.ndarray  # shape (2,) float, full-frame


@dataclass
class AssociationEngine:
    _next_id: int = 1
    tracks: dict[int, TrackState] = field(default_factory=dict)

    def _set_track(self, tid: int, blob: TrackedBlob, mask_crop: np.ndarray) -> None:
        self.tracks[tid] = TrackState(
            track_id=tid,
            bbox=(blob.bbox_x, blob.bbox_y, blob.bbox_w, blob.bbox_h),
            mask_crop=(mask_crop.astype(np.uint8) > 0).astype(np.uint8),
            centroid_xy=np.array([blob.centroid_x, blob.centroid_y], dtype=np.float64),
        )

    def associate(
        self,
        blobs: list[TrackedBlob],
        masks: list[np.ndarray],
        cost_max: float = 0.92,
        min_iou: float = 0.02,
    ) -> dict[int, int]:
        """
        Returns mapping blob_index -> track_id for each blob in this frame.
        """
        assert len(blobs) == len(masks)
        if not blobs:
            self.tracks.clear()
            return {}

        t_ids = list(self.tracks.keys())
        n_t, n_b = len(t_ids), len(blobs)
        matches: dict[int, int] = {}

        if n_t == 0:
            for i, b in enumerate(blobs):
                tid = self._next_id
                self._next_id += 1
                self._set_track(tid, b, masks[i])
                matches[i] = tid
            return matches

        cost = np.full((n_t, n_b), 1e6, dtype=np.float64)
        for ti, tid in enumerate(t_ids):
            ts = self.tracks[tid]
            for bi, blob in enumerate(blobs):
                iou = _mask_iou(ts.bbox, ts.mask_crop, blob, masks[bi])
                if iou >= min_iou:
                    cost[ti, bi] = 1.0 - float(iou)

        if linear_sum_assignment is not None:
            ri, ci = linear_sum_assignment(cost)
            used_b: set[int] = set()
            for ti, bi in zip(ri, ci):
                if cost[ti, bi] < cost_max:
                    matches[bi] = t_ids[ti]
                    used_b.add(bi)
        else:
            matches = _greedy_match(cost, t_ids)
            used_b = set(matches.keys())

        matched_tids = set()
        for bi, b in enumerate(blobs):
            if bi in matches:
                tid = matches[bi]
                self._set_track(tid, b, masks[bi])
                matched_tids.add(tid)
            else:
                tid = self._next_id
                self._next_id += 1
                self._set_track(tid, b, masks[bi])
                matches[bi] = tid
                matched_tids.add(tid)

        self.tracks = {k: v for k, v in self.tracks.items() if k in matched_tids}
        return matches


def _greedy_match(cost: np.ndarray, t_ids: list[int]) -> dict[int, int]:
    """Fallback greedy by smallest cost."""
    matches: dict[int, int] = {}
    used_rows: set[int] = set()
    used_cols: set[int] = set()
    flat = []
    for ti in range(cost.shape[0]):
        for bi in range(cost.shape[1]):
            flat.append((float(cost[ti, bi]), ti, bi))
    flat.sort(key=lambda x: x[0])
    for c, ti, bi in flat:
        if c >= 1e5:
            break
        if ti in used_rows or bi in used_cols:
            continue
        used_rows.add(ti)
        used_cols.add(bi)
        matches[bi] = t_ids[ti]
    return matches


def _mask_iou(
    box_a: tuple[int, int, int, int],
    mask_a: np.ndarray,
    blob: TrackedBlob,
    mask_b: np.ndarray,
) -> float:
    ax, ay, aw, ah = box_a
    bx, by, bw, bh = blob.bbox_x, blob.bbox_y, blob.bbox_w, blob.bbox_h
    ux = min(ax, bx)
    uy = min(ay, by)
    u2x = max(ax + aw, bx + bw)
    u2y = max(ay + ah, by + bh)
    uw, uh = u2x - ux, u2y - uy
    if uw <= 0 or uh <= 0:
        return 0.0
    canvas_a = np.zeros((uh, uw), dtype=np.uint8)
    canvas_b = np.zeros((uh, uw), dtype=np.uint8)
    _paste(canvas_a, mask_a, ax - ux, ay - uy)
    _paste(canvas_b, mask_b, bx - ux, by - uy)
    inter = np.logical_and(canvas_a > 0, canvas_b > 0).sum()
    union = np.logical_or(canvas_a > 0, canvas_b > 0).sum()
    if union == 0:
        return 0.0
    return float(inter) / float(union)


def _paste(canvas: np.ndarray, patch: np.ndarray, ox: int, oy: int) -> None:
    h, w = patch.shape[:2]
    H, W = canvas.shape[:2]
    y0 = max(0, oy)
    x0 = max(0, ox)
    y1 = min(H, oy + h)
    x1 = min(W, ox + w)
    if y0 >= y1 or x0 >= x1:
        return
    py0 = y0 - oy
    px0 = x0 - ox
    canvas[y0:y1, x0:x1] = np.maximum(
        canvas[y0:y1, x0:x1], patch[py0 : py0 + (y1 - y0), px0 : px0 + (x1 - x0)]
    )
