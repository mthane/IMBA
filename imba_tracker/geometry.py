"""
Spine / contour / metrics from a binary larva mask.

Uses skeleton-based interior path + width pairs (see ``spine_extract``).
"""

from __future__ import annotations

import cv2
import numpy as np

from imba_tracker.cv_angle import cv_angle_from_binary_mask
from imba_tracker.spine_extract import compute_spine_and_pairs
from imba_tracker.types import LarvaGeometry, TrackedBlob


def _largest_contour(mask_u8: np.ndarray) -> np.ndarray | None:
    contours, _ = cv2.findContours(mask_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return None
    return max(contours, key=cv2.contourArea)


def build_contour_csv_order(
    spine12: np.ndarray,
    pairs: list[tuple[np.ndarray, np.ndarray]],
    rev: bool,
) -> np.ndarray:
    """12×2 contour points in legacy ``csvLine`` order (ROI coords)."""
    if len(pairs) != 5:
        return np.zeros((12, 2), dtype=np.float32)
    if rev:
        pts = [spine12[0]]
        for p in pairs:
            pts.append(p[1])
        pts.append(spine12[11])
        for p in reversed(pairs):
            pts.append(p[0])
    else:
        pts = [spine12[11]]
        for p in reversed(pairs):
            pts.append(p[0])
        pts.append(spine12[0])
        for p in pairs:
            pts.append(p[1])
    return np.stack(pts, axis=0).astype(np.float32)


def compute_larva_geometry(
    mask_bool: np.ndarray,
    grey_roi: np.ndarray,
    blob: TrackedBlob,
    prev_centroid_xy: np.ndarray | None,
) -> LarvaGeometry:
    mask_u8 = (mask_bool.astype(np.uint8)) * 255
    cnt = _largest_contour(mask_u8)
    if cnt is None or len(cnt) < 3:
        cxy = np.array(
            [blob.centroid_x - blob.bbox_x, blob.centroid_y - blob.bbox_y], dtype=np.float32
        )
        z = np.tile(cxy, (12, 1))
        return LarvaGeometry(
            spine_xy=z,
            spine_pairs=[],
            contour_csv_xy=z.copy(),
            length_px=1.0,
            width_px=1.0,
            perimeter_px=1.0,
            grey_sum=float(np.sum(grey_roi[mask_bool])) if mask_bool.any() else 0.0,
            orientation_rad=0.0,
            round_flag=1,
            head_tail_confusion=True,
            csv_rev=True,
        )

    head_roi: np.ndarray | None = None
    if prev_centroid_xy is not None:
        head_roi = np.array(
            [
                prev_centroid_xy[0] - blob.bbox_x,
                prev_centroid_xy[1] - blob.bbox_y,
            ],
            dtype=np.float32,
        )

    spine12, pairs, (tail_tip, _head_tip), confusion = compute_spine_and_pairs(
        mask_bool, head_roi
    )

    if confusion or len(pairs) != 5:
        confusion = True

    # legacy ``rev``: spine[0] coincides with tail end (compare to tail_tip on contour)
    d0 = float(np.linalg.norm(spine12[0] - tail_tip))
    d1 = float(np.linalg.norm(spine12[11] - tail_tip))
    rev = d0 <= d1

    contour_csv = build_contour_csv_order(spine12, pairs, rev)

    seg_len = float(
        np.sum(np.linalg.norm(np.diff(spine12, axis=0), axis=1))
    )
    if pairs:
        widths = [
            float(np.linalg.norm(p[0].astype(np.float64) - p[1].astype(np.float64)))
            for p in pairs
        ]
        width_px = float(np.mean(widths))
    else:
        width_px = float(cv2.contourArea(cnt) / max(seg_len, 1.0))

    peri = float(cv2.arcLength(cnt, True))
    grey_sum = float(np.sum(grey_roi[mask_bool])) if mask_bool.any() else 0.0
    orient = cv_angle_from_binary_mask(mask_u8)

    return LarvaGeometry(
        spine_xy=spine12,
        spine_pairs=pairs,
        contour_csv_xy=contour_csv,
        length_px=seg_len,
        width_px=width_px,
        perimeter_px=peri,
        grey_sum=grey_sum,
        orientation_rad=orient,
        round_flag=1 if confusion else 0,
        head_tail_confusion=confusion,
        csv_rev=rev,
    )


def blob_mask_from_labels(
    labels: np.ndarray, grey: np.ndarray, blob: TrackedBlob
) -> tuple[np.ndarray, np.ndarray]:
    """Crop (mask_bool, grey_roi) aligned to blob bbox."""
    x, y, w, h = blob.bbox_x, blob.bbox_y, blob.bbox_w, blob.bbox_h
    lab = labels[y : y + h, x : x + w]
    g = grey[y : y + h, x : x + w]
    m = lab == blob.cc_index
    return m, g
