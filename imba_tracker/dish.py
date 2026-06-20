"""
Petri dish ROI from offline background, matching extract_background_offline (~3569–3651).

Legacy uses adaptiveThreshold + cvLabel + area filter + roundness; we use
connectedComponents + contour perimeter.
"""

from __future__ import annotations

import math

import cv2
import numpy as np

# lrvTrackOL.cpp ~3610: cvFilterByArea(dishBlob, 900400 ,6292529);
_AREA_MIN = 900_400
_AREA_MAX = 6_292_529


def _roundness(perimeter: float, area: float) -> float:
    return (perimeter * perimeter) / (2.0 * math.pi * max(area, 1.0))


def detect_petri_dish_circle(grey_bg: np.ndarray) -> tuple[float, float, float] | None:
    """
    Returns (cx, cy, radius) in pixel coordinates, or None if no candidate.

    Mirrors: greyBgFrame.copyTo(ctout); normalize; adaptiveThreshold block 355;
    label; area filter; best roundness; bbox center; radius = min(w,h)*0.48.
    """
    if grey_bg.ndim != 2:
        raise ValueError("grey_bg must be single-channel uint8")
    h, w = grey_bg.shape[:2]
    ctout = grey_bg.copy()
    cv2.normalize(ctout, ctout, 0, 255, cv2.NORM_MINMAX)

    block = min(355, max(3, (min(h, w) // 2) * 2 + 1))
    ctout = cv2.adaptiveThreshold(
        ctout,
        255,
        cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
        cv2.THRESH_BINARY_INV,
        block,
        0,
    )

    n, labels, stats, _ = cv2.connectedComponentsWithStats(ctout, connectivity=8)
    best_idx: int | None = None
    best_r = float("inf")
    for i in range(1, n):
        area = int(stats[i, cv2.CC_STAT_AREA])
        if area < _AREA_MIN or area > _AREA_MAX:
            continue
        mask = (labels == i).astype(np.uint8) * 255
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            continue
        cnt = max(contours, key=cv2.contourArea)
        peri = float(cv2.arcLength(cnt, True))
        r = _roundness(peri, area)
        if r < best_r:
            best_r = r
            best_idx = i

    if best_idx is None:
        return None

    x = int(stats[best_idx, cv2.CC_STAT_LEFT])
    y = int(stats[best_idx, cv2.CC_STAT_TOP])
    bw = int(stats[best_idx, cv2.CC_STAT_WIDTH])
    bh = int(stats[best_idx, cv2.CC_STAT_HEIGHT])
    cx = x + bw / 2.0
    cy = y + bh / 2.0
    # lrvTrackOL.cpp ~3633–3634: min span * 0.48 (C++ uses minx twice; use both axes)
    span = min(float(bw), float(bh))
    radius = 0.48 * span
    if int(radius) % 2 == 1:
        radius += 1.0
    return (cx, cy, float(radius))
