"""
Foreground segmentation matching lrvTrackOL.cpp process_frame (~2648–2774).

Steps:
  addWeighted (bg subtract) → normalize MINMAX → mask by dish circle ROI (or full frame)
  → Otsu threshold → dilate MORPH_CROSS 3x3 → processed = fgImage & binary_mask
"""

from __future__ import annotations

import cv2
import numpy as np

from imba_tracker.config import TrackerConfig


def process_frame(
    new_frame: np.ndarray,
    grey_bg: np.ndarray,
    cfg: TrackerConfig,
    color_frame: np.ndarray | None = None,
) -> np.ndarray:
    """
    new_frame, grey_bg: uint8 grey, same shape.
    color_frame: optional BGR; if provided and dish circle is set, draws green circle outline.
    Returns processedFrame (uint8 grey, larvae regions non-zero).
    """
    if new_frame.shape != grey_bg.shape:
        raise ValueError("new_frame and grey_bg must match shape")
    if cfg.extract_offline_bg:
        beta = cfg.bg_subtract_bg_weight_offline
    else:
        beta = cfg.bg_subtract_bg_weight

    fg_frame = cv2.addWeighted(new_frame, 1.0, grey_bg, beta, 0.0)
    cv2.normalize(fg_frame, fg_frame, 0, 255, cv2.NORM_MINMAX)

    h, w = fg_frame.shape[:2]
    fg_roi = np.zeros((h, w), dtype=np.uint8)

    if cfg.dish_circle_xyr is not None:
        cx, cy, r = cfg.dish_circle_xyr
        cv2.circle(fg_roi, (int(cx), int(cy)), int(r), 255, -1)
        if color_frame is not None and color_frame.shape[:2] == (h, w):
            cv2.circle(color_frame, (int(cx), int(cy)), int(r), (0, 255, 0), 1)
    else:
        fg_roi[:, :] = 255

    fg_image = cv2.bitwise_and(fg_frame, fg_roi)

    _, binary = cv2.threshold(fg_image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    element = cv2.getStructuringElement(cv2.MORPH_CROSS, (3, 3))
    binary = cv2.dilate(binary, element)
    processed = cv2.bitwise_and(fg_image, binary)
    return processed
