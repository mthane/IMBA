"""
Connected components + area filter replacing cvblob cvLabel + cvFilterByArea.

Legacy: lrvTrackOL.cpp main loop ~5717–5748 (cvLabel, cvFilterByArea, larvaToRing optional).

cvblob connectivity: we use OpenCV 8-connectivity (default) for components on mask > 0.
"""

from __future__ import annotations

import cv2
import numpy as np

from imba_tracker.config import TrackerConfig
from imba_tracker.types import TrackedBlob


def extract_blobs(processed_frame: np.ndarray, cfg: TrackerConfig) -> list[TrackedBlob]:
    """Labels connected components on processed_frame > 0."""
    mask = (processed_frame > 0).astype(np.uint8)
    n, labels, stats, centroids = cv2.connectedComponentsWithStats(
        mask, connectivity=8
    )
    out: list[TrackedBlob] = []
    next_label = 1
    for i in range(1, n):
        area = int(stats[i, cv2.CC_STAT_AREA])
        if area < cfg.min_obj_size or area > cfg.max_obj_size:
            continue
        x = int(stats[i, cv2.CC_STAT_LEFT])
        y = int(stats[i, cv2.CC_STAT_TOP])
        w = int(stats[i, cv2.CC_STAT_WIDTH])
        h = int(stats[i, cv2.CC_STAT_HEIGHT])
        cx, cy = centroids[i]
        out.append(
            TrackedBlob(
                label=next_label,
                area=area,
                centroid_x=float(cx),
                centroid_y=float(cy),
                bbox_x=x,
                bbox_y=y,
                bbox_w=w,
                bbox_h=h,
            )
        )
        next_label += 1
    return out
