"""Orientation angle compatible with cvblob-style central second moments (radians)."""

from __future__ import annotations

import math

import cv2
import numpy as np


def cv_angle_from_binary_mask(mask_u8: np.ndarray) -> float:
    """
    Principal axis angle (radians), consistent with common cvblob ``cvAngle`` implementations.
    """
    m = cv2.moments(mask_u8)
    if m["m00"] < 1e-9:
        return 0.0
    mu20 = m["mu20"]
    mu02 = m["mu02"]
    mu11 = m["mu11"]
    return 0.5 * math.atan2(2.0 * mu11, mu20 - mu02)
