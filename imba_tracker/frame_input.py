"""
Frame acquisition + intensity preprocessing matching lrvTrackOL.cpp get_next_frame (Mat overload).

Legacy: lrvTrackOL.cpp get_next_frame ~3373–3464 (BGR→gray, optional invert, brightnessContrastGamma,
per-pixel mean clamp, NORM_MINMAX, GRAY2BGR for color output).
"""

from __future__ import annotations

import cv2
import numpy as np

from imba_tracker.config import TrackerConfig


def brightness_contrast_gamma(
    o: np.ndarray, b: float, c: float, g: float
) -> None:
    """lrvTrackOL.cpp brightnessContrastGamma ~656–668. In-place on uint8 grey."""
    if c != 1.0 or b != 0.0:
        # addWeighted(o, 0, o, c, b - 255, o)
        o[:] = cv2.addWeighted(o, 0, o, c, b - 255)
    if g != 1.0:
        o_f = o.astype(np.float32)
        np.power(o_f, g, out=o_f)
        np.clip(o_f, 0, 255, out=o_f)
        o[:] = o_f.astype(np.uint8)


def preprocess_raw_frame(
    bgr_or_gray: np.ndarray,
    cfg: TrackerConfig,
    previous_raw: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """
    Returns:
        grey uint8 HxW
        color BGR uint8 HxWx3 (visualization; grey replicated to 3 channels like C++)
        updated previous_raw for corruption check (same shape/channels as capture.read)
    """
    raw = np.array(bgr_or_gray, copy=True)
    if previous_raw is not None and raw.shape == previous_raw.shape:
        diff = cv2.absdiff(bgr_or_gray, previous_raw)
        if np.linalg.norm(diff) > raw.shape[0] * raw.shape[1] * cfg.corruption_norm_threshold:
            raw = np.array(previous_raw, copy=True)

    if raw.ndim == 2:
        grey = raw
    else:
        grey = cv2.cvtColor(raw, cv2.COLOR_BGR2GRAY)

    if cfg.invert:
        grey = cv2.addWeighted(grey, 1, grey, -2, 255)
        grey = cv2.convertScaleAbs(grey.astype(np.float32))

    brightness_contrast_gamma(grey, cfg.brightness, cfg.contrast, cfg.gamma)

    mean_val = float(np.mean(grey))
    grey = np.maximum(grey, mean_val).astype(np.uint8)

    cv2.normalize(grey, grey, 0, 255, cv2.NORM_MINMAX)

    color = cv2.cvtColor(grey, cv2.COLOR_GRAY2BGR)

    return grey, color, raw
