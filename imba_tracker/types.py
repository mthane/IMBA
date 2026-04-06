from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class TrackedBlob:
    """One connected component after area filtering (per-frame label, not track ID)."""

    label: int
    cc_index: int
    area: int
    centroid_x: float
    centroid_y: float
    bbox_x: int
    bbox_y: int
    bbox_w: int
    bbox_h: int


@dataclass
class FrameDetections:
    frame_index: int
    blobs: list[TrackedBlob] = field(default_factory=list)


@dataclass
class LarvaGeometry:
    """
    Pixel-space geometry in ROI coordinates (relative to blob bbox origin).

    ``spine_pairs`` holds 5 (left, right) contour points for legacy CSV contour block.
    """

    spine_xy: np.ndarray  # (12, 2) float32, head → tail along interior path
    spine_pairs: list[tuple[np.ndarray, np.ndarray]]  # len 5, each (2,) xy
    contour_csv_xy: np.ndarray  # (12, 2) ordered for CSV (legacy spinePairs layout)
    length_px: float
    width_px: float
    perimeter_px: float
    grey_sum: float  # L1 sum in mask (matches C++ getGreyValue style)
    orientation_rad: float  # cvAngle / moments
    round_flag: int
    head_tail_confusion: bool
    csv_rev: bool  # legacy ``rev`` in larvaObject::csvLine


@dataclass
class TrackedLarvaView:
    """One larva in a frame after ID assignment."""

    track_id: int
    blob: TrackedBlob
    geometry: LarvaGeometry


@dataclass
class FrameTracking:
    frame_index: int
    color_bgr: np.ndarray
    grey: np.ndarray
    larvae: list[TrackedLarvaView] = field(default_factory=list)


@dataclass(frozen=True)
class CollisionEvent:
    """Merge/split/touch between frames (research heuristics)."""

    frame_index: int
    kind: str  # "merge" | "split" | "touch"
    track_ids: tuple[int, ...]
    detail: str = ""
