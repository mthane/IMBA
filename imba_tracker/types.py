from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class TrackedBlob:
    """One connected component after area filtering (per-frame label, not track ID)."""

    label: int
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
