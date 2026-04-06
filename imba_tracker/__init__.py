"""Minimal IMBA tracking core: video → preprocess → per-frame detections (no ID linking).

Call from Python or use ``python -m imba_tracker.cli`` (e.g. from R Shiny via ``system2`` / ``processx``).
"""

from imba_tracker.config import TrackerConfig
from imba_tracker.pipeline import (
    BackgroundMode,
    iter_frame_detections,
    run_experiment_pipeline,
    run_ui_style_pipeline,
    run_video_detections,
)
from imba_tracker.types import FrameDetections, TrackedBlob

__all__ = [
    "TrackerConfig",
    "BackgroundMode",
    "FrameDetections",
    "TrackedBlob",
    "iter_frame_detections",
    "run_video_detections",
    "run_ui_style_pipeline",
    "run_experiment_pipeline",
]
