"""IMBA tracking: video → preprocess → detections; optional multi-frame tracks + overlay MP4."""

from imba_tracker.config import TrackerConfig
from imba_tracker.pipeline import (
    BackgroundMode,
    iter_frame_detections,
    iter_frame_tracking,
    run_experiment_pipeline,
    run_tracked_experiment_with_overlay,
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
    "iter_frame_tracking",
    "run_video_detections",
    "run_ui_style_pipeline",
    "run_experiment_pipeline",
    "run_tracked_experiment_with_overlay",
]
