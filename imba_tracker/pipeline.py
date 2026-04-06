"""
Video → background estimate → per-frame process_frame → blob extraction.

This intentionally does NOT implement full extract_background() (Hough dish, cups,
bg_without_larvae, etc.). The CLI is the supported integration surface for R Shiny.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Generator

import cv2
import numpy as np

from imba_tracker.config import TrackerConfig
from imba_tracker.detect import extract_blobs
from imba_tracker.dish import detect_petri_dish_circle
from imba_tracker.frame_input import preprocess_raw_frame
from imba_tracker.preprocess import process_frame
from imba_tracker.export_ui import (
    UiExportPaths,
    prepare_ui_paths,
    write_detections_csv,
    write_metadata_stub,
    write_stderr_empty,
    write_stdout_log,
)
from imba_tracker.types import FrameDetections


class BackgroundMode(str, Enum):
    """How to obtain grey_bg for process_frame."""

    FIRST_FRAME = "first_frame"
    OFFLINE_MIN = "offline_min"


def _read_raw(cap: cv2.VideoCapture) -> np.ndarray | None:
    ok, frame = cap.read()
    if not ok or frame is None:
        return None
    return frame


def build_background_offline_min(
    path: Path | str,
    cfg: TrackerConfig,
) -> tuple[np.ndarray, tuple[int, int]]:
    """
    Matches extract_background_offline with LRVTRACK_EXTRACT_OFFLINEBG_MIN=true (~3519–3548):
    get_next_frame(..., step) reads one frame then advances by step (grab step-1 times).
    Up to offline_bg_max_frames samples, each at 0, step, 2*step, ... frame indices.
    """
    cap = cv2.VideoCapture(str(path))
    if not cap.isOpened():
        raise FileNotFoundError(f"Cannot open video: {path}")
    w = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    h = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    result = np.full((h, w), 255, dtype=np.uint8)
    sample_count = 0
    prev_for_corrupt: np.ndarray | None = None
    step = max(1, cfg.offline_bg_frame_step)
    while sample_count < cfg.offline_bg_max_frames:
        ok, raw = cap.read()
        if not ok:
            break
        grey, _, prev_for_corrupt = preprocess_raw_frame(raw, cfg, prev_for_corrupt)
        result = np.minimum(result, grey)
        sample_count += 1
        for _ in range(step - 1):
            if not cap.grab():
                break
    cap.release()
    return result, (w, h)


def build_background_first_frame(
    cap: cv2.VideoCapture, cfg: TrackerConfig
) -> tuple[np.ndarray, np.ndarray | None]:
    """First preprocessed grey frame as static background."""
    raw = _read_raw(cap)
    if raw is None:
        raise RuntimeError("Empty video")
    grey, _, prev = preprocess_raw_frame(raw, cfg, None)
    return grey, prev


def iter_frame_detections(
    path: Path | str,
    cfg: TrackerConfig | None = None,
    *,
    background_mode: BackgroundMode = BackgroundMode.FIRST_FRAME,
    grey_background: np.ndarray | None = None,
    max_frames: int | None = None,
) -> Generator[FrameDetections, None, None]:
    """
    Yields FrameDetections for each processed frame.

    If background_mode is FIRST_FRAME, the first decoded frame is used only as the
    grey background (matching C++: extract_background consumes one frame, then the
    main loop reads the rest).

    If grey_background is provided, background_mode is ignored for building bg.
    """
    cfg = cfg or TrackerConfig()
    cap = cv2.VideoCapture(str(path))
    if not cap.isOpened():
        raise FileNotFoundError(f"Cannot open video: {path}")

    prev_raw: np.ndarray | None = None

    if grey_background is not None:
        bg = grey_background.copy()
    elif background_mode == BackgroundMode.OFFLINE_MIN:
        cap.release()
        bg, _ = build_background_offline_min(path, cfg)
        cap = cv2.VideoCapture(str(path))
        if not cap.isOpened():
            raise FileNotFoundError(f"Cannot reopen video: {path}")
    else:
        bg, prev_raw = build_background_first_frame(cap, cfg)

    n_out = 0
    frame_index = 0
    while True:
        raw = _read_raw(cap)
        if raw is None:
            break
        if max_frames is not None and n_out >= max_frames:
            break
        grey, color, prev_raw = preprocess_raw_frame(raw, cfg, prev_raw)
        proc = process_frame(grey, bg, cfg, color_frame=color)
        blobs = extract_blobs(proc, cfg)
        yield FrameDetections(frame_index=frame_index, blobs=blobs)
        frame_index += 1
        n_out += 1

    cap.release()


def run_video_detections(
    path: Path | str,
    cfg: TrackerConfig | None = None,
    **kwargs,
) -> list[FrameDetections]:
    return list(iter_frame_detections(path, cfg, **kwargs))


def run_ui_style_pipeline(
    video_path: Path | str,
    output_dir: Path | str,
    cfg: TrackerConfig | None = None,
    *,
    max_frames: int | None = None,
    auto_dish: bool = True,
    odor_side: str = "left",
) -> UiExportPaths:
    """
    Write ``stdout.log``, ``stderr.log``, ``metadata.txt``, and
    ``<timestamp>-data/detections.csv`` under ``output_dir`` (for R Shiny / CLI).

    Uses offline min background, default ``TrackerConfig.ui_defaults()`` blob sizes,
    and optional dish ROI from the same adaptiveThreshold path as C++ offline extract.
    """
    video_path = Path(video_path)
    output_dir = Path(output_dir)
    cfg = cfg or TrackerConfig.ui_defaults()

    bg, _ = build_background_offline_min(video_path, cfg)
    dish: tuple[float, float, float] | None = None
    if auto_dish:
        dish = detect_petri_dish_circle(bg)
        if dish:
            cfg.dish_circle_xyr = dish

    frames = list(
        iter_frame_detections(
            video_path,
            cfg,
            grey_background=bg,
            max_frames=max_frames,
        )
    )

    paths = prepare_ui_paths(output_dir)
    n_obs = sum(len(f.blobs) for f in frames)
    write_detections_csv(paths.detections_csv, frames)
    write_stdout_log(
        paths.stdout_log,
        video_path=str(video_path.resolve()),
        n_frames=len(frames),
        n_blobs_total=n_obs,
        dish=dish,
        bg_samples=cfg.offline_bg_max_frames,
        bg_step=cfg.offline_bg_frame_step,
        data_dir=paths.data_dir,
    )
    write_metadata_stub(
        paths.metadata_txt, odor_side=odor_side, petri_dish_mm=cfg.petri_dish_mm
    )
    write_stderr_empty(paths.output_dir / "stderr.log")

    return paths


# Preferred name for callers (R Shiny, scripts); kept alias for older imports.
run_experiment_pipeline = run_ui_style_pipeline
