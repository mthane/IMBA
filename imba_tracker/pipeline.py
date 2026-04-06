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

from imba_tracker.association import AssociationEngine
from imba_tracker.config import TrackerConfig
from imba_tracker.detect import extract_blobs, extract_blobs_with_labels
from imba_tracker.dish import detect_petri_dish_circle
from imba_tracker.frame_input import preprocess_raw_frame
from imba_tracker.geometry import blob_mask_from_labels, compute_larva_geometry
from imba_tracker.legacy_csv import format_larva_csv_line
from imba_tracker.overlay_video import OverlayVideoWriter
from imba_tracker.preprocess import process_frame
from imba_tracker.export_ui import (
    UiExportPaths,
    prepare_ui_paths,
    write_detections_csv,
    write_metadata_stub,
    write_stderr_empty,
    write_stdout_log,
)
from imba_tracker.types import FrameDetections, FrameTracking


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


def dish_center_and_ppm(
    cfg: TrackerConfig,
    dish: tuple[float, float, float] | None,
    frame_hw: tuple[int, int],
) -> tuple[tuple[float, float], float]:
    """Petri dish center (px) and pixels-per-mm for legacy CSV scaling."""
    h, w = frame_hw
    if dish is not None:
        cx, cy, r = dish
        ppm = float(cfg.petri_dish_mm) / (2.0 * float(r))
        return (cx, cy), ppm
    return (w / 2.0, h / 2.0), float(cfg.petri_dish_mm) / float(max(min(w, h), 1))


def iter_frame_tracking(
    path: Path | str,
    cfg: TrackerConfig | None = None,
    *,
    background_mode: BackgroundMode = BackgroundMode.FIRST_FRAME,
    grey_background: np.ndarray | None = None,
    max_frames: int | None = None,
) -> Generator[FrameTracking, None, None]:
    """
    Per-frame blob detection + mask IoU association + spine/contour geometry.
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

    engine = AssociationEngine()
    prev_centroids: dict[int, np.ndarray] = {}
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
        blobs, labels = extract_blobs_with_labels(proc, cfg)
        masks: list[np.ndarray] = []
        for b in blobs:
            m, _ = blob_mask_from_labels(labels, grey, b)
            masks.append(m.astype(np.uint8))

        matches = engine.associate(blobs, masks)
        larvae: list[TrackedLarvaView] = []
        for bi, b in enumerate(blobs):
            tid = matches[bi]
            prev_xy = prev_centroids.get(tid)
            m, g_roi = blob_mask_from_labels(labels, grey, b)
            geom = compute_larva_geometry(m, g_roi, b, prev_xy)
            prev_centroids[tid] = np.array(
                [b.centroid_x, b.centroid_y], dtype=np.float64
            )
            larvae.append(
                TrackedLarvaView(track_id=tid, blob=b, geometry=geom)
            )

        yield FrameTracking(
            frame_index=frame_index,
            color_bgr=color,
            grey=grey,
            larvae=larvae,
        )
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


def run_tracked_experiment_with_overlay(
    video_path: Path | str,
    output_dir: Path | str,
    cfg: TrackerConfig | None = None,
    *,
    max_frames: int | None = None,
    auto_dish: bool = True,
    odor_side: str = "left",
    overlay_filename: str = "overlay_test.mp4",
) -> UiExportPaths:
    """
    Offline-min background, dish ROI, per-larva ``<id>.csv`` (legacy column order),
    aggregate ``detections.csv``, logs, and an MP4 with spine/contour overlays.
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

    cap = cv2.VideoCapture(str(video_path))
    if not cap.isOpened():
        raise FileNotFoundError(f"Cannot open video: {video_path}")
    w = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    h = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    fps = float(cap.get(cv2.CAP_PROP_FPS)) or 30.0
    cap.release()

    dish_center, ppm = dish_center_and_ppm(cfg, dish, (h, w))

    paths = prepare_ui_paths(output_dir)
    overlay_path = paths.data_dir / overlay_filename
    paths.overlay_mp4 = overlay_path

    fd_list: list[FrameDetections] = []
    csv_handles: dict[int, object] = {}
    try:
        vw = OverlayVideoWriter(overlay_path, (w, h), fps)
        for ft in iter_frame_tracking(
            video_path,
            cfg,
            grey_background=bg,
            max_frames=max_frames,
        ):
            fd_list.append(
                FrameDetections(
                    frame_index=ft.frame_index,
                    blobs=[lv.blob for lv in ft.larvae],
                )
            )
            vw.write_frame(ft)
            for lv in ft.larvae:
                tid = lv.track_id
                p = paths.data_dir / f"{tid}.csv"
                if tid not in csv_handles:
                    csv_handles[tid] = p.open("a", encoding="utf-8")
                csv_handles[tid].write(
                    format_larva_csv_line(
                        ft.frame_index,
                        lv.blob,
                        lv.geometry,
                        dish_center,
                        ppm,
                    )
                )
        vw.release()
    finally:
        for f in csv_handles.values():
            f.close()

    n_obs = sum(len(f.blobs) for f in fd_list)
    write_detections_csv(paths.detections_csv, fd_list)
    lines = [
        "IMBA Python tracker — tracked larvae with legacy-shaped per-ID CSV rows and overlay MP4.",
        f"Video: {video_path.resolve()}",
        f"Overlay: {overlay_path.resolve()}",
        f"Petri dish mm: {cfg.petri_dish_mm}",
    ]
    if dish:
        lines.append(
            f"Petri dish ROI (cx, cy, r): {dish[0]:.2f}, {dish[1]:.2f}, {dish[2]:.2f}"
        )
    lines.extend(
        [
            f"Frames processed: {len(fd_list)}",
            f"Total blob observations (detections rows): {n_obs}",
            f"Data directory: {paths.data_dir}",
            "Exit: OK",
        ]
    )
    paths.stdout_log.write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_metadata_stub(
        paths.metadata_txt, odor_side=odor_side, petri_dish_mm=cfg.petri_dish_mm
    )
    write_stderr_empty(paths.output_dir / "stderr.log")

    return paths


# Preferred name for callers (R Shiny, scripts); kept alias for older imports.
run_experiment_pipeline = run_ui_style_pipeline
