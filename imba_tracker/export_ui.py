"""
Structured experiment output for non-interactive use (CLI, R Shiny via system2/processx).

Writes stdout.log, stderr.log, metadata.txt, and <DATE>-data/detections.csv.
Full C++ lrvTrack also emits per-larva CSVs and CNTTRACK lines; MVP writes aggregate detections only.
"""

from __future__ import annotations

import csv
import time
from dataclasses import dataclass
from pathlib import Path

from imba_tracker.types import FrameDetections


def make_date_stamp() -> str:
    """Match lrvTrackOL.cpp strftime(cDate,16,\"%Y%m%d_%H%M%S\",...)."""
    return time.strftime("%Y%m%d_%H%M%S")


@dataclass
class UiExportPaths:
    output_dir: Path
    data_dir: Path
    detections_csv: Path
    stdout_log: Path
    metadata_txt: Path
    overlay_mp4: Path | None = None


def prepare_ui_paths(output_dir: Path, date_stamp: str | None = None) -> UiExportPaths:
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stamp = date_stamp or make_date_stamp()
    data_dir = output_dir / f"{stamp}-data"
    data_dir.mkdir(parents=True, exist_ok=True)
    return UiExportPaths(
        output_dir=output_dir,
        data_dir=data_dir,
        detections_csv=data_dir / "detections.csv",
        stdout_log=output_dir / "stdout.log",
        metadata_txt=output_dir / "metadata.txt",
    )


def write_metadata_stub(
    path: Path, *, odor_side: str = "left", petri_dish_mm: int = 138
) -> None:
    """INI fields referenced by createCollisionsCSVinExperiment.sh / columnize script."""
    path.write_text(
        "[Trial Data]\n"
        f"OdorA={odor_side}\n"
        "OdorALocation=0, 0\n"
        f"PetriDishSize_mm={petri_dish_mm}\n",
        encoding="utf-8",
    )


def write_detections_csv(
    path: Path,
    frames: list[FrameDetections],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "frame_index",
                "blob_index",
                "centroid_x",
                "centroid_y",
                "area",
                "bbox_x",
                "bbox_y",
                "bbox_w",
                "bbox_h",
            ]
        )
        for fd in frames:
            for b in fd.blobs:
                w.writerow(
                    [
                        fd.frame_index,
                        b.label,
                        f"{b.centroid_x:.4f}",
                        f"{b.centroid_y:.4f}",
                        b.area,
                        b.bbox_x,
                        b.bbox_y,
                        b.bbox_w,
                        b.bbox_h,
                    ]
                )


def write_stdout_log(
    path: Path,
    *,
    video_path: str,
    n_frames: int,
    n_blobs_total: int,
    dish: tuple[float, float, float] | None,
    bg_samples: int,
    bg_step: int,
    data_dir: Path,
) -> None:
    lines = [
        "IMBA Python tracker (MVP) — per-frame blob detections.",
        "NOTE: This is not the full C++ lrvTrack binary. Outputs do not include",
        "CNTTRACK lines, per-larva spine CSVs, or collision-resolution MODEL lines.",
        "Call the C++ lrvTrack build when you need full pipeline parity (CNTTRACK, per-larva CSVs).",
        "",
        f"Video: {video_path}",
        f"Offline background: min over up to {bg_samples} samples, frame step {bg_step}",
    ]
    if dish:
        lines.append(f"Petri dish ROI (cx, cy, r): {dish[0]:.2f}, {dish[1]:.2f}, {dish[2]:.2f}")
    else:
        lines.append("Petri dish ROI: not detected (full-frame foreground mask)")
    lines.extend(
        [
            f"Frames processed: {n_frames}",
            f"Total blob observations (rows): {n_blobs_total}",
            f"Data directory: {data_dir}",
            "Exit: OK",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_stderr_empty(path: Path) -> None:
    path.write_text("", encoding="utf-8")
