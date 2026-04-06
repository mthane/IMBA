"""Command-line entry for the Python tracker (intended to be called from R Shiny via system2/processx)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from imba_tracker.config import TrackerConfig
from imba_tracker.pipeline import BackgroundMode, iter_frame_detections, run_ui_style_pipeline


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="IMBA per-frame blob detection (Python MVP; not full C++ lrvTrack). "
        "Designed for non-interactive use (shell, CI, R Shiny)."
    )
    p.add_argument("video", type=Path, help="Input video path (e.g. example_video.mp4)")
    p.add_argument(
        "--experiment-output",
        "--ui-style",
        dest="experiment_output",
        action="store_true",
        help="Write stdout.log, stderr.log, metadata.txt, and <date>-data/detections.csv "
        "(stable folder layout for R / downstream tools; same shape as legacy batch runs).",
    )
    p.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help="With --experiment-output: directory for logs and data (default: same folder as the video).",
    )
    p.add_argument(
        "--small-dish",
        action="store_true",
        help="Use small petri dish metadata (84 mm); default matches big dish (138 mm).",
    )
    p.add_argument("--max-frames", type=int, default=None)
    p.add_argument("--no-auto-dish", action="store_true", help="Disable petri dish ROI detection")

    # JSONL mode (default when not --experiment-output)
    p.add_argument(
        "--jsonl",
        type=Path,
        default=None,
        help="Write one JSON object per line (default: stdout). Ignored with --experiment-output.",
    )
    p.add_argument(
        "--bg",
        choices=("first_frame", "offline_min"),
        default="first_frame",
        help="Background estimate when not using --experiment-output (default: first_frame).",
    )
    p.add_argument("--min-area", type=int, default=20)
    p.add_argument("--max-area", type=int, default=2000)
    p.add_argument(
        "--offline-bg",
        action="store_true",
        help="Use -1 bg weight in process_frame (with --bg offline_min).",
    )

    args = p.parse_args(argv)

    if not args.video.exists():
        print(f"Error: video not found: {args.video}", file=sys.stderr)
        return 2

    if args.experiment_output:
        out_dir = args.output_dir if args.output_dir is not None else args.video.parent
        cfg = TrackerConfig.ui_defaults()
        if args.small_dish:
            cfg.petri_dish_mm = 84
        run_ui_style_pipeline(
            args.video,
            out_dir,
            cfg,
            max_frames=args.max_frames,
            auto_dish=not args.no_auto_dish,
        )
        print(f"Wrote experiment output under: {out_dir.resolve()}", file=sys.stderr)
        return 0

    cfg = TrackerConfig(
        min_obj_size=args.min_area,
        max_obj_size=args.max_area,
        extract_offline_bg=args.offline_bg,
    )
    mode = (
        BackgroundMode.OFFLINE_MIN if args.bg == "offline_min" else BackgroundMode.FIRST_FRAME
    )

    out = open(args.jsonl, "w", encoding="utf-8") if args.jsonl else sys.stdout
    try:
        for fd in iter_frame_detections(
            args.video,
            cfg,
            background_mode=mode,
            max_frames=args.max_frames,
        ):
            row = {
                "frame_index": fd.frame_index,
                "n_blobs": len(fd.blobs),
                "blobs": [
                    {
                        "label": b.label,
                        "area": b.area,
                        "cx": b.centroid_x,
                        "cy": b.centroid_y,
                        "bbox": [b.bbox_x, b.bbox_y, b.bbox_w, b.bbox_h],
                    }
                    for b in fd.blobs
                ],
            }
            out.write(json.dumps(row) + "\n")
    finally:
        if args.jsonl:
            out.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
