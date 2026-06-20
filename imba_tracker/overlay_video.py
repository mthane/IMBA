"""Draw spine + contour overlays and write MP4."""

from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np

from imba_tracker.types import FrameTracking


def _bright_colors(n: int) -> list[tuple[int, int, int]]:
    base = [
        (0, 255, 0),
        (255, 128, 0),
        (0, 128, 255),
        (255, 0, 255),
        (0, 255, 255),
        (255, 255, 0),
        (128, 255, 0),
        (255, 0, 128),
    ]
    return [base[i % len(base)] for i in range(n)]


def draw_tracking_overlay(ft: FrameTracking) -> np.ndarray:
    """Return BGR image copy with polylines and track IDs."""
    img = ft.color_bgr.copy()
    colors = _bright_colors(max(8, len(ft.larvae)))
    for i, lv in enumerate(ft.larvae):
        col = colors[i % len(colors)]
        g = lv.geometry
        ox, oy = lv.blob.bbox_x, lv.blob.bbox_y
        spine = (g.spine_xy + np.array([ox, oy], dtype=np.float32)).astype(np.int32)
        cont = (g.contour_csv_xy + np.array([ox, oy], dtype=np.float32)).astype(np.int32)
        if len(spine) >= 2:
            cv2.polylines(img, [spine.reshape(-1, 1, 2)], False, col, 2, cv2.LINE_AA)
        if len(cont) >= 2:
            cv2.polylines(
                img,
                [np.vstack([cont, cont[:1]]).reshape(-1, 1, 2)],
                True,
                (int(col[0] // 2), int(col[1] // 2), 255),
                1,
                cv2.LINE_AA,
            )
        cx, cy = int(lv.blob.centroid_x), int(lv.blob.centroid_y)
        cv2.putText(
            img,
            str(lv.track_id),
            (cx + 4, cy - 4),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.6,
            (255, 255, 255),
            2,
            cv2.LINE_AA,
        )
        cv2.putText(
            img,
            str(lv.track_id),
            (cx + 4, cy - 4),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.6,
            col,
            1,
            cv2.LINE_AA,
        )
    return img


class OverlayVideoWriter:
    def __init__(self, path: Path, size: tuple[int, int], fps: float) -> None:
        self.path = Path(path)
        self._w, self._h = size
        self._fps = max(fps, 1.0)
        fourcc = cv2.VideoWriter_fourcc(*"mp4v")
        self._vw = cv2.VideoWriter(
            str(self.path), fourcc, self._fps, (self._w, self._h), True
        )
        if not self._vw.isOpened():
            raise RuntimeError(f"Could not open VideoWriter for {self.path}")

    def write_frame(self, ft: FrameTracking) -> None:
        vis = draw_tracking_overlay(ft)
        if vis.shape[1] != self._w or vis.shape[0] != self._h:
            vis = cv2.resize(vis, (self._w, self._h), interpolation=cv2.INTER_AREA)
        self._vw.write(vis)

    def release(self) -> None:
        self._vw.release()
