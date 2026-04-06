"""Serialize one CSV line compatible with legacy ``larvaObject::csvLine`` column order."""

from __future__ import annotations

from imba_tracker.types import LarvaGeometry, TrackedBlob


def _fmt(v: float) -> str:
    return f"{v:.8f}"


def format_larva_csv_line(
    frame_index: int,
    blob: TrackedBlob,
    geom: LarvaGeometry,
    dish_center_xy: tuple[float, float],
    ppm: float,
) -> str:
    """
    dish_center_xy: (cx, cy) petri dish center in pixels (full frame).
    ppm: pixels per mm.
    """
    ccx, ccy = dish_center_xy
    ox, oy = float(blob.bbox_x), float(blob.bbox_y)

    def mm_spine(sx: float, sy: float) -> tuple[str, str]:
        return _fmt((sx - ccx) * ppm), _fmt((-sy + ccy) * ppm)

    def mm_contour(cx: float, cy: float) -> tuple[str, str]:
        return _fmt((cx - ccx) * ppm), _fmt((-cy + ccy) * ppm)

    parts: list[str] = [str(int(frame_index))]

    if geom.head_tail_confusion:
        for _ in range(12):
            parts.extend(["na", "na"])
    else:
        d = geom.spine_xy
        rev = geom.csv_rev
        if rev:
            for i in range(12):
                sx = d[i, 0] + ox
                sy = d[i, 1] + oy
                parts.extend(mm_spine(sx, sy))
        else:
            for i in range(11, -1, -1):
                sx = d[i, 0] + ox
                sy = d[i, 1] + oy
                parts.extend(mm_spine(sx, sy))

    if geom.head_tail_confusion or len(geom.spine_pairs) != 5:
        for _ in range(12):
            parts.extend(["na", "na"])
    else:
        sp = geom.spine_pairs
        s = geom.spine_xy
        rev = geom.csv_rev
        if rev:
            parts.extend(mm_contour(s[0, 0] + ox, s[0, 1] + oy))
            for p in sp:
                parts.extend(mm_contour(p[1][0] + ox, p[1][1] + oy))
            parts.extend(mm_contour(s[11, 0] + ox, s[11, 1] + oy))
            for p in reversed(sp):
                parts.extend(mm_contour(p[0][0] + ox, p[0][1] + oy))
        else:
            parts.extend(mm_contour(s[11, 0] + ox, s[11, 1] + oy))
            for p in reversed(sp):
                parts.extend(mm_contour(p[0][0] + ox, p[0][1] + oy))
            parts.extend(mm_contour(s[0, 0] + ox, s[0, 1] + oy))
            for p in sp:
                parts.extend(mm_contour(p[1][0] + ox, p[1][1] + oy))

    cx = blob.centroid_x
    cy = blob.centroid_y
    parts.append(_fmt((cx - ccx) * ppm))
    parts.append(_fmt((cy - ccy) * ppm))
    parts.append(_fmt(geom.orientation_rad))
    parts.append(str(int(blob.area)))
    parts.append(_fmt(geom.grey_sum))
    parts.append(_fmt(geom.length_px * ppm))
    parts.append(_fmt(geom.width_px * ppm))
    parts.append(_fmt(geom.perimeter_px * ppm))
    parts.append(str(int(geom.round_flag)))
    parts.append(_fmt(ppm))
    return ",".join(parts) + "\n"
