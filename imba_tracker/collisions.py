"""Heuristic merge / split / touch events between consecutive frames."""

from __future__ import annotations

from imba_tracker.types import CollisionEvent, TrackedBlob, TrackedLarvaView


def detect_collision_events(
    prev: list[TrackedLarvaView],
    curr: list[TrackedLarvaView],
    frame_index: int,
    *,
    bbox_iou_gate: float = 0.05,
) -> list[CollisionEvent]:
    """
    - **merge**: two or more larvae in ``prev``, one blob in ``curr``, previous centroids
      fall inside the current blob bbox (expanded heuristic).
    - **split**: one larva in ``prev``, two+ blobs in ``curr``, both new centroids inside
      previous bbox.
    - **touch**: two larvae in ``curr`` with overlapping bboxes (IoU on bbox).
    """
    events: list[CollisionEvent] = []
    if not prev:
        return events

    def _bbox(b: TrackedBlob) -> tuple[int, int, int, int]:
        return (b.bbox_x, b.bbox_y, b.bbox_w, b.bbox_h)

    def _centroid_in_bbox(
        cx: float, cy: float, bb: tuple[int, int, int, int]
    ) -> bool:
        x, y, w, h = bb
        return x <= cx <= x + w and y <= cy <= y + h

    def _bbox_iou(a: tuple[int, int, int, int], b: tuple[int, int, int, int]) -> float:
        ax, ay, aw, ah = a
        bx, by, bw, bh = b
        x1 = max(ax, bx)
        y1 = max(ay, by)
        x2 = min(ax + aw, bx + bw)
        y2 = min(ay + ah, by + bh)
        iw = max(0, x2 - x1)
        ih = max(0, y2 - y1)
        inter = iw * ih
        area_a = aw * ah
        area_b = bw * bh
        union = area_a + area_b - inter
        if union <= 0:
            return 0.0
        return float(inter) / float(union)

    # merge
    if len(prev) >= 2 and len(curr) == 1:
        bb = _bbox(curr[0].blob)
        cents = [(lv.blob.centroid_x, lv.blob.centroid_y) for lv in prev]
        if all(_centroid_in_bbox(cx, cy, bb) for cx, cy in cents):
            events.append(
                CollisionEvent(
                    frame_index=frame_index,
                    kind="merge",
                    track_ids=tuple(sorted({lv.track_id for lv in prev})),
                    detail="multi_prev_to_one_blob",
                )
            )

    # split
    if len(prev) == 1 and len(curr) >= 2:
        bb = _bbox(prev[0].blob)
        if all(
            _centroid_in_bbox(lv.blob.centroid_x, lv.blob.centroid_y, bb)
            for lv in curr
        ):
            events.append(
                CollisionEvent(
                    frame_index=frame_index,
                    kind="split",
                    track_ids=tuple(sorted({lv.track_id for lv in curr})),
                    detail="one_prev_to_multi_blob",
                )
            )

    # touch (pairwise bbox IoU)
    n = len(curr)
    for i in range(n):
        for j in range(i + 1, n):
            iou = _bbox_iou(_bbox(curr[i].blob), _bbox(curr[j].blob))
            if iou >= bbox_iou_gate:
                events.append(
                    CollisionEvent(
                        frame_index=frame_index,
                        kind="touch",
                        track_ids=(curr[i].track_id, curr[j].track_id),
                        detail=f"bbox_iou={iou:.4f}",
                    )
                )

    return events
