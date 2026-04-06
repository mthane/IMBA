"""Defaults mirror global flags in legacy/tracker_cpp/lrvTrack.hpp (where applicable)."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class TrackerConfig:
    # lrvTrack.hpp
    min_obj_size: int = 20
    max_obj_size: int = 2000
    extract_offline_bg: bool = False
    # process_frame: addWeighted(new, 1, bg, -4, 0) vs -1 when offline
    bg_subtract_bg_weight: float = -4.0
    bg_subtract_bg_weight_offline: float = -1.0
    # get_next_frame corruption gate: norm(absdiff) > cols*rows*0.095
    corruption_norm_threshold: float = 0.095
    # Dish ROI: if None, use full frame (C++ fgROI = all 255 when circles.size()==0)
    dish_circle_xyr: tuple[float, float, float] | None = None
    # Frame input (lrvTrack.hpp defaults)
    invert: bool = False
    brightness: float = 0.0
    contrast: float = 1.0
    gamma: float = 1.0
    # Offline background: up to 100 frames, step 10 (extract_background_offline)
    offline_bg_max_frames: int = 100
    offline_bg_frame_step: int = 10
    # TrackExperiment.py: big dish 138 mm, small 84 (informational for metadata / future use)
    petri_dish_mm: int = 138
    # Multi-frame association (Hungarian on 1 - IoU)
    association_cost_max: float = 0.92
    association_min_iou: float = 0.02
    # Heuristic collision events (bbox overlap / centroid containment)
    collision_bbox_iou_gate: float = 0.05

    @classmethod
    def ui_defaults(cls) -> "TrackerConfig":
        """Defaults for ``--experiment-output`` / legacy batch (min area 110, offline BG, etc.)."""
        return cls(
            min_obj_size=110,
            max_obj_size=5000,
            extract_offline_bg=True,
            offline_bg_max_frames=100,
            offline_bg_frame_step=10,
            petri_dish_mm=138,
        )
