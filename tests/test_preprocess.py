import numpy as np

from imba_tracker.config import TrackerConfig
from imba_tracker.detect import extract_blobs
from imba_tracker.frame_input import preprocess_raw_frame
from imba_tracker.preprocess import process_frame


def test_process_frame_bright_blob_on_uniform_bg():
    """Synthetic: region brighter than bg; fg = new + (-4)*bg must stay > 0 in blob."""
    cfg = TrackerConfig()
    bg = np.full((64, 64), 30, dtype=np.uint8)
    new = np.full((64, 64), 30, dtype=np.uint8)
    new[20:40, 20:40] = 200
    proc = process_frame(new, bg, cfg, color_frame=None)
    assert proc.shape == (64, 64)
    blobs = extract_blobs(proc, cfg)
    assert len(blobs) >= 1
    b0 = max(blobs, key=lambda b: b.area)
    assert b0.area > 50


def test_preprocess_mean_clamp_and_normalize():
    cfg = TrackerConfig()
    raw = np.zeros((32, 32, 3), dtype=np.uint8)
    raw[5:15, 5:15] = 200
    g, color, prev = preprocess_raw_frame(raw, cfg, None)
    assert g.shape == (32, 32)
    assert color.shape == (32, 32, 3)
    assert prev.shape == raw.shape
    assert g.max() == 255


def test_extract_blobs_filters_area():
    cfg = TrackerConfig(min_obj_size=100, max_obj_size=200)
    proc = np.zeros((50, 50), dtype=np.uint8)
    proc[10:12, 10:12] = 200  # area 4 — too small
    assert extract_blobs(proc, cfg) == []
