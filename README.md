# IMBA 2.0: Individual Maggot Behavior Analyzer

Research codebase for larval tracking and analysis: **legacy** C++/Qt stack (preserved from IMBA 1.x), a **Python** tracker (`imba_tracker`), and an **R/Shiny** app (`shiny_app` package `imba`).

## Repository layout

| Path | Role |
|------|------|
| `legacy/tracker_cpp/` | Canonical C++ `lrvTrack` build (CMake, OpenCV, cvblob). |
| `legacy/lrv_track/` | Second copy of the C++ tree (historical duplicate; pick one for builds). |
| `legacy/qt_gui/` | PyQt5 UI (`TrackUI.py`, `.ui` files). |
| `legacy/scripts/` | Bash + R helpers (e.g. collision CSV pipeline). |
| `legacy/docker/` | Dockerfile for C++ tracker image (build from repo root). |
| `legacy/data/` | Example media (e.g. `example_video.mp4`). |
| `imba_tracker/` | Python package: `pip install -e .` from repository root. |
| `tests/` | Python tests (pytest). |
| `shiny_app/` | R Golem/Shiny visualizer; install with `devtools::install()` from this folder. |
| `shared/` | Placeholder for cross-language schemas and example payloads. |
| `evaluation/` | Placeholder for benchmarks vs legacy outputs. |
| `scripts/` | Helper scripts for the Python tracker (e.g. `run_python_tracker_example.*`). |

## Python tracker

From the repository root:

```bash
pip install -e .
pytest
```

Run a quick detection demo (uses `legacy/data/example_video.mp4` by default):

```bash
bash scripts/run_python_tracker_example.sh
```

On Windows PowerShell:

```powershell
.\scripts\run_python_tracker_example.ps1
```

## Legacy C++ tracker (Linux / WSL)

Build dependencies and compile (typical Ubuntu/WSL):

```bash
cd legacy/tracker_cpp
sudo bash setup.sh
```

Run `lrvTrack` (see `legacy/docs/tracker_command.txt` for example flags).

Docker (context = repository root):

```bash
docker build -f legacy/docker/Dockerfile -t imba-lrvtrack .
```

## Qt GUI

Legacy GUI sources live under `legacy/qt_gui/`. Install a Python environment with PyQt5 and OpenCV, then run `TrackUI.py` from that directory or extend `PYTHONPATH` accordingly.

## Shiny / R visualizer

```bash
cd shiny_app
sudo bash setup.sh
```

The setup script expects to be run from `shiny_app/` (package root). In R, `devtools::install()` then `imba::run_app()`.

## Legacy IMBA 1.0

The original monolithic layout is archived separately; this repository is structured for migration to Python + Shiny without mixing build artifacts into new code.
