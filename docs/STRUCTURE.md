# IMBA 2.0 layout

- **legacy/** — IMBA 1.x C++ tracker (two trees: `tracker_cpp` preferred for builds, `lrv_track` duplicate), Qt GUI, bash/R glue, Docker, sample data.
- **imba_tracker/** — installable Python package (`pip install -e .`).
- **shiny_app/** — R package `imba` (Golem Shiny UI).
- **shared/** — future JSON/YAML schemas and small fixtures for Python ↔ R contracts.
- **evaluation/** — future benchmarks comparing Python vs legacy outputs.
- **scripts/** — thin wrappers (e.g. run Python CLI on example video).

Build the Docker C++ image from the **repository root**: `docker build -f legacy/docker/Dockerfile .`
