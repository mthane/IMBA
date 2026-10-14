# IMBA 1.0 — Individual Maggot Behavior Analyzer

IMBA tracks individual *Drosophila* larvae in a group, including through collisions, and was used for Toshima et al., 2026, *eLife* (DOI: [10.7554/eLife.106148.3](https://doi.org/10.7554/eLife.106148.3)).

The tracker is the original C++ program `lrvTrack` in `IMBAtracker/IMBtracker_cpp`. This recovery runs that program, unmodified, in Docker. It does not implement the unfinished Python package `imba_tracker` that exists on `origin/main`. That package detects blobs per frame and does not link identities or resolve collisions.

There is no `LICENSE` file in this repository. One was not added here; that choice belongs to the rights holders. The build unpacks the bundled `cvblob` sources (`IMBAtracker/IMBtracker_cpp/cvblob.zip`), which carry their own notices.

## Prerequisites

- Docker Engine and Docker Compose v2.
- On Windows, [Docker Desktop](https://docs.docker.com/desktop/setup/install/windows-install/) with the WSL2 backend. The same `compose.yaml` works if you invoke `docker` from PowerShell or from a WSL distribution.
- A raw video on the host. This recovery was checked with `data/input/example_video.mp4` (80 frames, 2048×2048, 16 fps), taken from the repository history (`legacy/data/example_video.mp4` on `origin/main`). It is not stored in git.

The image is Ubuntu 18.04, the system named in the original install notes and in the historical Dockerfile. That release supplies the OpenCV 3.2 headers this source was written against (`CV_CAP_PROP_*`, `CV_FOURCC`, and related names).

| Package | Version in the image | Role |
|---|---|---|
| g++ / cmake | GCC 7.4, CMake 3.10 | Build |
| libopencv-dev | 3.2.0+dfsg-4ubuntu0.1 | Video I/O and image processing |
| libboost-all-dev | 1.65.1 | Program options, filesystem, log |
| libalglib-dev | 3.11.0 | Spine fitting |
| liblpsolve55-dev | 5.5.0.15 | Identity assignment |
| libtbb-dev | 2017 U7 | Parallel blob update |
| cvblob (bundled zip) | sources in `cvblob.zip` | Connected components |

`libtbb-dev` is required by `CMakeLists.txt` (`find_package(TBB REQUIRED)`). The original `setup.sh` does not install it. No MATLAB, GPU, or display is required. The batch command does not open a window. The PyQt program `IMBAtracker/TrackExperiment.py` is only a wrapper that shells out to `lrvTrack`.

## Build

From the repository root:

```bash
docker compose build
```

The image is tagged `imba-lrvtrack:1.0`. Videos are not copied into the image.

## Input and output mounts

| Host | Container | Access |
|---|---|---|
| `./data/input` | `/data/input` | read-only |
| `./data/output` | `/data/output` | read/write |
| `./config` | `/app/config` | read-only |

Put the video in `data/input`. `lrvTrack` writes its result files in the process working directory. The container entrypoint `cd`s to `/data/output` before exec, so the dated overlay and CSV directory appear on the host. The input mount stays read-only.

Large collections can stay outside the repository:

```bash
# Linux / WSL
export IMBA_INPUT_DIR=/data/larva/raw
export IMBA_OUTPUT_DIR=/data/larva/tracked
export IMBA_VIDEO=trial01.mp4
docker compose run --rm tracker
```

```powershell
# Windows PowerShell
$env:IMBA_INPUT_DIR = "D:\larva\raw"
$env:IMBA_OUTPUT_DIR = "D:\larva\tracked"
$env:IMBA_VIDEO = "trial01.mp4"
docker compose run --rm tracker
```

On Linux, files created in the output mount are owned by the container user (root unless you override it):

```bash
IMBA_UID="$(id -u)" IMBA_GID="$(id -g)" docker compose run --rm tracker
```

Docker Desktop on Windows already exposes those files to your Windows user.

`config/metadata.txt` is the sample trial file that was stored next to the example video (odor side and FPS). The default command does not read it. If you pass `--metadata-file`, copy that file into the writable output directory first: with an odor side of `left` or `right`, `lrvTrack` writes cup coordinates back into the metadata file, and `/app/config` is read-only.

## Run the example video

```bash
docker compose build
docker compose run --rm tracker
```

That is the command recorded in the original README, with container paths:

```text
lrvTrack -x -z 138 --min-obj-size 110 --max-obj-size 5000 -p --thread-count 9 -o -v 16 -d 13 -w 0.04 -u -t -r . -i /data/input/example_video.mp4
```

`-z` is the petri-dish diameter in millimetres. The original GUI also has a small-dish preset of 84 mm. On this example clip, `-z 138` scales every segmented animal above the program's 6 mm² cluster cutoff (`blob area × mm_per_pixel²`), so `lrvTrack` keeps four identities and draws them on the overlay, but it does not write spine CSV files. `-z 84` puts three of the four under that cutoff. Their lengths are then about 4.1–4.5 mm, which is a plausible larval body length. The fourth object stays above 6 mm² and is still treated as a cluster. To emit those trajectories:

```bash
docker compose run --rm tracker \
  -x -z 84 --min-obj-size 110 --max-obj-size 5000 \
  -p --thread-count 9 -o -v 16 -d 13 -w 0.04 -u -t -r . \
  -i /data/input/example_video.mp4
```

PowerShell does not use backslash continuation. Pass the same arguments on one line.

Dish diameter changes the millimetre scale. It does not change the tracking source.

## Expected output

In `data/output` (names use the run's local timestamp):

- `YYYYMMDD_HHMMSS_PROC.mp4` — input frames with the detected dish circle, frame index, and larva ids.
- `YYYYMMDD_HHMMSS-data/<id>.csv` — one file per animal that the original rules did not mark as a cluster.

Each CSV row is one frame. Fields, in order:

1. Frame index.
2. Spine and body-outline coordinates in millimetres, origin at the dish centre, y flipped so image-down is negative.
3. Centroid x, y (mm), ellipse angle, area (pixels), grey value, length, width, perimeter (mm), round flag, and millimetres per pixel.

A row can contain `na` pairs when head and tail cannot be ordered. Animals classified as clusters are omitted from this CSV; their ids are still drawn on the overlay when the tag pass runs.

Logs from a redirected run (`stdout.log`, `stderr.log`) are optional and are not written by `lrvTrack` itself.

## What was checked

On `example_video.mp4` the unmodified binary finished both passes (`FINISHED!!!`).

With `-z 138` it found one dish circle (centre 935, 1027 px, radius 634.56 px, 0.1087 mm/px) and kept identities 1–4 for all 80 frames. No spine CSV was written, because each blob was over 6 mm². That is the original cutoff, not a failed write.

With `-z 84` (0.0662 mm/px) it wrote `1.csv`, `2.csv`, and `3.csv`, 80 rows each. Centroids stay inside the 42 mm dish radius, lengths are about 4.1–4.5 mm, widths about 1 mm, and the recorded millimetres-per-pixel match `84 / (2 × 634.56)`. Identity 4 remained a cluster and has no CSV.

No historical result files were in this branch to compare against. A log saved on `origin/main` (`legacy/tracker_cpp/stdout.log`, 2024-05-13) stops at "Error setting up the capture device" after a Windows path (`C:/Users/.../example_video.mp4`) was given to the Linux binary. This recovery does not claim numerical equivalence to a published figure.

## Troubleshooting

- **Error setting up the capture device.** The file is missing from `/data/input`, or a Windows path was passed into the container. Use `-i /data/input/<name>.mp4`.
- **No CSV directory.** The run still tracked blobs, but every object was over the 6 mm² cluster cutoff for the chosen `-z`. Inspect the overlay, then rerun with the other original dish preset (84 or 138) if the physical dish matches it.
- **`H264` / `avc1` message on stderr.** Ubuntu 18.04's FFmpeg build of OpenCV rewrites the H264 fourcc to `avc1` and still writes the mp4. The process continues.
- **`set_int: Column 0 out of range`.** Printed by lp_solve during the original assignment step. The run used here still reported `LPSOLVE: Found 1 solutions` and `RES: OK`.
- **Output files owned by root.** Set `IMBA_UID` and `IMBA_GID` as shown above.
- **`/data/output` is not writable.** The output directory must exist on the host before `docker compose run`.

## Native build

`IMBAtracker/IMBtracker_cpp/setup.sh` is the original Debian/Ubuntu script: it installs the libraries above, unpacks `cvblob.zip`, and runs `cmake` and `make` in place. It calls `apt` without `-y`, so it is interactive. `IMBAtracker/setup.sh` then starts the PyQt GUI. Prefer the Docker image for a non-interactive run. The Shiny app under `IMBAvisualizer/` reads tracker tables; it is not required to produce tracks.
