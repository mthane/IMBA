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

`libtbb-dev` is required by `CMakeLists.txt` (`find_package(TBB REQUIRED)`). The original `setup.sh` does not install it. No MATLAB or GPU is required. The batch command does not open a window. The PyQt program `IMBAtracker/TrackExperiment.py` is a wrapper that shells out to `lrvTrack`; how to launch it is in [Track through the GUI](#track-through-the-gui). The Shiny app in `IMBAvisualizer/shiny_app` reads the resulting tables; it does not track. See [IMBAvisualizer](#imbavisualizer).

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

## Track through the GUI

`IMBAtracker/TrackExperiment.py` is the original PyQt5 window. It does not link larvae itself. It builds the same `lrvTrack` command and runs `IMBAtracker/IMBtracker_cpp/lrvTrack`. That path is a small launcher script: it `docker run`s the `imba-lrvtrack:1.0` binary, with the trial directory as the working directory, and passes the GUI's arguments through. Build the tracker image before opening the window (`docker compose build tracker`). A native `make` in `IMBtracker_cpp` overwrites the launcher with an ELF binary; restore the script from git if that happens.

The window was run from WSL (Ubuntu 22.04) with `python3-pyqt5`, `python3-opencv`, and Docker. WSLg supplies the display:

```bash
cd IMBAtracker
export DISPLAY=:0
export WAYLAND_DISPLAY=wayland-0
export XDG_RUNTIME_DIR=/mnt/wslg/runtime-dir
export QT_QPA_PLATFORM=xcb
python3 TrackExperiment.py
```

`XDG_RUNTIME_DIR` is owned by the normal WSL user. A Qt warning is printed if the process is root; the window still opens.

Point the file list at a video nested four directories deep. The folder names become Experiment, Group, Reciprocal, and Trial:

```text
data/gui_experiment/ExampleExp/GroupA/Rec1/Trial1/example_video.mp4
data/gui_experiment/ExampleExp/GroupA/Rec1/Trial1/metadata.txt
```

Choose an output directory, for example `data/gui_output`. The dish radios are the original presets: **15 cm** is 138 mm (the default), **9 cm (Normal)** is 84 mm. On `example_video.mp4`, 9 cm is the preset that writes spine CSVs. **Track!!** runs one video at a time.

After a zero exit code the GUI copies `metadata.txt` and the logs into `output/Group/Reciprocal/Trial/vidAndLogs`, moves `*-data/*.csv` up into the trial folder, and moves the overlay `2*.mp4` into `vidAndLogs`. That tree is what IMBAvisualizer reads. A checked run of this clip wrote `data/gui_output/GroupA/Rec1/Trial1/1.csv`, `2.csv`, and `3.csv`, plus the overlay and logs.

The GUI then runs `bash createCollisionsCSVinExperiment.sh` on the stdout log. Those shell scripts must be LF, not CRLF, or bash reports `unexpected end of file`. The script's last step calls `Rscript`, which is not part of the tracker image. Tracking has already finished by then; a missing `Rscript` does not change the CSVs.

## IMBAvisualizer

The app that runs is `IMBAvisualizer/shiny_app` (`ui.R`, `server.R`, `global.R`, and `r_files/`). `global.R` sources those scripts directly. The old golem package (`R/`, `DESCRIPTION`, `NAMESPACE`, and `imba::run_app()`) is kept under `IMBAvisualizer/legacy` and is not part of this pipeline. The app does not call `lrvTrack`. It reads a directory of digit-named CSVs (`1.csv`, `2.csv`, …) laid out as:

```text
<experiment>/<group>/<condition>/<trial>/<id>.csv
<experiment>/<group>/<condition>/<trial>/vidAndLogs/metadata.txt
```

The GUI output above matches this when the chosen directory is the parent of the group folder (`data/gui_output`, not the trial folder and not `data/` itself). The CLI writer does not build that tree. To plot a CLI run, copy `YYYYMMDD_HHMMSS-data/*.csv` into a `group/condition/trial/` folder and put `metadata.txt` in that trial's `vidAndLogs/`.

`metadata.txt` is INI. Cup coordinates are the sixth data row (`OdorALocation=x, y`), which `lrvTrack` appends when `--metadata-file` is set and the odor side is left or right. The sample file has no coordinates, so the odor position stays `(0, 0)`.

Each CSV row ends with millimetres-per-pixel after the collision flag. The reader keeps the columns it names (through `collision_flag`) and ignores that extra field. Contour columns 26–69 are the body outline; **Drop contours** leaves the 12 spine points.

Dish radius presets in the app are **42.5 mm** (the 84 mm / 9 cm dish) and **69 mm** (the 138 mm dish). Frames per second defaults to 16. **Video length** defaults to 180 seconds, which matches a 3-minute recording and is stored with the result. **Analyze!** keeps a track when it has more than `round(FPS / 0.5)` frames (32 frames at 16 fps). The 80-frame example passes. The later accepted-track view treats the threshold as seconds (default 90). Set that threshold below the clip length, or those 5-second tracks are counted as not accepted.

Build and open the app:

```bash
docker compose build visualizer
docker compose up visualizer
```

Then open `http://localhost:3838`. The base directory defaults to `/data/` inside the container (`IMBA_DATA_ROOT`). The whole `./data` tree is mounted there (`IMBA_DATA_DIR` overrides the host path). Choose `gui_output`, set the radius to 42.5 for the 9 cm example, and leave frames per second at 16. Analyzed CSV and RDS files are written under the selected parent, on the host mount.

The same function the button calls can be run without the window, against the GUI output of this clip:

```bash
docker compose run --rm --entrypoint Rscript visualizer -e \
  'source("global.R"); dir.create("/data/analysis/example", recursive=TRUE, showWarnings=FALSE); process_data_dt(input_dir="/data/gui_output", output_dir="/data/analysis/example", save_as="example", radius=42.5, drop_contours=TRUE, frame_rate=16, video_length=5)'
```

That runs in `/opt/imba/shiny_app`, the image working directory, so `global.R` finds `r_files/` and `global_variables.Rds`. `IMBAvisualizer/legacy/setup.sh` is the old interactive Ubuntu 18.04 + R 4.0 package install. The image is the non-interactive path for `shiny_app`. It uses R 4.5 so current CRAN builds of the libraries in `global.R` install; the analysis scripts are the ones in `shiny_app/r_files`.

## Troubleshooting

- **Error setting up the capture device.** The file is missing from `/data/input`, or a Windows path was passed into the container. Use `-i /data/input/<name>.mp4`.
- **No CSV directory.** The run still tracked blobs, but every object was over the 6 mm² cluster cutoff for the chosen `-z`. Inspect the overlay, then rerun with the other original dish preset (84 or 138) if the physical dish matches it.
- **`H264` / `avc1` message on stderr.** Ubuntu 18.04's FFmpeg build of OpenCV rewrites the H264 fourcc to `avc1` and still writes the mp4. The process continues.
- **`set_int: Column 0 out of range`.** Printed by lp_solve during the original assignment step. The run used here still reported `LPSOLVE: Found 1 solutions` and `RES: OK`.
- **Output files owned by root.** Set `IMBA_UID` and `IMBA_GID` as shown above.
- **`/data/output` is not writable.** The output directory must exist on the host before `docker compose run`.

## Native build

`IMBAtracker/IMBtracker_cpp/setup.sh` is the original Debian/Ubuntu script: it installs the libraries above, unpacks `cvblob.zip`, and runs `cmake` and `make` in place. It calls `apt` without `-y`, so it is interactive, and `make` replaces the GUI launcher script with the compiled binary. `IMBAtracker/setup.sh` then starts the PyQt GUI. Prefer the Docker image for a non-interactive run, and the [GUI](#track-through-the-gui) and [visualizer](#imbavisualizer) sections to drive the same outputs from the window and from Shiny.
