# IMBA 1.0 — Individual Maggot Behavior Analyzer

IMBA tracks individual *Drosophila* larvae in a group, including through collisions, and analyzes the resulting trajectories. 
Please cite this software using the original Paper: doi.org/10.1098/rsob.220308.
It was used for Toshima et al., 2026, *eLife* (DOI: [10.7554/eLife.106148.3](https://doi.org/10.7554/eLife.106148.3)).

Three programs make up a run:

1. **lrvTrack** tracks a video and writes one CSV per larva, plus an overlay video.
2. **TrackExperiment.py** is the PyQt window for a whole experiment. It calls the same tracker and arranges the files for analysis.
3. **IMBAvisualizer** is the Shiny app that reads those CSVs and computes behavioral measures.

This repository has no `LICENSE` file. The tracker build unpacks `IMBAtracker/IMBtracker_cpp/cvblob.zip`, which carries its own notices.

## What you need

- Docker Engine and Docker Compose v2.
- On Windows, [Docker Desktop](https://docs.docker.com/desktop/setup/install/windows-install/) with the WSL2 backend. Run the commands below from PowerShell or from WSL; both use the same `compose.yaml`.
- Your raw video. Put a copy at `data/input/example_video.mp4`, or set `IMBA_VIDEO` to another filename in `data/input`. Videos are not stored in git.

All commands below are run from the repository root.

## Track one video

Build the tracker image once:

```bash
docker compose build tracker
```

The image is `imba-lrvtrack:1.0` (Ubuntu 18.04, OpenCV 3.2). Then track the video in `data/input`:

```bash
docker compose run --rm tracker
```

That runs:

```text
lrvTrack -x -z 138 --min-obj-size 110 --max-obj-size 5000 -p --thread-count 9 -o -v 16 -d 13 -w 0.04 -u -t -r . -i /data/input/example_video.mp4
```

| Flag | Meaning |
|---|---|
| `-z` | Petri-dish diameter in millimetres. `138` is the 15 cm dish. `84` is the 9 cm dish. |
| `-v` | Frames per second. |
| `-i` | Video path inside the container. |
| `--min-obj-size`, `--max-obj-size` | Blob area limits, in pixels. |

Files are exchanged through these folders:

| On your computer | In the container | Access |
|---|---|---|
| `data/input` | `/data/input` | read-only |
| `data/output` | `/data/output` | read/write |
| `config` | `/app/config` | read-only |

`data/output` must exist before you run. Results appear there, named with the time of the run:

- `YYYYMMDD_HHMMSS_PROC.mp4` — the video with the dish circle, frame number, and larva ids.
- `YYYYMMDD_HHMMSS-data/<id>.csv` — one file per larva that was tracked as an individual.

Each CSV row is one frame: frame index, spine and body-outline coordinates in millimetres (origin at the dish centre), centroid, ellipse angle, area in pixels, grey value, length, width, perimeter, a round flag, and millimetres per pixel. A coordinate pair can be `na` when head and tail could not be ordered.

Larvae whose area is above 6 mm² are treated as clusters. They are drawn on the overlay and omitted from the CSV files. The millimetre scale comes from `-z`, so the diameter has to match the dish in the video. For a 9 cm dish, pass the small-dish preset:

```bash
docker compose run --rm tracker -x -z 84 --min-obj-size 110 --max-obj-size 5000 -p --thread-count 9 -o -v 16 -d 13 -w 0.04 -u -t -r . -i /data/input/example_video.mp4
```

PowerShell does not continue lines with a backslash. Put that command on one line.

To point at videos stored outside this repository:

```bash
export IMBA_INPUT_DIR=/data/larva/raw
export IMBA_OUTPUT_DIR=/data/larva/tracked
export IMBA_VIDEO=trial01.mp4
docker compose run --rm tracker
```

```powershell
$env:IMBA_INPUT_DIR = "D:\larva\raw"
$env:IMBA_OUTPUT_DIR = "D:\larva\tracked"
$env:IMBA_VIDEO = "trial01.mp4"
docker compose run --rm tracker
```

On Linux, output files are owned by root unless you set your user id:

```bash
IMBA_UID="$(id -u)" IMBA_GID="$(id -g)" docker compose run --rm tracker
```

`config/metadata.txt` is a sample trial file (odor side and frame rate). The command above does not read it. If you add `--metadata-file`, copy the file into `data/output` first: when the odor side is `left` or `right`, lrvTrack writes the cup coordinates back into that file.

## Track an experiment from the GUI

Use the GUI when you have the usual experiment tree and want files laid out for IMBAvisualizer. Build the tracker image first (`docker compose build tracker`). The window calls `IMBAtracker/IMBtracker_cpp/lrvTrack`, which starts that image. A native `make` in `IMBtracker_cpp` replaces this launcher with the compiled program; check the file back out from git if you need the Docker launcher again.

On Windows, start the window from WSL (Ubuntu 22.04) with Python, PyQt5, OpenCV, and Docker available (`python3-pyqt5`, `python3-opencv`). WSLg provides the display:

```bash
cd IMBAtracker
export DISPLAY=:0
export WAYLAND_DISPLAY=wayland-0
export XDG_RUNTIME_DIR=/mnt/wslg/runtime-dir
export QT_QPA_PLATFORM=xcb
python3 TrackExperiment.py
```

Place each video four folders deep. Those folder names become the experiment, group, reciprocal, and trial:

```text
data/gui_experiment/ExampleExp/GroupA/Rec1/Trial1/example_video.mp4
data/gui_experiment/ExampleExp/GroupA/Rec1/Trial1/metadata.txt
```

In the window, add that video, choose an output folder such as `data/gui_output`, and pick the dish. **15 cm** is 138 mm and is selected by default. **9 cm (Normal)** is 84 mm. Press **Track!!**. One video is tracked at a time.

When tracking finishes, the GUI writes:

```text
data/gui_output/GroupA/Rec1/Trial1/1.csv
data/gui_output/GroupA/Rec1/Trial1/2.csv
data/gui_output/GroupA/Rec1/Trial1/vidAndLogs/metadata.txt
data/gui_output/GroupA/Rec1/Trial1/vidAndLogs/<timestamp>_PROC.mp4
```

It then runs `createCollisionsCSVinExperiment.sh` on the tracker log. That step needs `Rscript` on the machine where the GUI is running. The CSVs and the overlay are already in place if that script stops.

## Analyze tracks

IMBAvisualizer reads the folder tree the GUI writes. Build it once and start it:

```bash
docker compose build visualizer
docker compose up visualizer
```

Open [http://localhost:3838](http://localhost:3838). The image is `imba-visualizer:1.0`. The app lives in `IMBAvisualizer/shiny_app`. An older packaged copy of the same scripts is kept in `IMBAvisualizer/legacy` and is not what this service starts.

In **Analysis Processing**:

1. Leave **Base directory** at `/data/` (that is the `data` folder of this repository, inside the container). `IMBA_DATA_DIR` changes which host folder is mounted there.
2. Choose `gui_output`. That directory is the parent of the group folders. The trial folder and `data` itself are the wrong level.
3. Set the dish radius to **42.5 mm** for a 9 cm dish, or **69 mm** for a 15 cm dish.
4. Set **Frames per second** to the rate used in tracking (16 for the bundled example).
5. Set **Video length** to the length of the recording in seconds. The default, 180, matches a 3-minute video.
6. Press **Analyze!**.

The app looks for:

```text
<experiment>/<group>/<condition>/<trial>/<id>.csv
<experiment>/<group>/<condition>/<trial>/vidAndLogs/metadata.txt
```

`metadata.txt` is an INI file. When tracking was given `--metadata-file` and the odor side is left or right, lrvTrack appends `OdorALocation=x, y`. Without that line the odor position used in the analysis is `(0, 0)`.

**Analyze!** keeps a track that is longer than half a second (`round(frames per second / 0.5)` frames). On the **Accepted tracks** tab the threshold is in seconds and starts at 90. Lower it when the recording is shorter than that, or the tracks are listed as not accepted.

A command-line run of the tracker does not create the group/condition/trial folders. Copy `YYYYMMDD_HHMMSS-data/*.csv` into a `group/condition/trial/` directory and put `metadata.txt` in that trial's `vidAndLogs` before opening the app on that folder.

Results are written next to the selected experiment as `<name>_analyzed.csv` and `<name>_analyzed.rds`, on the host under `data/`.

To run the same analysis without the browser, for the GUI output of a 5-second, 16 fps, 9 cm recording:

```bash
docker compose run --rm --entrypoint Rscript visualizer -e 'source("global.R"); dir.create("/data/analysis/example", recursive=TRUE, showWarnings=FALSE); process_data_dt(input_dir="/data/gui_output", output_dir="/data/analysis/example", save_as="example", radius=42.5, drop_contours=TRUE, frame_rate=16, video_length=5)'
```

## If something goes wrong

- **Error setting up the capture device.** The video is not in `/data/input`, or a Windows path such as `C:\...` was passed into the container. Use `-i /data/input/<name>.mp4`.
- **Overlay video, and no CSV folder.** Every detected object was above the 6 mm² cluster cutoff for the `-z` you chose. Check the overlay, then run again with the diameter of the dish you filmed (84 or 138).
- **`H264` / `avc1` on the tracker log.** OpenCV writes the mp4 with the `avc1` codec. The file is still produced.
- **`set_int: Column 0 out of range`.** Printed by the assignment solver while it searches. A finished run continues and prints `FINISHED!!!`.
- **Output files owned by root.** Set `IMBA_UID` and `IMBA_GID` as shown above.
- **`/data/output` is not writable.** Create `data/output` on the host before `docker compose run`.

## Install without Docker

`IMBAtracker/IMBtracker_cpp/setup.sh` installs the tracker libraries on Debian or Ubuntu, unpacks `cvblob.zip`, and compiles `lrvTrack` in place. It asks for confirmation at the `apt` prompts. The compile also needs `libtbb-dev`. Compiling in that directory replaces the GUI launcher, so use the compiled `lrvTrack` directly or restore the launcher from git afterwards.

`IMBAtracker/setup.sh` starts the PyQt GUI after that build. `IMBAvisualizer/legacy/setup.sh` is the older R-package install of the visualizer. The Docker commands above are the straightforward way to run both.
