# Install IMBA on Windows

This guide is for biology researchers and students using Windows 10 or Windows 11. You do not need previous experience with Linux, Docker, or the command line.

IMBA tracks individual *Drosophila* larvae in a group, including through collisions, and writes a trajectory file for each larva. The tracking program runs inside a Docker container so that different computers use the same system libraries. You do not install those libraries on Windows yourself.

You will install four things:

| Name | What it is |
|---|---|
| **WSL2** | A Windows feature that runs Linux. |
| **Ubuntu 24.04** | The Linux system you type commands in. |
| **Docker Engine** | The program that builds and runs containers. It is installed inside Ubuntu. |
| **Docker Compose** | A helper that reads IMBA's `compose.yaml` and starts the right container. It is installed together with Docker Engine. |

**Docker Desktop is not required.** This guide installs Docker Engine inside Ubuntu.

The Ubuntu you install (24.04) is only the place where you type commands. The tracker container itself is Ubuntu 18.04, because that is the system the IMBA tracker is built for. You do not install Ubuntu 18.04 as your Windows Linux.

## Which window to type in

Two different windows are used. The prompt (the text before the cursor) tells you which one you are in.

| Window | How you know | Prompt looks similar to |
|---|---|---|
| **Windows PowerShell** | Blue window, or a tab titled PowerShell | `PS C:\Users\YourName>` |
| **Ubuntu** | A window or tab titled Ubuntu | `yourname@computer:~$` |

Copy one command at a time. Paste it, then press Enter. When a command starts with `sudo`, Ubuntu asks for your Linux password. The cursor does not move while you type the password. Type it and press Enter.

## About the steps in this guide

The tracking commands in [Build and run IMBA](#6-build-and-run-imba) match `Dockerfile` and `compose.yaml` in this repository. They have been run with Docker Engine inside WSL.

The Windows and Ubuntu install steps are the published Microsoft and Docker instructions. They have not been repeated, as a full walkthrough, on a brand-new Windows PC while writing this page. Where a step can differ between computers, the page says so.

## 1. What this computer needs

These are requirements of Windows and of WSL2, from Microsoft's installation notes, plus what this repository itself states. This repository does not publish a minimum amount of RAM or disk space.

- **Windows 10 version 2004 or later (build 19041 or later), or Windows 11.** WSL2 does not install on older Windows 10 releases. Microsoft's current instructions are at [Install WSL](https://learn.microsoft.com/en-us/windows/wsl/install).
- **Virtualization turned on.** WSL2 uses a virtual machine. The setting is usually in the PC's BIOS or UEFI menu, under a name such as Intel VT-x, AMD-V, or SVM. Windows also needs the optional features "Virtual Machine Platform" and "Windows Subsystem for Linux". The `wsl --install` command turns those Windows features on when it can.
- **An internet connection** for the whole installation. Ubuntu, Docker, the tracker image, and (only if you later build the analysis app) a large set of R packages are downloaded.
- **A 64-bit PC.** Docker's Ubuntu packages support 64-bit Intel/AMD processors (`x86_64`).
- **No GPU and no MATLAB.** The tracker container does not use them. A screen is not required to track a video. The optional graphical window in the main README needs a Linux desktop session (WSLg), which current Windows 11 and updated Windows 10 provide with WSL2.

The tracker build compiles C++ code and downloads an Ubuntu 18.04 image. Leave several gigabytes free. An exact minimum is not stated in the repository.

## 2. Install Ubuntu with WSL2

Do this in **Windows PowerShell**, opened as Administrator.

1. Click Start, type `PowerShell`, right-click **Windows PowerShell**, and choose **Run as administrator**. If Windows asks "Do you want to allow this app to make changes?", choose Yes.
2. Install Ubuntu 24.04:

   ```powershell
   wsl --install -d Ubuntu-24.04
   ```

3. If Windows says a restart is required, restart the computer. After you sign in, open PowerShell again if Ubuntu does not continue on its own.
4. The first time Ubuntu starts, it asks you to create a **Linux username** and a **Linux password**. They can be different from your Windows account. Remember this password. You will type it whenever a command uses `sudo`.
5. Check that Ubuntu is version 2. In **PowerShell** (Administrator is no longer required):

   ```powershell
   wsl --list --verbose
   ```

   A successful result looks like this. The name and the star can vary. **VERSION must be 2.**

   ```text
     NAME            STATE           VERSION
   * Ubuntu-24.04    Running         2
   ```

Open Ubuntu later from the Start menu (**Ubuntu 24.04**). The rest of this guide, unless it says PowerShell, is typed in that Ubuntu window.

## 3. Install Docker Engine

These steps follow Docker's current instructions for Ubuntu: [Install Docker Engine on Ubuntu](https://docs.docker.com/engine/install/ubuntu/).

Use **Ubuntu's own packages only if you have not installed Docker yet.** Ubuntu also publishes packages named `docker.io` and `docker-compose-v2`. Those are not the packages from Docker's website. Mixing them with `docker-ce` causes the "conflicting packages" errors in the troubleshooting section.

### 3.1 Remove a previous unofficial Docker, if you installed one

On a new Ubuntu, skip this step.

If you already tried another Docker install, look first. In **Ubuntu**:

```bash
dpkg -l | grep -E 'docker|containerd|runc' || true
```

If that prints package names such as `docker.io` or `docker-compose-v2`, remove the conflicting packages. This does not delete folders in your Windows user directory. Docker's command is:

```bash
sudo apt remove $(dpkg --get-selections docker.io docker-compose docker-compose-v2 docker-doc docker-buildx podman-docker containerd runc | cut -f1)
```

`apt` may report that some of those names are not installed. That is fine.

### 3.2 Add Docker's repository

In **Ubuntu**, paste this whole block, then press Enter. It installs the download tools, saves Docker's signing key, and adds Docker's package list.

```bash
sudo apt update
sudo apt install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc
sudo tee /etc/apt/sources.list.d/docker.sources <<EOF
Types: deb
URIs: https://download.docker.com/linux/ubuntu
Suites: $(. /etc/os-release && echo "${UBUNTU_CODENAME:-$VERSION_CODENAME}")
Components: stable
Architectures: $(dpkg --print-architecture)
Signed-By: /etc/apt/keyrings/docker.asc
EOF
sudo apt update
```

### 3.3 Install Docker Engine, the CLI, containerd, Buildx, and Compose

In **Ubuntu**:

```bash
sudo apt install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

Confirm the two programs answer:

```bash
docker --version
docker compose version
```

Each command prints a version line. The exact numbers change as Docker publishes updates. The Compose command is `docker compose` with a space. The older `docker-compose` command is a different package and is not the one this guide installs.

### 3.4 Start Docker

In **Ubuntu**:

```bash
sudo systemctl enable --now docker
sudo systemctl status docker --no-pager
```

A running service includes the line `Active: active (running)`. Press `q` if the status view stays open.

Some WSL installations do not use systemd. If `systemctl` says the system has not been booted with systemd, start Docker this way instead:

```bash
sudo service docker start
sudo service docker status
```

### 3.5 Let your Linux user run Docker

Docker's post-install steps add your user to the `docker` group so you do not need `sudo` for every command. In **Ubuntu**, replace nothing in this command; `$USER` is your current name:

```bash
sudo usermod -aG docker $USER
```

Close the Ubuntu window completely and open **Ubuntu 24.04** again from the Start menu. Group membership is read when the window starts.

Check:

```bash
groups
```

The word `docker` should appear in the list.

### 3.6 Run Docker's test container

In **Ubuntu**:

```bash
docker run --rm hello-world
```

The first run downloads a tiny test image. A working install prints:

```text
Hello from Docker!
This message shows that your installation appears to be working correctly.
```

`--rm` deletes that test container when it exits. It does not delete IMBA.

## 4. Download IMBA

The repository is [https://github.com/mthane/IMBA](https://github.com/mthane/IMBA). The branch used for this guide is `main`. There are no git release tags in the repository.

Videos are not part of the download. After a clone, `data/input` and `data/output` exist as empty folders.

### Option A — Git

In **Ubuntu**, this saves IMBA in your Linux home directory:

```bash
cd ~
sudo apt install git
git clone https://github.com/mthane/IMBA.git
cd IMBA
```

`pwd` prints the folder you are in. It looks like `/home/yourname/IMBA`.

In Windows Explorer that folder is:

```text
\\wsl$\Ubuntu-24.04\home\yourname\IMBA
```

Replace `yourname` with the Linux username you created. If Explorer cannot open `\\wsl$`, the Ubuntu window has to be running.

### Option B — ZIP from GitHub

1. In a web browser, open [https://github.com/mthane/IMBA](https://github.com/mthane/IMBA).
2. Choose **Code**, then **Download ZIP**. The ZIP is a snapshot of the default branch.
3. Extract it under your Windows user folder, for example `C:\Users\YourName\Documents\IMBA`. Extracting onto the Desktop is fine as well. Avoid `C:\Program Files`, which is protected.
4. In **Ubuntu**, go to that folder. Windows drive C is `/mnt/c` inside Ubuntu. For the Documents example, with your real Windows user name:

   ```bash
   cd /mnt/c/Users/YourName/Documents/IMBA
   pwd
   ls
   ```

   `ls` should show `Dockerfile`, `compose.yaml`, `README.md`, and `data`.

Stay in this IMBA folder for every later command. If you open a new Ubuntu window, `cd` back to it first.

## 5. Build and run the tracker

Still in the IMBA folder, in **Ubuntu**.

### 5.1 Put a video where the tracker can read it

The default command reads `/data/input/example_video.mp4` inside the container. That path is the folder `data/input` in your IMBA copy. The file is not in the git or ZIP download.

1. Copy your recording into `data/input`.
2. Name it `example_video.mp4`, **or** keep its own filename and remember that name for the next section.

In Explorer, if you cloned with Git into the Linux home directory, open:

```text
\\wsl$\Ubuntu-24.04\home\yourname\IMBA\data\input
```

If you extracted the ZIP under Documents, open `C:\Users\YourName\Documents\IMBA\data\input` and copy the video in as usual.

Check from **Ubuntu**:

```bash
ls -lh data/input
```

You should see your `.mp4` file and a size, not an empty listing.

### 5.2 Build the tracker image

```bash
docker compose build tracker
```

The first build downloads Ubuntu 18.04 and compiles the tracker. It often takes several minutes. The image is named `imba-lrvtrack:1.0`. When it finishes, the command returns to the `$` prompt without an error line. This builds only the tracker. It does not build the optional analysis app.

### 5.3 Track the video

If the file is named `example_video.mp4`:

```bash
docker compose run --rm tracker
```

If the file has another name, still inside `data/input`:

```bash
export IMBA_VIDEO=my_recording.mp4
docker compose run --rm tracker
```

`export` lasts until you close the Ubuntu window.

That command is the one stored in `compose.yaml`:

```text
lrvTrack -x -z 138 --min-obj-size 110 --max-obj-size 5000 -p --thread-count 9 -o -v 16 -d 13 -w 0.04 -u -t -r . -i /data/input/example_video.mp4
```

`-z 138` is the 15 cm dish (138 mm). `-z 84` is the 9 cm dish. The diameter sets the millimetre scale. Objects larger than 6 mm² are treated as clusters: they are drawn on the overlay and left out of the CSV files. If you filmed a 9 cm dish, run:

```bash
docker compose run --rm tracker -x -z 84 --min-obj-size 110 --max-obj-size 5000 -p --thread-count 9 -o -v 16 -d 13 -w 0.04 -u -t -r . -i /data/input/example_video.mp4
```

Change the filename at the end if you used `IMBA_VIDEO`. A finished run prints:

```text
FINISHED!!!
```

### 5.4 Where the results are

In the IMBA folder:

```bash
ls -lh data/output
```

You should see files named with the date and time of the run:

- `YYYYMMDD_HHMMSS_PROC.mp4` — the video with the dish circle and larva numbers.
- `YYYYMMDD_HHMMSS-data/` — a folder of CSV files, one per larva that was tracked as an individual (`1.csv`, `2.csv`, …).

Open the same folder in Explorer:

- Git clone: `\\wsl$\Ubuntu-24.04\home\yourname\IMBA\data\output`
- ZIP under Documents: `C:\Users\YourName\Documents\IMBA\data\output`

### 5.5 Stop the container

`docker compose run --rm tracker` exits when tracking finishes, and `--rm` removes that container. You do not have a tracker server to shut down afterwards.

To stop a run that is still going, click the Ubuntu window and press **Ctrl+C**. A run stopped this way can leave an incomplete result.

The optional analysis app, if you start it later with `docker compose up visualizer`, keeps running. Stop it with Ctrl+C in that window, or in another Ubuntu window from the IMBA folder:

```bash
docker compose stop visualizer
```

## 6. Videos, folders, and Windows Explorer

Docker does not see your whole PC. A **bind mount** connects one folder on the computer to one folder inside the container. IMBA's `compose.yaml` sets these up:

| Folder in your IMBA copy | Folder inside the container | Who may write |
|---|---|---|
| `data/input` | `/data/input` | Read-only. The tracker cannot change or delete the raw videos. |
| `data/output` | `/data/output` | The tracker writes results here. |
| `config` | `/app/config` | Read-only. Sample trial notes live here. |

Putting a video in `data/input` on Windows is the same file the container reads as `/data/input/...`. Copy the video into that folder. Do not paste a Windows path such as `C:\Users\...` into the tracker's `-i` option. Inside the container that path does not exist.

The image build ignores videos on purpose (`.dockerignore` excludes `data` and `*.mp4`). A video copied into the image would make every rebuild huge, and you would have to rebuild to change a recording. Keep videos in `data/input`, or in an outside folder mounted the same way.

### Large collections outside the IMBA folder

Create two folders in Explorer, for example `D:\larva\raw` and `D:\larva\tracked`. Put videos only in `raw`. In **Ubuntu**, drive D is `/mnt/d`:

```bash
export IMBA_INPUT_DIR=/mnt/d/larva/raw
export IMBA_OUTPUT_DIR=/mnt/d/larva/tracked
export IMBA_VIDEO=trial01.mp4
docker compose run --rm tracker
```

Both folders must exist before you run the command. `/mnt/d/larva/raw` is still mounted read-only, because `compose.yaml` adds `:ro` to the input mount. Use a different folder for output so a later command cannot be pointed at the raw videos by mistake.

These `export` lines are for Docker Engine inside Ubuntu. Setting the same names in a Windows PowerShell window does not affect this Ubuntu Docker.

## 7. When something fails

Read the message in the terminal before changing anything. For a Docker service problem, look at the service log before reinstalling.

In **Ubuntu**:

```bash
sudo systemctl status docker --no-pager
journalctl -u docker -n 50 --no-pager
```

If `journalctl` is not available because systemd is off:

```bash
sudo service docker status
```

The log text is what tells you the next step. The cases below are the usual ones.

### `docker: command not found`

The Docker packages are not installed, or this window was opened before you installed them. Repeat [section 3.3](#33-install-docker-engine-the-cli-containerd-buildx-and-compose) in Ubuntu, then open a new Ubuntu window.

### `permission denied while trying to connect to the Docker API`

Your user is not in the `docker` group yet, or the Ubuntu window was already open when you ran `usermod`. Repeat [section 3.5](#35-let-your-linux-user-run-docker), close Ubuntu, open it again, and check `groups`. Until that works, `sudo docker run --rm hello-world` can test the engine. Day-to-day IMBA commands are meant to run without `sudo`.

### `Cannot connect to the Docker daemon`

The engine is not running. Start it with [section 3.4](#34-start-docker), then look at `systemctl status docker` or `service docker status`. Do not delete `/var/lib/docker`. That folder holds downloaded images.

### `docker.service is masked`

Masking is a switch that forbids the service to start. Check the status first. If the status line says `masked`, turn that switch off and start the service:

```bash
sudo systemctl unmask docker
sudo systemctl start docker
sudo systemctl status docker --no-pager
```

### Conflicting Docker Compose packages

If `docker compose version` fails and `docker-compose` exists, Ubuntu's older packages may be installed beside Docker's packages. Check with the `dpkg` command in [section 3.1](#31-remove-a-previous-unofficial-docker-if-you-installed-one), remove those unofficial names, and install `docker-compose-plugin` again. Use `docker compose` with a space afterwards.

### The build fails

Scroll to the first line that says `error:` or `E:`. Common causes:

- No internet, or a firewall blocked `archive.ubuntu.com`, `download.docker.com`, or `docker.io`.
- The command was not run from the IMBA folder (`ls` should show `Dockerfile` and `compose.yaml`).
- Disk space ran out. Free space on the Windows drive that holds Ubuntu, then run `docker compose build tracker` again. A failed build does not delete `data/input` or `data/output`.

### The video is missing

```text
Error setting up the capture device
```

`data/input` does not contain the filename you passed, or a Windows path was given to `-i`. From the IMBA folder, `ls data/input` must show the file. The path inside the container starts with `/data/input/`.

### You cannot open the result files

If the results were written as root, Ubuntu can read them but saving from another program may fail. For the next run, still in the IMBA folder:

```bash
export IMBA_UID="$(id -u)"
export IMBA_GID="$(id -g)"
docker compose run --rm tracker
```

To see who owns a file already written: `ls -l data/output`. Changing ownership of those files is `sudo chown -R "$USER:$USER" data/output`. That changes the owner. It does not delete the files.

### WSL or virtualization will not install

In **PowerShell**:

```powershell
wsl --status
wsl --list --verbose
```

- If **VERSION** is `1`, convert that distro: `wsl --set-version Ubuntu-24.04 2`.
- If Windows reports that virtualization is not available, enable it in the PC firmware, then in Windows optional features enable **Virtual Machine Platform** and **Windows Subsystem for Linux**, and restart. Microsoft describes this on the [WSL install page](https://learn.microsoft.com/en-us/windows/wsl/install).
- `wsl --install` needs an Administrator PowerShell and an internet connection.

## 8. Update IMBA without losing videos

Videos in `data/input` and results in `data/output` are not part of git. An update does not need those folders to be deleted.

If you cloned with Git, in **Ubuntu**, from the IMBA folder:

```bash
git pull
docker compose build tracker
```

If `git pull` stops and talks about local changes or conflicts, do not delete `data/` to clear the message. Leave the data folders as they are and ask someone who uses git before running any command that resets files.

If you used the ZIP, download a new ZIP, extract it to a **new** folder, and copy your `data/input` and `data/output` folders into the new copy. Do not extract the new ZIP on top of the old folder in a way that replaces `data`.

Rebuild after an update so the tracker image matches the new source. `docker compose build tracker` replaces the image, not your videos.

## 9. Record the version, and cite the work

There is no numbered release tag in this repository. The version of a run is the git commit, plus the image you built from it.

From the IMBA folder, in **Ubuntu**, record the output of:

```bash
git rev-parse HEAD
git status -sb
docker image inspect imba-lrvtrack:1.0 --format '{{.Id}}'
```

`git rev-parse HEAD` prints a long id such as `272d453...`. That id is the IMBA source you ran. `main` moves when people push new commits, so a lab notebook should store the full id, not only the word `main`.

This repository has no `LICENSE` file. Do not assume a licence. The tracker build unpacks `IMBAtracker/IMBtracker_cpp/cvblob.zip`, which has its own notices.

Cite the work as:

Toshima et al. (2026), *eLife*.  
DOI: [https://doi.org/10.7554/eLife.106148.3](https://doi.org/10.7554/eLife.106148.3)

The README also asks users to cite the software paper: [10.1098/rsob.220308](https://doi.org/10.1098/rsob.220308).

## Optional: the analysis app

Tracking does not require this app. IMBAvisualizer reads CSV files after tracking and opens in a web browser. Building it downloads a large R image and many R packages, which takes much longer than the tracker build.

```bash
docker compose build visualizer
docker compose up visualizer
```

In a browser on Windows, open [http://localhost:3838](http://localhost:3838). How to choose the experiment folder, dish radius, and frame rate is in the main [README](../README.md). Stop the app with Ctrl+C, or `docker compose stop visualizer`.

The PyQt tracking window is also described in the README. It is separate from `docker compose run tracker`, and it needs a graphical Ubuntu session.
