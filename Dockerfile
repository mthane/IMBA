# IMBA 1.0 tracker (lrvTrack).
# Ubuntu 18.04 is the environment documented in this repository and in the
# historical Dockerfile. Its OpenCV 3.2 packages still provide the CV_* APIs
# this source was written against. A newer OpenCV would require source changes.
FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC \
    LD_LIBRARY_PATH=/usr/local/lib

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        pkg-config \
        unzip \
        libalglib-dev \
        libboost-all-dev \
        liblpsolve55-dev \
        libopencv-dev \
        libtbb-dev \
    && rm -rf /var/lib/apt/lists/*

# libtbb-dev is required by CMakeLists.txt (find_package(TBB REQUIRED)) but is
# not listed in IMBAtracker/IMBtracker_cpp/setup.sh.

WORKDIR /opt/imba/IMBtracker_cpp
COPY IMBAtracker/IMBtracker_cpp/ ./

# Same build sequence as IMBAtracker/IMBtracker_cpp/setup.sh.
RUN unzip -o -q cvblob.zip \
    && cd cvblob \
    && cmake . \
    && make -j"$(nproc)" \
    && cd .. \
    && cmake . \
    && make -j"$(nproc)" \
    && install -m 755 lrvTrack /usr/local/bin/lrvTrack \
    && find cvblob -name 'libcvblob.so*' -exec cp -a {} /usr/local/lib/ ';' \
    && ldconfig

# Outputs (DATE-data/*.csv, overlay video, err.png) are written relative to
# the process working directory, not under -r. The published command uses
# -r . from the directory that should receive results.
RUN printf '%s\n' \
        '#!/bin/sh' \
        'set -eu' \
        'if [ ! -d /data/output ]; then' \
        '  echo "lrvTrack: /data/output is not mounted" >&2' \
        '  exit 1' \
        'fi' \
        'if [ ! -w /data/output ]; then' \
        '  echo "lrvTrack: /data/output is not writable" >&2' \
        '  exit 1' \
        'fi' \
        'cd /data/output' \
        'exec /usr/local/bin/lrvTrack "$@"' \
        > /usr/local/bin/imba-track \
    && chmod 755 /usr/local/bin/imba-track

WORKDIR /data/output
ENTRYPOINT ["/usr/local/bin/imba-track"]
# Published CLI from the repository README. The input path is a container
# path; compose.yaml supplies the video name.
CMD ["-x", "-z", "138", "--min-obj-size", "110", "--max-obj-size", "5000", "-p", "--thread-count", "9", "-o", "-v", "16", "-d", "13", "-w", "0.04", "-u", "-t", "-r", ".", "-i", "/data/input/example_video.mp4"]
