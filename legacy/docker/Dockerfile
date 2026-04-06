# Build from repository root:
#   docker build -f legacy/docker/Dockerfile .
FROM ubuntu:18.04
WORKDIR /app
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC
RUN apt-get update && \
    apt-get install -y --fix-missing build-essential zlib1g-dev \
    libncurses5-dev libgdbm-dev libnss3-dev \
    libssl-dev libreadline-dev libffi-dev curl software-properties-common \
    pyqt5-dev-tools \
    qttools5-dev-tools

RUN apt-get install -y --fix-missing \
    liblpsolve55-dev \
    libalglib-dev \
    libboost-all-dev \
    git unzip

RUN apt-get install -y --fix-missing libopencv-dev
RUN apt-get install -y cmake

COPY legacy/tracker_cpp /app/legacy/tracker_cpp
WORKDIR /app/legacy/tracker_cpp

RUN unzip -o cvblob.zip

WORKDIR /app/legacy/tracker_cpp/cvblob
RUN cmake .
RUN make

WORKDIR /app/legacy/tracker_cpp
RUN cmake .
RUN make

CMD ["./lrvTrack"]
