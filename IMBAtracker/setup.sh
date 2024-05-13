#!/bin/bash

# Update package lists and install necessary dependencies
apt-get update && \
apt-get install -y --fix-missing build-essential zlib1g-dev \
libncurses5-dev libgdbm-dev libnss3-dev \
libssl-dev libreadline-dev libffi-dev curl software-properties-common pyqt5-dev-tools \
qttools5-dev-tools libxkbcommon-x11-0 \
libx11-xcb1 \
liblpsolve55-dev \
libalglib-dev \
libboost-all-dev 

# Set noninteractive mode and timezone
export DEBIAN_FRONTEND=noninteractive
export TZ=UTC

# Update package lists and install wget and tar
apt-get update && \
apt-get install -yqq --no-install-recommends wget tar

# Download and extract Python 3.9 source code
wget https://www.python.org/ftp/python/3.9.0/Python-3.9.0.tar.xz && \
tar -xf Python-3.9.0.tar.xz && \
cd Python-3.9.0 && \
./configure && \
make && \
make install

# Install pip and setuptools
apt-get install -y python3-pip python3-setuptools python3-pyqt5

# Install additional Python dependencies
pip3 install --upgrade pip scikit-build opencv-contrib-python-headless 
pip install PyQt5

