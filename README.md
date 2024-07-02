# IMBA: Individual Maggot Behavior Analyzer


## Introduction <a name="introduction"></a>

Neuronally orchestrated muscular movement and locomotion are fundamental aspects of multicellular animals. The larva of the fruit fly Drosophila melanogaster provides a unique opportunity to study these processes due to its simple brain and genetic accessibility. However, traditional methods of studying larval locomotion often aggregate measurements across animals or test animals individually, limiting our ability to understand inter- and intra-individual variability in locomotion and its neurogenetic determinants.

To address this gap, we introduce the Individual Maggot Behavior Analyzer (IMBA), a software tool designed for analyzing the behavior of individual larvae within groups. IMBA enables researchers to reliably resolve individual larval identities across collisions, providing unprecedented insights into locomotion variability and its underlying genetic and neural mechanisms.

In this README, we provide an overview of IMBA's features, installation instructions, usage guidelines, and examples of its application in studying larval behavior in various biomedical research contexts. With IMBA, researchers can obtain a rich understanding of individual larval behavior, paving the way for deeper insights into neurogenetic pathways governing locomotion and its modulation.

## Table of Contents
1. [Introduction](#introduction)
2. [IMBAtracker](#imbatracker)
3. [IMBAvisualizer](#imbavisualizer)

---

# IMBAtracker <a name="imbatracker"></a>
The IMBAtracker consists of mainly two parts. The first part is the Tracking software that is written in c++. It is developed with the help of the OpenCV library and CVblob. In order to run the tracker this c++ software has to be installed. It can be used as a standalone tool via command line.
The second part is the GUI which is developed in Python 3.9 using the Qt library. The GUI is providing an easy use of the command line tool for tracking large amounts of data. First we will install the Tracker and then the GUI.

## Installation - Tracker

### Cloning Repository and using WSL
If you are using Ubuntu you can skip step 1 and 2.
1. Open command prompt and install WSL2 with Ubuntu 18.04 using:
```
wsl --install -d Ubuntu-18.04
```

This will open another command prompt with the Ubuntu-18.04 system. Here do the following:
2. Enter desired user name and password...
3. Install git and clone this repository:
```
sudo apt update
sudo apt install git
sudo git clone https://github.com/mthane/IMBA
```

4. Install dependencies

```
sudo apt-get update
sudo apt-get install -y --fix-missing build-essential zlib1g-dev \
    libncurses5-dev libgdbm-dev libnss3-dev \
    libssl-dev libreadline-dev libffi-dev curl software-properties-common \
    pyqt5-dev-tools qttools5-dev-tools

sudo apt-get install -y --fix-missing \
    liblpsolve55-dev \
    libalglib-dev \
    libboost-all-dev \
    git unzip

sudo apt-get install -y cmake
sudo apt-get install -y build-essential cmake git libgtk2.0-dev pkg-config libavcodec-dev libavformat-dev libswscale-dev
sudo apt-get install -y python3.6-dev python3-numpy libtbb2 libtbb-dev libjpeg-dev libpng-dev libtiff-dev libdc1394-22-dev
sudo apt-get install -y libavcodec-dev libavformat-dev libswscale-dev libv4l-dev
sudo apt-get install -y libxvidcore-dev libx264-dev
sudo apt-get install -y libgtk-3-dev
sudo apt-get install -y libatlas-base-dev gfortran
sudo apt-get install -y python3-dev

cd ~
git clone https://github.com/opencv/opencv.git
git clone https://github.com/opencv/opencv_contrib.git
cd opencv
mkdir build
cd build

cmake -D CMAKE_BUILD_TYPE=RELEASE \
      -D CMAKE_INSTALL_PREFIX=/usr/local \
      -D INSTALL_C_EXAMPLES=ON \
      -D INSTALL_PYTHON_EXAMPLES=ON \
      -D OPENCV_GENERATE_PKGCONFIG=ON \
      -D OPENCV_EXTRA_MODULES_PATH=~/opencv_contrib/modules \
      -D BUILD_EXAMPLES=ON ..

make -j8
sudo make install
sudo ldconfig

```

5. Change to tracker directory, unzip and build the cvblob package.
```

cd /IMBA/IMBAtracker/IMBAtracker_cpp
sudo unzip cvblob.zip
cd cvblob
sudo cmake .
sudo make

```

6. Change to tracker directory and build the tracker c++ code.

```
cd ..
sudo cmake .
sudo make

```

#### Run IMBAtracker
```
./lrvTrack -x -z 138 --min-obj-size 110 --max-obj-size 5000 -p --thread-count 9 -o -v 16 -d 13 -w 0.04 -u -t -r "." -i example_video.mp4

```

### Usage 


### Examples

## Installation - GUI

### Usage 


### Examples

# IMBAvisualizer <a name="imbavisualizer"></a>

Install Rtools 4.2 form https://cran.r-project.org/bin/windows/Rtools/
```
cd IMBA/IMBAvisualizer
sudo bash setup.sh
```


