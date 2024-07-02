#!/bin/sh

# install all dependencies
sudo apt install liblpsolve55-dev libalglib-dev libboost-all-dev libopencv-dev build-essential cmake git unzip

# git clone ...

# build cvblob
unzip cvblob.zip
cd cvblob
cmake .
make

# build lrvTrack
cd ..
cmake .
make
