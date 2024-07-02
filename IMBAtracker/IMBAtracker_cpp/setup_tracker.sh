sudo apt update
sudo apt install git
sudo git clone https://github.com/mthane/IMBA
cd IMBA/IMBAtracker/IMBAtracker_cpp

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


sudo unzip cvblob.zip
cd cvblob
sudo cmake .
sudo make

cd ..
sudo cmake .
sudo make