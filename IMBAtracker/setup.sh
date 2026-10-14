sudo apt update
sudo apt install python3
sudo apt --only-upgrade install python3
sudo apt install python3-pip
sudo apt install  -force-yes libxkbcommon-x11-0 libx11-xcb1 liblpsolve55-dev libalglib-dev libboost-all-dev libopencv-dev build-essential cmake git unzip
cd IMBtracker_cpp
setup.sh
cd ..
sudo apt install python-setuptools-scm
pip3 install --upgrade setuptools pip
pip3 install scikit-build
pip3 install opencv-python-headless

sudo apt-get install python3-pyqt5  
sudo apt-get install pyqt5-dev-tools
sudo apt-get install qttools5-dev-tools
sudo apt-get upgrade 

python3 TrackExperiment.py
