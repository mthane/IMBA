sudo apt update
sudo apt install python3
sudo apt --only-upgrade install python3
sudo apt install python3-pip
sudo apt-get install libxkbcommon-x11-0
sudo apt-get install libx11-xcb1
# install all dependencies
sudo apt install liblpsolve55-dev libalglib-dev libboost-all-dev libopencv-dev build-essential cmake git unzip
bash IMBtracker_cpp/setup.sh
#sudo apt install virtualenv
#virtualenv --python=python3 env
#source env/bin/activate

sudo apt install python-setuptools-scm
pip3 install --upgrade setuptools pip
pip3 install scikit-build
pip3 install opencv-python

pip3 install --user pyqt5  
sudo apt-get install python3-pyqt5  
sudo apt-get install pyqt5-dev-tools
sudo apt-get install qttools5-dev-tools

python3 TrackExperiment.py