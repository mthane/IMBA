# IMBA - Individual Maggot Behavior Analysis

## Installation (Ubuntu-18.04)
### Installation of Ubuntu on Windows

1. Open command prompt and install WSL2 with Ubuntu 18.04 using:
```
wsl --install -d Ubuntu-20.04
```

This will open another command prompt with the Ubuntu-18.04 system. Here do the following:
1. Enter desired user name and password...
3. Install git and clone this repository:
```
sudo apt update
sudo apt install git
sudo git clone https://github.com/mthane/IMBA
```
### Installation IMBAtracker

```
cd IMBA/IMBAtracker
bash setup.sh
```
### Run IMBAtracker via AppImage

```
sudo apt-get install fuse libfuse2
sudo apt-get install libfreetype6
sudo apt-get install libsm6
sudo apt-get install libgl1
sudo apt-get install libharfbuzz0b
sudo apt-get install libfontconfig1
sudo apt-get install libthai0
```

### Installation IMBAvisualizer

Install Rtools 4.2 form https://cran.r-project.org/bin/windows/Rtools/
```
cd IMBA/IMBAvisualizer
bash setup.sh
```
