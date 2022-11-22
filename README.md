# IMB - Individual Maggot Behavior Analysis

## Installation
### Windows

1. Open command prompt and install WSL2 with Ubuntu 18.04 using:
```
wsl --install -d Ubuntu-18.04
```
This will open another command prompt with the Ubuntu-18.04 system. Here do the following:

2. Install git and clone this repository:
```
sudo apt update
sudo apt install git
git clone https://github.com/mthane/IMB
```
3. Install the Tracking Software

```
cd IMB/IMBtracker
bash setup.sh
```
