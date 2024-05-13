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

## Installation 

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


```
cd IMBA/IMBAtracker
bash setup.sh
```

## Usage 


## Examples


# IMBAvisualizer <a name="imbavisualizer"></a>

Install Rtools 4.2 form https://cran.r-project.org/bin/windows/Rtools/
```
cd IMBA/IMBAvisualizer
bash setup.sh
```


