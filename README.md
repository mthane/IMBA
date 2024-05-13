# IMBA: Individual Maggot Behavior Analyzer

## Introduction

Neuronally orchestrated muscular movement and locomotion are fundamental aspects of multicellular animals. The larva of the fruit fly Drosophila melanogaster provides a unique opportunity to study these processes due to its simple brain and genetic accessibility. However, traditional methods of studying larval locomotion often aggregate measurements across animals or test animals individually, limiting our ability to understand inter- and intra-individual variability in locomotion and its neurogenetic determinants.

To address this gap, we introduce the Individual Maggot Behavior Analyzer (IMBA), a software tool designed for analyzing the behavior of individual larvae within groups. IMBA enables researchers to reliably resolve individual larval identities across collisions, providing unprecedented insights into locomotion variability and its underlying genetic and neural mechanisms.

In this README, we provide an overview of IMBA's features, installation instructions, usage guidelines, and examples of its application in studying larval behavior in various biomedical research contexts. With IMBA, researchers can obtain a rich understanding of individual larval behavior, paving the way for deeper insights into neurogenetic pathways governing locomotion and its modulation.

## Table of Contents
1. [Installation](#installation)
2. [Usage](#usage)
3. [Examples](#examples)
4. [Contributing](#contributing)
5. [License](#license)

## Installation (Ubuntu-18.04)
### Installation of Ubuntu on Windows

1. Open command prompt and install WSL2 with Ubuntu 18.04 using:
```
wsl --install -d Ubuntu-18.04
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

### Installation IMBAvisualizer

Install Rtools 4.2 form https://cran.r-project.org/bin/windows/Rtools/
```
cd IMBA/IMBAvisualizer
bash setup.sh
```
