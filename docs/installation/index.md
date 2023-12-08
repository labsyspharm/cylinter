---
layout: default-cylinter
title: Installation
nav_order: 2
has_children: false
---

# Installation

CyLinter is written in Python 3 and is compatible with MacOS, Windows, and Linux operating systems. The program can be installed via the cross-platform package manager, Conda.

## 1. Install Miniconda
The folllowing are examples of commands for quickly and quietly installing the latest version of the Miniconda installer for your operating system (MacOS, Windows, Linux). To install a different version or architecture of Miniconda for any platform, change the name of the .sh (or .exe) installer in the curl (or wget) command to the one appropriate for your machine. The latest Miniconda installers can be found [here](https://docs.conda.io/projects/miniconda/en/latest/index.html).

``` bash
# MacOS:
mkdir -p ~/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh

# Windows:
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -o miniconda.exe
start /wait "" miniconda.exe /S
del miniconda.exe

# Linux:
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

## 2. Set Libmamba as the default dependency solver
[Libmamba](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community) is a fast dependency solver that can greatly expidite the installation of Conda packages.

``` bash
conda update -n base conda
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

## 3. Install CyLinter
Install CyLinter into a dedicated virtual environment with the following command:  

``` bash
conda create -n cylinter -c conda-forge -c gjbaker -c labsyspharm python=3 cylinter=0.0.47
```
