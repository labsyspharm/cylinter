---
layout: default-cylinter
title: Installation
nav_order: 2
has_children: false
---

# Installation

CyLinter is written in Python 3 and is compatible with MacOS, Windows, and Linux operating systems. The program can be installed via the cross-platform package manager, Conda.

## 1. Install Miniconda

**NOTE**: If you already have Miniconda or Anaconda installed, [skip this section and jump to section 1B](#section-1b).

The folllowing are examples of commands for quickly and quietly installing the latest version of the Miniconda installer for your operating system (MacOS - M1 / Intel 64-bit, Windows, Linux - Intel 64-bit). For other platforms, [consult the Miniconda download page](https://docs.conda.io/projects/miniconda/en/latest/index.html).

### MacOS
Open Terminal and paste the following commands:
```bash
mkdir -p ~/miniconda3

# M1 chip
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh

# Intel 64-bit chip
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/miniconda3/miniconda.sh

bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -f ~/miniconda3/miniconda.sh
```

### Windows
Open a Command Prompt and paste the following commands:
```cmd
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -o miniconda.exe
start /wait "" miniconda.exe /S
del miniconda.exe
```

### Linux
Open a terminal window and paste the following commands:
```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -f ~/miniconda3/miniconda.sh
```

You should now [skip section 1B and go straight to section 2](#section-2).

## 1B. For existing conda installations, set libmamba as the default dependency solver
{: #section-1b}
CyLinter depends on a complex set of packages and older Conda installations will struggle with this. If you already have Miniconda or Anaconda installed, we suggest that you run the following commands to update Conda itself and enable the libmamba dependency solver. This will help ensure CyLinter can be installed efficiently.

``` bash
conda update -n base conda
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

## 2. Install CyLinter
{: #section-2}
Install CyLinter into a dedicated conda environment with the following command:  

``` bash
conda create -n cylinter -c conda-forge -c labsyspharm cylinter=0.0.53
```
