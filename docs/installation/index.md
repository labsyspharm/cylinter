---
layout: default-cylinter
title: Installation
nav_order: 2
has_children: false
---

# Installation

CyLinter is written in Python 3 and is compatible with MacOS, Linux, and Windows operating systems. The program can be installed via the cross-platform package manager, Mamba.

## 1. Install Mamba
Download the latest [Mambaforge installer](https://github.com/conda-forge/miniforge#mambaforge) for your operating system (e.g., Linux, MacOS, Windows).

MacOS and Linux users can install Mamba using the following command:
``` bash
sh <path_to_installer>
```
PC users can install Mamba by double-clicking the installer executable file and following the installation instructions.

## 2. Install CyLinter
Install CyLinter into a dedicated virtual environment with the following command:  

``` bash
mamba create -n cylinter -c conda-forge -c labsyspharm python=3.10 cylinter
```
