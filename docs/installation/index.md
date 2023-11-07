---
layout: default-cylinter
title: Installation
nav_order: 2
has_children: false
---

# Installation

CyLinter is written in Python 3 and is compatible with MacOS, Windows, and Linux operating systems. The program can be installed via the cross-platform package manager, Conda.

## 1. Install Miniconda
Download the latest [Miniconda installer](https://docs.conda.io/projects/miniconda/en/latest/index.html) for your operating system (e.g., MacOS, Windows, Linux) and set [Libmamba](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community) as the default dependency solver using the following commands:

``` bash
conda update -n base conda
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

## 2. Install CyLinter
Install CyLinter into a dedicated virtual environment with the following command:  

``` bash
conda create -n cylinter -c conda-forge -c gjbaker -c labsyspharm python=3 cylinter
```
