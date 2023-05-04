---
layout: default-cylinter
title: Installation
nav_order: 2
has_children: false
---

# Installation

CyLinter is written in Python 3 and is compatible with MacOS, Linux, and Windows operating systems. The program can be installed via the cross-platform package manager, Mamba.

## Step 1. Install Python 3

Check if Python 3 is installed.

``` bash
python3 --version
```

 Python 3 is installed if a version number is returned by the above command. If Python 3 is not installed, or a version less than 3.6.0 is installed, please install a version >= 3.6.0 using the official installer at [Python.org](https://www.python.org).

## 2. Install Mamba
Refer to the [Mamba Documentation](https://mamba.readthedocs.io/en/latest/installation.html) for platform-specific installation instructions.

## 3. Install CyLinter
To avoid conflicts between CyLinter and system-wide Python dependences, CyLinter is installed into a dedicated Mamba environment with the following command:  

``` bash
mamba create -n cylinter -c conda-forge -c labsyspharm cylinter
```
