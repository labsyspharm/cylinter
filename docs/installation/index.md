---
layout: default-cylinter
title: Installation
nav_order: 2
has_children: false
---

# Installation

CyLinter source code is written in Python 3.

## Install Python 3

Check if Python 3 is installed.

``` bash
python3 --version
```

 Python 3 is installed if a version number is returned by the above command. If Python 3 is not installed, or a version number less than 3.6.0 is installed, it is recommended to install a version >= 3.6.0 using the official installer at [Python.org](https://www.python.org).

## Create a Python virtual environment

To avoid conflicts between CyLinter and system-wide Python dependences, create a [Python virtual environment](https://docs.python.org/3/library/venv.html) dedicated to CyLinter analysis (conda environments are also acceptable). Create a virtual environment called "cylinter" in the home directory and activate it with the following commands:

``` bash
# for Mac
python3 -m venv ~/cylinter
source ~/cylinter/bin/activate  

# for PC:
python3 -m venv C:\Users\<username>\cylinter
cd C:\Users\<username>\cylinter\scripts
activate
```

## Install CyLinter
Use the `pip` command line tool to install CyLinter via the [Python Package Index (PyPI)](https://pypi.org/):

``` bash
pip install cylinter  
```
