# Installation

CyLinter is written in Python 3 and developed on macOS.

## Install Python 3

Check if Python 3 is already installed by running

``` bash
python3 --version
```

 If a software version number (e.g. Python 3.9.1) is returned by the above command, then Python 3 is already installed. Otherwise, install the latest version of Python3 using the official installer at [Python.org](https://www.python.org/downloads/mac-osx/).  

## Create a Python virtual environment

To avoid discrepancies between versions of CyLinter Python packages and those required by system-wide Python packages, create a dedicated [Python virtual environment](https://docs.python.org/3/library/venv.html) in which to install CyLinter and run projects.

``` bash
python3 -m venv $HOME/cylinter
```

## Install CyLinter

``` bash
source $HOME/cylinter/bin/activate  # activates virtual environment
pip install cylinter  # install CyLinter
```
