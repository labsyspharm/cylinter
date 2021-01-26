# Installation

CyLinter is developed on macOS and written in Python 3.

## Install Python 3

Check if Python 3 is already installed.

``` bash
python3 --version
```

 Python 3 is installed if a version number is returned by the above command. If not installed or a version other than 3.8.5, install 3.8.5 using the official installer at [Python.org](https://www.python.org/downloads/mac-osx/).

## Create a Python virtual environment

Create an isolated Python virtual environment [Python virtual environment](https://docs.python.org/3/library/venv.html) in which to run CyLinter projects. This will avoid potential compatibility issues between CyLinter and system-wide Python dependencies.

``` bash
# Create a virtual environment called "cylinter" in your home directory.
python3 -m venv $HOME/cylinter

# Activate the newly-created virtual environment.
source $HOME/cylinter/bin/activate  
```

## Install CyLinter

``` bash
pip install cylinter  
```
