# Installation

CyLinter is developed on macOS and written in Python 3.

## Install Python 3

Check if Python 3 is already installed.

``` bash
python3 --version
```

 Python 3 is installed if a software version number is returned by the above command. If not, install Python3 version 3.8.5 using the official installer at [Python.org](https://www.python.org/downloads/mac-osx/). Note, as of this writing CyLinter is incompatible with the latest version of Python 3 (3.9.1) ).  

## Create a Python virtual environment

Create an isolated Python virtual environment [Python virtual environment](https://docs.python.org/3/library/venv.html) for CyLinter projects to avoid incompatibles with system-wide Python packages.

``` bash
python3 -m venv $HOME/cylinter  # creates a virtual environment called "cylinter" in the user's home directory
source $HOME/cylinter/bin/activate  # activates the newly-created virtual environment
```

## Install CyLinter

``` bash
pip install cylinter  
```
