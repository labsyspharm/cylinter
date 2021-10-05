# Installation

CyLinter is developed on macOS and written in Python 3.

## Install Python 3

Check if Python 3 is installed.

``` bash
python3 --version
```

 Python 3 is installed if a version number is returned by the above command. If Python 3 is not installed, or a version number less than 3.6.0 is installed, it is recommended to install a version > 3.6.0 using the official installer at [Python.org](https://www.python.org/downloads/mac-osx/).

## Create a Python virtual environment

To avoid conflicts between CyLinter and system-wide Python dependences, create a [Python virtual environment](https://docs.python.org/3/library/venv.html) dedicated to CyLinter analysis.

``` bash
# Create a virtual environment called "cylinter" in the home directory.
python3 -m venv ~/cylinter

# Activate the newly-created virtual environment.
source ~/cylinter/bin/activate  
```

## Install CyLinter

``` bash
pip install cylinter  
```
