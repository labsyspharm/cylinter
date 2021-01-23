# Installation

CyLinter is written in Python 3 and developed on macOS.

## Install Python 3

Check if Python 3 is already installed by running
``` bash
python3 --version
```
 If a software version number (e.g. Python 3.8.5) is returned, then Python 3 is already installed. Otherwise, install Python 3 by downloading the official installer from [Python.org](https://www.python.org/downloads/mac-osx/). The macOS package manager [Homebrew](https://brew.sh/) is an acceptable alternative installation method.  

## Install virtualenv and virtualenvwrapper (recommended)

To avoid incompatibilities between Python packages and their versions, install CyLinter into a dedicated Python virtual environment (https://realpython.com/python-virtual-environments-a-primer/). We recommend using [virtualenv](https://virtualenv.pypa.io/en/latest/), a tool to create isolated Python environments, together with [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/), a set of virtualenv convenience extensions.

``` bash
# install virtualenv
python3 -m pip install --user virtualenv   
python3 -m virtualenv --help

# install virtualenvwrapper
python3 -m pip install virtualenvwrapper
echo 'export WORKON_HOME=$HOME/.virtualenvs' >> ~/.bashrc  # set virtual environments location
echo 'export PROJECT_HOME=$HOME/Devel' >> ~/.bashrc  # set development project directories location
source /usr/local/bin/virtualenvwrapper.sh  # set location of the script installed with virtualenvwrapper
source ~/.bashrc  # reload shell configuration
```

## Create a virtual environment for CyLinter projects
``` bash
mkvirtualenv -p python3 cylinter  # creates virtualenv and steps into it
```

## Install CyLinter
``` bash
pip install cylinter
```
