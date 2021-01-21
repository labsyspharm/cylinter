# Installation

## Install Python 3

If not already installed, install Python 3: https://www.python.org/downloads/mac-osx/

## Install virtualenv and virtualenvwrapper (recommended)

To isolate CyLinter and its dependencies from global Python settings, we recommend installing CyLinter into a dedicated Python virtual environment.

``` bash
python3 -m pip install --user virtualenv   
python3 -m virtualenv --help
python3 -m pip install virtualenvwrapper
echo 'export WORKON_HOME=$HOME/.virtualenvs' >> ~/.bashrc
echo 'export PROJECT_HOME=$HOME/Devel' >> ~/.bashrc
source /usr/local/bin/virtualenvwrapper.sh
source ~/.bashrc
```

## Create and step into CyLinter virtual environment
``` bash
mkvirtualenv -p python3 cylinter
```

## Install CyLinter and its dependencies
``` bash
pip install cylinter
```
