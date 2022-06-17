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

 Python 3 is installed if a version number is returned by the above command. If Python 3 is not installed, or a version number less than 3.6.0 or greater than 3.9.X is installed, install a version within this range using the official installer at [Python.org](https://www.python.org).

## Create a Python virtual environment

To avoid conflicts between CyLinter and system-wide Python dependences, create a [Python virtual environment](https://docs.python.org/3/library/venv.html) dedicated to CyLinter analysis (conda environments are also acceptable). Create a virtual environment called "cylinter" in the home directory:

``` bash
# for Mac
python3 -m venv ~/cylinter

# for PC:
python3 -m venv C:\Users\<username>\cylinter
```

Once the virtual environment has been created, it can be “activated” using a script in the virtual environment’s binary directory. The invocation of the script is platform-specific:
``` bash
| Platform | Shell           | Command to activate virtual environment                 |
|----------|-----------------|---------------------------------------------------------|
| POSIX    | bash/zsh        | $ source ~/cylinter/bin/activate                        |
|          | fish            | $ source ~/cylinter/bin/activate.fish                   |
|          | csh/tcsh        | $ source ~/cylinter/bin/activate.csh                    |
|          | PowerShell Core | $ ~/cylinter/bin/Activate.ps1                           |
| Windows  | cmd.exe         | C:\> \Users\<username>\cylinter\Scripts\activate.bat    |
|          | PowerShell      | PS C:\> \Users\<username>\cylinter\Scripts\Activate.ps1 |
```

## Install CyLinter
After activating the virtual environment, use the `pip` command line tool to install CyLinter into the virtual environment via the [Python Package Index (PyPI)](https://pypi.org/):

``` bash
pip install cylinter  
```
