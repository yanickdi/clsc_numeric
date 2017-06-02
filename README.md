# CLSC OnlineShop Model Numeric Solver

This project is part of my master's thesis. The goal is to create a full factorial analysis of some (new) mathematical  closed-loop supply chain models.

### Installation

You can download this project by clicking *Clone or download* and *Download as ZIP* or directly cloning it via [git](https://try.github.io/levels/1/challenges/).

##### Software requirements
This solver requires [Python3](https://www.python.org/) and [SciPy](https://www.scipy.org/) to run.

On Unix environments this may will work:

```sh
$ sudo apt-get install python3 python-scipy
```

Windows Users are encouraged to download the [Anaconda Installer](https://www.continuum.io/downloads) - This will automatically install Python3 and a lot of scientific stuff on your PC.

# Usage
documentation in progress.
If you want to check whether the solver is working - just open a terminal, change to your downloaded folder and start *generator.py* like:

```sh
# change to your Downloads folder and unzip the downloaded package:
$ cd ~/Downloads/
$ unzip clsc_numeric-master.zip
# now change to our unpacked project folder:
$ cd clsc_numeric-master
#
#
# and start my test cases:
$ python3 test_solver.py
#
#
# or print out some numeric results:
$ python3 generator.py --model 1 --output stdout
# or create a .csv file of some numeric results:
$ python3 generator.py --model 1 --output numeric_results.csv
```