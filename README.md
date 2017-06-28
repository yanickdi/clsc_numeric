# CLSC OnlineShop Model Numeric Solver

This project is part of my master's thesis. The goal is to create a full factorial analysis of some (new) mathematical  closed-loop supply chain models.

## Installation

You can download this project by clicking *Clone or download* and *Download as ZIP* or directly cloning it via [git](https://try.github.io/levels/1/challenges/1).

### Download via git (recommended)
First, you have to install git:
```sh
$ sudo apt-get install git
```
To download the repository, clone it to your pc using:
```sh
$ git clone https://github.com/yanickdi/clsc_numeric.git
```
A new folder called *clsc_numeric* was created.
Make sure to update your local copy by and by, using the git pull command:
```sh
$ cd clsc_numeric
$ git pull
```

### Download via direct zip download
Just download and unpack the zip folder, you can use your terminal:
```sh
# change to your Downloads folder:
$ cd ~/Downloads/
# download the zip file using wget
$ wget https://github.com/yanickdi/clsc_numeric/archive/master.zip
# and unzip the downloaded package:
$ unzip clsc_numeric-master.zip
$ mv clsc_numeric-master clsc_numeric
```

### Software requirements
This solver requires [Python3](https://www.python.org/) and [SciPy](https://www.scipy.org/) to run. To create html output you will also need [Jinja2 Template Engine](http://jinja.pocoo.org/).

On Unix environments this may will work:

```sh
$ sudo apt-get install python3 python3-scipy python3-pip
$ sudo pip3 install jinja2
```

# Usage
documentation in progress.
If you want to check whether the solver is working - just open a terminal, change to your downloaded folder and start *generator.py* like:

```sh
# change to your project folder:
$ cd clsc_numeric
#
# start all test cases:
$ python3 tests.py
# start test cases for model one:
$ python3 tests.py TestModelOneNumericalSolver
#
#

# CALC SOME NUMERIC RESULTS

$ python3 generator.py -model 1 -output stdout
# or create a .csv file of some numeric results:
$ python3 generator.py -model 1 -output numeric_results.csv

# csv results for model 2:
$ python3 generator.py -model 2 -output model_2.csv --german-comma

# html results for model 1:
$ python3 generator.py -model 1 -output model_2.html

# want to see only valid solutions:
$ python3 generator.py -model 1 -output model_2.html --only-valid-solutions

# PLOTTING

$ python3 plotter.py -plot profit-difference
or
$ python3 plotter.py -plot rho-difference
or
$ python3 plotter.py -plot cases-model-one
or
$ python3 plotter.py -plot cases-model-two
# if you are in a hurry
$ python3 plotter.py -plot any-plot-from-above --low-qual
# plot profit-difference and rho-difference not only be calculated with relative increase but with absolute:
$ python3 plotter.py -plot profit-difference --absolute
# to get a gray-scale colour map:
$ python3 plotter.py -plot cases-model-two --gray
```