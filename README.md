# Individual Based Simulator Project

## Project Description

Project creating individual based simulator for running epidemics.  Also investigating possible rate handling structures and optimisation.

## Basic Usage
For help on command line arguments, run `python .\simulator.py -h`.

Running the simulator always requires three main files:
1. Configuration file (e.g. config.ini)
1. Host position file (e.g. hosts.txt)
1. Host initial conditions file (e.g. hosts_init.txt)

The configuration file is a .ini key = value file giving many options for the simulations.  An example file with the minimal number of default parameters can be generated by running `python .\simulator.py -d`, or alternatively, with the `-k` flag the keyfile detailing all possible options is generated.

The first line of the two host files gives the number of hosts.  There is then one line per host with space delimited x,y positions or initial state string, as appropriate.

For a simple example see the example/ folder.

## Project Files and Directories

* README.md This document
* simulator.py Main module for simulator project.  Can be imported or run from command line
* code/ Directory containing majority of code for the simulator
* example/ Basic example of simulator usage
* test/ Testing code
* utilities/ Helpful modules for using the simulator

## To Do List
1. Documentation
  * Basic usage instructions within README
  * Example folder containing simple example with all necessary config files
  * Sphinx documentation
1. Improve intervention handling
  * Allow for intervention arguments e.g. budgets - not sure how best to implement this.
1. Testing
  * Fix calls to old simulator references
  * Implement unit tests
  * Make all tests part of unittest module


## Major Change Log
*(For more details see SVN log)*

### 13/06/2017
Log file now produced on output.  \_\_init\_\_ files so that functions as project.  __Note this will have broken some previous imports.__  Changed functionality of config parser so can take a string as well as a file.  Added utilities for reading output files and parameters.  BUGFIX in creating hosts - was not adding hostID.

### 12/06/2017
Host input changed significantly.  Now takes multiple files: a hosts.txt file detailing host positions, a hosts_init.txt file giving initial states, and optionally a hosts_regions.txt file giving the region for each host.  Also added some utilities for creating these files.

### 08/06/2017
Updated output functions - creates dataframe and then saves to csv if the key 'OutputFiles' is True.  Also handles file output stub properly.

### 07/06/2017
Added main function to simulator module - this is now what is called if invoked from the command line.  This way when invoked from within python, functionality is the same.  Have also implemented return of data when invoked from python.  Have used pandas data frames to hold data.
