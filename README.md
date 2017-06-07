# Individual Based Simulator Project

## Project Description

Project creating individual based simulator for running epidemics.  Also investigating possible rate handling structures and optimisation.

## Project Files and Directories

* README.md This document

## To Do List
1. Tidy up project structure
  * Rate structures contained within separate folder
  * All code within separate folder, so that just simulator.py in main folder
  * Documentation
2. ~~Main function so that can invoke simulator from python with similar interface to command line~~
3. Improve intervention handling
  * Allow for intervention arguments e.g. budgets - not sure how best to implement this.
4. Output should include log file or similar that details configuration options used, how simulator called, code version etc.
5. Improve output options
  * Turn on/off file output for each type? (host, event, summary)
  * Work out how best to do python data handling


## Log

### 07/06/2017
Added main function to simulator module - this is now what is called if invoked from the command line.  This way when invoked from within python, functionality is the same.
