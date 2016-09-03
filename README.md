# Sampling and Topological Data Analysis for real algebraic varieties


## Description 
Pipeline for performing topological data analysis to real algebraic varieties. 

## Usage

Before trying any examples, run the following with the root directory of the project as the working directory: 

`source source/environment_variables`

This script exports environment variables necessary for the Rtree library to function. After doing this, the usage syntax is as follows (with current working directory the project's root directory): 

`bin/python source/run_script.py DIRECTORY_NAME`


If you receive a DIPHA related error, you may need to remake the DIPHA executable. This can be done by changing your working directory to the dipha directory and running: 

`cmake .`
`make`

## Requirements

Matlab must be installed on your computer and runnable via the command "matlab". 