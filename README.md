# Group 5: Implementation of a Fast Fourier Transform

## Introduction

Implementation of mixed-radix, depth-first and breadth-first Fast Fourier transform.
The main algorithms are contained in ffts.hpp as function templates.

## Requirements
The external library [FFTW](https://www.fftw.org/download.html) is used.
Moreover, for the execution of the test scripts awk and [gnuplot](http://www.gnuplot.info/) are required.

## Compilation
A Makefile is provided that uses the shared Makefile provided during the course.
It has to be placed into ../shared/shared.mak.
The numerical tests for performance where compiled using the *fast* make target.

## Testing
The implementations can be tested using testit.cpp. 
Use ./testit -h to get help with the test options.

## Results
Results of the conducted experiments are found in the results/ folder. 
There also the scripts for the different tests can be found. 
