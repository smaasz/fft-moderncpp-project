# Implementation of a Fast Fourier Transform in Modern C++

## Introduction

Implementation of mixed-radix, depth-first and breadth-first Fast Fourier transform.
The main algorithms are contained in ffts.hpp as function templates.

## Requirements
The external library [FFTW](https://www.fftw.org/download.html) is used.
Moreover, for the execution of the test scripts awk and [gnuplot](http://www.gnuplot.info/) are required.

## Compilation
A Makefile is provided.
The numerical tests for performance where compiled using the *fast* make target.

## Testing
The implementations can be tested using testit.cpp. 
Use ./testit -h to get help with the test options.

## Results
Results of the conducted experiments are found in the results/ folder. 
There also the scripts for the different tests can be found. 

## Example
```
make fast
./testit -a 1 -g 1 -t 1 -p "Performance test: iterative, factors, powers-of-2"
```
Output (the timed numbers are hardware-dependent):
```
Performance test: iterative, factors, powers-of-2
  size   time (ms)   repetitions   average (ms)
    64      0.34          10        0.0337
   128      0.56          10        0.0561
   256      1.15          10        0.1154
   512      2.48          10        0.2479
  1024      5.60          10        0.5597
  2048     11.43          10        1.1428
  4096     25.45          10        2.5448
  8192     51.63          10        5.1628
 16384    109.67          10       10.9670
 32768    244.16          10       24.4157
```