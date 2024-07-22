#!/bin/bash

# The test script doesn't test all option combinations!

$1 -a 1 -g 1 -t 2 -p "Accuracy test: iterative, factors, powers-of-2"
$1 -a 2 -g 1 -t 2 -p "Accuracy test: recursive, factors, powers-of-2"
$1 -a 1 -g 3 -r 16 -t 2 -p "Accuracy test: iterative, thresholded (16), powers-of-2"