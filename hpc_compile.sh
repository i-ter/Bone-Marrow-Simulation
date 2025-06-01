#!/bin/bash

module load GCC OpenMPI

g++ -std=c++17 -fopenmp -O3 -march=icelake-server -Wall -Wextra main.cpp -o main