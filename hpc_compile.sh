#!/bin/bash

module load GCC OpenMPI

g++ -std=c++17 -fopenmp -Wall -Wextra -g main.cpp -o main
