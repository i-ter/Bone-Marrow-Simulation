# Bone Marrow Simulation

A parallelized particle-based simulation of cellular dynamics within the bone marrow.

## Overview

This project simulates the behavior and interactions of hematopoietic cells in the bone marrow of a mouse. Key features include:

- Cell motility: Cells move throughout in Euclidean space
- Volume exclusion: Inter-cell collisions modeled via Hookean (spring) potentials
- Full hematopoietic lineage: Realistic modeling of division along the complete hematopoietic lineage tree
- Cell exit and death: Cells can leave the marrow through blood vessels or die, reflecting physiological turnover

This simulation provides a computational framework for studying spatial organization and lineage dynamics of hematopoietic cells in a biologically relevant context. The main time-step of the simulation is 10 mins while pairwise interactions are considered every 1 min. 

## Build Instructions

The Makefile contains instructions for `clang` compiler

```bash

git clone https://github.com/i-ter/Bone-Marrow-Simulation.git
cd Bone-Marrow-Simulation
make all

# Run the simulation
./main [options] [--help for more information]
```

## Generate Visualization

After running the simulation, you can generate a movie with the following command:

```bash
./gen_movie_main [options]
```

Output frames are saved in the `frames` directory and movies in the `movies` directory. It uses the SFML library to render the frames.

## Design

The main classes are:

- `BoneMarrow`: Main class that simulates the bone marrow.
- `Cell`: Class that represents a cell.
- `BloodVessel`: Class that represents a blood vessel.
- `SpatialGrid`: Class that represents the spatial grid and keeps track of all the cells.

All these classes are defined in `main.cpp`.

`cell_config.h` contains the configuration for the simulation, including cell types, motilities, colors, and initial numbers.

`movie_generator.cpp` is used to generate the movie from the frames, and is called by `gen_movie_main.cpp`.

## Parallelization

The simulation is parallelized using OpenMP. The main parallelized parts are the cell movement updates and interaction calculations. It uses the `#pragma omp parallel for` directive.

There is a fallback to sequential execution if OpenMP is not available, which is supported out of the box.

## How to run on Imperial HPC cluster

```bash
ssh -X {username}@login.cx3.hpc.imperial.ac.uk
git clone https://github.com/i-ter/Bone-Marrow-Simulation.git
cd Bone-Marrow-Simulation
bash hpc_compile.sh
qsub run_sim.pbs
qstat -u {username}
```

