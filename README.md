# Bone Marrow Simulation

A particle simulation of cells inside the bone marrow. This simulation models the movement and interactions of cells within the bone marrow. It is parallelized using OpenMP.

## Build Instructions

```bash
# Clone the repository
git clone <repository-url>
cd <repository-name>

# Build the simulation
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

