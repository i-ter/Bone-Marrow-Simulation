# Bone Marrow Simulation

A particle simulation of cells inside the bone marrow. This simulation models the movement and interactions of cells within the bone marrow.

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

After running the simulation, you can generate a movie from the saved frames:

```bash
./gen_movie_main [options]
```

Output frames are saved in the `frames` directory and movies in the `movies` directory. 


## Design

The simulation is implemented in C++ and uses the SFML library for rendering. The main classes are:

- `BoneMarrow`: Main class that simulates the bone marrow.
- `Cell`: Class that represents a cell.
- `BloodVessel`: Class that represents a blood vessel.
- `SpatialGrid`: Class that represents the spatial grid and keeps track of all the cells.
They are all defined in `main.cpp`.

`cell_config.h` contains the configuration for the simulation. It includes the cell types, their motilities, colors, and initial numbers.

`movie_generator.cpp` is used to generate the movie from the frames. It is called by `gen_movie_main.cpp`.

