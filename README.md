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
./main [options]
```

## Generate Visualization

After running the simulation, you can generate a movie from the saved frames:

```bash
./gen_movie_main [options]
```

Output frames are saved in the `frames` directory and movies in the `movies` directory. 