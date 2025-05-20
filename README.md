# Bone Marrow Simulation

A particle simulation of cells inside the bone marrow. This simulation models the movement and interactions of cells within the bone marrow.

<details>
<summary><strong>Build Instructions</strong></summary>

```bash
# Clone the repository
git clone <repository-url>
cd <repository-name>

# Build the simulation
make all

# Run the simulation
./main [options] [--help for more information]
```

</details>

<details>
<summary><strong>Generate Visualization</strong></summary>

After running the simulation, you can generate a movie from the saved frames:

```bash
./gen_movie_main [options]
```

Output frames are saved in the `frames` directory and movies in the `movies` directory.

</details>

<details>
<summary><strong>Design</strong></summary>

The simulation is implemented in C++ and uses the SFML library for rendering. The main classes are:

- <code>BoneMarrow</code>: Main class that simulates the bone marrow.
- <code>Cell</code>: Class that represents a cell.
- <code>BloodVessel</code>: Class that represents a blood vessel.
- <code>SpatialGrid</code>: Class that represents the spatial grid and keeps track of all the cells.

All these classes are defined in <code>main.cpp</code>.

<code>cell_config.h</code> contains the configuration for the simulation, including cell types, motilities, colors, and initial numbers.

<code>movie_generator.cpp</code> is used to generate the movie from the frames, and is called by <code>gen_movie_main.cpp</code>.

</details>

<details>
<summary><strong>Parallelization</strong></summary>

The simulation leverages parallel computing to improve performance and scalability. Computationally intensive tasks, such as cell movement updates and interaction calculations, are parallelized using modern C++ techniques (e.g., OpenMP, std::thread, or other parallel libraries). This allows the simulation to efficiently utilize multiple CPU cores, significantly speeding up large-scale or high-resolution simulations.

<!--
If you use a specific library or approach (e.g., OpenMP, TBB, custom thread pools), you can add more details here:
- Which parts of the code are parallelized
- How to control the number of threads (e.g., via environment variables or config)
- Any caveats or tips for running in parallel
-->

</details>

