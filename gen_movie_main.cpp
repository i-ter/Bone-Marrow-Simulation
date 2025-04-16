#include "MovieGenerator.h"
#include <iostream>

int main(int argc, char* argv[]) {
    try {
        // Default simulation name
        std::string simName = "bm_sim";
        
        // If a simulation name is provided as command line argument, use it
        if (argc > 1) {
            simName = argv[1];
        }
        
        std::cout << "Generating movie for simulation: " << simName << std::endl;
        
        // Create a movie generator with the same dimensions as the simulation
        MovieGenerator movieGen(1500, 1000);
        
        // Generate the movie with 30 fps
        movieGen.generateMovie(simName, 30);
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
} 