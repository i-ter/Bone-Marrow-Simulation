#include "movie_generator.h"
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {
    try {
        // Default simulation name
        std::string simName = "bm_sim";
        
        // If a simulation name is provided as command line argument, use it
        if (argc > 1) {
            simName = argv[1];
        }
        
        std::cout << "Generating movie for simulation: " << simName << std::endl;

        std::string paramsFilename = "data/" + simName + "_params.csv";
        std::ifstream paramsFile(paramsFilename);
        if (!paramsFile.is_open()) {
            throw std::runtime_error("Failed to open parameters file: " + paramsFilename);
        }

        int width, height;
        std::string line;
        while (std::getline(paramsFile, line)) {
            if (line.find("width") != std::string::npos) {
                width = std::stoi(line.substr(line.find(",") + 1));
            }
            if (line.find("height") != std::string::npos) {
                height = std::stoi(line.substr(line.find(",") + 1));
            }
        }
        paramsFile.close();
        std::cout << "width: " << width << " height: " << height << std::endl;

        // Create a movie generator with the same dimensions as the simulation
        MovieGenerator movieGen(width, height);
        
        // Generate the movie with 30 fps
        movieGen.generateMovie(simName, 30);
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
} 
