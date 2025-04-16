#pragma once

#include <SFML/Graphics.hpp>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iomanip>
#include <cstdlib>
#include "cell_config.h"

namespace fs = std::filesystem;

class Cell;

class MovieGenerator {
private:
    int width, height;
    std::string framesDir;
    std::string dataDir;
    std::string videoDir;
    
    // Cell data structure for reading from CSV
    struct CellData {
        CellType type;
        float x, y;
        float radius;
        bool active;
    };
    
    // Vessel data structure for reading from CSV
    struct VesselData {
        int id;
        float start_x, start_y;
        float end_x, end_y;
        float radius;
    };
    
    // Map of step number to vector of cells for that step
    std::map<int, std::vector<CellData>> stepData;
    
    // Vector to store blood vessel data
    std::vector<VesselData> vesselData;
    
    void readCellDataFromFile(const std::string& filename);
    void readVesselDataFromFile(const std::string& vesselFilePath);
    void generateFrames();
    void createVideoFromFrames(const std::string& outputVideo, int fps);

public:
    MovieGenerator(int width, int height);
    ~MovieGenerator();
    
    // Main function to generate movie from a data file
    void generateMovie(const std::string& simName, int fps = 30);
}; 