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

namespace fs = std::filesystem;

class Cell;

class MovieGenerator {
private:
    int width, height;
    std::string framesDir;
    std::string dataDir;
    std::string videoDir;
    std::map<std::string, sf::Color> cellColors;
    std::map<std::string, float> cellRadii;
    
    // Cell data structure for reading from CSV
    struct CellData {
        std::string type;
        float x, y;
        float radius;
        bool active;
    };
    
    // Map of step number to vector of cells for that step
    std::map<int, std::vector<CellData>> stepData;
    
    void createFramesDirectory();
    void loadCellColors();
    void loadCellRadii();
    void readCellDataFromFile(const std::string& filename);
    void generateFrames();
    void createVideoFromFrames(const std::string& outputVideo, int fps);

public:
    MovieGenerator(int width, int height);
    ~MovieGenerator();
    
    // Main function to generate movie from a data file
    void generateMovie(const std::string& simName, int fps = 30);
}; 