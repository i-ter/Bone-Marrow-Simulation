#include "MovieGenerator.h"

MovieGenerator::MovieGenerator(int width, int height) 
    : width(width), height(height), framesDir("frames"), dataDir("data"), videoDir("movies") {

    // create the movies directory if it doesn't exist
    if (!fs::exists(videoDir)) {
        fs::create_directory(videoDir);
    }
    if (!fs::exists(dataDir)) {
        fs::create_directory(dataDir);
    }
    if (!fs::exists(framesDir)) {
        fs::create_directory(framesDir);
    }
    loadCellColors();
    loadCellRadii();
}

MovieGenerator::~MovieGenerator() {
    // Nothing specific to clean up
}

void MovieGenerator::loadCellColors() {
    // Color mapping for each cell type
    cellColors["HSC"] = sf::Color(255, 0, 0);           // Red
    cellColors["MPP"] = sf::Color(255, 165, 0);         // Orange
    cellColors["CMP"] = sf::Color(0, 255, 0);           // Green
    cellColors["CLP"] = sf::Color(0, 255, 0);           // Green
    cellColors["Granulocyte"] = sf::Color(0, 255, 255); // Cyan
    cellColors["Erythrocyte"] = sf::Color(0, 255, 255); // Cyan
    cellColors["Lymphocyte"] = sf::Color(0, 255, 255);  // Cyan
}

void MovieGenerator::loadCellRadii() {
    // Map for cell radii
    cellRadii["HSC"] = 5.0;
    cellRadii["MPP"] = 4.0;
    cellRadii["CMP"] = 4.0;
    cellRadii["CLP"] = 4.0;
    cellRadii["Granulocyte"] = 3.0;
    cellRadii["Erythrocyte"] = 3.0;
    cellRadii["Lymphocyte"] = 3.0;
}

void MovieGenerator::generateMovie(const std::string& simName, int fps) {
    std::string dataFilePath = dataDir + "/" + simName + "_all_steps.csv";
    std::string outputVideo = videoDir + "/" + simName + ".mp4";
    
    std::cout << "Generating movie from data: " << dataFilePath << std::endl;

    readCellDataFromFile(dataFilePath);
    
    generateFrames();
    
    createVideoFromFrames(outputVideo, fps);
    
    std::cout << "Movie generation complete: " << outputVideo << std::endl;
}


void MovieGenerator::readCellDataFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    std::string line;
    // Skip header line
    std::getline(file, line);
    
    // Read data lines
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        
        // Parse step number
        std::getline(ss, token, ',');
        int step = std::stoi(token);
        
        CellData cell;
        
        // Parse cell type
        std::getline(ss, cell.type, ',');
        
        // Parse x coordinate
        std::getline(ss, token, ',');
        cell.x = std::stof(token);
        
        // Parse y coordinate
        std::getline(ss, token, ',');
        cell.y = std::stof(token);
        
        // Skip dx and dy
        std::getline(ss, token, ',');  // dx
        std::getline(ss, token, ',');  // dy
        
        // Parse status
        std::getline(ss, token, ',');
        cell.active = (token == "0"); // Assuming 0 means active
        
        // Get radius based on cell type
        cell.radius = cellRadii[cell.type];
        
        // Add to step data
        stepData[step].push_back(cell);
    }
    
    std::cout << "Read data for " << stepData.size() << " simulation steps" << std::endl;
}

void MovieGenerator::generateFrames() {
    if (stepData.empty()) {
        throw std::runtime_error("No data to generate frames from");
    }
    
    // Setup render texture for drawing frames using constructor
    sf::RenderTexture renderTexture(sf::Vector2u(width, height));
    
    // Generate frames for each step
    for (const auto& [step, cells] : stepData) {
        renderTexture.clear(sf::Color::Black);
        
        for (const auto& cell : cells) {
            if (cell.active) {
                // Create a circle shape for each cell
                sf::CircleShape cellShape(cell.radius);
                cellShape.setOrigin(sf::Vector2f(cell.radius, cell.radius));
                cellShape.setPosition(sf::Vector2f(cell.x, cell.y));
                
                sf::Color cellColor = cellColors[cell.type];
                cellShape.setFillColor(cellColor);
                
                renderTexture.draw(cellShape);
            }
        }
        
        renderTexture.display();
        
        // Save the current frame as an image
        sf::Image screenshot = renderTexture.getTexture().copyToImage();
        std::stringstream ss;
        ss << framesDir << "/frame_" << std::setw(4) << std::setfill('0') << step << ".png";
        if (!screenshot.saveToFile(ss.str())) {
            throw std::runtime_error("Failed to save frame " + std::to_string(step));
        }
        
        if (step % 100 == 0) {
            std::cout << "Generated frame " << step << std::endl;
        }
    }
}

void MovieGenerator::createVideoFromFrames(const std::string& outputVideo, int fps) {
    std::string command = "ffmpeg -y -framerate " + std::to_string(fps) +
                     " -i " + framesDir + "/frame_%04d.png" +
                     " -c:v libx264 -pix_fmt yuv420p " + outputVideo;

    std::cout << "Creating video with command: " << command << std::endl;
    int result = system(command.c_str());

    if (result == 0) {
        std::cout << "Video created successfully: " << outputVideo << std::endl;
    } else {
        std::cerr << "Failed to create video. Make sure FFmpeg is installed." << std::endl;
    }
} 