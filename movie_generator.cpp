#include "movie_generator.h"

MovieGenerator::MovieGenerator(int width, int height) 
    : width(width), height(height), framesDir("frames"), dataDir("data"), videoDir("movies") {

    // create the movies directory if it doesn't exist
    if (!fs::exists(videoDir)) {
        fs::create_directory(videoDir);
    }
    if (!fs::exists(dataDir)) {
        fs::create_directory(dataDir);
    }
    
    if (fs::exists(framesDir)) {
        fs::remove_all(framesDir);
    }
    fs::create_directory(framesDir);
}

MovieGenerator::~MovieGenerator() {
    // Nothing specific to clean up
}


void MovieGenerator::generateMovie(const std::string& simName, int fps) {
    std::string dataFilePath = dataDir + "/" + simName + "_all_steps.csv";
    std::string vesselFilePath = dataDir + "/" + simName + "_vessels.csv";

    std::string outputVideo = videoDir + "/" + simName + ".mp4";
    
    std::cout << "Generating movie from data: " << dataFilePath << std::endl;

    readCellDataFromFile(dataFilePath);
    readVesselDataFromFile(vesselFilePath);  // Read vessel data
    
    generateFrames();
    
    createVideoFromFrames(outputVideo, fps);
    
    std::cout << "Movie generation complete: " << outputVideo << std::endl;
}

void MovieGenerator::readVesselDataFromFile(const std::string& vesselFilePath) {
    std::ifstream file(vesselFilePath);
    if (!file.is_open()) {
        std::cerr << "Could not open vessel file: " << vesselFilePath << std::endl;
        return;  // If file doesn't exist, we just won't draw vessels
    }
    
    std::string line;
    // Skip header line
    std::getline(file, line);
    
    // Read data lines
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        
        VesselData vessel;
        
        // Parse vessel id
        std::getline(ss, token, ',');
        vessel.id = std::stoi(token);
        
        // Parse start_x
        std::getline(ss, token, ',');
        vessel.start_x = std::stof(token);
        
        // Parse start_y
        std::getline(ss, token, ',');
        vessel.start_y = std::stof(token);
        
        // Parse end_x
        std::getline(ss, token, ',');
        vessel.end_x = std::stof(token);
        
        // Parse end_y
        std::getline(ss, token, ',');
        vessel.end_y = std::stof(token);
        
        // Parse radius
        std::getline(ss, token, ',');
        vessel.radius = std::stof(token);
        
        // Add to vessel data
        vesselData.push_back(vessel);
    }
    
    std::cout << "Read data for " << vesselData.size() << " blood vessels" << std::endl;
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
        std::getline(ss, token, ',');
        cell.type = static_cast<CellType>(std::stoi(token));
        
        // Parse x coordinate
        std::getline(ss, token, ',');
        cell.x = std::stof(token);
        
        // Parse y coordinate
        std::getline(ss, token, ',');
        cell.y = std::stof(token);
        
        // Skip dx and dy
        std::getline(ss, token, ',');  // dx
        std::getline(ss, token, ',');  // dy

        // parse clone_id
        std::getline(ss, token, ',');
        cell.clone_id = std::stoi(token);

        // parse vessel_neighbourhood
        std::getline(ss, token, ',');
        cell.vessel_neighbourhood = std::stoi(token);
        
        // Parse status
        std::getline(ss, token, ',');
        cell.active = (token == "0"); // Assuming 0 means active
        
        // Get radius based on cell type
        cell.radius = CELL_RADII.at(cell.type);
        
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
        
        // Draw blood vessels (only once per frame)
        for (const auto& vessel : vesselData) {
            // Calculate the direction vector of the vessel
            float dx = vessel.end_x - vessel.start_x;
            float dy = vessel.end_y - vessel.start_y;
            float length = sqrt(dx * dx + dy * dy);
            
            // Apply scaling factor for blood vessels
            float scaledRadius = vessel.radius;
            
            // Create a convex shape for the vessel (rectangle with rounded ends)
            sf::ConvexShape vesselShape;
            vesselShape.setPointCount(4);
            
            // Normal vector (perpendicular to vessel direction)
            float nx = -dy / length;
            float ny = dx / length;
            
            // Set the points of the rectangle using the scaled radius
            vesselShape.setPoint(0, sf::Vector2f(vessel.start_x + nx * scaledRadius, vessel.start_y + ny * scaledRadius));
            vesselShape.setPoint(1, sf::Vector2f(vessel.start_x - nx * scaledRadius, vessel.start_y - ny * scaledRadius));
            vesselShape.setPoint(2, sf::Vector2f(vessel.end_x - nx * scaledRadius, vessel.end_y - ny * scaledRadius));
            vesselShape.setPoint(3, sf::Vector2f(vessel.end_x + nx * scaledRadius, vessel.end_y + ny * scaledRadius));
            
            // Set vessel color (dark red for blood vessels)
            sf::Color vesselColor(120, 0, 0);
            vesselShape.setFillColor(vesselColor);
            
            // Draw the vessel
            renderTexture.draw(vesselShape);
            
            // Draw rounded caps at both ends of the vessel with scaled radius
            sf::CircleShape startCap(scaledRadius);
            startCap.setPosition(sf::Vector2f(vessel.start_x - scaledRadius, vessel.start_y - scaledRadius));
            startCap.setFillColor(vesselColor);
            renderTexture.draw(startCap);
            
            sf::CircleShape endCap(scaledRadius);
            endCap.setPosition(sf::Vector2f(vessel.end_x - scaledRadius, vessel.end_y - scaledRadius));
            endCap.setFillColor(vesselColor);
            renderTexture.draw(endCap);
        }
        
        // Draw cells over the vessels
        for (const auto& cell : cells) {
            if (cell.active) {
                size_t pointCount = 30;
                sf::CircleShape cellShape(cell.radius, pointCount);
                // cellShape.setOrigin(sf::Vector2f(cell.radius, cell.radius));
                cellShape.setPosition(sf::Vector2f(cell.x, cell.y));
                
                const auto& color = CELL_COLORS.at(cell.type);
                sf::Color cellColor = sf::Color(std::get<0>(color), std::get<1>(color), std::get<2>(color));
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