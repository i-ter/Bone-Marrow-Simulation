#include <iostream>
#include <vector>
#include <random>
#include <map>
#include <memory>
#include <string>
#include <filesystem>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>  // Add for file operations
#include "cell_config.h"  // Import the cell configuration

using namespace std;
namespace fs = std::filesystem;


// Random number generator
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> unif_01(0, 1);
uniform_real_distribution<float> angle_dist(0, 2 * M_PI);
uniform_real_distribution<float> angle_change(-0.5, 0.5);
// uniform_real_distribution<float> vessel_length_dist(100, 600);  // Blood vessel length distribution
// uniform_real_distribution<float> vessel_radius_dist(20.0, 50);    // Blood vessel radius distribution


class BloodVessel
{
public:
    // Start and end points of the blood vessel
    float start_x, start_y;
    float end_x, end_y;
    float radius;  // Radius of the blood vessel
    float length;  // Length of the blood vessel
    float angle;   // Angle of the blood vessel
    
    // Constructor for a blood vessel
    BloodVessel(float start_x, float start_y, float length, float radius, float angle) 
        : start_x(start_x), start_y(start_y), radius(radius), length(length), angle(angle) {
        // Calculate end point based on start point, length, and angle
        end_x = start_x + length * cos(angle);
        end_y = start_y + length * sin(angle);
    }
    
    // Check if a point is inside this blood vessel
    bool contains(float x, float y) const {
        // Vector from start point to test point
        float v1x = x - start_x;
        float v1y = y - start_y;
        
        // Vector representing the vessel direction
        float v2x = end_x - start_x;
        float v2y = end_y - start_y;
        
        // Length of vessel vector
        float vessel_length_squared = v2x * v2x + v2y * v2y;
        
        // Calculate dot product
        float dot_product = v1x * v2x + v1y * v2y;
        
        // Calculate projection ratio (how far along the vessel line)
        float projection_ratio = dot_product / vessel_length_squared;
        
        // If projection is outside the vessel segment, point is not in vessel
        if (projection_ratio < 0.0 || projection_ratio > 1.0) {
            return false;
        }
        
        // Find closest point on vessel line to test point
        float closest_x = start_x + projection_ratio * v2x;
        float closest_y = start_y + projection_ratio * v2y;
        
        // Calculate distance from test point to closest point on vessel line
        float distance_squared = (x - closest_x) * (x - closest_x) + (y - closest_y) * (y - closest_y);
        
        // Check if distance is less than vessel radius
        return distance_squared <= radius * radius;
    }
    
    // Calculate distance from point to vessel centerline
    float distanceFrom(float x, float y) const {
        // Vector from start point to test point
        float v1x = x - start_x;
        float v1y = y - start_y;
        
        // Vector representing the vessel direction
        float v2x = end_x - start_x;
        float v2y = end_y - start_y;
        
        // Length of vessel vector
        float vessel_length_squared = v2x * v2x + v2y * v2y;
        
        // Calculate dot product
        float dot_product = v1x * v2x + v1y * v2y;
        
        // Calculate projection ratio (how far along the vessel line)
        float projection_ratio = dot_product / vessel_length_squared;
        
        // Clamp projection to vessel segment
        projection_ratio = max(0.0f, min(1.0f, projection_ratio));
        
        // Find closest point on vessel line to test point
        float closest_x = start_x + projection_ratio * v2x;
        float closest_y = start_y + projection_ratio * v2y;
        
        // Calculate distance from test point to closest point on vessel line
        return sqrt((x - closest_x) * (x - closest_x) + (y - closest_y) * (y - closest_y));
    }
    
    // Calculate vessel-to-vessel intersection
    bool intersectsWith(const BloodVessel& other) const {
        // Check if any endpoint of either vessel is inside the other vessel
        if (distanceFrom(other.start_x, other.start_y) <= other.radius + radius) return true;
        if (distanceFrom(other.end_x, other.end_y) <= other.radius + radius) return true;
        if (other.distanceFrom(start_x, start_y) <= other.radius + radius) return true;
        if (other.distanceFrom(end_x, end_y) <= other.radius + radius) return true;
        
        // Calculate minimum distance between centerlines
        // Using the formula for distance between two line segments
        float ux = end_x - start_x;
        float uy = end_y - start_y;
        float vx = other.end_x - other.start_x;
        float vy = other.end_y - other.start_y;
        float wx = start_x - other.start_x;
        float wy = start_y - other.start_y;
        
        float a = ux * ux + uy * uy;        // Squared length of u
        float b = ux * vx + uy * vy;        // Dot product of u and v
        float c = vx * vx + vy * vy;        // Squared length of v
        float d = ux * wx + uy * wy;        // Dot product of u and w
        float e = vx * wx + vy * wy;        // Dot product of v and w
        
        float D = a * c - b * b;            // Denominator of the equations
        float sc, sN, sD = D;               // s = sN / sD, default sD = D >= 0
        float tc, tN, tD = D;               // t = tN / tD, default tD = D >= 0
        
        // Compute parameters and clamp them to respective intervals
        if (D < 0.00001f) {
            sN = 0.0;  sD = 1.0;  tN = e;  tD = c;
        } else {
            sN = (b * e - c * d);
            tN = (a * e - b * d);
            if (sN < 0.0) {
                sN = 0.0;  tN = e;  tD = c;
            } else if (sN > sD) {
                sN = sD;  tN = e + b;  tD = c;
            }
        }
        
        if (tN < 0.0) {
            tN = 0.0;
            if (-d < 0.0)      sN = 0.0;
            else if (-d > a)   sN = sD;
            else {             sN = -d;  sD = a; }
        } else if (tN > tD) {
            tN = tD;
            if ((-d + b) < 0.0)  sN = 0;
            else if ((-d + b) > a)  sN = sD;
            else {                  sN = (-d + b);  sD = a; }
        }
        
        sc = (abs(sN) < 0.00001f) ? 0.0 : sN / sD;
        tc = (abs(tN) < 0.00001f) ? 0.0 : tN / tD;
        
        // Calculate minimum distance point
        float minimum_distance_x = start_x + sc * ux - (other.start_x + tc * vx);
        float minimum_distance_y = start_y + sc * uy - (other.start_y + tc * vy);
        
        float minimum_distance = sqrt(minimum_distance_x * minimum_distance_x + minimum_distance_y * minimum_distance_y);
        
        // Check if minimum distance is less than the sum of radii
        return minimum_distance <= (radius + other.radius);
    }
};


class Cell
{
protected:
    CellType cell_type;

public:
    float x, y;  // position
    float radius;
    float dx, dy;  // velocity
    bool is_dead = false;
    bool is_leaving = false;
    bool in_vessel_neighbourhood = false;  // Flag to track if cell is currently inside a blood vessel
    
    Cell(float x, float y, float radius, CellType type) : cell_type(type), x(x), y(y), radius(radius) {
        // Initialize random movement direction
        float angle = angle_dist(gen);
        float speed = unif_01(gen) * MOTILITY.at(cell_type); // Adjust speed by motility
        dx = cos(angle) * speed;
        dy = sin(angle) * speed;
    }

    ~Cell() {}

    const CellType& getType() const {
        return cell_type;
    }

    void move(const float& width, const float& height) {
        // If cell is not moving, randomly change direction
        if (dx == 0 && dy == 0) {
            if (unif_01(gen) < 0.3) {
                float angle = angle_dist(gen);
                float speed = unif_01(gen) * MOTILITY.at(cell_type); // Adjust speed by motility
                dx = cos(angle) * speed;
                dy = sin(angle) * speed;
            }
        }

        x += dx;
        y += dy;
        
        handleBoundaryCollision(width, height);
        
        if (unif_01(gen) < 0.2) {  // 10% chance to change direction each step
            float angle = atan2(dy, dx);
            angle += angle_change(gen);
            float speed = sqrt(dx*dx + dy*dy);
            dx = cos(angle) * speed;
            dy = sin(angle) * speed;
        }
    }

    void handleBoundaryCollision(const float& width, const float& height) {
        if (x - radius < 0) {
            x = radius;
            dx = -dx;
        } 
        else if (x + radius > width) {
            x = width - radius;
            dx = -dx;
        }
        if (y - radius < 0) {
            y = radius;
            dy = -dy;
        } else if (y + radius > height) {
            y = height - radius;
            dy = -dy;
        }
    }

    bool collidesWith(const Cell& other) const {
        float dx = x - other.x;
        float dy = y - other.y;
        float distance_squared = dx*dx + dy*dy;
        float min_distance = radius + other.radius;
        return distance_squared < min_distance * min_distance;
    }
    
    // Check if cell collides with a blood vessel
    bool inVesselNeighbour(const BloodVessel& vessel, const float& distance_threshold) const {
        return vessel.distanceFrom(x, y) <= distance_threshold;
    }
    
    // Handle collision with a blood vessel
    void resolveVesselCollision(const BloodVessel& vessel) {
        // Calculate the closest point on the vessel centerline
        float v1x = x - vessel.start_x;
        float v1y = y - vessel.start_y;
        float v2x = vessel.end_x - vessel.start_x;
        float v2y = vessel.end_y - vessel.start_y;
        float vessel_length_squared = v2x * v2x + v2y * v2y;
        float dot_product = v1x * v2x + v1y * v2y;
        float projection_ratio = max(0.0f, min(1.0f, dot_product / vessel_length_squared));
        
        float closest_x = vessel.start_x + projection_ratio * v2x;
        float closest_y = vessel.start_y + projection_ratio * v2y;
        
        // Calculate distance and normal vector from cell to vessel
        float dx = x - closest_x;
        float dy = y - closest_y;
        float distance = sqrt(dx * dx + dy * dy);
        
        // Normalize normal vector
        float nx = dx / distance;
        float ny = dy / distance;
        
        // Calculate penetration depth
        float penetration = radius - distance;
        
        if (distance < 0.1) {
            // Generate random normal if too close to centerline
            float angle = angle_dist(gen);
            nx = cos(angle);
            ny = sin(angle);
            distance = 0.1;
        }
        
        if (penetration > 0) {
            // Move cell out of the vessel wall if it's penetrating
            if (!in_vessel_neighbourhood) {
                x += nx * penetration;
                y += ny * penetration;
            }
            
            // Reflect velocity along the normal direction
            float dot = dx * nx + dy * ny;
            dx -= 2 * dot * nx;
            dy -= 2 * dot * ny;
            
            // Add some randomization to the bounce
            float speed = sqrt(dx * dx + dy * dy);
            float angle = atan2(dy, dx) + (unif_01(gen) - 0.5) * 0.2;
            dx = cos(angle) * speed;
            dy = sin(angle) * speed;
        }
    }
    
    void resolveCollision(Cell& other) {
        // Calculate displacement vector
        float x_diff = x - other.x;
        float y_diff = y - other.y;
        float distance = sqrt(x_diff * x_diff + y_diff * y_diff);

        // Ensure cells are not exactly on top of each other. leads to numerical instability.
        if (distance < 0.1) {
            distance = adjustMinimalSeparation(x_diff, y_diff);
        }

        // Calculate minimum separation distance
        float min_distance = radius + other.radius;

        // If cells are overlapping
        if (distance < min_distance) {
            float overlap = min_distance - distance;
            separateCells(other, x_diff, y_diff, distance, overlap);
            exchangeMomentum(other, x_diff, y_diff);
            capVelocities();
            other.capVelocities();
        }
    }

    float adjustMinimalSeparation(float& x_diff, float& y_diff) {
        float angle = angle_dist(gen);
        x_diff = cos(angle);
        y_diff = sin(angle);
        return 0.1;
    }

    void separateCells(Cell& other, float x_diff, float y_diff, float distance, float overlap) {
        float nx = x_diff / distance;
        float ny = y_diff / distance;
        float total_radius = radius + other.radius;
        float ratio1 = radius / total_radius;
        float ratio2 = other.radius / total_radius;

        // Move this cell
        x += nx * overlap * ratio2;
        y += ny * overlap * ratio2;

        // Move other cell
        other.x -= nx * overlap * ratio1;
        other.y -= ny * overlap * ratio1;
    }

    void exchangeMomentum(Cell& other, float nx, float ny) {
        float dot1 = this->dx * nx + this->dy * ny;
        float dot2 = other.dx * nx + other.dy * ny;

        // Update velocities
        this->dx += (dot2 - dot1) * nx;
        this->dy += (dot2 - dot1) * ny;
        other.dx += (dot1 - dot2) * nx;
        other.dy += (dot1 - dot2) * ny;
    }

    void capVelocities() {
        float speed = sqrt(dx * dx + dy * dy);
        if (speed > MAX_SPEED) {
            dx = (dx / speed) * MAX_SPEED;
            dy = (dy / speed) * MAX_SPEED;
        }
    }
    
    void setVelocity(float dx, float dy) {
        this->dx = dx;
        this->dy = dy;
    }

    bool should_divide() {
        return unif_01(gen) < DIVISION_PROB.at(cell_type);
    }

    // Check if cell should die (all cell types)
    bool should_die() {
        return unif_01(gen) < CELL_DEATH_PROB.at(cell_type);
    }
    
    // Check if cell should leave (based on probability from LEAVE_PROB map)
    bool should_leave() {
        float multiplier = in_vessel_neighbourhood ? VESSEL_LEAVING_MULTIPLIER : 1.0;
        return unif_01(gen) < LEAVE_PROB.at(cell_type) * multiplier;
    }
};

class BoneMarrow
{
public:
    float width, height;

    vector<unique_ptr<Cell>> cells;
    vector<BloodVessel> blood_vessels;  // Added blood vessels collection
    int num_vessels;  // Number of blood vessels
    
    string sim_name;
    string framesDir = "frames";
    string dataDir = "data";  // Directory for cell data files
    ofstream consolidatedDataFile;  // File stream for consolidated data
    ofstream vesselDataFile;  // File for blood vessel data
    ofstream paramsFile;      // File for simulation parameters
    
    // Statistics tracking
    struct Stats {
        int total_deaths = 0;
        int total_leaving = 0;
        map<CellType, int> deaths_by_type;
        map<CellType, int> leaving_by_type;
    } stats;
    
    // Spatial partitioning for collision detection optimization
    struct SpatialGrid {
        float block_size;
        int grid_width, grid_height;
        vector<vector<vector<Cell*>>> grid;  //grid[x][y] returns a vector of cells in block (x,y)
        
        SpatialGrid(float width, float height, float block_size) : block_size(block_size) {
            grid_width = ceil(width / block_size);
            grid_height = ceil(height / block_size);
            grid.resize(grid_width, vector<vector<Cell*>>(grid_height));
        }
        
        void clear() {
            for (auto& col : grid) {
                for (auto& block : col) {
                    block.clear();
                }
            }
        }
        
        void insert(Cell* cell) {
            int grid_x = min(grid_width - 1, max(0, static_cast<int>(cell->x / block_size)));
            int grid_y = min(grid_height - 1, max(0, static_cast<int>(cell->y / block_size)));
            grid[grid_x][grid_y].push_back(cell);
        }
        
        vector<Cell*> getPotentialCollisions(const Cell* cell) {
            vector<Cell*> result;
            int grid_x = static_cast<int>(cell->x / block_size);
            int grid_y = static_cast<int>(cell->y / block_size);
            
            // Check cell's grid cell and neighboring cells
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    int nx = grid_x + dx;
                    int ny = grid_y + dy;
                    
                    if (nx >= 0 && nx < grid_width && ny >= 0 && ny < grid_height) {
                        for (Cell* other : grid[nx][ny]) {
                            if (other != cell) {
                                result.push_back(other);
                            }
                        }
                    }
                }
            }
            
            return result;
        }
    };
    
    SpatialGrid spatial_grid;

    // Generate a random blood vessel
    BloodVessel generateRandomVessel() {
        float start_x = unif_01(gen) * width;
        float start_y = unif_01(gen) * height;
        float length = 100 + unif_01(gen) * 500;
        float radius = 20 + unif_01(gen) * 30;
        float angle = angle_dist(gen);
        
        return BloodVessel(start_x, start_y, length, radius, angle);
    }
    
    // Check if a vessel is valid (within bounds and not intersecting other vessels too much)
    bool isVesselValid(const BloodVessel& vessel) {
        // Check if vessel endpoints are within bounds with some margin
        float margin = vessel.radius;
        if (vessel.start_x < margin || vessel.start_x > width - margin ||
            vessel.start_y < margin || vessel.start_y > height - margin ||
            vessel.end_x < margin || vessel.end_x > width - margin ||
            vessel.end_y < margin || vessel.end_y > height - margin) {
            return false;
        }
        
        // Check intersection with other vessels
        // for (const auto& existing_vessel : blood_vessels) {
        //     if (vessel.intersectsWith(existing_vessel)) {
        //         return false;
        //     }
        // }
        
        return true;
    }
    
    // Generate all blood vessels
    void generateBloodVessels(int num_vessels) {
        int max_attempts = num_vessels * 10; // Maximum generation attempts
        int attempts = 0;
        
        while (blood_vessels.size() < static_cast<size_t>(num_vessels) && attempts < max_attempts) {
            BloodVessel vessel = generateRandomVessel();
            if (isVesselValid(vessel)) {
                blood_vessels.push_back(vessel);
            }
            attempts++;
        }
        
        // Write vessel data to file
        string vesselFilename = dataDir + "/" + sim_name + "_vessels.csv";
        vesselDataFile.open(vesselFilename);
        vesselDataFile << "id,start_x,start_y,end_x,end_y,radius" << endl;
        
        for (size_t i = 0; i < blood_vessels.size(); i++) {
            const auto& vessel = blood_vessels[i];
            vesselDataFile << i << ","
                          << fixed << setprecision(4) << vessel.start_x << ","
                          << fixed << setprecision(4) << vessel.start_y << ","
                          << fixed << setprecision(4) << vessel.end_x << ","
                          << fixed << setprecision(4) << vessel.end_y << ","
                          << fixed << setprecision(4) << vessel.radius << endl;
        }
        
        vesselDataFile.close();
        
        cout << "Generated " << blood_vessels.size() << " blood vessels out of " 
             << num_vessels << " requested (after " << attempts << " attempts)" << endl;
    }

    // init
    BoneMarrow(float width, float height, int initial_cells, string sim_name="bm", int num_vessels=3) 
        : width(width), height(height), num_vessels(num_vessels), sim_name(sim_name), spatial_grid(width, height, 12.)
    {
        if (!fs::exists(dataDir))
        {
            fs::create_directory(dataDir);
        }
        
        // Save simulation parameters in CSV format
        string paramsFilename = dataDir + "/" + sim_name + "_params.csv";
        paramsFile.open(paramsFilename);
        
        // Write global parameters
        paramsFile << "parameter,value" << endl
                  << "width," << width << endl
                  << "height," << height << endl
                  << "initial_cells," << initial_cells << endl
                  << "num_vessels," << num_vessels << endl
                  << "spatial_grid_block_size," << 12.0 << endl
                  << "max_speed," << MAX_SPEED << endl
                  << "vessel_distance_threshold," << VESSEL_DISTANCE_THRESHOLD << endl
                  << "vessel_leaving_multiplier," << VESSEL_LEAVING_MULTIPLIER << endl;
        paramsFile.close();
        
        // Save cell type specific parameters in a separate CSV
        string cellParamsFilename = dataDir + "/" + sim_name + "_cell_params.csv";
        paramsFile.open(cellParamsFilename);
        
        // Write header
        paramsFile << "cell_type,radius,division_probability,death_probability,leave_probability,motility" << endl;
        
        // Write cell type specific parameters
        for (const auto& [type, radius] : CELL_RADII) {
            paramsFile << getCellTypeName(type) << ","
                      << radius << ","
                      << DIVISION_PROB.at(type) << ","
                      << CELL_DEATH_PROB.at(type) << ","
                      << LEAVE_PROB.at(type) << ","
                      << MOTILITY.at(type) << endl;
        }
        paramsFile.close();
        
        // Initialize consolidated data file
        string consolidatedFilename = dataDir + "/" + sim_name + "_all_steps.csv";
        consolidatedDataFile.open(consolidatedFilename);
        consolidatedDataFile << "step,cell_type,x,y,dx,dy,status" << endl;
        
        for (int i = 0; i < initial_cells; ++i)
        {
            float x = static_cast<float>(unif_01(gen) * width);
            float y = static_cast<float>(unif_01(gen) * height);
            cells.push_back(make_unique<Cell>(x, y, CELL_RADII.at(HSC), HSC));
        }

        generateBloodVessels(num_vessels);
    }
    ~BoneMarrow() {
        // destructor to close the consolidated data file
        consolidatedDataFile.close();
    }

    void step()
    {
        vector<unique_ptr<Cell>> new_cells;
        
        // 1. Move cells
        for (auto& cell : cells) {
            cell->move(width, height);
            
        }
        
        // 2. Handle blood vessel collisions 
        for (auto& cell : cells) {
            cell->in_vessel_neighbourhood = false;
            for (const auto& vessel : blood_vessels) {
                if (cell->inVesselNeighbour(vessel, VESSEL_DISTANCE_THRESHOLD)) {
                    cell->in_vessel_neighbourhood = true;
                    break;
                }
            }
        }
        
        // 3. Handle cell-cell collisions using spatial partitioning O(n log n)
        spatial_grid.clear();
        for (auto& cell : cells) {
            spatial_grid.insert(cell.get());
        }
        
        // Check for collisions and resolve them
        for (auto& cell : cells) {
            auto potential_collisions = spatial_grid.getPotentialCollisions(cell.get());
            for (Cell* other : potential_collisions) {
                if (cell->collidesWith(*other)) {
                    cell->resolveCollision(*other);
                }
            }
        }
        
        // 4. Cell division
        for (auto& cell : cells) {
            if (cell->should_divide()) {
                if (cell->is_dead && cell->is_leaving ) {throw runtime_error("Cell is dead or leaving. something is wrong.");}
                
                auto& possible_types = LINEAGE_TREE.at(cell->getType());
                CellType new_type = possible_types[rand() % possible_types.size()];
                
                float offset = cell->radius;
                float new_x = cell->x + (unif_01(gen) * 2. - 1.) * offset;
                float new_y = cell->y + (unif_01(gen) * 2. - 1.) * offset;
                
                new_cells.push_back(make_unique<Cell>(new_x, new_y, CELL_RADII.at(new_type), new_type));

                cell->setVelocity(0, 0);
            }
        }
        
        // 5. Check for cell death and leaving
        for (auto& cell : cells) {
            CellType cellType = cell->getType();
            
            // Check for cell death (all cell types)
            if (cell->should_die()) {
                cell->is_dead = true;
                stats.total_deaths++;
                stats.deaths_by_type[cellType]++;
            }
            // Check for cell leaving (only terminal cells)
            else if (cell->should_leave()) {
                cell->is_leaving = true;
                stats.total_leaving++;
                stats.leaving_by_type[cellType]++;
            }
        }
        
        // Remove dead and leaving cells from the simulation
        auto isCellInactive = [](const unique_ptr<Cell>& cell) {
            return cell->is_dead || cell->is_leaving;
        };
        auto newEnd = remove_if(cells.begin(), cells.end(), isCellInactive);
        cells.erase(newEnd, cells.end());
        
        // Add new cells
        for (auto& new_cell : new_cells) {
            cells.push_back(std::move(new_cell));
        }
    }

    // Write cell data to file
    void writeCellDataToFile(int step_num)
    {        
        // Write data for each cell
        for (const auto& cell : cells)
        {
            bool status = cell->is_dead || cell->is_leaving;
            
            // Write to consolidated file
            consolidatedDataFile << step_num << ","
                     << cell->getType() << ","
                     << fixed << setprecision(4) << cell->x << ","
                     << fixed << setprecision(4) << cell->y << ","
                     << fixed << setprecision(4) << cell->dx << ","
                     << fixed << setprecision(4) << cell->dy << ","
                     << cell->in_vessel_neighbourhood << ","
                     << status << endl;
        }
        
    }

    // Print detailed statistics
    void printStatistics() {
        cout << "\n===== Simulation Statistics =====\n";
        cout << "Total cells remaining: " << cells.size() << endl;
        cout << "Total deaths: " << stats.total_deaths << endl;
        cout << "Total cells that left: " << stats.total_leaving << endl;
        cout << "Total blood vessels: " << blood_vessels.size() << endl;
        
        cout << "\nDeaths by cell type:\n";
        for (const auto& [type, count] : stats.deaths_by_type) {
            cout << "  " << type << ": " << count << endl;
        }
        
        cout << "\nLeaving by cell type:\n";
        for (const auto& [type, count] : stats.leaving_by_type) {
            cout << "  " << type << ": " << count << endl;
        }
        
        cout << "\nRemaining cells by type:\n";
        map<CellType, int> remaining_by_type;
        for (const auto& cell : cells) {
            remaining_by_type[cell->getType()]++;
        }
        for (const auto& [type, count] : remaining_by_type) {
            cout << "  " << type << ": " << count << endl;
        }
        cout << "==============================\n";
    }

    void run(int steps)
    {
        for (int current_step = 0; current_step < steps; ++current_step)
        {
            // Perform simulation step
            step();
            if (current_step % 50 == 0) {
                cout << "Step " << current_step << ": " << cells.size() << " cells";
                cout << " | Deaths: " << stats.total_deaths << " | Leaving: " << stats.total_leaving << endl;
            }

            // Write cell data to file for the current step
            writeCellDataToFile(current_step);
        }
        printStatistics();
    }
};

int main(int argc, char* argv[])
{
    // Default values
    float width = 1500.0;
    float height = 1000.0;
    int initial_cells = 20;
    int steps = 200;
    string sim_name = "bm_sim";
    int num_vessels = 10;  // Default number of blood vessels
    
    // Parse command-line arguments
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            string arg = argv[i];
            if ((arg == "--width" || arg == "-w") && i + 1 < argc) {
                width = stof(argv[++i]);
            } else if ((arg == "--height" || arg == "-h") && i + 1 < argc) {
                height = stof(argv[++i]);
            } else if ((arg == "--cells" || arg == "-c") && i + 1 < argc) {
                initial_cells = stoi(argv[++i]);
            } else if ((arg == "--steps" || arg == "-s") && i + 1 < argc) {
                steps = stoi(argv[++i]);
            } else if ((arg == "--name" || arg == "-n") && i + 1 < argc) {
                sim_name = argv[++i];
            } else if ((arg == "--vessels" || arg == "-v") && i + 1 < argc) {
                num_vessels = stoi(argv[++i]);
            } else if (arg == "--help") {
                cout << "Bone Marrow Simulation\n"
                     << "Usage: " << argv[0] << " [options]\n"
                     << "Options:\n"
                     << "  --width WIDTH      Set simulation width (default: 1500.0)\n"
                     << "  --height HEIGHT    Set simulation height (default: 1000.0)\n"
                     << "  --cells NUM        Set initial cell count (default: 40)\n"
                     << "  --steps NUM        Set simulation steps (default: 200)\n"
                     << "  --name NAME        Set simulation name (default: bm_sim)\n"
                     << "  --vessels NUM      Set number of blood vessels (default: 10)\n"
                     << "  --help             Display this help message\n";
                return 0;
            }
        }
    }
    
    cout << "Starting bone marrow simulation...\n";
    cout << "Configuration:\n";
    cout << "  Width: " << width << "\n";
    cout << "  Height: " << height << "\n";
    cout << "  Initial cells: " << initial_cells << "\n";
    cout << "  Blood vessels: " << num_vessels << "\n";
    cout << "  Steps: " << steps << "\n";
    cout << "  Simulation name: " << sim_name << "\n\n";
    
    BoneMarrow model(width, height, initial_cells, sim_name, num_vessels);
    model.run(steps);
    
    cout << "Simulation completed.\n";
    return 0;
}
