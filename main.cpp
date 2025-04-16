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


// Base Cell Class
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
    
    Cell(float x, float y, float radius, CellType type) : cell_type(type), x(x), y(y), radius(radius) {
        // Initialize random movement direction
        float angle = angle_dist(gen);
        float speed = (0.5 + unif_01(gen) * 0.5) * MOTILITY.at(cell_type); // Adjust speed by motility
        dx = cos(angle) * speed;
        dy = sin(angle) * speed;
    }

    virtual ~Cell() {}

    const CellType& getType() const {
        return cell_type;
    }

    void move(const float&  width, const float& height) {
        x += dx;
        y += dy;
        
        handleBoundaryCollision(width, height);

        // // Occasionally change direction slightly
        // if (unif_01(gen) < 0.5) {
        //     float angle = atan2(dy, dx);
        //     angle += angle_change(gen);
        //     float speed = sqrt(dx*dx + dy*dy);
        //     dx = cos(angle) * speed;
        //     dy = sin(angle) * speed;
        // }
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

    virtual bool should_divide() {
        return unif_01(gen) < DIVISION_PROB.at(cell_type);
    }

    // Check if cell should die (all cell types)
    virtual bool should_die() {
        return unif_01(gen) < CELL_DEATH_PROB.at(cell_type); // 1% chance to die for all cells
    }
    
    // Check if cell should leave (based on probability from LEAVE_PROB map)
    virtual bool should_leave() {
        return unif_01(gen) < LEAVE_PROB.at(cell_type);
    }
};

class BoneMarrow
{
public:
    float width, height;

    string sim_name;
    vector<unique_ptr<Cell>> cells;
    string framesDir = "frames";
    string dataDir = "data";  // Directory for cell data files
    ofstream consolidatedDataFile;  // File stream for consolidated data
    
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

    // init
    BoneMarrow(float width, float height, int initial_cells, string sim_name="bm") 
        : width(width), height(height), spatial_grid(width, height, 12.) // Grid block size of 20 units
    {
        if (!fs::exists(dataDir))
        {
            fs::create_directory(dataDir);
        }
        
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
        
        // 2. Handle collisions using spatial partitioning O(n log n)
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
        
        // 3. Cell division
        for (auto& cell : cells) {
            if (cell->should_divide()) {
                if (cell->is_dead && cell->is_leaving ) {throw runtime_error("Cell is dead or leaving. something is wrong.");}
                CellType new_type = LINEAGE_TREE.at(cell->getType())[rand() % LINEAGE_TREE.at(cell->getType()).size()];
                // Create new cell at a slightly offset position. overlap handling will deal w this.
                float offset = cell->radius;
                float new_x = cell->x + (unif_01(gen) * 2. - 1.) * offset;
                float new_y = cell->y + (unif_01(gen) * 2. - 1.) * offset;
                
                new_cells.push_back(make_unique<Cell>(new_x, new_y, CELL_RADII.at(new_type), new_type));
            }
        }
        
        // 4. Check for cell death and leaving
        for (auto& cell : cells) {
            if (!cell->is_dead && !cell->is_leaving) {
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
                     << status << endl;
        }
        
    }

    // Print detailed statistics
    void printStatistics() {
        cout << "\n===== Simulation Statistics =====\n";
        cout << "Total cells remaining: " << cells.size() << endl;
        cout << "Total deaths: " << stats.total_deaths << endl;
        cout << "Total cells that left: " << stats.total_leaving << endl;
        
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
            } else if (arg == "--help") {
                cout << "Bone Marrow Simulation\n"
                     << "Usage: " << argv[0] << " [options]\n"
                     << "Options:\n"
                     << "  --width WIDTH      Set simulation width (default: 1500.0)\n"
                     << "  --height HEIGHT    Set simulation height (default: 1000.0)\n"
                     << "  --cells NUM        Set initial cell count (default: 40)\n"
                     << "  --steps NUM        Set simulation steps (default: 200)\n"
                     << "  --name NAME        Set simulation name (default: bm_sim)\n"
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
    cout << "  Steps: " << steps << "\n";
    cout << "  Simulation name: " << sim_name << "\n\n";
    
    BoneMarrow model(width, height, initial_cells, sim_name);
    model.run(steps);
    
    cout << "Simulation completed.\n";
    return 0;
}
