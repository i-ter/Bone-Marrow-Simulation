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
#include <cassert>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <mutex>
#include <unordered_set>
#include "cell_config.h"



// OpenMP support with fallback. fallback runs single threaded, OpenMP not needed to be linked.
#ifdef _OPENMP
    #include <omp.h>
#else  
    inline int omp_get_max_threads() { return 1; }
    inline int omp_get_thread_num() { return 0; }
#endif

// Global seed for reproducibility
unsigned int GLOBAL_SEED = 10; 

// Function to get thread-safe random generators
std::mt19937& get_thread_rng() {
    #ifdef _OPENMP
        thread_local std::mt19937 tl_gen(GLOBAL_SEED + omp_get_thread_num());
    #else
        static std::mt19937 tl_gen(GLOBAL_SEED);
    #endif
        return tl_gen;
}

// Thread-local distributions
#ifdef _OPENMP
    thread_local std::uniform_real_distribution<> tl_unif_01(0, 1);
    thread_local std::uniform_real_distribution<float> tl_angle_dist(0, 2 * M_PI);
    thread_local std::normal_distribution<float> tl_normal_dist(0, 1);
#else  
    std::uniform_real_distribution<> tl_unif_01(0, 1);
    std::uniform_real_distribution<float> tl_angle_dist(0, 2 * M_PI);
    std::normal_distribution<float> tl_normal_dist(0, 1);
#endif



using namespace std;


// Wrapper functions for easier migration
inline float get_random_01() { return tl_unif_01(get_thread_rng()); }
inline float get_random_angle() { return tl_angle_dist(get_thread_rng()); }
inline float get_random_normal() { return tl_normal_dist(get_thread_rng()); }

// Helper function for random integer in range [0, max_val)
inline int get_random_int(int max_val) {
    std::uniform_int_distribution<int> int_dist(0, max_val - 1);
    return int_dist(get_thread_rng());
}

class BloodVessel
{
public:
    // Start and end points of the blood vessel
    float start_x, start_y;
    float end_x, end_y;
    float radius; // Radius of the blood vessel
    float length; // Length of the blood vessel
    float angle;  // Angle of the blood vessel

    // Constructor for a blood vessel
    BloodVessel(float start_x, float start_y, float length, float radius, float angle)
        : start_x(start_x), start_y(start_y), radius(radius), length(length), angle(angle)
    {
        // Calculate end point based on start point, length, and angle
        end_x = start_x + length * cos(angle);
        end_y = start_y + length * sin(angle);
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
};

class Cell 
{
private:
    unsigned int id; // unique cell id
    static unsigned int next_id; // static counter for unique IDs

protected:
    CellType cell_type;

public:
    float x, y; // position
    float radius;
    float dx, dy; // velocity
    bool is_dead = false;
    bool is_leaving = false;
    bool in_vessel_neighbourhood = false; // Flag to track if cell is currently inside a blood vessel
    int clone_id=-100;
    int mass=1;
    int grid_x, grid_y;
    float swap_motility_prob=0.0;  // probability of swapping with random neighbour
    mutable std::mutex cell_mutex;

    void updateType(CellType type) {
        cell_type = type;
        radius = get_with_default(CELL_RADII, type, DEFAULT_CELL_RADII);
    }

    Cell(float x, float y, float radius, CellType type, int clone_id, int mass=1)
        : cell_type(type), x(x), y(y), radius(radius), clone_id(clone_id), mass(mass) {
        id = next_id++;
        swap_motility_prob = get_with_zero(SWAP_MOTILITY, type);
    }

    ~Cell() {}

    int getId() const { return id; }

    CellType getType() const {
        return cell_type;
    }

    void move(const float &width, const float &height, const float &multiplier=1.0f) {
        // Define diffusion coefficient
        float dt = 1.0; // Define a time step, can be adjusted
        float D = get_with_default(MOTILITY, cell_type, DEFAULT_CELL_MOTILITY) * multiplier;
        // Update position with spatial diffusion

        dx = D * sqrt(dt) * get_random_normal();
        dy = D * sqrt(dt) * get_random_normal();
        x += dx;
        y += dy;

        handleBoundaryCollision(width, height);
    }

    void handleBoundaryCollision(const float &width, const float &height) {
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
        }
        else if (y + radius > height) {
            y = height - radius;
            dy = -dy;
        }
    }

    // bool collidesWith(const Cell &other) const {
    //     float dx = x - other.x;
    //     float dy = y - other.y;
    //     float distance_squared = dx * dx + dy * dy;
    //     float min_distance = radius + other.radius;
    //     return distance_squared < min_distance * min_distance;
    // }

    // Check if cell collides with a blood vessel
    bool inVesselNeighbour(const BloodVessel &vessel, const float &distance_threshold) const {
        return vessel.distanceFrom(x, y) <= distance_threshold;
    }

    // void resolveCollision(Cell &other, const float &k, const float &dt) {
    //     // TODO: This is technically Gauss Seidel method.
    //     // Calculates the force of the collision. F_ij = k * overlap * n_ij.
    //     // k is the spring constant. assumes unit mass.

    //     // Calculate displacement vector
    //     float x_diff = x - other.x;
    //     float y_diff = y - other.y;
    //     float distance = sqrt(x_diff * x_diff + y_diff * y_diff);

    //     // Ensure cells are not exactly on top of each other. leads to numerical instability.
    //     if (distance < 0.1) {
    //         distance = 0.1;
    //     }

    //     float overlap = radius + other.radius - distance;

    //     // If cells are overlapping
    //     if (overlap < 0){
    //         return;
    //     }

    //     // unit vector in the direction of the collision
    //     float nx = x_diff / distance;
    //     float ny = y_diff / distance;

    //     // calculate the displacement of the cells (Hookean spring)
    //     float displacement = k * overlap * dt;
    //     // calculate the mass of the cells
    //     float total_mass = mass + other.mass;
    //     float displacement_this = displacement * other.mass / total_mass;
    //     float displacement_other = displacement * mass / total_mass;

    //     // separate the cells
    //     if (this->getType() != STROMA) {
    //         x += nx * displacement_this;
    //         y += ny * displacement_this;
    //     }
    //     if (other.getType() != STROMA) {
    //         other.x -= nx * displacement_other;
    //         other.y -= ny * displacement_other;
    //     }

    // }

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

    float distanceFrom(const Cell &other) const {
        float dx = x - other.x;
        float dy = y - other.y;
        return sqrt(dx * dx + dy * dy);
    }

    void setSwapMotility(float prob) {
        swap_motility_prob = prob;
    }

    bool should_divide() {
        return get_random_01() < get_with_default(DIVISION_PROB, cell_type, DEFAULT_DIVISION_PROB);
    }

    // Check if cell should die (all cell types)
    bool should_die() {
        return get_random_01() < get_with_default(CELL_DEATH_PROB, cell_type, DEFAULT_CELL_DEATH_PROB);
    }

    // Check if cell should leave (based on probability from LEAVE_PROB map)
    bool should_leave() {
        float multiplier = in_vessel_neighbourhood ? VESSEL_LEAVING_MULTIPLIER : 1.0;
        return get_random_01() < get_with_zero(LEAVE_PROB, cell_type) * multiplier;
    }
    
    // if true, cell will swap with a random neighbour. 
    bool should_swap() {
        return get_random_01() < swap_motility_prob;
    }
};

unsigned int Cell::next_id = 0;

// Spatial partitioning for collision detection optimization
struct SpatialGrid
{
    float block_size;
    int grid_width, grid_height;
    vector<vector<vector<Cell *>>> grid; // grid[x][y] returns a vector of cells in block (x,y)
    vector<vector<int>> grid_block_cell_counts;
    float grid_block_density_limit;

    SpatialGrid(int width, int height, float block_size) : block_size(block_size)
    {
        grid_width = ceil(width / block_size);
        grid_height = ceil(height / block_size);
        grid.resize(grid_width, vector<vector<Cell *>>(grid_height));
        grid_block_cell_counts.resize(grid_width, vector<int>(grid_height, 0));

        // upper bound. 
        grid_block_density_limit = (block_size * block_size) / (M_PI * 5 * 5);
    }


    void clear()
    {
        for (auto &col : grid) {
            for (auto &block : col) {
                block.clear();
            }
        }
    }

    void insert(Cell *cell)
    {
        int grid_x = min(grid_width - 1, max(0, static_cast<int>(cell->x / block_size)));
        int grid_y = min(grid_height - 1, max(0, static_cast<int>(cell->y / block_size)));
        grid[grid_x][grid_y].push_back(cell);
        cell->grid_x = grid_x;
        cell->grid_y = grid_y;
    }

    void calculateGridBlockCellCounts() {
        for (int i = 0; i < grid_width; i++) {
            for (int j = 0; j < grid_height; j++) {
                grid_block_cell_counts[i][j] = grid[i][j].size();
            }
        }
    }

    bool isGridBlockOvercrowded(const int &x, const int &y) {
        return grid_block_cell_counts[x][y] > grid_block_density_limit;
    }

    bool isGridBlockNeighbourOvercrowded(const int &x, const int &y) {
        int valid_blocks = 0;
        int overcrowded_blocks = 0;
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int nx = x + dx;
                int ny = y + dy;
                if (nx >= 0 && nx < grid_width && ny >= 0 && ny < grid_height) {
                    ++valid_blocks;
                    if (isGridBlockOvercrowded(nx, ny)) {
                        ++overcrowded_blocks;
                    }
                }
            }
        }
        return overcrowded_blocks > valid_blocks / 2;
    }
    
    int countHSPCNeighbours(const Cell *cell) {
        int count = 0;
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                int nx = cell->grid_x + dx;
                int ny = cell->grid_y + dy;
                if (nx >= 0 && nx < grid_width && ny >= 0 && ny < grid_height) {
                    for (Cell *other : grid[nx][ny]) {
                        // check if other is close by HSPC
                        if (
                            other != cell &&
                            other->getType() < Erythroblast1 &&
                            cell->distanceFrom(*other) <= CELL_SENSE_RADIUS
                        ) {
                            count++;
                        }
                    }
                }
            }
        }
        return count;
    }

    vector<Cell *> getPotentialCollisions(const Cell *cell) {
        vector<Cell *> result;
        // Check cell's grid cell and neighboring cells
        for (int dx = -1; dx <= 1; dx++)
        {
            for (int dy = -1; dy <= 1; dy++)
            {
                int nx = cell->grid_x + dx;
                int ny = cell->grid_y + dy;

                if (nx >= 0 && nx < grid_width && ny >= 0 && ny < grid_height)
                {
                    for (Cell *other : grid[nx][ny])
                    {
                        if (other != cell)
                        {
                            result.push_back(other);
                        }
                    }
                }
            }
        }

        return result;
    }
    vector<Cell *> getValidNeighbours(const Cell *cell, int search_radius) {
        vector<Cell *> result;
        auto potential_collisions = getPotentialCollisions(cell);
        for (auto &other : potential_collisions) {
            // TODO: Find a better way to handle this.
            if (other->radius > cell->radius) {
                continue;
            }

            if (cell->distanceFrom(*other) <= search_radius) {
                result.push_back(other);
            }
        }
        return result;
    }

    void remove(Cell* cell_to_remove) {
        // Use the same logic as insert to determine the grid cell, including clamping.
        int grid_x_coord = std::min(grid_width - 1, std::max(0, static_cast<int>(cell_to_remove->x / block_size)));
        int grid_y_coord = std::min(grid_height - 1, std::max(0, static_cast<int>(cell_to_remove->y / block_size)));

        auto& bucket = grid[grid_x_coord][grid_y_coord];
        
        // Efficiently remove the cell pointer from the bucket.
        // std::remove shifts matching elements to the end and returns an iterator to the new 'logical' end.
        // vector::erase then removes the elements from this new logical end to the physical end.
        auto it = std::remove(bucket.begin(), bucket.end(), cell_to_remove);
        bucket.erase(it, bucket.end());
    }
};

enum CellNumState {
    NORMAL,
    LOW,
    HIGH,
};

class BoneMarrow {
public:
    int width, height;

    vector<unique_ptr<Cell>> cells;
    vector<BloodVessel> blood_vessels; // Added blood vessels collection
    int num_vessels;                   // Number of blood vessels

    string sim_name;
    string framesDir = "frames";
    string dataDir = "data";       // Directory for cell data files
    ofstream consolidatedDataFile; // File stream for consolidated data
    ofstream vesselDataFile;       // File for blood vessel data
    ofstream paramsFile;           // File for simulation parameters
    bool cold_start;
    SpatialGrid spatial_grid;
    int stop_motility_at_step=-1;
    int immotile_hsc_clone_id=-1;
    int data_write_freq = 1;
    std::map<std::string, std::chrono::duration<double, std::milli>> latest_step_timings;
    CellNumState cell_num_state = NORMAL;

    // Statistics tracking
    struct Stats {
        int total_deaths = 0;
        int total_leaving = 0;
        map<CellType, int> deaths_by_type;
        map<CellType, int> leaving_by_type;
    } stats;

    int equilibriate_steps;
     
    // init
    BoneMarrow(
        int width,
        int height, 
        int initial_cells, 
        string sim_name = "bm",
        int num_vessels = 3, 
        bool cold_start = false, 
        int stop_motility_at_step = -1,
        int data_write_freq = 1)
        : width(width),
          height(height),
          num_vessels(num_vessels),
          sim_name(sim_name),
          cold_start(cold_start),
          spatial_grid(width, height, SPATIAL_GRID_BLOCK_SIZE),
          stop_motility_at_step(stop_motility_at_step),
          data_write_freq(data_write_freq)
    {
        if (!std::filesystem::exists(dataDir)) {
            std::filesystem::create_directory(dataDir);
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
                   << "vessel_leaving_multiplier," << VESSEL_LEAVING_MULTIPLIER << endl
                   << "random_seed," << GLOBAL_SEED << endl
                   << "cold_start," << cold_start << endl
                   << "stop_motility_at_step," << stop_motility_at_step << endl;
        paramsFile.close();

        // Save cell type specific parameters in a separate CSV
        string cellParamsFilename = dataDir + "/" + sim_name + "_cell_params.csv";
        paramsFile.open(cellParamsFilename);

        // Write header
        paramsFile << "cell_type,radius,division_probability,death_probability,leave_probability,motility,initial_number" << endl;

        // Write cell type specific parameters
        for (const auto &[type, _] : LINEAGE_TREE) {
            paramsFile << getCellTypeName(type) << ","
                       << get_with_default(CELL_RADII, type, DEFAULT_CELL_RADII) << ","
                       << get_with_default(DIVISION_PROB, type, DEFAULT_DIVISION_PROB) << ","
                       << get_with_default(CELL_DEATH_PROB, type, DEFAULT_CELL_DEATH_PROB) << ","
                       << get_with_zero(LEAVE_PROB, type) << ","
                       << get_with_default(MOTILITY, type, DEFAULT_CELL_MOTILITY) << ","
                       << get_with_zero(INITIAL_CELL_NUMBERS, type) << endl;
        }
        paramsFile.close();

        // Initialize consolidated data file
        string consolidatedFilename = dataDir + "/" + sim_name + "_all_steps.csv";
        consolidatedDataFile.open(consolidatedFilename);
        consolidatedDataFile << "step,cell_type,x,y,dx,dy,clone_id,vessel_neighbourhood,status,cell_id" << endl;

        if (cold_start) {
            cout << "--- COLD STARTING THE SIMULATION ---" << endl;
            for (int i = 0; i < initial_cells; ++i)
            {
                float x = static_cast<float>(get_random_01() * width);
                float y = static_cast<float>(get_random_01() * height);
                cells.push_back(make_unique<Cell>(x, y, get_with_default(CELL_RADII, HSC, DEFAULT_CELL_RADII), HSC, i));
            }
        } else {

            int starting_cells = 0;
            int HSPC_num = 0; // tag HSPC cells with clone_id

            // initial cell_numbers are configured for 500x500 grid

            float cell_num_multiplier = width/500.0f * height/500.0f;
            cout << "initial cell numbers multiplier: " << cell_num_multiplier << endl;

            for (auto &[type, num] : INITIAL_CELL_NUMBERS) {
                num *= cell_num_multiplier;
                for (int i = 0; i < num; ++i) {
                    float x = static_cast<float>(get_random_01() * width);
                    float y = static_cast<float>(get_random_01() * height);
                    float radius = get_with_default(CELL_RADII, type, DEFAULT_CELL_RADII);

                    int clone_id = -1;
                    if (type < 11) 
                    {
                        clone_id = HSPC_num;
                        HSPC_num++;
                    }

                    cells.push_back(make_unique<Cell>(x, y, radius, type, clone_id));
                }
                starting_cells += num;
            }
            cout << "--- WARM STARTING THE SIMULATION WITH " << starting_cells << " CELLS ---" << endl;
        }

        float totalArea = width * height;
        float nStromaCells = totalArea / 10000 * CXCL_DENSITY_PER_100_AREA;
        float stroma_radius = get_with_default(CELL_RADII, STROMA, DEFAULT_CELL_RADII);

        for (int i = 0; i < nStromaCells; ++i)
        {
            float x = static_cast<float>(get_random_01() * width);
            float y = static_cast<float>(get_random_01() * height);
            cells.push_back(make_unique<Cell>(x, y, stroma_radius, STROMA, -1, 10));
        }

        cout << "Generated " << nStromaCells << " stroma cells" << endl;

        generateBloodVessels();
    }

    void generateFixedVessels()
    {
        float start_x = width / 4;
        float start_y = height / 2;
        float length = width / 2;
        float radius = 20;
        float angle = 0;

        blood_vessels.push_back(BloodVessel(start_x, start_y, length, radius, angle));
        blood_vessels.push_back(BloodVessel(start_x, height / 4, length, radius, angle));
        blood_vessels.push_back(BloodVessel(start_x, 3 * height / 4, length, radius, angle));

        return;
    }

    // Check if a vessel is valid (within bounds and not intersecting other vessels too much)
    bool isVesselValid(const BloodVessel &vessel)
    {
        // Check if vessel endpoints are within bounds with some margin
        float margin = vessel.radius;
        if (vessel.start_x < margin || vessel.start_x > width - margin ||
            vessel.start_y < margin || vessel.start_y > height - margin ||
            vessel.end_x < margin || vessel.end_x > width - margin ||
            vessel.end_y < margin || vessel.end_y > height - margin)
        {
            return false;
        }

        return true;
    }

    // Generate all blood vessels
    void generateBloodVessels() {

        generateFixedVessels();

        // Write vessel data to file
        string vesselFilename = dataDir + "/" + sim_name + "_vessels.csv";
        vesselDataFile.open(vesselFilename);
        vesselDataFile << "id,start_x,start_y,end_x,end_y,radius" << endl;

        for (size_t i = 0; i < blood_vessels.size(); i++)
        {
            const auto &vessel = blood_vessels[i];
            vesselDataFile << i << ","
                           << fixed << setprecision(4) << vessel.start_x << ","
                           << fixed << setprecision(4) << vessel.start_y << ","
                           << fixed << setprecision(4) << vessel.end_x << ","
                           << fixed << setprecision(4) << vessel.end_y << ","
                           << fixed << setprecision(4) << vessel.radius << endl;
        }

        vesselDataFile.close();

        cout << "Generated " << blood_vessels.size() << " blood vessels" << endl;
    }

    ~BoneMarrow() {
        // destructor to close the consolidated data file
        consolidatedDataFile.close();
        vesselDataFile.close();
        paramsFile.close();
    }

    void swapCells(Cell *cell, Cell *neighbour) {
        float temp_x = cell->x;
        float temp_y = cell->y;
        cell->x = neighbour->x;
        cell->y = neighbour->y;
        neighbour->x = temp_x;
        neighbour->y = temp_y;
    }

    CellType sampleRandomType(const CellType &type) {
        const auto &possible_types = LINEAGE_TREE.at(type);
        
        if (possible_types.size() == 1) {
            return possible_types[0].first;
        }

        float r = get_random_01();
        float cumulative_prob = 0.0f;
        
        for (const auto &[type, prob]: possible_types) {
            cumulative_prob += prob;
            if (r < cumulative_prob) {
                return type;
            }
        }
        throw std::runtime_error("Failed to sample random type");
    }

    void updateSpatialGrid() {
        spatial_grid.clear();
        for (auto &cell : cells) {
            spatial_grid.insert(cell.get());
        }
        spatial_grid.calculateGridBlockCellCounts();
    }

    void moveCells() {
        #pragma omp parallel for schedule(static)
        for (size_t idx = 0; idx < cells.size(); ++idx) {
            Cell* cell = cells[idx].get();
            // swapping cells with neighbours
            if (cell->should_swap()) {
                // find a random neighbour
                vector<Cell *> neighbours = spatial_grid.getValidNeighbours(cell, 10);
                if (neighbours.size() > 0) {
                    Cell *random_neighbour = neighbours[get_random_int(neighbours.size())];

                    #pragma omp critical
                    {
                        // Remove both cells from their current grid positions BEFORE their coordinates are changed.
                        // This ensures they are removed from the correct buckets based on their pre-swap locations.
                        spatial_grid.remove(cell);
                        spatial_grid.remove(random_neighbour);

                        swapCells(cell, random_neighbour);

                        spatial_grid.insert(cell);
                        spatial_grid.insert(random_neighbour);
                    }
                }
            }
            // float local_density = spatial_grid.grid_block_cell_counts[cell->grid_x][cell->grid_y];
            // float density_factor = 1.0 + (local_density / spatial_grid.grid_block_density_limit);

            cell->move(width, height);
        }
    }

    void handleVesselCollisions() {
        #pragma omp parallel for schedule(static)
        for (size_t idx = 0; idx < cells.size(); ++idx) {
            // Cell* cell = cells[idx].get();
            cells[idx]->in_vessel_neighbourhood = false;
            for (const auto &vessel : blood_vessels){
                if (cells[idx]->inVesselNeighbour(vessel, VESSEL_DISTANCE_THRESHOLD)) {
                    cells[idx]->in_vessel_neighbourhood = true;
                    break;
                }
            }
        }   
    }

    float resolveCollisions(const float &dt) {
        // Store displacement vectors for each cell (Jacobi iteration)
        std::vector<std::pair<float, float>> displacements(cells.size(), {0.0f, 0.0f});
        
        // Create a map from cell pointer to index for O(1) lookup
        std::unordered_map<Cell*, size_t> cell_to_index;
        for (size_t i = 0; i < cells.size(); ++i) {
            cell_to_index[cells[i].get()] = i;
        }
        
        // Calculate all displacements
        #pragma omp parallel for schedule(static)
        for (size_t idx = 0; idx < cells.size(); ++idx) {
            Cell* c1 = cells[idx].get();

            auto potential_collisions = spatial_grid.getPotentialCollisions(c1);

            for (Cell *c2 : potential_collisions) {
                // Only process if c2 comes after c1 in the array to avoid double-processing
                auto c2_iter = cell_to_index.find(c2);
                if (c2_iter == cell_to_index.end() || c2_iter->second <= idx) {
                    continue;
                }
                size_t c2_idx = c2_iter->second;
                
                float x_diff = c1->x - c2->x;
                float y_diff = c1->y - c2->y;
                float distance = sqrt(x_diff * x_diff + y_diff * y_diff);

                // Ensure cells are not exactly on top of each other
                if (distance < 0.1) {
                    distance = 0.1;
                }

                float overlap = c1->radius + c2->radius - distance;
                
                // Define repulsion parameters
                const float repulsion_multiplier = 1.5f; // How much larger than contact distance for repulsion
                float repulsion_distance = (c1->radius + c2->radius) * repulsion_multiplier;
                
                // Check if cells are within repulsion range
                if (distance < repulsion_distance) {
                    float nx = x_diff / distance;
                    float ny = y_diff / distance;
                    float total_mass = c1->mass + c2->mass;
                    
                    float displacement = 0.0f;
                    
                    const float k = 0.8f;
                    displacement = k * (repulsion_distance - distance) * dt;
                    
                    float displacement_c1 = displacement * c2->mass / total_mass;
                    float displacement_c2 = displacement * c1->mass / total_mass;

                    std::scoped_lock lock(c1->cell_mutex, c2->cell_mutex);
                    
                    if (c1->getType() != STROMA) {
                        displacements[idx].first += nx * displacement_c1;
                        displacements[idx].second += ny * displacement_c1;
                    }
                    if (c2->getType() != STROMA) {
                        displacements[c2_idx].first -= nx * displacement_c2;
                        displacements[c2_idx].second -= ny * displacement_c2;
                    }
                }
            }
        }
        
        double total_displacement_x = 0.0;
        double total_displacement_y = 0.0;
        
        // Apply all displacements simultaneously (Jacobi step)
        // reduction creates private copies of the variables
        #pragma omp parallel for reduction(+ : total_displacement_x, total_displacement_y)
        for (size_t idx = 0; idx < cells.size(); ++idx) {
            Cell* cell = cells[idx].get();
            cell->x += displacements[idx].first;
            cell->y += displacements[idx].second;
            cell->handleBoundaryCollision(width, height);
            
            total_displacement_x += std::abs(displacements[idx].first);
            total_displacement_y += std::abs(displacements[idx].second);
        }
        return total_displacement_x;
    }

    void handleCellDivision(std::vector<std::unique_ptr<Cell>> &new_cells) {
        #pragma omp parallel for schedule(static)
        for (size_t idx = 0; idx < cells.size(); ++idx) {
            auto &cell = cells[idx];

            if (cell->getType() == HSC) {
                if (spatial_grid.countHSPCNeighbours(cell.get()) > MAX_HSPC_NEIGHBOURS) {
                    continue;
                }
            }

            if (cell->should_divide() && !spatial_grid.isGridBlockNeighbourOvercrowded(cell->grid_x, cell->grid_y)) 
            {                
                CellType new_type = sampleRandomType(cell->getType());
                
                float new_radius = get_with_default(CELL_RADII, new_type, DEFAULT_CELL_RADII);
                float offset = cell->radius + new_radius;
                float angle = get_random_angle();

                float new_x = cell->x + offset * cos(angle);
                float new_y = cell->y + offset * sin(angle);

                auto new_cell = std::make_unique<Cell>(new_x, new_y, new_radius, new_type, cell->clone_id);
                if (new_cell->clone_id == immotile_hsc_clone_id) {
                    new_cell->setSwapMotility(0.0);
                }

                #pragma omp critical 
                {
                    new_cells.push_back(std::move(new_cell));
                }

                if (cell->getType() != HSC) {
                    // two daughter cells, each w own fate decision. modify original cell.
                    cell->updateType(sampleRandomType(cell->getType()));
                }
            }
        }
    }
    
    void removeCells() {
        auto isCellInactive = [](const std::unique_ptr<Cell> &cell) {
            return cell->is_dead || cell->is_leaving;
        };
        auto newEnd = std::remove_if(cells.begin(), cells.end(), isCellInactive);
        cells.erase(newEnd, cells.end());
    }


    void handleCellDeathAndLeaving() {
        int num_threads = omp_get_max_threads();
        std::vector<BoneMarrow::Stats> thread_stats(num_threads);

        #pragma omp parallel for schedule(static)
        for (size_t idx = 0; idx < cells.size(); ++idx) {
            int tid = omp_get_thread_num();
            Cell* cell = cells[idx].get();

            if (cell->should_die()) {
                cell->is_dead = true;
                thread_stats[tid].total_deaths++;
                thread_stats[tid].deaths_by_type[cell->getType()]++;
            } else if (cell->should_leave()) {
                cell->is_leaving = true;
                thread_stats[tid].total_leaving++;
                thread_stats[tid].leaving_by_type[cell->getType()]++;
            }
        }

        // After the parallel region, merge thread_stats into stats
        for (const auto& tstat : thread_stats) {
            stats.total_deaths += tstat.total_deaths;
            stats.total_leaving += tstat.total_leaving;
            for (const auto& [type, count] : tstat.deaths_by_type) {
                stats.deaths_by_type[type] += count;
            }
            for (const auto& [type, count] : tstat.leaving_by_type) {
                stats.leaving_by_type[type] += count;
            }
        }
    }

    void step() {
        auto step_total_start = std::chrono::high_resolution_clock::now();
        std::map<std::string, std::chrono::duration<double, std::milli>> timings;


        // 0. setup spatial grid for efficient handling of cell neighbourhood
        auto grid_setup_start = std::chrono::high_resolution_clock::now();
        updateSpatialGrid();
        auto grid_setup_end = std::chrono::high_resolution_clock::now();
        timings["grid_setup"] = grid_setup_end - grid_setup_start;

        // 1. Move cells
        auto move_cells_start = std::chrono::high_resolution_clock::now();
        moveCells();
        auto move_cells_end = std::chrono::high_resolution_clock::now();
        timings["move_cells"] = move_cells_end - move_cells_start;

        // 2. Handle blood vessel collisions
        auto vessel_collision_start = std::chrono::high_resolution_clock::now();
        handleVesselCollisions();
        auto vessel_collision_end = std::chrono::high_resolution_clock::now();
        timings["vessel_collision"] = vessel_collision_end - vessel_collision_start;

        // 3. Check for collisions and resolve them
        auto resolve_collision_start = std::chrono::high_resolution_clock::now();

        float small_time_step = 0.1; // 1 min
        int num_inner_steps = 100; // big loop is 10 mins

        float disp;
        int i;
        for (i = 0; i < num_inner_steps; ++i) {
            disp = resolveCollisions(small_time_step);
            if ( i > 5 && disp < 50) {
                // cout << "Stopping at step " << i << endl;
                break;
            }
        }
        equilibriate_steps = i;
        // cout << "Equilibriated in " << equilibriate_steps << " steps, with disp " << disp << endl;

        auto resolve_collision_end = std::chrono::high_resolution_clock::now();
        timings["resolve_collision"] = resolve_collision_end - resolve_collision_start;

        // 4. Cell division
        auto cell_division_start = std::chrono::high_resolution_clock::now();
        std::vector<std::unique_ptr<Cell>> new_cells;

        handleCellDivision(new_cells);

        auto cell_division_end = std::chrono::high_resolution_clock::now();
        timings["cell_division"] = cell_division_end - cell_division_start;

        // 5. Check for cell death and leaving
        auto death_leaving_start = std::chrono::high_resolution_clock::now();
        handleCellDeathAndLeaving();

        auto death_leaving_end = std::chrono::high_resolution_clock::now();
        timings["death_leaving"] = death_leaving_end - death_leaving_start;

        // Remove dead and leaving cells from the simulation
        auto remove_inactive_start = std::chrono::high_resolution_clock::now();
        removeCells();
        auto remove_inactive_end = std::chrono::high_resolution_clock::now();
        timings["remove_inactive"] = remove_inactive_end - remove_inactive_start;

        // Add new cells
        auto add_new_cells_start = std::chrono::high_resolution_clock::now();
        cells.insert(cells.end(), 
                    std::make_move_iterator(new_cells.begin()), // move iterators to instead of copying
                    std::make_move_iterator(new_cells.end()));
        auto add_new_cells_end = std::chrono::high_resolution_clock::now();
        timings["add_new_cells"] = add_new_cells_end - add_new_cells_start;

        auto step_total_end = std::chrono::high_resolution_clock::now();
        timings["step_total"] = step_total_end - step_total_start;
        latest_step_timings = timings;
    }

    // Write cell data to file
    void writeCellDataToFile(int step_num) {
        // Write data for each cell
        for (const auto &cell : cells) {
            bool status = cell->is_dead || cell->is_leaving;

            // Write to consolidated file
            consolidatedDataFile << step_num << ","
                                 << cell->getType() << ","
                                 << fixed << setprecision(4) << cell->x << ","
                                 << fixed << setprecision(4) << cell->y << ","
                                 << fixed << setprecision(4) << cell->dx << ","
                                 << fixed << setprecision(4) << cell->dy << ","
                                 << cell->clone_id << ","
                                 << cell->in_vessel_neighbourhood << ","
                                 << status << ","
                                 << cell->getId() << endl;
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
        for (const auto &[type, count] : stats.deaths_by_type) {
            cout << "  " << getCellTypeName(type) << ": " << count << endl;
        }

        cout << "\nLeaving by cell type:\n";
        for (const auto &[type, count] : stats.leaving_by_type) {
            cout << "  " << getCellTypeName(type) << ": " << count << endl;
        }

        cout << "\nRemaining cells by type:\n";
        map<CellType, int> remaining_by_type;
        for (const auto &cell : cells) {
            remaining_by_type[cell->getType()]++;
        }
        for (const auto &[type, count] : remaining_by_type) {
            cout << "  " << getCellTypeName(type) << ": " << count << endl;
        }
        cout << "==============================\n";
    }

    // adjust division probability based on cell number
    void adjustDivisionProb() {
        float max_cells = MAX_CELLS*width/500*height/500;
        
        if (cells.size() > max_cells && cell_num_state == NORMAL) {
            cell_num_state = HIGH;
            DEFAULT_DIVISION_PROB /= 2;
            for (auto &[type, prob] : DIVISION_PROB) {
                prob /= 2;
            }
            DEFAULT_CELL_DEATH_PROB *= 2;
            for (auto &[type, prob] : CELL_DEATH_PROB) {
                prob *= 2;
            }
            for (auto &[type, prob] : LEAVE_PROB) {
                prob *= 2;
            }
        }
        else if (cells.size() < max_cells && cell_num_state == HIGH) {
            cell_num_state = NORMAL;
            DEFAULT_DIVISION_PROB *= 2;
            for (auto& [type, prob] : DIVISION_PROB) {
                prob *= 2;
            }
            DEFAULT_CELL_DEATH_PROB /= 2;
            for (auto& [type, prob] : CELL_DEATH_PROB) {
                prob /= 2;
            }   
            for (auto &[type, prob] : LEAVE_PROB) {
                prob /= 2;
            }
        }
        else if (cells.size() < max_cells*MAX_CELLS_LOW_MULTIPLIER && cell_num_state == NORMAL) {
            cell_num_state = LOW;
            DEFAULT_DIVISION_PROB *= 2;
            for (auto& [type, prob] : DIVISION_PROB) {
                prob *= 2;
            }
        }
        else if (cells.size() > max_cells*MAX_CELLS_LOW_MULTIPLIER && cell_num_state == LOW) {
            cell_num_state = NORMAL;
            DEFAULT_DIVISION_PROB /= 2;
            for (auto& [type, prob] : DIVISION_PROB) {
                prob /= 2;
            }
        }
        if (DEFAULT_DIVISION_PROB > 0.5f){
             throw std::runtime_error("Division probability is too high. Somehting wrong with the logic.");
        }
    }

    // stop motility mechanism
    void stopMotility(const float &step_num) {
        if (step_num == stop_motility_at_step) {
            cout << "Stopped motility at step " << step_num << endl;

            std::vector<int> hsc_clone_ids; 
            
            for (auto &cell : cells) {
                if (cell->getType() == HSC) {
                    hsc_clone_ids.push_back(cell->clone_id);
                }
            }
            
            if (!hsc_clone_ids.empty()) {
                // Find the HSC with the lowest clone_id
                std::sort(hsc_clone_ids.begin(), hsc_clone_ids.end());
                immotile_hsc_clone_id = hsc_clone_ids[0];
                
                cout << "HSC with clone_id " << immotile_hsc_clone_id << " will lose motility" << endl;
                
                // Set motility to 0 for HSC
                for (auto &cell : cells) {
                    if (cell->getType() == HSC && cell->clone_id == immotile_hsc_clone_id) {
                        cell->setSwapMotility(0.0);
                        cout << "It's position: " << cell->x << ", " << cell->y << endl;
                    }
                }
            }
        }
    }
    
    void run(int steps) {
        auto start = std::chrono::system_clock::now();
        auto step_start_time = start;

        writeCellDataToFile(0);
        
        for (int current_step = 1; current_step < steps+1; ++current_step) {
            
            stopMotility(current_step);
            
            step();

            if (current_step % 1000 == 0) {
                adjustDivisionProb();
            }

            if (current_step % 100 == 0) {
                auto time_now = std::chrono::system_clock::now();   
                std::chrono::seconds time_elapsed_total = std::chrono::duration_cast<std::chrono::seconds>(time_now - step_start_time);
                int minutes = time_elapsed_total.count() / 60;
                int seconds = time_elapsed_total.count() % 60;

                cout << "Step " << current_step << ": " << cells.size() << " cells";
                cout << " | Deaths: " << stats.total_deaths << " | Leaving: " << stats.total_leaving;
                cout << " | Iter time: " << minutes << "m " << seconds << "s";
                cout << " | Equil in " << equilibriate_steps << " steps";
                cout << endl;
                
                if (false) {
                    if (!latest_step_timings.empty()) {
                        cout << "Timings (ms):" << endl;
                        for(const auto& pair : latest_step_timings) {
                            cout << pair.first << ": " << fixed << setprecision(3) << pair.second.count() << endl;
                        }
                    }
                    cout << "\n\n" << endl;
                }
                step_start_time = std::chrono::system_clock::now();
            }

            // Write cell data to file for the current step
            if (current_step % data_write_freq == 0){
                writeCellDataToFile(current_step);
            }
        }

        auto end = std::chrono::system_clock::now();
        std::chrono::seconds duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        int minutes = duration.count() / 60;
        int seconds = duration.count() % 60;

        printStatistics();

        cout << "\nTotal simulation time: " << minutes << " minutes " << seconds << " seconds\n" << endl;

    }
    
    void runTest() {
        cells.clear();
        
        float radius = get_with_default(CELL_RADII, Myeloblast1, DEFAULT_CELL_RADII);

        cells.push_back(std::make_unique<Cell>(20, 20, radius, Myeloblast1, 0));
        // cells.push_back(std::make_unique<Cell>(20, 30, radius, Myeloblast1, 0));

        cells.push_back(std::make_unique<Cell>(21, 21, radius, Myeloblast1, 0));
        cells.push_back(std::make_unique<Cell>(19, 21, radius, Myeloblast1, 0));
        cells.push_back(std::make_unique<Cell>(19, 19, radius, Myeloblast1, 0));
        cells.push_back(std::make_unique<Cell>(21, 19, radius, Myeloblast1, 0));

        cells.push_back(std::make_unique<Cell>(22, 18, radius, Myeloblast1, 0));
        cells.push_back(std::make_unique<Cell>(22, 22, radius, Myeloblast1, 0));
        cells.push_back(std::make_unique<Cell>(18, 22, radius, Myeloblast1, 0));
        cells.push_back(std::make_unique<Cell>(18, 18, radius, Myeloblast1, 0));


        writeCellDataToFile(0);

        for (int i = 0; i < 100; i++) {
            updateSpatialGrid();
            resolveCollisions(0.1);
            writeCellDataToFile(i+1);
        }
        
    }
};

int main(int argc, char *argv[])
{
    // Default values
    int width = 500;
    int height = 500;
    int initial_cells = 2;
    int steps = 100;
    string sim_name = "dbg";
    int num_vessels = 10; // Default number of blood vessels
    bool cold_start = true;
    int stop_motility_at_step = -999;
    int data_write_freq = 1;
    bool smc = false;

    // Parse command-line arguments
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            string arg = argv[i];
            if ((arg == "--width" || arg == "-w") && i + 1 < argc) {
                width = stof(argv[++i]);
            }
            else if ((arg == "--height" || arg == "-h") && i + 1 < argc) {
                height = stof(argv[++i]);
            }
            else if ((arg == "--cells" || arg == "-c") && i + 1 < argc) {
                initial_cells = stoi(argv[++i]);
            }
            else if ((arg == "--steps" || arg == "-s") && i + 1 < argc) {
                steps = stoi(argv[++i]);
            }
            else if ((arg == "--name" || arg == "-n") && i + 1 < argc) {
                sim_name = argv[++i];
            }
            else if ((arg == "--vessels" || arg == "-v") && i + 1 < argc) {
                num_vessels = stoi(argv[++i]);
            }
            else if ((arg == "--seed" || arg == "-seed") && i + 1 < argc) {
                GLOBAL_SEED = stoi(argv[++i]);
            }
            else if ((arg == "--cold_start" || arg == "-cs") && i + 1 < argc) {
                string cold_start_str = argv[++i];
                // transform(cold_start_str.begin(), cold_start_str.end(), cold_start_str.begin(), ::tolower);
                cold_start = (cold_start_str == "true" || cold_start_str == "1" || cold_start_str == "yes" || cold_start_str == "y");
            }
            else if ((arg == "--stop_motility" || arg == "-sm") && i + 1 < argc) {
                stop_motility_at_step = stoi(argv[++i]);
            }
            else if ((arg == "--stop_motility_completely" || arg == "-smc") && i + 1 < argc) {
                string smc_str = argv[++i];
                smc = (smc_str == "true" || smc_str == "1" || smc_str == "yes" || smc_str == "y");
            }
            else if ((arg == "--data_write_freq" || arg == "-dwf") && i + 1 < argc) {
                data_write_freq = stoi(argv[++i]);
            }
            else if (arg == "--help") {
                cout << "Bone Marrow Simulation\n"
                     << "Usage: " << argv[0] << " [options]\n"
                     << "Options:\n"
                     << "  --width WIDTH      Set simulation width (default: 1500.0)\n"
                     << "  --height HEIGHT    Set simulation height (default: 1000.0)\n"
                     << "  --cells NUM        Set initial cell count (default: 40)\n"
                     << "  --steps NUM        Set simulation steps (default: 200)\n"
                     << "  --name NAME        Set simulation name (default: bm_sim)\n"
                     << "  --vessels NUM      Set number of blood vessels (default: 10)\n"
                     << "  --seed SEED        Set random seed for reproducibility (default: 42)\n"
                     << "  --help             Display this help message\n"
                     << "  --cold_start       Cold start the simulation (default: false)\n"
                     << "  --stop_motility   Stop motility at step (default: -1)\n"
                     << "  --data_write_freq  Write data to file every n steps (default: 1)\n";
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
    cout << "  Simulation name: " << sim_name << "\n";
    cout << "  Random seed: " << GLOBAL_SEED << "\n";
    cout << "  Cold start: " << cold_start << "\n";
    cout << "  Stop motility: " << stop_motility_at_step << "\n";
    cout << "  Stop motility completely: " << smc << "\n";
    #ifdef _OPENMP
        cout << "OpenMP is supported" << " | " << omp_get_max_threads() << " threads" << endl;
    #else
        cout << "OpenMP is not supported" << endl;
    #endif

    if (smc) {
        for (auto& [type, prob] : SWAP_MOTILITY) {
            prob = 0.0;
        }
    }

    BoneMarrow model(width, height, initial_cells, sim_name, num_vessels, cold_start, stop_motility_at_step, data_write_freq);
    model.run(steps);
    // model.runTest();

    cout << "Simulation completed.\n";
    return 0;
}
