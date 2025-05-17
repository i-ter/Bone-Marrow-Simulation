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
#include <fstream>       // Add for file operations
#include <algorithm>     // For std::remove, std::min, std::max
#include "cell_config.h" // Import the cell configuration

using namespace std;
namespace fs = std::filesystem;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> unif_01(0, 1);
std::uniform_real_distribution<float> angle_dist(0, 2 * M_PI);
// std::uniform_real_distribution<float> angle_change(-0.5, 0.5);
std::normal_distribution<float> normal_dist(0, 1);
// uniform_real_distribution<float> vessel_length_dist(100, 600);  // Blood vessel length distribution
// uniform_real_distribution<float> vessel_radius_dist(20.0, 50);    // Blood vessel radius distribution

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
};

class Cell 
{
protected:
    CellType cell_type;

public:
    float x, y; // position
    float radius;
    float dx, dy; // velocity
    bool is_dead = false;
    bool is_leaving = false;
    bool in_vessel_neighbourhood = false; // Flag to track if cell is currently inside a blood vessel
    int clone_id;

    void updateType(CellType type) {
        cell_type = type;
        radius = get_with_default(CELL_RADII, type, DEFAULT_CELL_RADII);
    }

    Cell(float x, float y, float radius, CellType type, int clone_id)
        : cell_type(type), x(x), y(y), radius(radius), clone_id(clone_id) {}

    ~Cell() {}

    const CellType &getType() const {
        return cell_type;
    }

    void move(const float &width, const float &height) {
        // Define diffusion coefficient
        float dt = 1.0; // Define a time step, can be adjusted
        float D = get_with_default(MOTILITY, cell_type, DEFAULT_CELL_MOTILITY);
        // Update position with spatial diffusion

        dx = D * sqrt(dt) * normal_dist(gen);
        dy = D * sqrt(dt) * normal_dist(gen);
        x += dx;
        y += dy;

        handleBoundaryCollision(width, height);
    }

    void handleBoundaryCollision(const float &width, const float &height) {
        if (x - radius < 0)
        {
            x = radius;
            dx = -dx;
        }
        else if (x + radius > width)
        {
            x = width - radius;
            dx = -dx;
        }
        if (y - radius < 0)
        {
            y = radius;
            dy = -dy;
        }
        else if (y + radius > height)
        {
            y = height - radius;
            dy = -dy;
        }
    }

    bool collidesWith(const Cell &other) const {
        float dx = x - other.x;
        float dy = y - other.y;
        float distance_squared = dx * dx + dy * dy;
        float min_distance = radius + other.radius;
        return distance_squared < min_distance * min_distance;
    }

    // Check if cell collides with a blood vessel
    bool inVesselNeighbour(const BloodVessel &vessel, const float &distance_threshold) const {
        return vessel.distanceFrom(x, y) <= distance_threshold;
    }

    void resolveCollision(Cell &other, const float &k, const float &dt) {
        // Calculates the force of the collision. F_ij = k * overlap * n_ij.
        // k is the spring constant. assumes unit mass.

        // Calculate displacement vector
        float x_diff = x - other.x;
        float y_diff = y - other.y;
        float distance = sqrt(x_diff * x_diff + y_diff * y_diff);

        // Ensure cells are not exactly on top of each other. leads to numerical instability.
        if (distance < 0.1) {
            distance = 0.1;
        }

        float overlap = radius + other.radius - distance;

        // If cells are overlapping
        if (overlap < 0){
            return;
        }

        // unit vector in the direction of the collision
        float nx = x_diff / distance;
        float ny = y_diff / distance;

        // calculate the displacement of the cells (Hookean spring)
        // TODO: consider damping this
        float displacement = k * overlap / 2  * dt;  // equal split of the overlap scaled by dt and k

        // separate the cells
        x += nx * displacement;
        y += ny * displacement;
        other.x -= nx * displacement;
        other.y -= ny * displacement;
    }

    // void exchangeMomentum(Cell &other, float nx, float ny) {
    //     float dot1 = this->dx * nx + this->dy * ny;
    //     float dot2 = other.dx * nx + other.dy * ny;

    //     // Update velocities
    //     this->dx += (dot2 - dot1) * nx;
    //     this->dy += (dot2 - dot1) * ny;
    //     other.dx += (dot1 - dot2) * nx;
    //     other.dy += (dot1 - dot2) * ny;
    // }

    void capVelocities()
    {
        float speed = sqrt(dx * dx + dy * dy);
        if (speed > MAX_SPEED)
        {
            dx = (dx / speed) * MAX_SPEED;
            dy = (dy / speed) * MAX_SPEED;
        }
    }

    void setVelocity(float dx, float dy)
    {
        this->dx = dx;
        this->dy = dy;
    }

    float distanceFrom(const Cell &other) const {
        float dx = x - other.x;
        float dy = y - other.y;
        return sqrt(dx * dx + dy * dy);
    }

    bool should_divide()
    {
        return unif_01(gen) < get_with_default(DIVISION_PROB, cell_type, DEFAULT_DIVISION_PROB);
    }

    // Check if cell should die (all cell types)
    bool should_die()
    {
        return unif_01(gen) < get_with_zero(CELL_DEATH_PROB, cell_type);
    }

    // Check if cell should leave (based on probability from LEAVE_PROB map)
    bool should_leave()
    {
        float multiplier = in_vessel_neighbourhood ? VESSEL_LEAVING_MULTIPLIER : 1.0;
        return unif_01(gen) < get_with_zero(LEAVE_PROB, cell_type) * multiplier;
    }
    
    // if true, cell will swap with a random neighbour. 
    bool should_swap() {
        
        return unif_01(gen) < get_with_default(MOTILITY, cell_type, DEFAULT_CELL_MOTILITY);
    }
};

// Spatial partitioning for collision detection optimization
struct SpatialGrid
{
    float block_size;
    int grid_width, grid_height;
    vector<vector<vector<Cell *>>> grid; // grid[x][y] returns a vector of cells in block (x,y)

    SpatialGrid(float width, float height, float block_size) : block_size(block_size)
    {
        grid_width = ceil(width / block_size);
        grid_height = ceil(height / block_size);
        grid.resize(grid_width, vector<vector<Cell *>>(grid_height));
    }

    void clear()
    {
        for (auto &col : grid)
        {
            for (auto &block : col)
            {
                block.clear();
            }
        }
    }

    void insert(Cell *cell)
    {
        int grid_x = min(grid_width - 1, max(0, static_cast<int>(cell->x / block_size)));
        int grid_y = min(grid_height - 1, max(0, static_cast<int>(cell->y / block_size)));
        grid[grid_x][grid_y].push_back(cell);
    }

    vector<Cell *> getPotentialCollisions(const Cell *cell) {
        vector<Cell *> result;
        int grid_x = static_cast<int>(cell->x / block_size);
        int grid_y = static_cast<int>(cell->y / block_size);

        // Check cell's grid cell and neighboring cells
        for (int dx = -1; dx <= 1; dx++)
        {
            for (int dy = -1; dy <= 1; dy++)
            {
                int nx = grid_x + dx;
                int ny = grid_y + dy;

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
    vector<Cell *> getValidNeighbours(const Cell *cell, int radius) {
        vector<Cell *> result;
        auto potential_collisions = getPotentialCollisions(cell);
        for (auto &other : potential_collisions) {
            // TODO: Find a better way to handle this.
            if (other->radius > cell->radius) {
                continue;
            }

            if (cell->distanceFrom(*other) <= radius) {
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

class BoneMarrow {
public:
    float width, height;

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

    // Statistics tracking
    struct Stats {
        int total_deaths = 0;
        int total_leaving = 0;
        map<CellType, int> deaths_by_type;
        map<CellType, int> leaving_by_type;
    } stats;
     
    // init
    BoneMarrow(float width, float height, int initial_cells, string sim_name = "bm", int num_vessels = 3, bool cold_start = false)
        : width(width),
          height(height),
          num_vessels(num_vessels),
          sim_name(sim_name),
          cold_start(cold_start),
          spatial_grid(width, height, 8.0)
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
        paramsFile << "cell_type,radius,division_probability,death_probability,leave_probability,motility,initial_number" << endl;

        // Write cell type specific parameters
        for (const auto &[type, _] : LINEAGE_TREE)
        {
            paramsFile << getCellTypeName(type) << ","
                       << get_with_default(CELL_RADII, type, DEFAULT_CELL_RADII) << ","
                       << get_with_default(DIVISION_PROB, type, DEFAULT_DIVISION_PROB) << ","
                       << get_with_zero(CELL_DEATH_PROB, type) << ","
                       << get_with_zero(LEAVE_PROB, type) << ","
                       << get_with_default(MOTILITY, type, DEFAULT_CELL_MOTILITY) << ","
                       << get_with_zero(INITIAL_CELL_NUMBERS, type) << endl;
        }
        paramsFile.close();

        // Initialize consolidated data file
        string consolidatedFilename = dataDir + "/" + sim_name + "_all_steps.csv";
        consolidatedDataFile.open(consolidatedFilename);
        consolidatedDataFile << "step,cell_type,x,y,dx,dy,clone_id,vessel_neighbourhood,status" << endl;

        if (cold_start)
        {
            cout << "--- COLD STARTING THE SIMULATION ---" << endl;
            for (int i = 0; i < initial_cells; ++i)
            {
                float x = static_cast<float>(unif_01(gen) * width);
                float y = static_cast<float>(unif_01(gen) * height);
                cells.push_back(make_unique<Cell>(x, y, get_with_default(CELL_RADII, HSC, DEFAULT_CELL_RADII), HSC, i));
            }
        }
        else
        {

            int starting_cells = 0;
            int HSPC_num = 0; // tag HSPC cells with clone_id

            for (auto &[type, num] : INITIAL_CELL_NUMBERS)
            {
                num /= 3;
                for (int i = 0; i < num; ++i)
                {
                    float x = static_cast<float>(unif_01(gen) * width);
                    float y = static_cast<float>(unif_01(gen) * height);
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
            float x = static_cast<float>(unif_01(gen) * width);
            float y = static_cast<float>(unif_01(gen) * height);
            cells.push_back(make_unique<Cell>(x, y, stroma_radius, STROMA, -1));
        }

        cout << "Generated " << nStromaCells << " stroma cells" << endl;

        generateBloodVessels(num_vessels);
    }

    // Generate a random blood vessel
    BloodVessel generateRandomVessel()
    {
        float start_x = unif_01(gen) * width;
        float start_y = unif_01(gen) * height;
        float length = 100 + unif_01(gen) * 500;
        float radius = 20 + unif_01(gen) * 30;
        float angle = angle_dist(gen);

        return BloodVessel(start_x, start_y, length, radius, angle);
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
    void generateBloodVessels(int num_vessels) {

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

    void step() {
        vector<unique_ptr<Cell>> new_cells;

        // 0. setup spatial grid for efficient handling of cell neighbourhood
        spatial_grid.clear();
        for (auto &cell : cells) {
            spatial_grid.insert(cell.get());
        }

        // 1. Move cells
        for (auto &cell : cells) {
            // diffusion movement handled by cell
            // cell->move(width, height);

            // swapping cells with neighbours
            if (cell->should_swap()) {
                // find a random neighbour
                vector<Cell *> neighbours = spatial_grid.getValidNeighbours(cell.get(), 10);
                if (neighbours.size() > 0) {
                    Cell *current_cell_ptr = cell.get(); // Get raw pointer from unique_ptr for the current cell
                    Cell *random_neighbour = neighbours[rand() % neighbours.size()];

                    // Remove both cells from their current grid positions BEFORE their coordinates are changed.
                    // This ensures they are removed from the correct buckets based on their pre-swap locations.
                    spatial_grid.remove(current_cell_ptr);
                    spatial_grid.remove(random_neighbour);

                    swapCells(current_cell_ptr, random_neighbour);

                    spatial_grid.insert(current_cell_ptr);
                    spatial_grid.insert(random_neighbour);
                }
            }
        }

        // 2. Handle blood vessel collisions
        for (auto &cell : cells) {
            cell->in_vessel_neighbourhood = false;
            for (const auto &vessel : blood_vessels){
                if (cell->inVesselNeighbour(vessel, VESSEL_DISTANCE_THRESHOLD)) {
                    cell->in_vessel_neighbourhood = true;
                    break;
                }
            }
        }

        float small_time_step = 0.1; // 1 min
        int num_inner_steps = 10; // big loop is 10 mins

        // 3. Check for collisions and resolve them
        for (int i = 0; i < num_inner_steps; ++i) {

            for (auto &cell : cells) {
                auto potential_collisions = spatial_grid.getPotentialCollisions(cell.get());
                for (Cell *other : potential_collisions) {
                    if (cell->collidesWith(*other)) {
                        cell->resolveCollision(*other, 1, small_time_step);
                    }
                }
            }
        }

        // 4. Cell division
        for (auto &cell : cells) {
            if (cell->should_divide()) {
                assert(!(cell->is_dead || cell->is_leaving) && "Cell is dead or leaving. Something is wrong.");

                auto &possible_types = LINEAGE_TREE.at(cell->getType());
                CellType new_type = possible_types[std::rand() % possible_types.size()];

                float new_radius = get_with_default(CELL_RADII, new_type, DEFAULT_CELL_RADII);
                float offset = cell->radius + new_radius;
                float angle = angle_dist(gen);

                // Calculate new position using polar coordinates
                float new_x = cell->x + offset * cos(angle);
                float new_y = cell->y + offset * sin(angle);

                // set velocity to 0
                new_cells.push_back(make_unique<Cell>(new_x, new_y, new_radius, new_type, cell->clone_id));

                if (cell->getType() != HSC) {
                    // two daughter cells, each w own fate decision. modify original cell.
                    if (possible_types.size() == 2) {
                        cell->updateType(possible_types[std::rand() % possible_types.size()]);
                    } else {
                        cell->updateType(new_type);
                    }
                }
            }
        }

        // 5. Check for cell death and leaving
        for (auto &cell : cells) {
            // Check for cell death (all cell types)
            if (cell->should_die()) {
                cell->is_dead = true;
                stats.total_deaths++;
                stats.deaths_by_type[cell->getType()]++;
            }
            // Check for cell leaving (only terminal cells)
            else if (cell->should_leave()) {
                cell->is_leaving = true;
                stats.total_leaving++;
                stats.leaving_by_type[cell->getType()]++;
            }
        }

        // Remove dead and leaving cells from the simulation
        auto isCellInactive = [](const unique_ptr<Cell> &cell) {
            return cell->is_dead || cell->is_leaving;
        };
        auto newEnd = remove_if(cells.begin(), cells.end(), isCellInactive);
        cells.erase(newEnd, cells.end());

        // Add new cells
        for (auto &new_cell : new_cells) {
            cells.push_back(std::move(new_cell));
        }
    }

    // Write cell data to file
    void writeCellDataToFile(int step_num)
    {
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

    void run(int steps) {
        for (int current_step = 0; current_step < steps; ++current_step) {
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

int main(int argc, char *argv[])
{
    // Default values
    float width = 50.0;
    float height = 50.0;
    int initial_cells = 20;
    int steps = 100;
    string sim_name = "bm_sim";
    int num_vessels = 10; // Default number of blood vessels
    bool cold_start = true;

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
            else if ((arg == "--cold_start" || arg == "-cs") && i + 1 < argc) {
                string cold_start_str = argv[++i];
                // transform(cold_start_str.begin(), cold_start_str.end(), cold_start_str.begin(), ::tolower);
                cold_start = (cold_start_str == "true" || cold_start_str == "1" || cold_start_str == "yes" || cold_start_str == "y");
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
                     << "  --help             Display this help message\n"
                     << "  --cold_start       Cold start the simulation (default: false)\n";
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
    cout << "  Cold start: " << cold_start << "\n\n";

    BoneMarrow model(width, height, initial_cells, sim_name, num_vessels, cold_start);
    model.run(steps);

    cout << "Simulation completed.\n";
    return 0;
}
