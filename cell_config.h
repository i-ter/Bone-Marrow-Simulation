#ifndef CELL_CONFIG_H
#define CELL_CONFIG_H

#include <vector>
#include <map>
#include <cmath>
#include <tuple>

constexpr float MAX_SPEED = 1.0f;
constexpr float VESSEL_DISTANCE_THRESHOLD = 20.0f;
constexpr float VESSEL_LEAVING_MULTIPLIER = 5.0f;
constexpr float CXCL_DENSITY_PER_100_AREA = 6.0f;

enum CellType {
    HSC,
    MPP1,
    MPP2,
    MPP3,
    CMP, // RPP1
    CLP, // RRP1
    MEP, // RPP2
    GMP, // RRP2
    Erythrocyte,
    Granulocyte,
    Lymphocyte, 
    STROMAL,
};



// Cell type names for output
inline const char* getCellTypeName(CellType type) {
    static const char* names[] = {
        "HSC", "MPP1", "MPP2", "MPP3", "CMP", "CLP", "MEP", "GMP", 
        "Erythrocyte", "Granulocyte", "Lymphocyte", "Stromal"
    };
    return names[static_cast<int>(type)];
}

// Define cell types in lineage tree
inline const std::map<CellType, std::vector<CellType>> LINEAGE_TREE = {
    {HSC, {MPP1}},
    {MPP1, {MPP2}},
    {MPP2, {MPP3}},
    {MPP3, {CMP, CLP}},
    {CMP, {MEP, GMP}},
    {MEP, {Erythrocyte}},
    {GMP, {Granulocyte}},
    {CLP, {Lymphocyte}},
    {STROMAL, {}}
};

// Map for cell radii
inline const std::map<CellType, float> CELL_RADII = {
    {HSC, 4.0f},
    {MPP1, 4.0f},
    {MPP2, 4.0f},
    {MPP3, 4.0f},
    {CMP, 4.0f},
    {CLP, 4.0f},
    {MEP, 4.0f},
    {GMP, 4.0f},
    {Granulocyte, 3.0f},
    {Erythrocyte, 3.0f},
    {Lymphocyte, 3.0f},
    {STROMAL, 6.0f}
};

// Map for cell division probabilities
inline const std::map<CellType, float> DIVISION_PROB = {
    {HSC, 0.01f},
    {MPP1, 0.05f},
    {MPP2, 0.05f},
    {MPP3, 0.05f},
    {CMP, 0.1f},
    {CLP, 0.1f},
    {MEP, 0.1f},
    {GMP, 0.1f},
    {Granulocyte, 0.0f},
    {Erythrocyte, 0.0f},
    {Lymphocyte, 0.0f},
    {STROMAL, 0.0f}
};

// Map for cell leaving probabilities
constexpr float FINAL_LEAVING_PROB = 0.05f;

inline const std::map<CellType, float> LEAVE_PROB = {
    {HSC, 0.0f},
    {MPP1, 0.0f},
    {MPP2, 0.0f},
    {MPP3, 0.0f},
    {CMP, 0.0f},
    {CLP, 0.0f},
    {MEP, 0.0f},
    {GMP, 0.0f},
    {Granulocyte, FINAL_LEAVING_PROB},
    {Erythrocyte, FINAL_LEAVING_PROB},
    {Lymphocyte, FINAL_LEAVING_PROB},
    {STROMAL, 0.0f}
};

// Map for cell death probabilities
inline const std::map<CellType, float> CELL_DEATH_PROB = {
    {HSC, 0.000f},
    {MPP1, 0.000f},
    {MPP2, 0.000f},
    {MPP3, 0.000f},
    {CMP, 0.000f},
    {CLP, 0.000f},
    {MEP, 0.000f},
    {GMP, 0.000f},
    {Granulocyte, 0.000f},
    {Erythrocyte, 0.000f},
    {Lymphocyte, 0.000f},
    {STROMAL, 0.000f}
};

// Map for cell motility
inline const std::map<CellType, float> MOTILITY = {
    {HSC, 1.0f},
    {MPP1, 0.8f},
    {MPP2, 0.7f},
    {MPP3, 0.6f},
    {CMP, 0.5f},
    {CLP, 0.5f},
    {MEP, 0.5f},
    {GMP, 0.5f},
    {Granulocyte, 0.3f},
    {Erythrocyte, 0.3f},
    {Lymphocyte, 0.3f},
    {STROMAL, 0.f}
}; 

// Map for cell motility
// inline const std::map<CellType, float> MOTILITY = {
//     {HSC, 5.0f},
//     {MPP1, 4.0f},
//     {MPP2, 3.0f},
//     {MPP3, 3.0f},
//     {CMP, 3.0f},
//     {CLP, 3.0f},
//     {MEP, 3.0f},
//     {GMP, 3.0f},
//     {Granulocyte, 3.0f},
//     {Erythrocyte, 3.0f},
//     {Lymphocyte, 3.0f}
// }; 

inline const std::map<CellType, std::tuple<int, int, int>> CELL_COLORS = {
    {HSC, {255, 0, 0}},           // Red
    {MPP1, {255, 165, 0}},        // Orange
    {MPP2, {255, 165, 0}},        // Orange
    {MPP3, {255, 165, 0}},        // Orange
    {CMP, {0, 255, 0}},           // Green
    {CLP, {0, 255, 0}},           // Green
    {MEP, {0, 255, 0}},           // Green
    {GMP, {0, 255, 0}},           // Green
    {Granulocyte, {150, 200, 255}}, // Light Blue
    {Erythrocyte, {0, 150, 255}},     // Blue
    {Lymphocyte, {0, 255, 255}},     // Cyan
    {STROMAL, {255, 255, 255}}     // White
};


inline const std::map<CellType, float> INITIAL_CELL_NUMBERS = {
    {HSC, 1},
    {MPP1, 5},
    {MPP2, 50},
    {MPP3, 100},
    {CMP, 100},
    {CLP, 100},
    {MEP, 100},
    {GMP, 200},
    {Granulocyte, 2000},
    {Erythrocyte, 2000},
    {Lymphocyte, 1000},
    {STROMAL, -1}
};



#endif