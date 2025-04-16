#ifndef CELL_CONFIG_H
#define CELL_CONFIG_H

#include <vector>
#include <map>
#include <cmath>
#include <tuple>

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
    Lymphocyte
};

// Cell type names for output
inline const char* getCellTypeName(CellType type) {
    static const char* names[] = {
        "HSC", "MPP1", "MPP2", "MPP3", "CMP", "CLP", "MEP", "GMP", 
        "Erythrocyte", "Granulocyte", "Lymphocyte"
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
    {CLP, {Lymphocyte}}
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
    {Lymphocyte, 3.0f}
};

// Map for cell division probabilities
inline const std::map<CellType, float> DIVISION_PROB = {
    {HSC, 0.01f},
    {MPP1, 0.02f},
    {MPP2, 0.03f},
    {MPP3, 0.04f},
    {CMP, 0.05f},
    {CLP, 0.05f},
    {MEP, 0.05f},
    {GMP, 0.05f},
    {Granulocyte, 0.0f},
    {Erythrocyte, 0.0f},
    {Lymphocyte, 0.0f}
};

// Map for cell leaving probabilities
inline const std::map<CellType, float> LEAVE_PROB = {
    {HSC, 0.0f},
    {MPP1, 0.0f},
    {MPP2, 0.0f},
    {MPP3, 0.0f},
    {CMP, 0.0f},
    {CLP, 0.0f},
    {MEP, 0.0f},
    {GMP, 0.0f},
    {Granulocyte, 0.1f},
    {Erythrocyte, 0.1f},
    {Lymphocyte, 0.1f}
};

// Map for cell death probabilities
inline const std::map<CellType, float> CELL_DEATH_PROB = {
    {HSC, 0.000f},
    {MPP1, 0.01f},
    {MPP2, 0.01f},
    {MPP3, 0.01f},
    {CMP, 0.01f},
    {CLP, 0.01f},
    {MEP, 0.01f},
    {GMP, 0.01f},
    {Granulocyte, 0.05f},
    {Erythrocyte, 0.05f},
    {Lymphocyte, 0.05f}
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
    {Lymphocyte, 0.3f}
}; 

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
    {Lymphocyte, {0, 255, 255}}     // Cyan
};

constexpr float MAX_SPEED = 1.0f;

#endif