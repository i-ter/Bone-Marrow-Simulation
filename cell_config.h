#ifndef CELL_CONFIG_H
#define CELL_CONFIG_H

#include <vector>
#include <map>
#include <cmath>
#include <tuple>
#include "magic_enum.hpp"
#include <iostream>

constexpr float MAX_SPEED = 1.0f;
constexpr float VESSEL_DISTANCE_THRESHOLD = 20.0f;
constexpr float VESSEL_LEAVING_MULTIPLIER = 10.0f;
constexpr float CXCL_DENSITY_PER_100_AREA = 6.0f;

enum CellType {
    HSC,
    MPP1,
    MPP2,
    MPP3,
    MPP4, 
    MPP5,
    CMP, // RPP1
    CLP, // RRP1
    MEP, // RPP2
    GMP, // RRP2
    Erythroblast1,
    Erythroblast2,
    Erythroblast3,
    Erythroblast4,
    Erythroblast5,
    Erythroblast6,
    Erythroblast7,
    Erythroblast8,
    RBC,
    Myeloblast1,
    Myeloblast2,
    Myeloblast3,
    Myeloblast4,
    Myeloblast5,
    Myeloblast6,
    Myeloblast7,
    Myeloid,
    Megakaryocyte,
    Platelet,
    Lymphocyte1,
    Lymphocyte2,
    Lymphocyte3,
    Lymphocyte4,
    Lymphocyte5,
    Lymphocyte6,
    Lymphocyte7,
    Bcell,
    STROMA,
};


inline const char* getCellTypeName(CellType type) {
    return magic_enum::enum_name(type).data();  // returns a string_view's c_str
}


// Define cell types in lineage tree
inline const std::map<CellType, std::vector<CellType>> LINEAGE_TREE = {
    {HSC, {MPP1}},
    {MPP1, {MPP2}},
    {MPP2, {MPP3}},
    {MPP3, {MPP4}},
    {MPP4, {MPP5}},
    {MPP5, {CMP, CLP}},
    {CMP, {MEP, GMP}},
    {MEP, {Megakaryocyte}},
    {Megakaryocyte, {Platelet}},
    {MEP, {Erythroblast1}},
    {Erythroblast1, {Erythroblast2}},
    {Erythroblast2, {Erythroblast3}},
    {Erythroblast3, {Erythroblast4}},
    {Erythroblast4, {Erythroblast5}},
    {Erythroblast5, {Erythroblast6}},
    {Erythroblast6, {Erythroblast7}},
    {Erythroblast7, {Erythroblast8}},
    {Erythroblast8, {RBC}},
    {GMP, {Myeloblast1}},
    {Myeloblast1, {Myeloblast2}},
    {Myeloblast2, {Myeloblast3}},
    {Myeloblast3, {Myeloblast4}},
    {Myeloblast4, {Myeloblast5}},
    {Myeloblast5, {Myeloblast6}},
    {Myeloblast6, {Myeloblast7}},
    {Myeloblast7, {Myeloid}},
    {CLP, {Lymphocyte1}},
    {Lymphocyte1, {Lymphocyte2}},
    {Lymphocyte2, {Lymphocyte3}},
    {Lymphocyte3, {Lymphocyte4}},
    {Lymphocyte4, {Lymphocyte5}},
    {Lymphocyte5, {Lymphocyte6}},
    {Lymphocyte6, {Lymphocyte7}},
    {Lymphocyte7, {Bcell}},
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
    {stroma, 6.0f}
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
    {stroma, 0.0f}
};

// Map for cell leaving probabilities
constexpr float FINAL_LEAVING_PROB = 0.001f;

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
    {stroma, 0.0f}
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
    {stroma, 0.000f}
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
    {stroma, 0.f}
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
    {stroma, {255, 255, 255}}     // White
};


inline const std::map<CellType, float> INITIAL_CELL_NUMBERS = {
    {HSC, 1},
    {MPP1, 5},
    {MPP2, 50},
    {MPP3, 100},
    {CMP, 200},
    {CLP, 200},
    {MEP, 200},
    {GMP, 200},
    {Granulocyte, 2000},
    {Erythrocyte, 2000},
    {Lymphocyte, 1000},
    {stroma, -1}
};


#endif