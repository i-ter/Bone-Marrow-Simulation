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
inline const std::unordered_map<CellType, std::vector<CellType>> LINEAGE_TREE = {
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

constexpr float DEFAULT_CELL_RADII = 4.0f;
// unordered_map for cell radii
inline const std::unordered_map<CellType, float> CELL_RADII = {
    {STROMA, 6.0f},
    {Megakaryocyte, 10.0f},
    {Platelet, 1.0f},
    {RBC, 1.0f}
};


constexpr float DEFAULT_DIVISION_PROB = 0.01f;
// unordered_map for cell division probabilities
inline const std::unordered_map<CellType, float> DIVISION_PROB = {
    {HSC, 0.01f},
    {MPP1, 0.05f},
    {MPP2, 0.05f},
    {MPP3, 0.05f},
    {MPP4, 0.05f}, 
    {MPP5, 0.05f},
    {CMP, 0.1f}, // RPP1
    {CLP, 0.1f}, // RRP1
    {MEP, 0.1f}, // RPP2
    {GMP, 0.1f}, // RRP2
    {Erythroblast8, 0.1f},
    {Myeloblast7, 0.1f},
    {Megakaryocyte, 0.1f},
    {Lymphocyte7, 0.1f},
    {Bcell, 0.0f},
    {Myeloid, 0.0f},
    {RBC, 0.0f},
    {Platelet, 0.0f},
    {STROMA, 0.0f}
};

// unordered_map for cell leaving probabilities
constexpr float FINAL_LEAVING_PROB = 0.01f;

inline const std::unordered_map<CellType, float> LEAVE_PROB = {
    {RBC, FINAL_LEAVING_PROB},
    {Platelet, FINAL_LEAVING_PROB},
    {Bcell, FINAL_LEAVING_PROB},
    {Myeloid, FINAL_LEAVING_PROB},
    {Megakaryocyte, FINAL_LEAVING_PROB},
};

// unordered_map for cell death probabilities
inline const std::unordered_map<CellType, float> CELL_DEATH_PROB = {
};


constexpr float DEFAULT_CELL_MOTILITY = 0.0f;

// unordered_map for cell motility
inline const std::unordered_map<CellType, float> MOTILITY = {
//     {HSC, 0.5f},
//     {MPP1, 0.5f},
//     {MPP2, 0.5f},
//     {MPP3, 0.5f},
//     {CMP, 0.5f},
//     {CLP, 0.5f},
//     {MEP, 0.5f},
//     {GMP, 0.5f},
};

inline const std::unordered_map<CellType, std::tuple<int, int, int>> CELL_COLORS = {
    // HSC/MPP lineage (Red-Orange)
    {HSC, {220, 20, 60}},         // Crimson
    {MPP1, {255, 69, 0}},         // OrangeRed
    {MPP2, {255, 99, 71}},        // Tomato
    {MPP3, {255, 140, 0}},        // DarkOrange
    {MPP4, {255, 165, 0}},        // Orange
    {MPP5, {255, 200, 0}},        // Lighter Orange

    // CMP/MEP/Erythroblast/RBC (Blue family)
    {CMP, {30, 144, 255}},        // DodgerBlue
    {MEP, {65, 105, 225}},        // RoyalBlue
    {Erythroblast1, {100, 149, 237}}, // CornflowerBlue
    {Erythroblast2, {70, 130, 180}},  // SteelBlue
    {Erythroblast3, {135, 206, 250}}, // LightSkyBlue
    {Erythroblast4, {176, 224, 230}}, // PowderBlue
    {Erythroblast5, {173, 216, 230}}, // LightBlue
    {Erythroblast6, {135, 206, 235}}, // SkyBlue
    {Erythroblast7, {0, 191, 255}},   // DeepSkyBlue
    {Erythroblast8, {0, 0, 255}},     // Blue
    {RBC, {25, 25, 112}},             // MidnightBlue

    // GMP/Myeloblast/Myeloid (Green family)
    {GMP, {34, 139, 34}},         // ForestGreen
    {Myeloblast1, {60, 179, 113}},// MediumSeaGreen
    {Myeloblast2, {46, 139, 87}}, // SeaGreen
    {Myeloblast3, {0, 128, 0}},   // Green
    {Myeloblast4, {0, 255, 127}}, // SpringGreen
    {Myeloblast5, {144, 238, 144}}, // LightGreen
    {Myeloblast6, {152, 251, 152}}, // PaleGreen
    {Myeloblast7, {50, 205, 50}}, // LimeGreen
    {Myeloid, {0, 100, 0}},       // DarkGreen

    // CLP/Lymphocyte/Bcell (Purple family)
    {CLP, {138, 43, 226}},        // BlueViolet
    {Lymphocyte1, {148, 0, 211}}, // DarkViolet
    {Lymphocyte2, {186, 85, 211}},// MediumOrchid
    {Lymphocyte3, {221, 160, 221}},// Plum
    {Lymphocyte4, {218, 112, 214}},// Orchid
    {Lymphocyte5, {199, 21, 133}}, // MediumVioletRed
    {Lymphocyte6, {153, 50, 204}}, // DarkOrchid
    {Lymphocyte7, {128, 0, 128}},  // Purple
    {Bcell, {75, 0, 130}},         // Indigo

    // Megakaryocyte/Platelet (Yellow/Gold)
    {Megakaryocyte, {255, 215, 0}}, // Gold
    {Platelet, {255, 255, 102}},    // Light Yellow

    // Stroma
    {STROMA, {169, 169, 169}},      // DarkGray
};


inline const std::unordered_map<CellType, float> INITIAL_CELL_NUMBERS = {
    {HSC, 1},
    {MPP1, 5},
    {MPP2, 5},
    {MPP3, 5},
    {MPP4, 5},  
    {MPP5, 5},
    {CMP, 30},
    {CLP, 30},
    {MEP, 30},
    {GMP, 30},
    {Erythroblast1, 30},
    {Erythroblast2, 60},
    {Erythroblast3, 120},
    {Erythroblast4, 180},
    {Erythroblast5, 320},
    {Erythroblast6, 340},
    {Erythroblast7, 200},
    {Erythroblast8, 100},
    {RBC, 100},
    {Myeloblast1, 30},
    {Myeloblast2, 60},
    {Myeloblast3, 120},
    {Myeloblast4, 180},
    {Myeloblast5, 320},
    {Myeloblast6, 240},
    {Myeloblast7, 100},
    {Myeloid, 500},
    {Megakaryocyte, 30},
    {Lymphocyte1, 30},
    {Lymphocyte2, 60},
    {Lymphocyte3, 120},
    {Lymphocyte4, 180},
    {Lymphocyte5, 320},
    {Lymphocyte6, 240},
    {Lymphocyte7, 100},
    {Bcell, 100},
};



template<typename Value> // basically int or float
inline Value get_with_default(const std::unordered_map<CellType, Value>& m, CellType k, Value def) {
    try {
        return m.at(k);
    } catch (const std::out_of_range&) {
        return def;
    }
}

template<typename Value>
inline Value get_with_zero(const std::unordered_map<CellType, Value>& m, CellType k) {
    return get_with_default(m, k, Value{0});
}


#endif