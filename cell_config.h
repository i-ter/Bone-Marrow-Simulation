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

// consdier the time step to be 10 mins
constexpr float TIME_UNITS_PER_DAY = 144.0f; 


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
    Lymphocyte8,
    Lymphocyte9,
    Lymphocyte10,
    Lymphocyte11,
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
    {MEP, {Megakaryocyte, Erythroblast1}},
    // {Megakaryocyte, {Platelet}}, // megas just release platelets into blood stream
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
    {Lymphocyte7, {Lymphocyte8}},
    {Lymphocyte8, {Lymphocyte9}},
    {Lymphocyte9, {Lymphocyte10}},
    {Lymphocyte10, {Lymphocyte11}},
    {Lymphocyte11, {Bcell}},
};

constexpr float DEFAULT_CELL_RADII = 4.0f;
// unordered_map for cell radii
inline const std::unordered_map<CellType, float> CELL_RADII = {
    {STROMA, 6.0f},
    {Megakaryocyte, 8.0f},
    // {Platelet, 2.0f},
    // {RBC, 2.0f}
};

// per day division probability divided by time step
constexpr float DEFAULT_DIVISION_PROB = .5f / TIME_UNITS_PER_DAY;
// unordered_map for cell division probabilities
inline const std::unordered_map<CellType, float> DIVISION_PROB = {
    {HSC, 0.3f / TIME_UNITS_PER_DAY},
    {MPP5, 2.0f / TIME_UNITS_PER_DAY},
    {CMP, 2.0f / TIME_UNITS_PER_DAY}, // RPP1
    {MEP, 2.0f / TIME_UNITS_PER_DAY}, // RPP2
    {Bcell, 0.0f},
    {Myeloid, 0.0f},
    {RBC, 0.0f},
    {Megakaryocyte, 0.0f},
    {STROMA, 0.0f}
};

// unordered_map for cell leaving probabilities
inline const std::unordered_map<CellType, float> LEAVE_PROB = {
    {RBC, 0.76f/TIME_UNITS_PER_DAY},
    {Platelet, 0.0f},
    {Bcell, 0.32f/TIME_UNITS_PER_DAY},
    {Myeloid, 0.92f/TIME_UNITS_PER_DAY},
    {Megakaryocyte, 0.24f/TIME_UNITS_PER_DAY},
};

// unordered_map for cell death probabilities
inline const std::unordered_map<CellType, float> CELL_DEATH_PROB = {
};


constexpr float DEFAULT_CELL_MOTILITY = 0.0f;

// unordered_map for cell motility
inline const std::unordered_map<CellType, float> MOTILITY = {
    {STROMA, 0.0f},
    {HSC, 0.5f},
    {MPP1, 0.5f},
    {MPP2, 0.5f},
    {MPP3, 0.5f},
    {CMP, 0.5f},
    {CLP, 0.5f},
    {MEP, 0.5f},
    {GMP, 0.5f},
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
    {Lymphocyte8, {128, 0, 128}},  // Purple
    {Lymphocyte9, {128, 0, 128}},  // Purple
    {Lymphocyte10, {128, 0, 128}},  // Purple
    {Lymphocyte11, {128, 0, 128}},  // Purple
    {Bcell, {75, 0, 130}},         // Indigo

    // Megakaryocyte/Platelet (Yellow/Gold)
    {Megakaryocyte, {255, 215, 0}}, // Gold
    {Platelet, {255, 255, 102}},    // Light Yellow

    // Stroma
    {STROMA, {169, 169, 169}},      // DarkGray
};

// 500 x 500 micron grid
inline std::unordered_map<CellType, float> INITIAL_CELL_NUMBERS = {
    {HSC, 1},
    {MPP1, 2},
    {MPP2, 4},
    {MPP3, 8},
    {MPP4, 16},  
    {MPP5, 32},
    {CMP, 10},
    {CLP, 5},
    {MEP, 10},
    {GMP, 10},
    {Erythroblast1, 15},
    {Erythroblast2, 30},
    {Erythroblast3, 60},
    {Erythroblast4, 120},
    {Erythroblast5, 240},
    {Erythroblast6, 480},
    {Erythroblast7, 960},
    {Erythroblast8, 1920},
    {RBC, 100},
    {Myeloblast1, 10},
    {Myeloblast2, 20},
    {Myeloblast3, 40},
    {Myeloblast4, 80},
    {Myeloblast5, 160},
    {Myeloblast6, 320},
    {Myeloblast7, 640},
    {Myeloid, 1500},
    {Megakaryocyte, 10},
    {Lymphocyte1, 6},
    {Lymphocyte2, 12},
    {Lymphocyte3, 24},
    {Lymphocyte4, 48},
    {Lymphocyte5, 96},
    {Lymphocyte6, 192},
    {Lymphocyte7, 384},
    {Lymphocyte8, 700},
    {Lymphocyte9, 800},
    {Lymphocyte10, 900},
    {Lymphocyte11, 1000},
    {Bcell, 768},
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