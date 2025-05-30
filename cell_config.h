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
constexpr float VESSEL_LEAVING_MULTIPLIER = 20.0f;
constexpr float CXCL_DENSITY_PER_100_AREA = 6.0f;
constexpr int MAX_CELLS = 3500;
constexpr int SPATIAL_GRID_BLOCK_SIZE = 15;

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
    Bcell,
    STROMA,
};


inline const char* getCellTypeName(CellType type) {
    return magic_enum::enum_name(type).data();  // returns a string_view's c_str
}


// Define cell types in lineage tree with probabilities (sums to 1)
inline const std::unordered_map<CellType, std::vector<std::pair<CellType, float>>> LINEAGE_TREE = {
    {HSC, {{MPP1, 1.0f}}},
    {MPP1, {{MPP2, 1.0f}}},
    {MPP2, {{MPP3, 1.0f}}},
    {MPP3, {{MPP4, 1.0f}}},
    {MPP4, {{MPP5, 1.0f}}},
    {MPP5, {{CMP, 0.95f}, {CLP, 0.05f}}},
    {CMP, {{MEP, 0.6f}, {GMP, 0.4f}}},
    {MEP, {{Megakaryocyte, 0.1f}, {Erythroblast1, 0.9f}}},
    // {Megakaryocyte, {Platelet}}, // megas just release platelets into blood stream
    {Erythroblast1, {{Erythroblast2, 1.0f}}},
    {Erythroblast2, {{Erythroblast3, 1.0f}}},
    {Erythroblast3, {{Erythroblast4, 1.0f}}},
    {Erythroblast4, {{Erythroblast5, 1.0f}}},
    {Erythroblast5, {{Erythroblast6, 1.0f}}},
    {Erythroblast6, {{Erythroblast7, 1.0f}}},
    {Erythroblast7, {{Erythroblast8, 1.0f}}},
    {Erythroblast8, {{RBC, 1.0f}}},
    {GMP, {{Myeloblast1, 1.0f}}},
    {Myeloblast1, {{Myeloblast2, 1.0f}}},
    {Myeloblast2, {{Myeloblast3, 1.0f}}},
    {Myeloblast3, {{Myeloblast4, 1.0f}}},
    {Myeloblast4, {{Myeloblast5, 1.0f}}},
    {Myeloblast5, {{Myeloblast6, 1.0f}}},
    {Myeloblast6, {{Myeloblast7, 1.0f}}},
    {Myeloblast7, {{Myeloid, 1.0f}}},
    {CLP, {{Lymphocyte1, 1.0f}}},
    {Lymphocyte1, {{Lymphocyte2, 1.0f}}},
    {Lymphocyte2, {{Lymphocyte3, 1.0f}}},
    {Lymphocyte3, {{Lymphocyte4, 1.0f}}},
    {Lymphocyte4, {{Lymphocyte5, 1.0f}}},
    {Lymphocyte5, {{Lymphocyte6, 1.0f}}},
    {Lymphocyte6, {{Lymphocyte7, 1.0f}}},
    {Lymphocyte7, {{Lymphocyte8, 1.0f}}},
    {Lymphocyte8, {{Lymphocyte9, 1.0f}}},
    {Lymphocyte9, {{Lymphocyte10, 1.0f}}},
    {Lymphocyte10, {{Bcell, 1.0f}}},
    {STROMA, {{STROMA, 1.0f}}},
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
constexpr float DEFAULT_DIVISION_PROB = 1.0f / TIME_UNITS_PER_DAY;
// unordered_map for cell division probabilities
inline const std::unordered_map<CellType, float> DIVISION_PROB = {
    {HSC, 0.3f / TIME_UNITS_PER_DAY},
    {MPP1, 0.4f / TIME_UNITS_PER_DAY},
    {MPP2, 0.4f / TIME_UNITS_PER_DAY},
    {MPP3, 0.4f / TIME_UNITS_PER_DAY},
    {MPP4, 0.4f / TIME_UNITS_PER_DAY},
    {MPP5, 0.5f / TIME_UNITS_PER_DAY},
    {CMP, 1.5f / TIME_UNITS_PER_DAY},
    {MEP, 1.5f / TIME_UNITS_PER_DAY},
    {GMP, 1.5f / TIME_UNITS_PER_DAY},
    {CLP, 1.5f / TIME_UNITS_PER_DAY},
    {Bcell, 0.0f},
    {Myeloid, 0.0f},
    {RBC, 0.0f},
    {Megakaryocyte, 0.0f},
    {STROMA, 0.01f/TIME_UNITS_PER_DAY}
};

// unordered_map for cell leaving probabilities
inline const std::unordered_map<CellType, float> LEAVE_PROB = {
    {RBC, 1.17f/TIME_UNITS_PER_DAY},
    {Platelet, 0.0f},
    {Bcell, 0.65f/TIME_UNITS_PER_DAY},
    {Myeloid, 0.8f/TIME_UNITS_PER_DAY},
    {Megakaryocyte, 0.24f/TIME_UNITS_PER_DAY},
};

constexpr float DEFAULT_CELL_DEATH_PROB = 0.05f/TIME_UNITS_PER_DAY;

inline const std::unordered_map<CellType, float> CELL_DEATH_PROB = {
    {HSC, 0.0f},
    {MPP1, 0.01f/TIME_UNITS_PER_DAY},
    {MPP2, 0.01f/TIME_UNITS_PER_DAY},
    {MPP3, 0.01f/TIME_UNITS_PER_DAY},
    {MPP4, 0.01f/TIME_UNITS_PER_DAY},
    {MPP5, 0.01f/TIME_UNITS_PER_DAY},
    {CMP, 0.02f/TIME_UNITS_PER_DAY},
    {MEP, 0.02f/TIME_UNITS_PER_DAY},
    {GMP, 0.02f/TIME_UNITS_PER_DAY},
    {CLP, 0.02f/TIME_UNITS_PER_DAY},
    {STROMA, 0.01f/TIME_UNITS_PER_DAY},
};


constexpr float DEFAULT_CELL_MOTILITY = 0.2f;

inline const std::unordered_map<CellType, float> MOTILITY = {
    {STROMA, 0.0f},
};

inline std::unordered_map<CellType, float> SWAP_MOTILITY = {
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
    {MPP1, 1},
    {MPP2, 1},
    {MPP3, 1},
    {MPP4, 1},  
    {MPP5, 1},
    {CMP, 3},
    {CLP, 1},
    {MEP, 4},
    {GMP, 3},
    {Erythroblast1, 5},
    {Erythroblast2, 10},
    {Erythroblast3, 20},
    {Erythroblast4, 40},
    {Erythroblast5, 80},
    {Erythroblast6, 160},
    {Erythroblast7, 320},
    {Erythroblast8, 640},
    {RBC, 300},
    {Myeloblast1, 3},
    {Myeloblast2, 6},
    {Myeloblast3, 12},
    {Myeloblast4, 24},
    {Myeloblast5, 48},
    {Myeloblast6, 96},
    {Myeloblast7, 192},
    {Myeloid, 130},
    {Megakaryocyte, 3},
    {Lymphocyte1, 1},
    {Lymphocyte2, 1},
    {Lymphocyte3, 1},
    {Lymphocyte4, 3},
    {Lymphocyte5, 5},
    {Lymphocyte6, 10},
    {Lymphocyte7, 20},
    {Lymphocyte8, 40},
    {Lymphocyte9, 80},
    {Lymphocyte10, 160},
    {Bcell, 130},
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