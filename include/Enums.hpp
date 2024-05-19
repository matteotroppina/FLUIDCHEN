#pragma once

// If no geometry file is provided in the input file, lid driven cavity case
// will run by default. In the Grid.cpp, geometry will be created following
// PGM convention
namespace LidDrivenCavity {
const int moving_wall_id = 8;
const int fixed_wall_id = 4;
const double wall_velocity = 1.0;
} // namespace LidDrivenCavity

// # Legend:
// # 0  Fluid
// # 1  Inflow
// # 2  Outflow
// # 3  Wall/Obstacle (adiabatic)
// # 4  Wall/Obstacle (hot)
// # 5  Wall/Obstacle (cold)

namespace GeometryIDs {
const int fluid = 0;
const int fixed_velocity = 1;
const int zero_gradient = 2;
const int fixed_wall = 3;
const int moving_wall = 8;
//const int hot_wall = 4;
//const int cold_wall = 5;
} // namespace GeometryIDs


enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};

enum class cell_type {
    FLUID,
    FIXED_WALL,
    MOVING_WALL,
    FIXED_VELOCITY,
    ZERO_GRADIENT,
    DEFAULT,
//  HOT_WALL,
//  COLD_WALL,
};
