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

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};


namespace GeometryIDs {
const int FLUID = 0;
const int FIXED_VELOCITY = 1;
const int ZERO_GRADIENT = 2;
const int FIXED_WALL = 3;
const int HOT_WALL = 4;
const int COLD_WALL = 5;
const int INNER_OBSTACLE = 6;
const int DEFAULT = 7;
const int MOVING_WALL = 8;
} // namespace GeometryIDs

enum class cell_type {
    FLUID,
    FIXED_VELOCITY,
    ZERO_GRADIENT,
    FIXED_WALL,
    HOT_WALL,
    COLD_WALL,
    INNER_OBSTACLE,
    DEFAULT,
    MOVING_WALL
};

enum DIRECTIONS {DOWN, UP, LEFT, RIGHT };