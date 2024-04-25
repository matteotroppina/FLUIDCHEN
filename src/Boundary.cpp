#include "Boundary.hpp"

Boundary::Boundary(std::vector<Cell *> cells) : _cells(cells) {}

void Boundary::applyFlux(Fields &field) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::applyVelocity(Fields &field) {
    for (Cell* cell : _cells){
        field.u(cell->i(), cell->j()) = 0;
        field.v(cell->i(), cell->j()) = 0;
    };
}

void FixedWallBoundary::applyPressure(Fields &field) {}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : Boundary(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::applyVelocity(Fields &field) {
    for (Cell* cell : _cells){
        field.u(cell->i(), cell->j()) = 1;
        field.v(cell->i(), cell->j()) = 0;
    };
}

void MovingWallBoundary::applyPressure(Fields &field) {}
