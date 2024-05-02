#include "Boundary.hpp"

Boundary::Boundary(std::vector<Cell *> cells) : _cells(cells) {}

void Boundary::applyFlux(Fields &field) {
    for (auto cell: _cells){
        int i = cell -> i();
        int j = cell -> j();
        if (cell -> is_border(border_position::RIGHT)) {
            field.f(i,j) = field.u(i,j); //why not directly 0?
        }
        if (cell -> is_border(border_position::LEFT)) {
            field.f(i,j) = field.u(i,j);
        }
        if (cell -> is_border(border_position::TOP)) {
            field.g(i,j) = field.v(i,j);
        } 
        if (cell -> is_border(border_position::BOTTOM)) {
            field.g(i,j) = field.v(i,j);
        }        
    }
}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}



std::string cellTypeToString(cell_type type) {
    std::map<cell_type, std::string> cellTypeMap = {
        {cell_type::FLUID, "FLUID"},
        {cell_type::FIXED_WALL, "FIXED_WALL"},
        {cell_type::MOVING_WALL, "MOVING_WALL"},
        {cell_type::DEFAULT, "DEFAULT"}
    };

    return cellTypeMap[type];
}

void FixedWallBoundary::applyVelocity(Fields &field) {
    for (auto cell: _cells){
        int i = cell -> i();
        int j = cell -> j();

        if (cell -> is_border(border_position::RIGHT)) {
            field.u(i,j) = 0;
            field.v(i,j) = - field.v(i+1,j);
        }

        if (cell -> is_border(border_position::LEFT)) {
            field.u(i,j) = 0;
            field.v(i,j) = - field.v(i-1,j);
        }

        if (cell -> is_border(border_position::TOP)) {
            field.v(i,j) = 0;
            field.u(i,j) = - field.u(i,j+1);
        }

        if (cell -> is_border(border_position::BOTTOM)) {
            field.v(i,j) = 0;
            field.u(i,j) = - field.u(i,j-1);
        }
    }
}


void FixedWallBoundary::applyPressure(Fields &field) {}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : Boundary(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::applyVelocity(Fields &field) {

        for (auto cell: _cells){
            int i = cell -> i();
            int j = cell -> j();

            if (cell -> is_border(border_position::BOTTOM)) {
                field.v(i,j) = 0;
                field.u(i,j) = 1;
            }
        }

}


void MovingWallBoundary::applyPressure(Fields &field) {}
