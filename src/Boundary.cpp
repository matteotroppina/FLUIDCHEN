#include "Boundary.hpp"

Boundary::Boundary(std::vector<Cell *> cells) : _cells(cells) {}

void Boundary::applyFlux(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();
        if (cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j); // why not directly 0?
        }
        if (cell->is_border(border_position::LEFT)) {
            field.f(i - 1, j) = field.u(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.g(i, j) = field.v(i, j);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.g(i, j - 1) = field.v(i, j - 1);
        }
    }
}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::applyVelocity(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j) = -field.v(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
            field.v(i, j) = -field.v(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = 0;
            field.u(i, j) = -field.u(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j - 1) = 0;
            field.u(i, j) = -field.u(i, j - 1);
        }
    }
    //probably we should go through all cases ex. (border_position::RIGHT && border_position::TOP) etc...
}

void FixedWallBoundary::applyPressure(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1);
        }
    }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : Boundary(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::applyVelocity(Fields &field) {

    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j - 1) = 0;
            field.u(i, j) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(i, j - 1);
        }
        if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = 0;
            field.u(i, j) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(i, j + 1);
        }

        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] - field.v(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
            field.v(i, j) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] - field.v(i - 1, j);
        }

    }
}

void MovingWallBoundary::applyPressure(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1);
        }
    }
}

FixedVelocity::FixedVelocity(std::vector<Cell *> cells, std::map<int, double> inflow_u_velocity, std::map<int, double> inflow_v_velocity,
                                       std::map<int, double> wall_temperature)
    : Boundary(cells), _inflow_u_velocity(inflow_u_velocity), _inflow_v_velocity(inflow_v_velocity), _wall_temperature(wall_temperature) {}

void FixedVelocity::applyVelocity(Fields &field) {

    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        //field.u(i,j) = _inflow_u_velocity;
        //field.v(i,j) = _inflow_v_velocity;
    }
}

void FixedVelocity::applyPressure(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1);
        }
    }
    //Neumann condition, can we just leave it like this?
}

ZeroGradient::ZeroGradient(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void ZeroGradient::applyVelocity(Fields &field) {

    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        //field.u(i,j) = _outflow_u_velocity;
        //field.v(i,j) = _outflow_v_velocity;  //how are the outflow velocities??
        // It should be zero gradient fluid speed for outflow - Daniels
    }
}

void ZeroGradient::applyPressure(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1);
        }
    }
    //Neumann condition, can we just leave it like this?
}
