#include "Boundary.hpp"

Boundary::Boundary(std::vector<Cell *> cells) : _cells(cells) {}

void Boundary::applyFlux(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();
        if (cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j);
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

        // B_NW cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
            field.f(i - 1, j) = field.u(i - 1, j);
            field.g(i, j) = field.v(i, j);
        }
        // B_SE cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j);
            field.g(i, j - 1) = field.v(i, j - 1);
        }
        // B_NE cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j);
            field.g(i, j) = field.v(i, j);
        }
        // B_SW cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
            field.f(i, j) = field.u(i - 1, j);
            field.g(i, j - 1) = field.v(i, j - 1);
        }
    }
}

void ZeroGradientBoundary::applyFlux(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();
    }
}


InnerObstacle::InnerObstacle(std::vector<Cell *> cells) : Boundary(cells) {}

// TODO: there's a bug currently where the velocity is changing in the obstacle cells
//  currently I am setting the velocity to zero in the obstacle cells.
//  This works because the boundaries are set in order,
//  and the obstacle cells are set before the fixed wall cells
void InnerObstacle::applyVelocity(Fields &field) {} // do nothing
void InnerObstacle::applyFlux(Fields &field) {}
void InnerObstacle::applyPressure(Fields &field) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : Boundary(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::applyVelocity(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        // B_E cell
        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j) = -field.v(i + 1, j);
        }
        // B_W cell
        if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
            field.v(i, j) = -field.v(i - 1, j);
        }
        // B_N cell
        if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = 0;
            field.u(i, j) = -field.u(i, j + 1);
        }
        // B_S cell
        if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j - 1) = 0;
            field.u(i, j) = -field.u(i, j - 1);
        }

        // // B_NW cell
        // if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
        //     field.u(i,j) = -field.u(i,j+1);
        //     field.v(i,j) = 0;
        //     field.u(i-1,j) = 0;
        //     field.v(i,j-1) = -field.v(i-1,j-1);
        // }
        // // B_SE cell
        // if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
        //     field.u(i,j) = 0;
        //     field.v(i,j-1) = 0;
        //     field.u(i-1,j) = -field.u(i-1,j-1);
        //     field.v(i,j) = -field.v(i+1,j);
        // }
        // // B_NE cell
        // if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
        //     field.u(i,j) = 0;
        //     field.v(i,j) = 0;
        //     field.u(i-1,j) = -field.u(i-1,j+1);
        //     field.v(i,j-1) = -field.v(i+1,j-1);
        // }
        // // B_SW cell
        // if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
        //     field.u(i,j) = -field.u(i,j-1);
        //     field.v(i,j) = -field.v(i-1,j);
        //     field.u(i-1,j) = 0;
        //     field.v(i,j-1) = 0;
        // }

        // forbidden cells with two opposite borders or three boundaries
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::TOP)) {
            std::cout << "there are forbidden cells with two opposite borders or three boundaries \n";
        }
        if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::RIGHT)) {
            std::cout << "there are forbidden cells with two opposite borders or three boundaries \n";
        }
        // forbidden cells with obstacles only consisting of 1 cell
        if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::RIGHT) &&
            cell->is_border(border_position::TOP) && cell->is_border(border_position::BOTTOM)) {
            std::cout << "there are forbidden cells with four boundaries \n";
        }
    }
    // TODO: probably we should go through all cases ex. (border_position::RIGHT && border_position::TOP) etc... is it
    // maybe useful to change the order
    //  and do else if and else statements? Does the order of if clauses make a difference anyway??
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

        // B_NW cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
            field.p(i, j) = (field.p(i, j + 1) + field.p(i - 1, j)) / 2;
        }
        // B_SE cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = (field.p(i + 1, j) + field.p(i, j - 1)) / 2;
        }
        // B_SW cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
            field.p(i, j) = (field.p(i - 1, j) + field.p(i, j - 1)) / 2;
        }
        // B_NE cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = (field.p(i, j + 1) + field.p(i + 1, j)) / 2;
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
            field.u(i, j) = 2 * _wall_velocity[GeometryIDs::moving_wall] - field.u(i, j - 1);
            field.v(i, j - 1) = 0;
        }
        if (cell->is_border(border_position::TOP)) {
            field.u(i, j) = 2 * _wall_velocity[GeometryIDs::moving_wall] - field.u(i, j + 1);
            field.v(i, j) = 0;
        }

        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j) = 2 * _wall_velocity[GeometryIDs::moving_wall] - field.v(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
            field.v(i, j) = 2 * _wall_velocity[GeometryIDs::moving_wall] - field.v(i - 1, j);
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

FixedVelocityBoundary::FixedVelocityBoundary(std::vector<Cell *> cells, double inflow_u_velocity,
                                             double inflow_v_velocity)
    : Boundary(cells) {
    _inflow_u_velocity.insert(std::pair<int, double>(GeometryIDs::fixed_velocity, inflow_u_velocity));
    _inflow_v_velocity.insert(std::pair<int, double>(GeometryIDs::fixed_velocity, inflow_v_velocity));
}

FixedVelocityBoundary::FixedVelocityBoundary(std::vector<Cell *> cells, std::map<int, double> inflow_u_velocity,
                                             std::map<int, double> inflow_v_velocity,
                                             std::map<int, double> wall_temperature)
    : Boundary(cells), _inflow_u_velocity(inflow_u_velocity), _inflow_v_velocity(inflow_v_velocity),
      _wall_temperature(wall_temperature) {}

void FixedVelocityBoundary::applyVelocity(Fields &field) {

    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::BOTTOM)) {
            field.u(i, j) = 2 * _inflow_u_velocity[GeometryIDs::fixed_velocity] - field.u(i, j - 1);
            field.v(i, j - 1) = _inflow_v_velocity[GeometryIDs::fixed_velocity];
        }
        if (cell->is_border(border_position::TOP)) {
            field.u(i, j) = 2 * _inflow_u_velocity[GeometryIDs::fixed_velocity] - field.u(i, j + 1);
            field.v(i, j) = _inflow_v_velocity[GeometryIDs::fixed_velocity];
        }

        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = _inflow_u_velocity[GeometryIDs::fixed_velocity];
            field.v(i, j) = 2 * _inflow_v_velocity[GeometryIDs::fixed_velocity] - field.v(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = _inflow_u_velocity[GeometryIDs::fixed_velocity];
            field.v(i, j) = 2 * _inflow_v_velocity[GeometryIDs::fixed_velocity] - field.v(i - 1, j);
        }
    }
}

void FixedVelocityBoundary::applyPressure(Fields &field) {
    // Neumann condition
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
    // Neumann condition, can we just leave it like this?
}
ZeroGradientBoundary::ZeroGradientBoundary(std::vector<Cell *> cells) : Boundary(cells) {}
ZeroGradientBoundary::ZeroGradientBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void ZeroGradientBoundary::applyVelocity(Fields &field) {
    // Neumann condition !!! CHange needed
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = field.u(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.u(i, j) = field.u(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = field.v(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j) = field.v(i, j - 1);
            // field.u(i, j + 1) = field.u(i, j);
        }
        //        //field.u(i,j) = _outflow_u_velocity;
        //        //field.v(i,j) = _outflow_v_velocity;  //how are the outflow velocities??
        //        // It should be zero gradient fluid speed for outflow - Daniels
    }
}

void ZeroGradientBoundary::applyPressure(Fields &field) {
    // Dirichlet condition pressure on boundary = 0
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = -field.p(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = -field.p(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = -field.p(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = -field.p(i, j - 1);
        }
    }
}