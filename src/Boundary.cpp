#include "Boundary.hpp"

Boundary::Boundary(std::vector<Cell *> cells) : _cells(cells) {}
void Boundary::applyFlux(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        // B_NW cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
            field.f(i - 1, j) = field.u(i - 1, j);
            field.g(i, j) = field.v(i, j);
            continue;
        }
        // B_SE cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j);
            field.g(i, j - 1) = field.v(i, j - 1);
            continue;
        }
        // B_NE cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j);
            field.g(i, j) = field.v(i, j);
            continue;
        }
        // B_SW cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
            field.f(i, j) = field.u(i - 1, j);
            field.g(i, j - 1) = field.v(i, j - 1);
            continue;
        }

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
    }
}
void Boundary::applyTemperature(Fields &field) {}

InnerObstacle::InnerObstacle(std::vector<Cell *> cells) : Boundary(cells) {}
void InnerObstacle::applyVelocity(Fields &field) {
    for (auto cell: _cells) {
        int i = cell->i();
        int j = cell->j();

        field.u(i,j) = 0;
        field.v(i,j) = 0;
    }
} // do nothing
void InnerObstacle::applyFlux(Fields &field) {}
void InnerObstacle::applyPressure(Fields &field) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : Boundary(cells) {}


FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, double wall_temperature)
    : Boundary(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::applyVelocity(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        // B_NW cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
            field.v(i, j) = 0;
            field.u(i, j) = -field.u(i, j+1);
            field.v(i, j-1) = - field.v(i-1, j-1);

            continue;
        }
        // B_SE cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j-1) = 0;
            field.u(i-1, j) = -field.u(i-1, j-1);
            field.v(i, j) = -field.v(i+1, j);

            continue;
        }
        // B_NE cell
        if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
            field.v(i, j) = 0;
            field.u(i-1, j) = -field.u(i-1, j+1);
            field.v(i, j-1) = - field.v(i+1, j-1);

            continue;
        }
        // B_SW cell
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
            field.u(i-1, j) = 0;
            field.v(i, j-1) = 0;
            field.u(i, j) = -field.u(i, j-1);
            field.v(i,j) = -field.v(i-1, j);

            continue;
        }

         // B_E cell
         if (cell->is_border(border_position::RIGHT)) {
             field.v(i, j) = -field.v(i + 1, j);
             field.u(i, j) = 0;
         }
         // B_W cell
         if (cell->is_border(border_position::LEFT)) {
             field.v(i, j) = -field.v(i - 1, j);
             field.u(i - 1, j) = 0;
         }
         // B_N cell
         if (cell->is_border(border_position::TOP)) {
             field.u(i, j) = -field.u(i, j + 1);
             field.v(i, j) = 0;
         }
         // B_S cell
         if (cell->is_border(border_position::BOTTOM)) {
             field.u(i, j) = -field.u(i, j - 1);
             field.v(i, j - 1) = 0;
         }

        // forbidden cells with two opposite borders or three boundaries
        if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::TOP)) {
            std::cerr << "there are forbidden cells with two opposite borders or three boundaries \n";
        }
        if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::RIGHT)) {
            std::cerr << "there are forbidden cells with two opposite borders or three boundaries \n";
        }
        // forbidden cells with obstacles only consisting of 1 cell
        if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::RIGHT) &&
            cell->is_border(border_position::TOP) && cell->is_border(border_position::BOTTOM)) {
            std::cerr << "there are forbidden cells with four boundaries \n";
        }
    }
}

void FixedWallBoundary::applyTemperature(Fields &field) {

    if (_wall_temperature == -1) { // Neumann
        for (auto cell : _cells) {
            int i = cell->i();
            int j = cell->j();

            if (cell->is_border(border_position::RIGHT)) {
                field.T(i, j) = field.T(i + 1, j);
            }

            if (cell->is_border(border_position::LEFT)) {
                field.T(i, j) = field.T(i - 1, j);
            }

            if (cell->is_border(border_position::TOP)) {
                field.T(i, j) = field.T(i, j + 1);
            }

            if (cell->is_border(border_position::BOTTOM)) {
                field.T(i, j) = field.T(i, j - 1);
            }
        }
    } else { // Dirichlet

        for (auto cell : _cells) {
            int i = cell->i();
            int j = cell->j();

            if (cell->is_border(border_position::BOTTOM)) {
                field.T(i, j) = 2 * _wall_temperature - field.T(i, j - 1);
            }
            if (cell->is_border(border_position::TOP)) {
                field.T(i, j) = 2 * _wall_temperature - field.T(i, j + 1);
            }

            if (cell->is_border(border_position::RIGHT)) {
                field.T(i, j) = 2 * _wall_temperature - field.T(i + 1, j);
            }

            if (cell->is_border(border_position::LEFT)) {
                field.T(i, j) = 2 * _wall_temperature - field.T(i - 1, j);
            }
        }
    }
}

void FixedWallBoundary::applyPressure(Fields &field) {
    for (auto cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        if (cell->is_border(border_position::RIGHT)) {
            field.p(i, j) = field.p(i + 1, j); // i = 0
        }

        if (cell->is_border(border_position::LEFT)) {
            field.p(i, j) = field.p(i - 1, j); // i = imax + 1
        }

        if (cell->is_border(border_position::TOP)) {
            field.p(i, j) = field.p(i, j + 1); // j = 0
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.p(i, j) = field.p(i, j - 1); // j = jmax + 1
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
            field.v(i, j) = field.v(i + 1, j);
        }

        if (cell->is_border(border_position::LEFT)) {
            field.u(i, j) = field.u(i - 2, j);
            field.u(i - 1, j) = field.u(i, j);
            field.v(i, j) = field.v(i - 1, j);
        }

        if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = field.v(i, j + 1);
            field.u(i, j) = field.u(i, j + 1);
        }

        if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j) = field.v(i, j - 2);
            field.v(i, j -1) = field.v(i, j);
            field.u(i, j) = field.u(i, j -1);
        }
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
            field.v(i, j) = 2 * _wall_velocity[GeometryIDs::moving_wall] - field.v(i + 1, j);
            field.u(i, j) = 0;
        }

        if (cell->is_border(border_position::LEFT)) {
            field.v(i, j) = 2 * _wall_velocity[GeometryIDs::moving_wall] - field.v(i - 1, j);
            field.u(i - 1, j) = 0;

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